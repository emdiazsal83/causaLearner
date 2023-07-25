# demo of causality package

remove(list=ls())
#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag5-ME32-Add_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v3_cc"
hs_cmfm_ob_version <- "v2_cc"


#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)
print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

library(FNN) #  (in learn.qhsic, myGpCreate)
# FNN sometimes has a problem: cant find C functions perhaps because of some
# unloading of dll libraries: see https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16446


#########################################################################################################################################3
# Load/simulate data on which we will evaluate performance of causal learners
#########################################################################################################################################

count.graphs(type="all-dags", 5)
allDags5 <- genAllDAGS(5)
# pm <- proc.time()
# sizeME <- apply(allDags5, 3, function(dagMat) dim(getMarkovEquivClass(dagMat))[3])
# proc.time() - pm # 5 mins
# save(sizeME, file="./pkg_causaLearner/experiments/CMEM/results/sizeME_allDags5.RData")
load(file="./pkg_causaLearner/experiments/CMEM/results/sizeME_allDags5.RData")
table(sizeME)
indxME32 <- which(sizeME == 32)[1]

q <- 100
set.seed(4)
p <- 5
ps <- rep(p, q)
ns <- rep(100, q)
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)
dagMat <- allDags5[,,indxME32]
rownames(dagMat) <- colnames(dagMat) <- names(nodes)
dagMatME <- getMarkovEquivClass(dagMat)
dim(dagMatME)


# defAttrs <- getDefaultAttrs()
# defAttrs$node$fontsize <- "8"
# plot(as(dagMat, "graphNEL"), attrs=defAttrs)
# defAttrs$node$fontsize <- "6"
# par(mfrow=c(6,6))
# for(i in 1:dim(dagMatME)[3]) plot(as(dagMatME[,,i], "graphNEL"), attrs=defAttrs)
# plot(as(dagMatME[,,2], "graphNEL"))



 # set.seed(5)
 # dataList <- simRandSEMs(q, ps, ns, nodess, sigma=1, dagMat=dagMat, markovEquiv=T)
 # pairs(dataList$xs[[3]])
 # ggpairs(as.data.frame(dataList$xs[[13]]),aes(alpha = 0.4))
 # save(dataList, file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))



block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#block <- 1

pm <- proc.time()
load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs

# weights for each data set
ws <- rep(1, q)
ws <- ws/sum(ws)
names(ws) <- dataList$names

#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm") #, "stdr"
approxDagSetMethodList <- c("oracleME1")
#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")

hypScorerList <- c(cmem_hypScorer_pack_none_1, "rnd1") #[c(5,6)]

#hypScorerList <- cmem_hypScorer_pack_lambda_1
#hypScorerList <- cmem_hypScorer_pack_kernParsX_1
#hypScorerList <- cmem_hypScorer_pack_kernParsXY_1
       



recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)



#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=6, plot=TRUE)
proc.time()-pm 

# 

# seems like its best to parallelize on applyLearnRecipes level and on parameters but not on folsds

save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName,"_",hs_cmem_ob_version,"_", block, ".RData", sep=""))

print("finish")
q("no")

blockFiles <- dir("./pkg_causaLearner/experiments/blocks/")
indx <- grep(experimentName, blockFiles)
blockFiles <- blockFiles[indx]
blockNumber <- as.numeric(sapply(strsplit(sapply(strsplit(blockFiles, "_"), function(el) el[[5]]), "\\."), function(el) el[[1]]))
(blockFiles <- blockFiles[order(blockNumber)])

scores <- lapply(blockFiles, function(el){
  # el <- blockFiles[1]
  load(file=paste("./pkg_causaLearner/experiments/blocks/", el, sep=""))
  res <- scores
})
names(scores) <- dataList$names


# For running all on my pc

pm <- proc.time()
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=1, plot=FALSE)
proc.time() - pm 

# 12.4 mins for none_1      with 3-1-2 cores on my pc with 1-1 in CV
# for kernParsX_1 with 1-1-6 cores on my pc with 1-1 in CV

save(scores, file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName,"_", hs_cmem_ob_version, ".RData", sep=""))

unlink(paste("./pkg_causaLearner/experiments/blocks/",blockFiles, sep=""))





#########################################################################################################################################*
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################*

load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName,"_",hs_cmem_ob_version, ".RData", sep=""))


# clean all scores
# first lets see how many Inf, -Inf, NA or NULL or NaN
aux <- deliverFuncToScore_recipeTree(scores, func ="printErrors", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))

# actually we only have a problem when there is a NA, nan, or both all hypothesis are non-finite with the same sign 
# lets check first nas and nans

aux <- deliverFuncToScore_recipeTree(scores, func ="printErrors2", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))



# lets see which scores are the main problem as far as na's and nan's

aux <- deliverFuncToScore_recipeTree(scores, func ="tableErrors2")

errors <- unwrapScoreFunc(scoresFunc=aux, dataList=dataList, ws)

cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)/length(x)*100/2)

# what about non decision problem -> all hypothesis scores equal

aux <- deliverFuncToScore_recipeTree(scores, func ="printNoDecision", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))

aux <- deliverFuncToScore_recipeTree(scores, func ="tableNoDecision")

noDecisions <- unwrapScoreFunc(scoresFunc=aux, dataList=dataList, ws)


cast(noDecisions, score~. ,value="value", fun.aggregate=function(x) sum(x, na.rm=T)/2)



#########################################################################################################################################3
# Obtain and visualize performance measures 
#########################################################################################################################################


###############################################################
# Obtain dag-classification measures
###############################################################

msrs <- deliverFuncToScore_recipeTree(scores, func="measuresDagDistEdge", dags=dataList$dags, ppTab=1,
                                      getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                      getDataLevelPars=getHypots)


# analyze these measures

measures <- unwrapScoreFunc(scoresFunc=msrs, dataList=dataList, ws)


msrss <- c("tss","edgeAUC","ccr")
scorses <-   c("KCDC","rnd")
indx <- with(measures, which(measure %in% msrss & score %in% scorses))
p <- ggplot(measures[indx,])
p <- p + geom_boxplot(aes(x=recipeFull, y=value, color=recipeFull))
p <- p + facet_wrap(~measure, scales="free")
p

# Total data sets per recipe and p
dts <- c("norm")
msrss <- c("ccr","edgeAUC","totEdgeD") #"edgeD","fnr","fpr","msr","nonEdgeD","npp", "npp", "ppp","sens", "spec", "totEdgeD", "tss"
scorses <-   c("KCDC", "KCDCrel", "rnd")
indx <- with(measures, which(measure %in% msrss  & dt %in%dts))
cast(measures[indx,], recipeFull~measure, value="value", margins="grand_col")
# mean stats
tab <- cast(measures[indx,], recipeFull~measure, value="value", fun.aggregate="mean", margins="grand_col")
tab
tab[1:13,]
tab[14:26,]
summary(tab)

tab[order(tab$ccr, decreasing=TRUE),]


#############################################
# first choose parameters for aggregating 
# scores of sets of dags into two classes
#############################################


groupFuncs <- c("agg_truevsRest")
groupFuncPars <- list(list())

#############################################
# Aggregate scores according to previous 
# parameters
#############################################

# first aggregate by min score then rank by rank2func associated to the score

aggScores <- deliverFuncToScore_recipeTree(scores, func="aggregateScores", dags=dataList$dags, 
                                           getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                           getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1)


#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")




n.boot <- 500

statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
statsDB <- statsArrayToDB(statsArr)


msrss <- c("ccr")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
tab

tab[order(tab$KCDC, decreasing=T),]
tab[order(tab$KCMC, decreasing=T),]
tab[order(tab$KCNSC, decreasing=T),]
tab[order(tab$KCRDC, decreasing=T),]
tab[order(tab$KCSC, decreasing=T),]

msrss <- c("auc")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
tab

msrss <- c("ccr")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("tss")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p




