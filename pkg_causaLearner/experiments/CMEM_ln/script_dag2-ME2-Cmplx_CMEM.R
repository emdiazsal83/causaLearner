# demo of causality package

remove(list=ls())
#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2-ME2-Cmplx_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v3_cc"
hs_cmfm_ob_version <- "v2_cc"

#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"

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


q <- 100
set.seed(4)
p <- 2
ps <- rep(p, q)
ns <- rep(100, q)
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)

dagMat <- matrix(c(0,0,1,0), 2, 2)
rownames(dagMat) <- colnames(dagMat) <- names(nodes)


# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
# pm <- proc.time()
# set.seed(5)
# dataList <- simRandSEMs(q, ps, ns, nodess, sigma=1, sigmaErr=1 ,markovEquiv=T, dagMat=dagMat, geU=function(y, nois, scale, constant) y)
# proc.time() - pm # 6.2 secs
# plotPairsList(dataList, sizePts=0.2, alpha=0.1)
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

hypScorerList <- cmem_hypScorer_pack_none_1
#hypScorerList <- cmem_hypScorer_pack_lambda_1
#hypScorerList <- cmem_hypScorer_pack_kernParsX_1
#hypScorerList <- cmem_hypScorer_pack_kernParsXY_1
       



recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)




#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE)
proc.time()-pm 

# 14 secs for none_1 on my pc


save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName, "_", block, ".RData", sep=""))

print("finish")
q("no")

blockFiles <- dir("./pkg_causaLearner/experiments/blocks/")
indx <- grep(experimentName, blockFiles)
blockFiles <- blockFiles[indx]
blockNumber <- as.numeric(sapply(strsplit(sapply(strsplit(blockFiles, "_"), function(el) el[[3]]), "\\."), function(el) el[[1]]))
(blockFiles <- blockFiles[order(blockNumber)])

scores <- lapply(blockFiles, function(el){
  # el <- blockFiles[1]
  load(file=paste("./pkg_causaLearner/experiments/blocks/", el, sep=""))
  res <- scores
})
names(scores) <- dataList$names

# For running all on my pc

pm <- proc.time()
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=6, plot=FALSE)
proc.time() - pm 

# 8.6 mins for none_1      with 1-1-6 cores on my pc with 1-1 in CV
# for kernParsX_1 with 1-1-6 cores on my pc with 1-1 in CV

save("scores", file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, "_", hs_cmem_ob_version,  ".RData", sep=""))

unlink(paste("./pkg_causaLearner/experiments/blocks/", blockFiles, sep=""))




#########################################################################################################################################3
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################

 

load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, "_", hs_cmem_ob_version,  ".RData", sep=""))

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


dagIds <- getHypID(sapply(dataList$dags, function(el) el, simplify="array"))

names(scores[[1]][["norm"]][["oracleME1"]])
scores[[2]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]]

rownames(scores[[1]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none_1_1_1"]]) %in% dagIds$id

i <- 5
dagIds$id[i]
scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]]
aggScores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]]
rnkDecs[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]]

which.max(scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none_1_1_1"]][,"KCMC"])



#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=dataList, ws)
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

i <- 5
dagIds$id[i]
which.min(scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]][,"KCMC"])
names(which.min(scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]][,"KCMC"]))



pred <- sapply(1:100, function(i) names(which.min(scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_laplace1_cn_L2_none"]][,"KCMC"])))
obs <- as.character(dagIds$id)
pred[i]; obs[i]
tab <- table(pred, obs)
(tab[1,1]+tab[2,2])/sum(tab)


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








