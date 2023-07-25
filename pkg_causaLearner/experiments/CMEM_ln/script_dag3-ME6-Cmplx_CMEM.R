# demo of causality package

remove(list=ls())
#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag3-ME6-Cmplx_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v3_cc"
hs_cmfm_ob_version <- "v1_ncnc"


#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"

dir(repos)
setwd(repos)
print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)



#########################################################################################################################################3
# Load/simulate data on which we will evaluate performance of causal learners
#########################################################################################################################################


q <- 100
set.seed(4)
p <- 3
(ps <- rep(p, q))
(ns <- rep(100, q))
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)
#dagMat <- matrix(c(0,0,0,1,0,0,0,1,0), 3, 3)
dagMat <- matrix(c(0,0,0,1,0,0,1,1,0), 3, 3)
rownames(dagMat) <- colnames(dagMat) <- names(nodes)
dagMatME <- getMarkovEquivClass(dagMat)
dim(dagMatME)



# set.seed(5)
# dataList <- simRandSEMs(q, ps, ns, nodess, sigma=1, sigmaErr=1, dagMat=dagMat, markovEquiv=T, geU=function(y, nois, scale, constant) y)
# ggpairs(as.data.frame(dataList$xs[[1]]),aes(alpha = 0.4))
# save(dataList, file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
# plotDataList(dataList)



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

hypScorerList <- c(cmem_hypScorer_pack_none_1, "rnd1")
#hypScorerList <- c(cmem_hypScorer_pack_none_1[c(5,6)], "rnd1")
#hypScorerList <- cmem_hypScorer_pack_lambda_1
#hypScorerList <- c(cmem_hypScorer_pack_kernParsX_1, "rnd1")
#hypScorerList <- cmem_hypScorer_pack_kernParsXY_1
       



recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)



#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE)
proc.time()-pm 

# 

# seems like its best to parallelize on applyLearnRecipes level and on parameters but not on folsds

save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName,"_",hs_cmem_ob_version, "_", block, ".RData", sep=""))

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
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=6, plot=FALSE)
proc.time() - pm 

#  45 mins for none_1      with 1-1-6 cores on my pc with 1-1 in CV
# for kernParsX_1 with 1-1-6 cores on my pc with 1-1 in CV

# scores nn - no centering, KCDC, KCDCrel, log, quad, rbf
# scores_v1 cc- both kernel matrices centered, KCDC, KCMC, etc, quad, laplace, poly, rbf, etc
# scores_v2 cc- both kernel matrices centered, KCDC, KCMC, etc, log, quad, rbf, etc

scores_v2 <- scores
scores_v3 <- scores

save(scores, file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName,"_",hs_cmem_ob_version,".RData", sep=""))


unlink(paste("./pkg_causaLearner/experiments/blocks/",blockFiles, sep=""))





#########################################################################################################################################*
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################*


load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName,"_",hs_cmem_ob_version,".RData", sep=""))

length(scores)
length(scores_v2)
length(scores_v3)
names(scores[[1]])
names(scores_v2[[1]])
names(scores_v3[[1]])
names(scores[[1]][["norm"]][["oracleME1"]])
names(scores_v2[[1]][["norm"]][["oracleME1"]])
names(scores_v3[[1]][["norm"]][["oracleME1"]])

scores[[1]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_quad_nn_L2_none"]]
scores_v2[[1]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_quad_cc_L2_none"]]
scores_v3[[1]][["norm"]][["oracleME1"]][["hs_cmem_kern_quadANDkern_quad_cc_L2_none"]]
                                          

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


cast(noDecisions, score~hypSc ,value="value", fun.aggregate=function(x) sum(x, na.rm=T)/2)



#########################################################################################################################################3
# Obtain and visualize performance measures 
#########################################################################################################################################


names(scores[[1]][[1]][[1]])
scores <- lapply(1:100, function(i){
  # i <- 1
  # names(scores[[i]])
  #print(i)
  scrs <- lapply(1:length(scores[[i]]), function(j){
    # j <- 1
    # names(scores[[i]][[j]])
    scrs <- lapply(1:length(scores[[i]][[j]]), function(k){
      # k <- 1
      # names(scores[[i]][[j]][[k]])
      # get rid of rand1
      scrs <- lapply(1:(length(scores[[i]][[j]][[k]])-1), function(l){
        # l <- 1
        scores[[i]][[j]][[k]][[l]] #[,c("KCDC","KCRDC","KCMC")]
      })
      names(scrs) <- names(scores[[i]][[j]][[k]])[1:(length(scores[[i]][[j]][[k]])-1)]
      return(scrs)
    })
    names(scrs) <- names(scores[[i]][[j]])
    return(scrs)
  })
  names(scrs) <- names(scores[[i]])  
  return(scrs)
})
names(scores) <- as.character(1:100)

class(scores[[1]][[1]][[1]][[1]])
scores[[1]][[1]][[1]][[1]]
###############################################################
# Obtain dag-classification measures
###############################################################

msrs <- deliverFuncToScore_recipeTree(scores, func="measuresDagDistEdge", dags=dataList$dags, ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1],
                                      getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                      getDataLevelPars=getHypots)


# analyze these measures

measures <- unwrapScoreFunc(scoresFunc=msrs, dataList=dataList, ws)

msrss <- c("SHD","edgeAUC","ccr")
indx <- with(measures, which(measure %in% msrss))
p <- ggplot(measures[indx,])
p <- p + geom_boxplot(aes(x=recipeFull, y=value, color=recipeFull))
p <- p + facet_wrap(~measure, scales="free")
p

# Total data sets per recipe and p
dts <- c("stdr","norm")
msrss <- c("ccr","SHD", "totEdgeD") #,"edgeD","fnr","fpr","msr","nonEdgeD","npp", "npp", "ppp","sens", "spec", "totEdgeD", "tss"
indx <- with(measures, which(measure %in% msrss & dt %in% dts))
cast(measures[indx,], recipeFull~measure, value="value", margins="grand_col")
# mean stats
tab <- cast(measures[indx,], recipeFull~measure, value="value", fun.aggregate="mean", margins="grand_col")
tab





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

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1])


#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")




n.boot <- 50

statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
statsDB <- statsArrayToDB(statsArr)


msrss <- c("ccr")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value", fun.aggregate="length")
tab[1:10,]
tab[11:20,]
tab[21:30,]

tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
tab[1:10,]
tab[11:20,]
tab[21:30,]

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




