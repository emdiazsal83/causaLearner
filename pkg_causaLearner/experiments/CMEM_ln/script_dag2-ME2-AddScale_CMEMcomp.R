# demo of causality package

remove(list=ls())
#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2-ME2-CmplxMixt_cmem-comp-4-l2"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v5_comp_ord"
hs_cmfm_ob_version <- "v5_comp"

#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
getwd()
#repos <- paste("/media/disk/users/emiliano/causaLearner/", sep="")
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/blocks/learners/"

repos <- paste("/home/emiliano/causaLearner/", sep="")
folderSave <- "/home/emiliano/causaLearner/pkg_causaLearner/experiments/blocks/learners/"

dir(folderSave)


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
n <- 20
ns <- rep(n, q)
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)

dagMat <- matrix(c(0,0,1,0), 2, 2)
rownames(dagMat) <- colnames(dagMat) <- names(nodes)


# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
 pm <- proc.time()
 set.seed(5)
 dataLists <- lapply(1:5, function(i){
   # i <- 1
  dataList <- simRandSEMs(q, ps, ns, nodess, sigma=1 ,markovEquiv=F, dagMat=dagMat)
  dataList$xs <- lapply(dataList$xs, function(el) apply(el, 2, norml))
  names(dataList$xs) <- dataList$names
  # plotPairsList(dataList)
  return(dataList)
})

 xs <- lapply(1:q, function(i){
   # i <- 1
   xs <- do.call(rbind, lapply(1:length(dataLists), function(j) j*10+dataLists[[j]]$xs[[i]]))
   return(xs)
 })
 noiss <- lapply(1:q, function(i){
   # i <- 1
   xs <- do.call(rbind, lapply(1:length(dataLists), function(j) dataLists[[j]]$noiss[[i]]))
   return(xs)
 })
 dags <- dataLists[[1]]$dags
 names <- dataLists[[1]]$names

 dataList <- list(dags=dags, xs=xs, noiss=noiss, names=names)
 plotPairsList(dataList)
 pm <- proc.time()
 i <- 6
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2])
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2], xlim=c(9,11), ylim=c(9,11))
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2], xlim=c(19,21), ylim=c(19,21))
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2], xlim=c(29,31), ylim=c(29,31))
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2], xlim=c(39,41), ylim=c(39,41))
 plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2], xlim=c(49,51), ylim=c(49,51))
save(dataList, file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))



block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#block <- 1

pm <- proc.time()
load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name

# weights for each data set
q <- length(dataList$names)
ws <- rep(1, q)
ws <- ws/sum(ws)
names(ws) <- dataList$names

#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm") #, "stdr"
approxDagSetMethodList <- c("oracleME1")
#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")


hyp_scorers_cmem <- c(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1,
                      #cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1,
cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)
#cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)


hypScorerList <- hyp_scorers_cmem

hypScorerList <- cmem_hypScorer_pack_none_1[2]

recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)

#recipe <- recipe[1:4,]

recipe <- recipe[c(5,8),]

dim(recipe)


#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE,  folderSave=folderSave)
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
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderSave)
proc.time() - pm 

scores
dagAux <- dataList$dags[[block]]
dim(dagAux) <- c(dim(dagAux), 1)
getHypID(dagAux)

# 8.6 mins for none_1      with 1-1-6 cores on my pc with 1-1 in CV
# for kernParsX_1 with 1-1-6 cores on my pc with 1-1 in CV

version <- "vCorr"#"vCLASS"
version <- "vNone"
save("scores", file=paste("./pkg_causaLearner/experiments/CMEM/results/", 
                          experimentName, "_", hs_cmem_ob_version, "_",
                          version  ,".RData", sep=""))

unlink(paste("./pkg_causaLearner/experiments/blocks/", blockFiles, sep=""))




#########################################################################################################################################3
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################

 

load(file=paste("./pkg_causaLearner/experiments/CMFM/results/", experimentName, "_", hs_cmfm_ob_version,"_numF_",num_f,  ".RData", sep=""))

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

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1])


dagIds <- getHypID(sapply(dataList$dags, function(el) el, simplify="array"))




#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")


n.boot <- 50
statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
statsDB <- statsArrayToDB(statsArr)

writeRepos <- "./pkg_causaLearner/experiments/CMEM/results/"
msrss <- c("ccr")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
#cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value", fun.aggregate=length)
tabCCR <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
write.csv(tabCCR, file=paste(writeRepos, experimentName, "_CCR.csv", sep=""))
msrss <- c("auc")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tabAUC <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
write.csv(tabAUC, file=paste(writeRepos, experimentName ,"_AUC.csv", sep=""))




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







