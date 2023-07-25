# demo of causality package

remove(list=ls())
server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"


experimentName <- "dag2-ME2-emulator_cmem-comp-nn-4-l2"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v5_comp"
hs_cmfm_ob_version <- "v5_comp"

#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
getwd()
#repos <- paste("/media/disk/users/emiliano/causaLearner/", sep="")

repos <- paste("/home/emiliano/causaLearner/", sep="")
folderSave <- "/home/emiliano/causaLearner/pkg_causaLearner/experiments/NeuralNetworks/results/prosailEmulation/experimentsnew/"
folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/NeuralNetworks/results/prosailEmulation/experimentsnew/"

dir(folderSave)
dir("/home/emiliano/causaLearner/pkg_causaLearner/experiments/NeuralNetworks/results/")
dir("/home/emiliano/causaLearner/pkg_causaLearner/experiments/NeuralNetworks/results/prosailEmulation/experimentsnew/")

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


# dataRepos <- "./data/DATA_CAUSAL/emulator_problems_v6/"
# n <- 100
# pm <- proc.time()
# dataList <- createEmulatedList(dataRepos, n)
# proc.time() - pm #3 mins
# save(dataList, file=paste("./data/DATA_CAUSAL/", dataName, "_sims.RData", sep=""))
# plotPairsList(dataList)

block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#block <- 3

pm <- proc.time()
load(file=paste("./data/DATA_CAUSAL/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dat$name

indxAp <- 1:91
dataList$xs <- dataList$xs[indxAp]
dataList$dags <- dataList$dags[indxAp]
dataList$noiss <- dataList$noiss[indxAp]
dataList$names <- dataList$names[indxAp]
lapply(dataList, length)  


# weights for each data set
q <- length(dataList$xs)
ws <- rep(1, q)
ws <- ws/sum(ws)
names(ws) <- dataList$names


#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm")
approxDagSetMethodList <- c("oracleME1")
#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")

hypScorerList <- "cmem_hypScorer_comp_nn1"


recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)


dim(recipe)

#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################

print("start")

pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE,  folderSave=folderSave)
proc.time()-pm 


# 5 mins with 1-1 cores on ERC with 10-1 in CV
# 48 secs with 1-1 cores on ERC with 10-1 in CV

print("done")

# seems like you can parallelize on data, scorers, parameters but not folds

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

i <- 1
fallas <- sapply(1:length(scores), function(i) all(is.na(scores[[i]][[1]][[1]][[1]])))
sum(fallas)

# 6 for 10_12
# 5 for 40_30
# for 60_50
# mins for kernParsX on my pc with 1-1-6 cores and 1-1 on CV
#  for none on my pc with 1-1-6 cores and 1-1 on CV

version <- "NNs_60_50"
save("scores", file=paste("./pkg_causaLearner/experiments/CMEM/results/", 
                          experimentName, "_", hs_cmem_ob_version,"_numF_",num_f[1], "_",
                          version  ,".RData", sep=""))

unlink(paste("./pkg_causaLearner/experiments/blocks/", blockFiles[indx], sep=""))



 #########################################################################################################################################3
 # Pre- processing of scores (inc. creation of new scores by combining old ones) 
 #########################################################################################################################################
 
load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, ".RData", sep=""))


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

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1, 
                                         hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1])


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

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p




