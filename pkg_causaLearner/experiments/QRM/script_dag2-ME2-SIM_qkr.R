# demo of causality package

remove(list=ls())
server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2-ME2-SIM-1000_qkr"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"



#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
getwd()
#repos <- paste("/media/disk/users/emiliano/causaLearner/", sep="")
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/blocks/learners/"
#repos <- paste("/home/emiliano/causaLearner/", sep="")                         
folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/QRM/learners/SIM_v3/"
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/NeuralNetworks/results/TCEP/experimentsnew/"


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


# Data downloaded from:
# https://webdav.tuebingen.mpg.de/cause-effect/
############################################
# SIM
############################################
# folder <- "data/TCEPs/14-518-appendix2/SIM/"
# fileName <- "pairmeta.txt"
# fileFolder <- paste(folder, fileName, sep="")
# meta <- read.csv(fileFolder, sep="", header=F)
# # From README file
# colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
# head(meta)
# # which one is cause and which one effect?
# indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
# length(indxUni)
# pairs <- meta$pairNumber[indxUni]
# summary(meta[indxUni,])
# # how many data points do they have?
# (nPerPair <- numDataPairs(pairs, folder))
# # get range for x and y for each pair
# summaryDataPairs(meta, pairs, folder)
# # create data test
# dataList1 <- createTCEPList2(pairs, meta, folder)
# plotPairsList(dataList1)
# ############################################
# # SIM-c
# ############################################
# folder <- "data/TCEPs/14-518-appendix2/SIM-c/"
# fileName <- "pairmeta.txt"
# fileFolder <- paste(folder, fileName, sep="")
# meta <- read.csv(fileFolder, sep="", header=F)
# # From README file
# colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
# head(meta)
# # which one is cause and which one effect?
# indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
# length(indxUni)
# pairs <- meta$pairNumber[indxUni]
# summary(meta[indxUni,])
# # how many data points do they have?
# (nPerPair <- numDataPairs(pairs, folder))
# # get range for x and y for each pair
# summaryDataPairs(meta, pairs, folder)
# # create data test
# dataList2 <- createTCEPList2(pairs, meta, folder)
# plotPairsList(dataList2)
# ############################################
# # SIM-G
# ############################################
# folder <- "data/TCEPs/14-518-appendix2/SIM-G/"
# fileName <- "pairmeta.txt"
# fileFolder <- paste(folder, fileName, sep="")
# meta <- read.csv(fileFolder, sep="", header=F)
# # From README file
# colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
# head(meta)
# # which one is cause and which one effect?
# indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
# length(indxUni)
# pairs <- meta$pairNumber[indxUni]
# summary(meta[indxUni,])
# # how many data points do they have?
# (nPerPair <- numDataPairs(pairs, folder))
# # get range for x and y for each pair
# summaryDataPairs(meta, pairs, folder)
# # create data test
# dataList3 <- createTCEPList2(pairs, meta, folder)
# plotPairsList(dataList3)
# ############################################
# # SIM-ln
# ############################################
# folder <- "data/TCEPs/14-518-appendix2/SIM-ln/"
# fileName <- "pairmeta.txt"
# fileFolder <- paste(folder, fileName, sep="")
# meta <- read.csv(fileFolder, sep="", header=F)
# # From README file
# colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
# head(meta)
# # which one is cause and which one effect?
# indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
# length(indxUni)
# pairs <- meta$pairNumber[indxUni]
# summary(meta[indxUni,])
# # how many data points do they have?
# (nPerPair <- numDataPairs(pairs, folder))
# # get range for x and y for each pair
# summaryDataPairs(meta, pairs, folder)
# # create data test
# dataList4 <- createTCEPList2(pairs, meta, folder)
# plotPairsList(dataList4)
# ############################################
# # ALL
# ############################################
# dataListList <- list(SIM=dataList1, SIMc=dataList2, SIMG=dataList3,SIMln=dataList4)
# dataList <- dataJoin(dataListList)
# 
# # weights for each data set
# ws <- rep(1, length(dataList$xs))
# ws <- ws/sum(ws)
# names(ws) <- dataList$names
# save(list=c("dataList","ws"), file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))

#block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
block <- 1

pm <- proc.time()
load(file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name

plot(dat$x[,1],dat$x[,2])
plot(dat$x[,2],dat$x[,1])
#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("stdr")
approxDagSetMethodList <- c("oracleME1")
#hypScorerList <- c("kqr2_re_red","kqr2_ho_red")
#hypScorerList <- c("cqr1_re_red","fqr1_re_red","nnqr1_re_red")
hypScorerList <- c("cqr1_re_red","fqr1_re_red","nnqr1_re_red")


recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)



#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################

print("start")
pm <- proc.time()
#data=dat;numCoresDt=1;numCoresHypSc=1;plot=TRUE;folderSave=folderSave
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE,  folderSave=folderSave)
proc.time()-pm 


# 14 secs for none_1 on my pc


save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName, "_", block, ".RData", sep=""))

print("finish")
q("no")

for(block in 1:length(dataList$xs)){
  # block <- 161
  print(paste("block: ", block))
  pm <- proc.time()
  load(file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))
  dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
  plot(dat$x[,1], dat$x[,2], main=paste("block: ", block))
  proc.time() - pm # 0.058 secs
  dataNm <- dat$name
  scores <- try(applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=TRUE,  folderSave=folderSave))
  if(class(scores)!= "try-error"){
    print("scores")
    print(scores)
    save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName, "_", block, ".RData", sep=""))
  } else{
    print("error")
  }
}



blockFiles <- dir("./pkg_causaLearner/experiments/blocks/")
indx <- grep(experimentName, blockFiles)
blockFiles <- blockFiles[indx]
blockNumber <- as.numeric(sapply(strsplit(sapply(strsplit(blockFiles, "_"), function(el) el[[3]]), "\\."), function(el) el[[1]]))
(blockFiles <- blockFiles[order(blockNumber)])

scores <- lapply(blockFiles, function(el){
  # el <- blockFiles[93]
  load(file=paste("./pkg_causaLearner/experiments/blocks/", el, sep=""))
  #scores2 <- lapply(scores, function(el) lapply(el, function(el) el[-c(4,7)]))
  #res <- scores2
  res <- scores
  return(res)
})
names(scores) <- dataList$names

# add dimnames
scores <- lapply(scores, function(scr) lapply(scr, function(scr) lapply(scr, function(scr) lapply(scr, function(scr){
  res <- scr
  names(dimnames(res)) <- c("dag","score")
  return(res)
}))))

# For running all on my pc

pm <- proc.time()
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderSave)
proc.time() - pm 

# 67 mins for kernParsX on my pc with 1-1-6 cores and 1-1 on CV
# 9.8 for none on my pc with 1-1-6 cores and 1-1 on CV

#version <- "vNC_CLASS"
#version <- "vNC_REG"
version <- "v1" #"v_none" #"vNC_CORSPEARMAN"
save("scores", file=paste("./pkg_causaLearner/experiments/QRM/results/", 
                          experimentName, "_",
                          version  ,".RData", sep=""))

unlink(paste("./pkg_causaLearner/experiments/blocks/", blockFiles[indx], sep=""))



 #########################################################################################################################################3
 # Pre- processing of scores (inc. creation of new scores by combining old ones) 
 #########################################################################################################################################
 
version <- "vNC_REG"
version <- "vNC_CLASS"
load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, "_",
                hs_cmem_ob_version,"_",version,".RData", sep=""))

names(scores)

setsGeo <- c("1","2","3","4",
              "20", "21",
              "42","43","44","45","46",
              "49","50","51",
              "72","73","77","78","79","80",
              "81","82","83","87",
              "89","90","91","92","93","94")

indxGeo <- which(sapply(dataList$xs, function(el) sum(apply(el,2, function(col) length(unique(col))))/length(el))>0.75)


indxGeo <- which(names(scores) %in% setsGeo)


dataList <- list(dags=dataList$dags[indxGeo], xs=dataList$xs[indxGeo], noiss=dataList$noiss[indxGeo], names=dataList$names[indxGeo])
scores <- scores[indxGeo]
ws <- rep(1, length(scores))
ws <- ws/sum(ws)
names(ws) <- dataList$names

plotPairsList(dataList)

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

# for vanilla_corspearman 93 doesnt work (98 in name)
# for vanilla_corspearman 67,93 doesnt work (72,98 in name)

aggScores <- deliverFuncToScore_recipeTree(scores, func="aggregateScores", dags=dataList$dags, 
                                           getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                           getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1])


#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")


n.boot <- 50
exps <- c("SIM","SIMc","SIMG", "SIMln")
res <- sapply(exps, function(exp){
  # exp <- "SIM"
  print("***************")
  print(exp)
  
  if(exp=="SIM"){
    indx <- 1:100
  } else{
    indx <- which( substr(dimnames(rnkDecsArr)$dataset, 1, nchar(exp))==exp)  
  }
  
  
  statsArr <- measures2DagClass(rankDecArray=rnkDecsArr[indx,,], ws, n.boot)
  statsDB <- statsArrayToDB(statsArr)
  length(unique(statsDB$hypSc))
  
  msrss <- c("ccr")
  stats <- c("orig")
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
  #cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value", fun.aggregate=length)
  tabCCR <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
  msrss <- c("auc")
  stats <- c("orig")
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
  tabAUC <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
  c.nms <- colnames(tabAUC)
  r.nms <- paste(tabAUC$dt, tabAUC$hypSc, sep=".")
  res <- abind(tabCCR[,3:ncol(tabCCR)],tabAUC[,3:ncol(tabCCR)], along=3)
  dim(res)
  dimnames(res) <- list(hypSc=r.nms, kd.msr=c.nms[3:length(c.nms)], a.msr=c("ccr", "auc"))               
  return(res)
}, simplify="array")
names(dimnames(res))[4] <- "data"

df <- melt(res)
writeRepos <- "./pkg_causaLearner/experiments/QRM/results/"
a.msrs <- "ccr"
indx <- which(df$a.msr == a.msrs)

tabCCR <- cast(df[indx,], data+hypSc~kd.msr, value="value")
write.csv(tabCCR, file=paste(writeRepos, experimentName,"_",version ,"_CCR.csv", sep=""))

a.msrs <- "auc"
indx <- which(df$a.msr == a.msrs)
tabAUC <- cast(df[indx,], data+hypSc~kd.msr, value="value")
write.csv(tabAUC, file=paste(writeRepos, experimentName ,"_",version,"_AUC.csv", sep=""))


linIndex <- sapply(dataList$xs, function(xs){
  # xs <- dataList$xs[[1]]
  df <- as.data.frame(xs)
  mod <- lm(y~x,df)
  #plot(df$x, residuals(mod))
  res <- dhsic.test(X=df[,"x",drop=F],Y=matrix(residuals(mod),ncol=1))$p.value
  res <- c(res, ks.test(x=residuals(mod), y="pnorm")$p.value)
  return(res)
})

# hsic: null independent -> linear bad
# ks  : null dist normal ->        bad

# so in both cases we favor rejections or low p.values
linIndex <- t(linIndex)
linIndex <- apply(linIndex, 1, function(col) col[1]*col[2])
ord <- order(linIndex, decreasing=F)
linIndex[ord]

dataList2 <- dataList
dataList2$xs <- dataList$xs[ord]
plotPairsList(dataList2)

i <- 99
plot(dataList2$xs[[i]][,1],dataList2$xs[[i]][,2])

rnkDecsArr2 <- rnkDecsArr[ord,,]
ws2 <- ws[ord]

ini <- 1
seqFin <- round(seq(10,102, length.out=20)) 

n.boot <- 50
msrss <- c("auc")
stats <- c("orig")

res <- numeric()

for(fin in seqFin){
  #fin <- seqFin[20] 
  print(fin)
  wsAux <- (ws2[ini:fin])/sum(ws2[ini:fin])
  statsArr <- measures2DagClass(rankDecArray=rnkDecsArr2[ini:fin,,], wsAux, n.boot)
  statsDB <- statsArrayToDB(statsArr)
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
  tabAUC <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
  res <- c(res, tabAUC[1,"KCMC_reglr_comp"])
  print(res)
}

plot(seqFin, res)
indxMax <- which.max(res)
lines(seqFin[indxMax], max(res), col="red", cex=2, type="p")
seqFin[(indxMax-1):indxMax]
res[(indxMax-1):indxMax]
linIndex[ord][seqFin[(indxMax-1)]:seqFin[indxMax]]

dataList3 <- dataList2
dataList3$xs <- dataList2$xs[ini:seqFin[which.max(res)]]
rnkDecsArr3 <- rnkDecsArr2[ini:seqFin[which.max(res)],,]
ws3 <- ws[ini:seqFin[which.max(res)]]
linIndex[ord][seqFin[which.max(res)]]

indxDisc <- sapply(dataList3$xs, function(el) sum(apply(el,2, function(col) length(unique(col))))/length(el))

ord2 <- order(indxDisc, decreasing=T)
indxDisc[ord2]

dataList4 <- dataList3
dataList4$xs <- dataList3$xs[ord2]
plotPairsList(dataList4)

i <- 80
plot(dataList4$xs[[i]][,1],dataList4$xs[[i]][,2])

rnkDecsArr4 <- rnkDecsArr3[ord2,,]
ws4 <- ws3[ord2]

ini <- 1
seqFin <- round(seq(10,83, length.out=20)) 

n.boot <- 50
msrss <- c("auc")
stats <- c("orig")

res2 <- numeric()

for(fin in seqFin){
  #fin <- seqFin[1] 
  print(fin)
  statsArr <- measures2DagClass(rankDecArray=rnkDecsArr4[ini:fin,,], ws4[ini:fin], n.boot)
  statsDB <- statsArrayToDB(statsArr)
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
  tabAUC <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
  res2 <- c(res2, tabAUC[1,"KCMC_reglr_comp"])
}

plot(seqFin, res2)
