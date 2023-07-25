# demo of causality package

remove(list=ls())
server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "tsXYZ-tanhPoly2_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"

#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
getwd()
#repos <- paste("/media/disk/users/emiliano/causaLearner/", sep="")
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/blocks/learners/"
repos <- paste("/home/emiliano/causaLearner/", sep="")
folderSave <- "/home/emiliano/causaLearner/pkg_causaLearner/experiments/blocks/learners/"
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/CMEM_ln_ts/results/learners/learnersNCE_SIMdens_f_getRbf4mod2_200/"
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/CMEM_ln_ts/results/learners/dummy/"
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/CMEM_ln_ts/results/learners/learnersNCE_CMPLXdens_f_getRbf4mod2_corrRKHS/"

dir.exists(folderSave)
dir(folderSave)[1]
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



block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#block <- 12

pm <- proc.time()
load(file=paste("./data/timeSeries/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name

plot(as.ts(dat$x))
ggpairs(as.data.frame(dat$x))

q <- length(dataList$xs)
# weights for each data set
ws <- rep(1, q)
ws <- ws/sum(ws)
names(ws) <- dataList$names

#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm")
approxDagSetMethodList <- c("timeSeriesXYZ_ts1")

hyp_scorers_cmem <- c(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
                      #cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1[1],
                      #cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1)#,
#cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[1:3])
#cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)


hypScorerList <- hyp_scorers_cmem #"rnd1"

#hypScorerList <- cmem_hypScorer_pack_none_1[2]

#hypScorerList <- "cmem_hypScorer_comp_nn1"

recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)

#recipe <- recipe[c(5,8),]

dim(recipe)


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

count <- 0
scrsTest <- sapply(scores, function(score){
  #score <- scores[[10]]
  count <<- count + 1
  print(count)
  if(!is.na(score[[1]][[1]][[4]][1])){
    res <- score[[1]][[1]][[4]][,"KCMC_reglr_comp"]
  } else{
    res <- c(NA,NA)
  }
})
scrsTest <- t(scrsTest)
sum(scrsTest[,1]<scrsTest[,2], na.rm=T)/sum(!is.na(scrsTest[,2]))


hit <- scrsTest[1:100,2]<scrsTest[1:100,1]
names(hit) <- dataList$names
dataList2 <- dataList
names(dataList2$xs) <- dataList$names
#plotPairsList(dataList=dataList2)
plotPairsList(dataList=dataList2, hit=hit)

plot(dataList$desc[,"lin"],dataList$desc[,"add"], col=c("red","blue")[hit*1+1])
plot(dataList$desc[,"lin"],dataList$desc[,"add"], col=c("red","blue")[hit*1+1], xlim=c(-5,0), ylim=c(-5,0))
df <- cbind(dataList$desc, hit=hit)
mod <- glm(hit~lin*add, data=as.data.frame(df))
summary(mod)


scores <- lapply(scores, function(score_dt){
  res <- lapply(score_dt, function(score_dt_apprx){ 
    res <- lapply(score_dt_apprx, function(score_dt_apprx2){
      # score_dt <- scores[[1]]; score_dt_apprx <- score_dt[[1]]; score_dt_apprx2 <- score_dt_apprx[[1]]
      res <- score_dt_apprx2
      names(res) <- names(scores[[400]][["norm"]][["oracleME1"]])
      return(res)
    })
    return(res)
  })
  return(res)
})


# For running all on my pc

pm <- proc.time()
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderSave)
proc.time() - pm 

# 67 mins for kernParsX on my pc with 1-1-6 cores and 1-1 on CV
# 9.8 for none on my pc with 1-1-6 cores and 1-1 on CV

#version <- "vNC_CLASS"
#version <- "vNC_REG"
version <- "v1" #"v_none" #"vNC_CORSPEARMAN"
save("scores", file=paste("./pkg_causaLearner/experiments/CMEM_ln_ts/results/", 
                          experimentName, "_", hs_cmem_ob_version, "_",
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

plot(as.ts(dataList$xs[[10]]))

[-c(10,24,40,65,66,72,78,80,82,87,100)]

aggScores <- deliverFuncToScore_recipeTree(scores, func="aggregateScores", dags=dataList$dags, 
                                           getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                           getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")



rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["timeSeriesXYZ_ts1"]])[1])


dagIds <- getHypID(sapply(dataList$dags, function(el) el, simplify="array"))




#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")

n.boot <- 50
statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
statsDB <- statsArrayToDB(statsArr)

msrss <- c("ccr")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
#cast(statsDB, dt+hypSc~cmplxFunc, value="value", fun.aggregate=length)
tabCCR <- cast(statsDB, dt+hypSc~cmplxFunc, value="value")
msrss <- c("auc")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tabAUC <- cast(statsDB, dt+hypSc~cmplxFunc, value="value")

writeRepos <- "./pkg_causaLearner/experiments/CMEM/results/"
a.msrs <- "ccr"
indx <- which(df$a.msr == a.msrs)
tabCCR <- cast(df[indx,], data+hypSc~kd.msr, value="value")
write.csv(tabCCR, file=paste(writeRepos, experimentName,"_",version ,"_CCR.csv", sep=""))

a.msrs <- "auc"
indx <- which(df$a.msr == a.msrs)
tabAUC <- cast(df[indx,], data+hypSc~kd.msr, value="value")
write.csv(tabAUC, file=paste(writeRepos, experimentName ,"_",version,"_AUC.csv", sep=""))


count <- 0
linIndex <- sapply(dataList$xs, function(xs){
  count <<- count + 1
  # count <- 1;  xs <- dataList$xs[[count]]
  print("*****************")
  print(dataList$names[count])
  # plot(xs[,"x"], xs[,"y"])
  # hist(xs[,"x"])
  # hist(xs[,"y"])
  df <- as.data.frame(xs)
  mod <- lm(y~x,df)
  #plot(df$x, residuals(mod))
  #hist(residuals(mod))
  xs2 <- apply(xs, 2, norml)
  #plot(xs2[,"x"], xs2[,"y"])
  trDat <- constructData(x=as.matrix(xs2[,"x"]), y=as.matrix(xs2[,"y"]))
  krrAux <- setParams(krr1, trDat)
  krrAux <- krrAux$learn(krrAux)
  predKrr <- pred.CV(krrAux, trDat)$test
  #o <- order(trDat$x); plot(trDat$x, trDat$y); lines(trDat$x[o], predKrr$gyh_class[o], col="red")
  residsLin <- residuals(mod)
  # plot(xs[,"x"], residsLin)
  residsNonLin <- krrAux$resids(krrAux, predKrr)[,"resid"]
  #plot(predKrr$x_class, residsNonLin)
  #plot(o, residsNonLin[o])
  res1 <- dhsic.test(X=df[,"x",drop=F],Y=matrix(residsLin,ncol=1))$p.value
  res2 <- dhsic.test(X=predKrr$x_class,Y=residsNonLin)$p.value
  #o <- order(predKrr$x_class)
  #res2b <- dhsic.test(X=o,Y=residsNonLin[o])$p.value
  res3 <- ks.test(x=residsLin, y="pnorm", mean=mean(residsLin), sd=sd(residsLin))$p.value
  res4 <- ks.test(x=residsNonLin, y="pnorm", mean=mean(residsNonLin), sd=sd(residsNonLin))$p.value
  res5 <- ks.test(x=xs2[,"x"], y="pnorm", mean=mean(xs2[,"x"]), sd=sd(xs2[,"x"]))$p.value
  res6 <- ks.test(x=xs2[,"x"], y="punif", min=min(xs2[,"x"]), max=max(xs2[,"x"]))$p.value
  res7 <- Shannon_KDP(xs2[,"x"])
  
  modDens <- kepdf(xs2[,"x"], eval.points = xs2[,"x"], kernel = "gaussian", bwtype = "adaptive")
  #hist(xs2[,"x"], prob=T); o <- order(modDens@eval.points); lines(modDens@eval.points[o], modDens@estimate[o], col="red")
  res8 <- max(modDens@estimate)-min(modDens@estimate)
  res <- c(res1, res2, res3, res4, res5, res6, res7, res8)
  names(res) <- c("lm_indep","add_indep","lm_gauss","add_gauss",
                  "cause_gauss","cause_unif", "cause_ent","cause_rngPdf")
  print(res)
  return(res)
})
linIndex <- t(linIndex)


types <- c("SIM","SIMc","SIMG","SIMln")[rep(c(1,2,3,4), rep(100, 4))]
dfSeg <- cbind(types, as.data.frame(linIndex))  

dim(rnkDecsArr)
aux <- strsplit(dimnames(rnkDecsArr)[[3]], "\\.")
sapply(aux, length)
table(sapply(aux, function(el) el[[3]]))
table(sapply(aux, function(el) el[[4]]))
lrnr <- sapply(aux, function(el) el[[3]])=="hs_cmem_kern_rbfANDkern_rbf_DEL1_negCE_makeCME_cond_logRegfeatsK_x_logRegInt3_DEL2_nc_lambda_kernParsXY_L2_lambda_kernParsX_L2"
msr <- sapply(aux, function(el) el[[4]])=="KCMC_reglr_comp"
indxScr <- which(lrnr & msr)
dfSeg <- cbind(dfSeg, rnkDecsArr[,,indxScr])
dir("./pkg_causaLearner/experiments/CMEM/results/segment_datasets/")
save(dfSeg, file="./pkg_causaLearner/experiments/CMEM/results/segment_datasets/SIM_seg.RData")





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
