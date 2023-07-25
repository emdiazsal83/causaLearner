# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)

hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)


# p <- 2
# nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
# nodes <- rep(list(nodesX),p)
# names(nodes) <- c("x","y")
# dag <- matrix(c(0,0,1,0),2,2)
# colnames(dag) <- rownames(dag) <- c("x","y")
# nTr <- 100
# nTe <- 110
# # random function
# set.seed(4)
# simTest_tr <- simRandSEM(p, nTr, nodes, sigma=2, sigmaErr=0, dagMat=dag)
# simTest_te <- simRandSEM(p, nTe, nodes, sigma=2, sigmaErr=0, dagMat=dag)
# #plot(getGraph(simTest_tr$dag))
# Xtr <- simTest_tr$x
# Xtr <- apply(Xtr, 2, norml)
# Xte <- simTest_te$x
# Xte <- apply(Xte, 2, norml)
# xtr <- Xtr[,"x"]
# ytr <- Xtr[,"y"]
# xte <- Xte[,"x"]
# yte <- Xte[,"y"]

#experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmem-comp-ord_nc"
experimentName <- "dag2-ME2-SIMclust_cmem-comp-4-l2"
dataName <- strsplit(experimentName, "_")[[1]][1]
#dataName <- "dag2-ME2-Cmplx-SinPlus-SinTimes"

block <- 20 # 4 (fails comp msr, hit NCE - logRegCME)
pm <- proc.time()
load(file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))
#load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name
dat$dag

Xtr <- dat$x
Xtr <- apply(Xtr, 2, norml)
Xte <- dat$x
Xte <- apply(Xte, 2, norml)


apply(Xtr, 2, min)
apply(Xtr, 2, max)
apply(Xte, 2, min)
apply(Xte, 2, max)


indx_Xv <- c(1,2)[match(paste(as.numeric(dat$dag), collapse="."),  
                      c(paste(c(0,0,1,0), collapse="."),paste(c(0,1,0,0), collapse=".")))]
indx_Yv <- setdiff(c(1,2), indx_Xv)

xtr <- Xtr[,indx_Xv]
ytr <- Xtr[,indx_Yv]
xte <- Xte[,indx_Xv]
yte <- Xte[,indx_Yv]


plot(xtr,ytr)
plot(ytr, xtr)

trainDataXY <- constructData(x=as.matrix(xtr), y=ytr)
trainDataYX <- constructData(x=as.matrix(ytr), y=xtr)

testDataXY <- constructData(x=as.matrix(xte), y=yte)
testDataYX <- constructData(x=as.matrix(yte), y=xte)


# cmfm_learner_pack_none_1
# cmfm_learner_pack_lambda_1
# cmfm_learner_pack_kernParsX_1
# cmfm_learner_pack_kernParsXY_1

hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

cmfm_learner_pack_none_1
cmfm_learner_pack_kernParsXY_1
indx_lrn <- 8
                
cmfm_learner <- cmem_learner_pack_none_1[2]#cmem_learner_pack_gkernParsXY_1[indx_lrn]
cmfm_learner <- eval(parse(text=cmfm_learner))



# train hyperparameters
pm <- proc.time()
cmfm_learner_xy <- setParams(learner=cmfm_learner, trainData=trainDataXY, plot=T)
cmfm_learner_yx <- setParams(learner=cmfm_learner, trainData=trainDataYX, plot=T)
proc.time() - pm 
# 11.4 mins for 27 gridpoints for cmfm (kappa=10)
# 3.2 mins for 27 gridpoints for cmem (kappa=10)
# 90 secs for 27 gridpoints for cmem (kappa=1)
# 8.8 mins secs for 27*4=108 gridpoints for cmfm (kappa=1)



expType <- toupper(strsplit(strsplit(experimentName, "_")[[1]][2], "-")[[1]][1])
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", sep=""))
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, sep=""))
folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, "/",sep="")

fileXY <- paste(dataNm, "_",cmem_learner_pack_kernParsXY_1[indx_lrn], "_1_",paste(indx_Yv, indx_Xv,sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_",cmem_learner_pack_kernParsXY_1[indx_lrn], "_1_",paste(indx_Xv, indx_Yv,sep="on"), ".RData", sep="")
dir(folderLearners)[13:18]
fileXY  <- paste(substr(fileXY, 1, as.numeric(gregexpr("DEL1", fileXY))-1), substr(fileXY, as.numeric(gregexpr("DEL2", fileXY))+5, nchar(fileXY)), sep="")
fileYX  <- paste(substr(fileYX, 1, as.numeric(gregexpr("DEL1", fileYX))-1), substr(fileYX, as.numeric(gregexpr("DEL2", fileYX))+5, nchar(fileYX)), sep="")
fileXY %in% dir(folderLearners); fileYX %in% dir(folderLearners)
load(file=paste(folderLearners, fileXY, sep=""))
cmfm_learner_xy <- cmemLearnerAux
load(file=paste(folderLearners, fileYX, sep=""))
cmfm_learner_yx <- cmemLearnerAux


grid_xy <- cmfm_learner_xy$hyperParams$data$grid
grid_yx <- cmfm_learner_yx$hyperParams$data$grid



df_xy <- melt(grid_xy)
df_xy$dir <- "xy"
df_yx <- melt(grid_yx)
df_yx$dir <- "yx"
df <- rbind(df_xy, df_yx)

type <- "kernelY"
ini <- regexpr(type,df$params)
fin <- regexpr("lambda",df$params)
df[[type]] <- substr(df$params, ini+nchar(type)+1, fin-2)

ini <- regexpr("lambda",df$params)
fin <- nchar(as.character(df$params))
df$lambda <- as.numeric(sapply(strsplit(substr(df$params, ini+nchar("lambda="), fin)," "), function(el) el[1]))
 

ini <- regexpr(".Y=", df$params, fixed=T)
fin <- nchar(as.character(df$params))
df$parY <- as.numeric(substr(df$params, ini+nchar(".Y="), fin))
df$RKHS_Y <- as.numeric(as.factor(paste(df[,type], df$parY, sep=".")))
table(df$RKHS_Y)

df2 <- cast(df, trainTest+params+RKHS_Y+dir+lambda~var, value="value")
colnames(df2) <- sapply(strsplit(colnames(df2), "\\."), function(el) el[1])
indx <-  which(df2$trainTest == "test" )
p <- ggplot(df2[indx,])
p <- p + geom_point(aes(x=log(lambda, 10), y=negCE, colour=dir), size=0.5, alpha=0.5)
p <- p + facet_wrap(RKHS_Y~., scales="free")
p

cmfm_learner$learnParams$learnBlambda
cmfm_learner$hyperParams$data$non_optimizable$lambda$val <- 0 

cmfm_learner_xy$hyperParams$data$optimizable$gamma <- 1e-4
cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.X <- 0.1
cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y <- 100
cmfm_learner_xy <- cmfm_learner$learn(learner=cmfm_learner_xy, forLoss=T)
pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy, data=trainDataXY, forLoss=T)
plotPredCMEM(cmfm_learner_xy, pred=pred_xy, var="k",indx=1:9)

plot(cmfm_learner_xy$learnParams$Lx, pred_xy$gy_k)
plot(pred_xy$gy_k, pred_xy$gyh_k)
pred_xy$gyh_k
apply(cmfm_learner_xy$learnParams$Ky,2,mean)


cmfm_learner_xy <- cmfm_learner$learn(learner=cmfm_learner_xy, forLoss=T)
cmfm_learner_yx <- cmfm_learner$learn(learner=cmfm_learner_yx, forLoss=T)



pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy, data=testDataXY, forLoss=T)
pred_yx <- cmfm_learner_yx$pred(learner=cmfm_learner_yx, data=testDataYX, forLoss=T)

plotPredCMEM(cmfm_learner_xy, pred=pred_xy, var="f",indx=1:9)
plotPredCMEM(cmfm_learner_xy, pred=pred_xy, var="k",indx=1:9)
plotPredCMEM(cmfm_learner_yx, pred_yx, var="f",indx=1:9)
plotPredCMEM(cmfm_learner_yx, pred_yx, var="k",indx=1:9)


log(getHyperPar(cmfm_learner_xy, "lambda"), 10)
log(getHyperPar(cmfm_learner_yx, "lambda"), 10)

log(getHyperPar(cmfm_learner_xy, "sigma.rbf.Y"), 10)
log(getHyperPar(cmfm_learner_yx, "sigma.rbf.Y"), 10)

# xy
krnX_short_xy <- paste(strsplit(getHyperPar(cmfm_learner_xy, "kernelX"),"_")[[1]][2],"X", sep=".") 
indx_X <- grep(krnX_short_xy, names(cmfm_learner_xy$hyperParams$data$optimizable))
parX1_xy <- names(cmfm_learner_xy$hyperParams$data$optimizable)[indx_X]
indx_X <- grep(krnX_short_xy, names(cmfm_learner_xy$hyperParams$data$non_optimizable))
parX2_xy <- names(cmfm_learner_xy$hyperParams$data$non_optimizable)[indx_X]
parX_xy <- c(parX1_xy, parX2_xy)
featKrn_xy <- c("kernelY","featureY")[(!is.null(getHyperPar(cmfm_learner_xy, "kernelY")))*1+(!is.null(getHyperPar(cmfm_learner_xy, "featureY")))*2]
krnY_short_xy <- paste(strsplit(getHyperPar(cmfm_learner_xy, featKrn_xy),"_")[[1]][2],"Y", sep=".") 
indx_Y <- grep(krnY_short_xy, names(cmfm_learner_xy$hyperParams$data$optimizable))
parY1_xy <- names(cmfm_learner_xy$hyperParams$data$optimizable)[indx_Y]
indx_Y <- grep(krnY_short_xy, names(cmfm_learner_xy$hyperParams$data$non_optimizable))
parY2_xy <- names(cmfm_learner_xy$hyperParams$data$non_optimizable)[indx_Y]
parY_xy <- c(parY1_xy, parY2_xy)
# yx
krnX_short_yx <- paste(strsplit(getHyperPar(cmfm_learner_yx, "kernelX"),"_")[[1]][2],"X", sep=".") 
indx_X <- grep(krnX_short_yx, names(cmfm_learner_yx$hyperParams$data$optimizable))
parX1_yx <- names(cmfm_learner_yx$hyperParams$data$optimizable)[indx_X]
indx_X <- grep(krnX_short_yx, names(cmfm_learner_yx$hyperParams$data$non_optimizable))
parX2_yx <- names(cmfm_learner_yx$hyperParams$data$non_optimizable)[indx_X]
parX_yx <- c(parX1_yx, parX2_yx)
featKrn_yx <- c("kernelY","featureY")[(!is.null(getHyperPar(cmfm_learner_yx, "kernelY")))*1+(!is.null(getHyperPar(cmfm_learner_yx, "featureY")))*2]
krnY_short_yx <- paste(strsplit(getHyperPar(cmfm_learner_yx, featKrn_yx),"_")[[1]][2],"Y", sep=".") 
indx_Y <- grep(krnY_short_yx, names(cmfm_learner_yx$hyperParams$data$optimizable))
parY1_yx <- names(cmfm_learner_yx$hyperParams$data$optimizable)[indx_Y]
indx_Y <- grep(krnY_short_yx, names(cmfm_learner_yx$hyperParams$data$non_optimizable))
parY2_yx <- names(cmfm_learner_yx$hyperParams$data$non_optimizable)[indx_Y]
parY_yx <- c(parY1_yx, parY2_yx)




log(getHyperPar(cmfm_learner_xy, parX_xy), 10)
log(getHyperPar(cmfm_learner_yx, parX_yx), 10)
log(getHyperPar(cmfm_learner_xy, parY_xy), 10)
log(getHyperPar(cmfm_learner_yx, parY_yx), 10)


# obtain loss out of sample
# loss <- "cmem_L2_k"
# this applies loss directly on predicted vector, no folds or CV
losses_xy <- unlist(sapply(names(cmfm_learner_xy$optimizeParams$losses), function(loss) do.call(loss, list(learner=cmfm_learner_xy, pred=pred_xy))))
losses_yx <- unlist(sapply(names(cmfm_learner_yx$optimizeParams$losses), function(loss) do.call(loss, list(learner=cmfm_learner_yx, pred=pred_yx))))
data.frame(losses_xy, losses_yx)



# calculate measures
pm <- proc.time()
msrs_xy <- cmfm_learner_xy$calcMsrs(cmemLearner=cmfm_learner_xy) # learner=cmfm_learner_xy
msrs_yx <- cmfm_learner_yx$calcMsrs(cmemLearner=cmfm_learner_yx)
proc.time() - pm # 0.415 secs
msrs_xy
msrs_yx


newKernelY_XY <- cmfm_learner_xy$hyperParams$data$optimizable[[featKrn_xy]]$val
newKernelY_YX <- cmfm_learner_yx$hyperParams$data$optimizable[[featKrn_yx]]$val

parY_xy
parY_yx

valsY_xy <- lapply(parY_xy, function(par){
  # par <- parY_xy[1]
  optimNonoptim <- c("optimizable","non_optimizable")[(par %in% names(cmfm_learner_xy$hyperParams$data$optimizable))*1+(par %in% names(cmfm_learner_xy$hyperParams$data$non_optimizable))*2]
  cmfm_learner_xy$hyperParams$data[[optimNonoptim]][[par]]$val
})
names(valsY_xy) <- parY_xy

valsY_yx <- lapply(parY_yx, function(par){
  # par <- parY_yx[1]
  optimNonoptim <- c("optimizable","non_optimizable")[(par %in% names(cmfm_learner_yx$hyperParams$data$optimizable))*1+(par %in% names(cmfm_learner_yx$hyperParams$data$non_optimizable))*2]
  cmfm_learner_yx$hyperParams$data[[optimNonoptim]][[par]]$val
})
names(valsY_yx) <- parY_yx
# now lets set the sigma parameter of xy and yx learners
# as the direct sum of the two optimal RKHS H_y  and compare
# measures 

newKer <- unique(c(newKernelY_XY, newKernelY_YX))

length(cmfm_learner_pack_kernParsX_1)
length(cmfm_learner_pack_kernParsXY_1)

indx_lrn <- 10
cmfm_learner_pack_kernParsX_1[indx_lrn]
cmfm_learner_pack_kernParsXY_1[indx_lrn]

cmfm_learner <- cmem_learner_pack_kernParsX_1[indx_lrn]
cmfm_learner <- eval(parse(text=cmfm_learner))


cmfm_learner$hyperParams$data$non_optimizable[[featKrn_xy]]$val <- newKer
names(cmfm_learner$hyperParams$data$non_optimizable)

parms <- unique(c(parY_xy, parY_yx))
for(parm in parms){
  # parm <- unique(c(parY_xy, parY_yx))[2]
  print(cmfm_learner$hyperParams$data$non_optimizable[[parm]])
  cmfm_learner$hyperParams$data$non_optimizable[[parm]]$val <- NULL
}
names(cmfm_learner$hyperParams$data$non_optimizable)

parms <- c(parY_xy, parY_yx)
vals <- c(valsY_xy, valsY_yx)
for(parm in parms){
  # parm <- parms[[1]]
  print("********************************")
  print(parm)
  print(names(cmfm_learner$hyperParams$data$non_optimizable))
  # parm <- unique(c(parY_xy, parY_yx))[2]
  indx <- match(parm, names(vals))
  oldVal <- cmfm_learner$hyperParams$data$non_optimizable[[parm]]$val
  cmfm_learner$hyperParams$data$non_optimizable[[parm]]$val <- unique(c(oldVal, vals[[indx]]))
}

names(cmfm_learner$hyperParams$data$non_optimizable)


pm <- proc.time()
cmfm_learner_xy_2 <- setParams(learner=cmfm_learner, trainData=trainDataXY, plot=F)
cmfm_learner_yx_2 <- setParams(learner=cmfm_learner, trainData=trainDataYX, plot=F)
proc.time() - pm # 3.8 mins

dir(folderLearners)
fileXY <- paste(dataNm, "_" ,cmem_learner_pack_kernParsX_1[indx_lrn],"_2_", paste(indx_Yv, indx_Xv,sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_" ,cmem_learner_pack_kernParsX_1[indx_lrn],"_2_", paste(indx_Xv, indx_Yv,sep="on"), ".RData", sep="")
dir(folderLearners)[19:22]
fileXY
fileXY %in% dir(folderLearners)
fileYX %in% dir(folderLearners)
load(file=paste(folderLearners, fileXY, sep=""))
cmfm_learner_xy_2 <- cmemLearnerAux
load(file=paste(folderLearners, fileYX, sep=""))
cmfm_learner_yx_2 <- cmemLearnerAux




grid_xy <- cmfm_learner_xy_2$hyperParams$data$grid
grid_yx <- cmfm_learner_yx_2$hyperParams$data$grid


df_xy <- melt(grid_xy)
df_xy$dir <- "xy"
df_yx <- melt(grid_yx)
df_yx$dir <- "yx"
df <- rbind(df_xy, df_yx)


# L2


df2 <- cast(df, trainTest+params+dir~var, value="value")
colnames(df2) <- sapply(strsplit(colnames(df2), "\\."), function(el) el[1])

indx <-  which(df2$trainTest == "test" )
p <- ggplot(df2[indx,])
p <- p + geom_point(aes(x=log(KCMC, 10), y=cmem_L2_f, colour=dir), size=0.5, alpha=0.5)
p

head(df2)
min(df2$negCE[which(df2$trainTest=="test")])

# now massage data so that we can do KCDC vs entropy loss

names(cmfm_learner_xy$optimizeParams$losses)
loss <- "cmem_L2_f"
reglrs <- c("KCDC","KCMC","KCRDC","KCSC","KCNSC","KCCC_ent","KCCC_pca_ent")
reglrs <- c("KCDC","KCMC","KCRDC","KCSC","KCNSC")
reglrs <- c("KCDC","KCMC")
sapply(reglrs, function(reglr){
  # reglr <- reglrs[1]
  #print(reglr)
  logReglr <- switch(reglr, KCCC_ent=F, KCCC_pca_ent=F,T)
  curveXY <- getLossReglrCurve(learner=cmfm_learner_xy_2, loss=loss, reglr=reglr, logReglr=logReglr)
  curveYX <- getLossReglrCurve(learner=cmfm_learner_yx_2, loss=loss, reglr=reglr, logReglr=logReglr)
  curveList <- list(curveXY, curveYX)
  names(curveList) <- c("xy","yx")
  compReglrVals <- getComparableReglrVal(curveList)
  plotCurves2(compReglrVals, curveList, xlim=c(1,1.5))
  #plotCurves(compReglrVals, curveList)
  return(compReglrVals$reglrs)
})


# train learn parameters
cmfm_learner_xy_2 <- cmfm_learner$learn(learner=cmfm_learner_xy_2, forLoss=T)
cmfm_learner_yx_2 <- cmfm_learner$learn(cmfm_learner_yx_2, forLoss=T)

# bin
#getHyperPar(cmfm_learner_xy_2, "num.bin.X")
#log(getHyperPar(cmfm_learner_yx_2, "sigma.rbf.Y"), 10)

#L2
log(getHyperPar(cmfm_learner_xy_2, "lambda"), 10)
log(getHyperPar(cmfm_learner_yx_2, "lambda"), 10)

parX_xy <- names(cmfm_learner_xy_2$hyperParams$data$optimizable)[grep(paste(strsplit(getHyperPar(cmfm_learner_xy_2, "kernelX"),"_")[[1]][2],"X",sep="."), names(cmfm_learner_xy_2$hyperParams$data$optimizable))]
parX_yx <- names(cmfm_learner_yx_2$hyperParams$data$optimizable)[grep(paste(strsplit(getHyperPar(cmfm_learner_yx_2, "kernelX"),"_")[[1]][2],"X",sep="."), names(cmfm_learner_yx_2$hyperParams$data$optimizable))]


log(getHyperPar(cmfm_learner_xy, parX_xy), 10)
log(getHyperPar(cmfm_learner_yx, parX_yx), 10)



# obtain pred out of sample
pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy_2, data=testDataXY, forLoss=T)
pred_yx <- cmfm_learner_yx$pred(learner=cmfm_learner_yx_2, data=testDataYX, forLoss=T)


plotPredCMEM(cmfm_learner_xy_2, pred=pred_xy, var="k",indx=1:9)
plotPredCMEM(cmfm_learner_yx_2, pred_yx, var="k",indx=1:9)



# obtain loss out of sample
# loss <- "cmem_L2_k"
sapply(cmfm_learner_xy$optimizeParams$losses, function(loss) do.call(loss, list(learner=cmfm_learner_xy_2, pred=pred_xy)))
sapply(cmfm_learner_yx$optimizeParams$losses, function(loss) do.call(loss, list(learner=cmfm_learner_yx_2, pred=pred_yx)))

pred_xy$gyh_class
pred_yx$gyh_class

# calculate measures
pm <- proc.time()
msrs_xy <- cmfm_learner_xy$calcMsrs(cmemLearner=cmfm_learner_xy_2) # learner=cmfm_learner_xy
msrs_yx <- cmfm_learner_yx$calcMsrs(cmemLearner=cmfm_learner_yx_2)
proc.time() - pm # 0.415 secs

msrs_xy
msrs_yx

# now lets build a hypothesis scorer with

# getLossRegCurve
# getComparableRegVal
# plotCurves

# step 1 - regressions : get regression for each node in dags
# as usual with CV on NCE for lambda, parsKernX, parsKernY

# step 2 - choose Hy : only regressions with same number of 
# regressors will be compared so we choose sum of Hys for
# same size regressions. Actually best do it by node?
# one Hy per node

# step 3 - regression: get regressions for each node in dags
# with CV on L2-loss for lambda and parsKernX

# step 4 - compare regressions: compare regressions of same size
# (there should be at least two otherwise all dags in hypothesis
# set have the same regression and we don't have to assess it)
# get envelope curves for each set of regressions, get measures
# and assign scores back to the regression.

# step 5 - assemble dags with scores as before

# simulate some data for dag with two 2-variable regressions
# and two 1-variable regressions so that I can test the 
# hypothesis scorer properly 
p <- 5
n <- 50
nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
nodes <- rep(list(nodesX),p)
names(nodes) <- seq(p)
dag2 <- matrix(c(0,0,1,0),2,2)
dag5 <- matrix(c(0,1,0,1,0, 0,0,0,1,0, 1,1,0,0,0, 0,0,0,0,1, 0,0,0,0,0),5,5)
dag <- switch(c("two","five")[match(p, c(2,5))], two=dag2, five=dag5)
colnames(dag) <- rownames(dag) <- seq(p)
#plot(getGraph(dag))
set.seed(4)
simTest <- simRandSEM(p, n, nodes, sigma=2, sigmaErr=0, dagMat=dag)
library(GGally)
ggpairs(as.data.frame(simTest$x))

#plot(getGraph(simTest$dag))
x <- simTest$x
x <- apply(x, 2, norml)

source("./pkg_causaLearner/approxDagSetMethods/func_approxDagSetMethods_v1.R", echo=FALSE)

hypArray <- oracleME(trueDAG=dag)
dim(hypArray)

source("./pkg_learner/func_learners_v4.R", echo=FALSE)
source("./pkg_learner/ob_cmf_learner_v3.R", echo=FALSE)

cmfm_learner_pack_compCombos_1
indx_lrn <- 1
cmemLearner <- cmfm_learner_pack_compCombos_1[indx_lrn]

source("./pkg_causaLearner/hypothesisScorers/cmem_hypothesisScorer/func_CMEM_hypothesisScorer_v2.R", echo=FALSE)
source("./pkg_causaLearner/hypothesisScorers/pkg_ANM_hypothesisScorer/func_ANM_hypothesisScorer_v2.R", echo=FALSE)
  

scores <- cmem_hypScorer_comp(x, hypArray, cmemLearner)

dagsChosen <- apply(scores, 2, which.min)
dagsChosen <- rownames(scores)[dagsChosen]
dagsChosen <- lapply(dagsChosen, function(dag) dagIDtoMatrix(as.numeric(dag),p))

# get structural hamming distance of chosen dag 
shds <- sapply(dagsChosen, function(dagChosen) SHD(dagChosen, dagTrue=dag))
names(shds) <- colnames(scores)
shds

# lets run whole recipe in preparation for running it on server
hs_cmem_ob_version <- "v3"
hs_cmfm_ob_version <- "v3"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

data <- dat
dataTreatmentList <- c("norm") #, "stdr"
approxDagSetMethodList <- c("oracleME1")
indx_lrn <- 1
hypScorerList <- cmfm_hypScorer_pack_compCombos_1[indx_lrn]
recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)


scrs <- applyLearnRecipesOne(recipe, data,  numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderLearners)

