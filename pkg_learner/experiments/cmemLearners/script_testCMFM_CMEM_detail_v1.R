# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)

hs_cmem_ob_version <- "v3_comp"
hs_cmfm_ob_version <- "v3_comp"
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

experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegDists"
dataName <- strsplit(experimentName, "_")[[1]][1]
#dataName <- "dag2-ME2-Cmplx-SinPlus-SinTimes"

block <- 250 # 4 (fails comp msr, hit NCE - logRegCME)
pm <- proc.time()
load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
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


indx_X <- c(1,2)[match(paste(as.numeric(dat$dag), collapse="."),  
                      c(paste(c(0,0,1,0), collapse="."),paste(c(0,1,0,0), collapse=".")))]
indx_Y <- setdiff(c(1,2), indx_X)

xtr <- Xtr[,indx_X]
ytr <- Xtr[,indx_Y]
xte <- Xte[,indx_X]
yte <- Xte[,indx_Y]


plot(xtr,ytr)
plot(ytr, xtr)

trainDataXY <- constructData(as.matrix(xtr), ytr)
trainDataYX <- constructData(as.matrix(ytr), xtr)

testDataXY <- constructData(as.matrix(xte), yte)
testDataYX <- constructData(as.matrix(yte), xte)


# cmfm_learner_pack_none_1
# cmfm_learner_pack_lambda_1
# cmfm_learner_pack_kernParsX_1
# cmfm_learner_pack_kernParsXY_1

hs_cmem_ob_version <- "v3_comp"
hs_cmfm_ob_version <- "v3_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

cmfm_learner_pack_none_1
indx_lrn <- 1
cmfm_learner <- cmfm_learner_pack_kernParsXY_1[indx_lrn]
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
fileXY <- paste(dataNm, "_1_", paste(indx_Y, indx_X,sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_1_", paste(indx_X, indx_Y,sep="on"), ".RData", sep="")
fileXY <- paste(dataNm, "_1_", paste("x", "y",sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_1_", paste("y", "x",sep="on"), ".RData", sep="")

dir(folderLearners)
load(file=paste(folderLearners, fileXY, sep=""))
cmfm_learner_xy <- cmemLearnerAux
load(file=paste(folderLearners, fileYX, sep=""))
cmfm_learner_yx <- cmemLearnerAux


grid_xy <- cmfm_learner_xy$hyperParams$data$grid
grid_yx <- cmfm_learner_yx$hyperParams$data$grid


# L2
dimnames(grid_xy)$sigma.rbf.Y <- seq(dim(grid_xy)["sigma.rbf.Y"])
dimnames(grid_xy)$sigma.rbf.X <- seq(dim(grid_xy)["sigma.rbf.X"])
dimnames(grid_yx)$sigma.rbf.Y <- seq(dim(grid_yx)["sigma.rbf.Y"])
dimnames(grid_yx)$sigma.rbf.X <- seq(dim(grid_yx)["sigma.rbf.X"])


df_xy <- melt(grid_xy)
df_xy$dir <- "xy"
df_yx <- melt(grid_yx)
df_yx$dir <- "yx"
df <- rbind(df_xy, df_yx)

# bins
# indx <-  which(df$trainTest == "test" & df$var %in% c("cmem_L2_f", "negCE", "MisCR"))
# p <- ggplot(df[indx,])
# p <- p + geom_line(aes(x=log(sigma.rbf.Y, 10), y=value, colour=dir))
# p <- p + geom_point(aes(x=log(sigma.rbf.Y, 10), y=value, colour=dir))
# p <- p + facet_grid(var~num.bin.X, scales="free")
# p 


# L2
indx <-  which(df$trainTest == "test" & df$var %in% c("MisCR"))
p <- ggplot(df[indx,])
p <- p + geom_line(aes(x=log(lambda, 10), y=value, colour=dir))
p <- p + geom_point(aes(x=log(lambda, 10), y=value, colour=dir))
p <- p + facet_grid(sigma.rbf.X~sigma.rbf.Y, scales="free")
p 

df2 <- cast(df, trainTest+lambda+sigma.rbf.X+sigma.rbf.Y+dir~var, value="value")
indx <-  which(df2$trainTest == "test" )
p <- ggplot(df2[indx,])
p <- p + geom_point(aes(x=log(KCMC, 10), y=negCE, colour=dir), size=0.5, alpha=0.5)
p <- p + facet_wrap(sigma.rbf.Y~., scales="free")
p



cmfm_learner_xy <- cmfm_learner$learn(learner=cmfm_learner_xy, forLoss=T)
cmfm_learner_yx <- cmfm_learner$learn(cmfm_learner_yx, forLoss=T)


dfBeta_xy <- as.data.frame(cbind(x=cmfm_learner_xy$hyperParams$trainData$x[,1], y=cmfm_learner_xy$hyperParams$trainData$y, beta=cmfm_learner_xy$learnParams$classifier$learnParams$model$beta[,1]))
summary(dfBeta_xy)
p <- ggplot(dfBeta_xy)
p <- p + geom_point(aes(x=x, y=y, col=abs(beta)))
p
dfBeta_yx <- as.data.frame(cbind(x=cmfm_learner_yx$hyperParams$trainData$x[,1], y=cmfm_learner_yx$hyperParams$trainData$y, beta=cmfm_learner_yx$learnParams$classifier$learnParams$model$beta[,1]))
summary(dfBeta_yx)
p <- ggplot(dfBeta_yx)
p <- p + geom_point(aes(x=x, y=y, col=abs(beta)))
p

pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy, data=testDataXY, forLoss=T)
pred_yx <- cmfm_learner_yx$pred(learner=cmfm_learner_yx, data=testDataYX, forLoss=T)

plotPredCMEM(cmfm_learner_xy, pred=pred_xy, var="k",indx=1:9)
plotPredCMEM(cmfm_learner_yx, pred_yx, var="k",indx=1:9)

log(getHyperPar(cmfm_learner_xy, "lambda"), 10)
log(getHyperPar(cmfm_learner_yx, "lambda"), 10)
log(getHyperPar(cmfm_learner_xy, "sigma.rbf.X"), 10)
log(getHyperPar(cmfm_learner_yx, "sigma.rbf.X"), 10)
log(getHyperPar(cmfm_learner_xy, "sigma.rbf.Y"), 10)
log(getHyperPar(cmfm_learner_yx, "sigma.rbf.Y"), 10)


# obtain loss out of sample
# loss <- "cmem_L2_k"
# this applies loss directly on predicted vector, no folds or CV
sapply(cmfm_learner_xy$optimizeParams$losses, function(loss) do.call(loss, list(learner=cmfm_learner_xy, pred=pred_xy)))
sapply(cmfm_learner_yx$optimizeParams$losses, function(loss) do.call(loss, list(learner=cmfm_learner_yx, pred=pred_yx)))


# calculate measures
pm <- proc.time()
msrs_xy <- cmfm_learner_xy$calcMsrs(cmemLearner=cmfm_learner_xy) # learner=cmfm_learner_xy
msrs_yx <- cmfm_learner_yx$calcMsrs(cmemLearner=cmfm_learner_yx)
proc.time() - pm # 0.415 secs
msrs_xy
msrs_yx


newSigmaXY <- cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y$val
newSigmaYX <- cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.Y$val

# now lets set the sigma parameter of xy and yx learners
# as the direct sum of the two optimal RKHS H_y  and compare
# measures 

newSigma <- unique(c(newSigmaXY, newSigmaYX))
log(newSigma, 10)

cmfm_learner <- cmfm_learner_pack_kernParsX_1[indx_lrn]
cmfm_learner <- eval(parse(text=cmfm_learner))
cmfm_learner$hyperParams$data$non_optimizable$sigma.rbf.Y$val <- newSigma

pm <- proc.time()
cmfm_learner_xy_2 <- setParams(learner=cmfm_learner, trainData=trainDataXY, plot=F)
cmfm_learner_yx_2 <- setParams(learner=cmfm_learner, trainData=trainDataYX, plot=F)
proc.time() - pm # 3.8 mins

fileXY <- paste(dataNm, "_2_1_", paste(indx_Y, indx_X,sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_2_1_", paste(indx_X, indx_Y,sep="on"), ".RData", sep="")
fileXY <- paste(dataNm, "_2_1_", paste("y", "x",sep="on"), ".RData", sep="")
fileYX <- paste(dataNm, "_2_1_", paste("x", "y",sep="on"), ".RData", sep="")

load(file=paste(folderLearners, fileXY, sep=""))
cmfm_learner_xy_2 <- cmemLearnerSumAux
load(file=paste(folderLearners, fileYX, sep=""))
cmfm_learner_yx_2 <- cmemLearnerSumAux




grid_xy <- cmfm_learner_xy_2$hyperParams$data$grid
grid_yx <- cmfm_learner_yx_2$hyperParams$data$grid

#L2
dimnames(grid_xy)$sigma.rbf.X <- seq(dim(grid_xy)["sigma.rbf.X"])
dimnames(grid_yx)$sigma.rbf.X <- seq(dim(grid_yx)["sigma.rbf.X"])


df_xy <- melt(grid_xy)
df_xy$dir <- "xy"
df_yx <- melt(grid_yx)
df_yx$dir <- "yx"
df <- rbind(df_xy, df_yx)

# bins
# indx <-  which(df$trainTest == "test" & df$var %in% c("cmem_L2_f", "negCE", "MisCR"))
# p <- ggplot(df[indx,])
# p <- p + geom_line(aes(x=num.bin.X, y=value, colour=dir))
# p <- p + geom_point(aes(x=num.bin.X, y=value, colour=dir))
# p <- p + facet_grid(var~., scales="free")
# p 


# L2

indx <-  which(df$trainTest == "test" & df$var %in% c("cmem_L2_f","MisCR","negCE"))
p <- ggplot(df[indx,])
p <- p + geom_line(aes(x=log(lambda, 10), y=value, colour=dir))
p <- p + geom_point(aes(x=log(lambda, 10), y=value, colour=dir))
p <- p + facet_grid(var~sigma.rbf.X, scales="free")
p 

df2 <- cast(df, trainTest+lambda+sigma.rbf.X+dir~var, value="value")
indx <-  which(df2$trainTest == "test" )
p <- ggplot(df2[indx,])
p <- p + geom_point(aes(x=log(KCMC, 10), y=cmem_L2_f, colour=dir), size=0.5, alpha=0.5)
p

head(df2)
min(df2$negCE[which(df2$trainTest=="test")])

# now massage data so that we can do KCDC vs entropy loss

names(cmfm_learner_xy$optimizeParams$losses)
loss <- "cmem_L2_f"
reglrs <- c("KCDC","KCMC","KCRDC","KCSC","KCNSC","KCCC_ent")
reglrs <- c("KCDC","KCMC","KCRDC","KCSC","KCNSC")
sapply(reglrs, function(reglr){
  # reglr <- reglrs[1]
  logReglr <- switch(reglr, KCCC_ent=F, T)
  curveXY <- getLossReglrCurve(learner=cmfm_learner_xy_2, loss=loss, reglr=reglr, logReglr=logReglr)
  curveYX <- getLossReglrCurve(learner=cmfm_learner_yx_2, loss=loss, reglr=reglr, logReglr=logReglr)
  curveList <- list(curveXY, curveYX)
  names(curveList) <- c("xy","yx")
  compReglrVals <- getComparableReglrVal(curveList)
  plotCurves(compReglrVals, curveList)
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
log(getHyperPar(cmfm_learner_xy_2, "sigma.rbf.X"), 10)
log(getHyperPar(cmfm_learner_yx_2, "sigma.rbf.X"), 10)
log(getHyperPar(cmfm_learner_xy_2, "sigma.rbf.Y"), 10)
log(getHyperPar(cmfm_learner_yx_2, "sigma.rbf.Y"), 10)

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

