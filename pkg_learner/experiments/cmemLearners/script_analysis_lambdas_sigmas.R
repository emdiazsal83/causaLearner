# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)

hs_cmem_ob_version <- "v3_comp"
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

experimentNames <- c("dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegDists",
                     "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegCME",
                     "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegCMEnorm",
                     "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegCMEodds",
                     "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-CMEprob")

indxs <- sapply(experimentNames, function(experimentName){
  print("********************")
  print(experimentName)
  #experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-logRegDists"
  dataName <- strsplit(experimentName, "_")[[1]][1]
  #dataName <- "dag2-ME2-Cmplx-SinPlus-SinTimes"

  indxs <- sapply(1:300, function(block){
    #block <- 250 # 4 (fails comp msr, hit NCE - logRegCME)
  pm <- proc.time()
  load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
  dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
  proc.time() - pm # 0.058 secs
  dataNm <- dat$name
  dat$dag

  indx_X <- c(1,2)[match(paste(as.numeric(dat$dag), collapse="."),  
                      c(paste(c(0,0,1,0), collapse="."),paste(c(0,1,0,0), collapse=".")))]
  indx_Y <- setdiff(c(1,2), indx_X)


  expType <- toupper(strsplit(strsplit(experimentName, "_")[[1]][2], "-")[[1]][1])
  dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", sep=""))
  dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, sep=""))
  folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, "/",sep="")
  fileXY_rand <- paste(dataNm, "_1_", paste(indx_Y, indx_X,sep="on"), ".RData", sep="")
  fileYX_rand <- paste(dataNm, "_1_", paste(indx_X, indx_Y,sep="on"), ".RData", sep="")
  fileXY_other <- paste(dataNm, "_1_", paste("y", "x",sep="on"), ".RData", sep="")
  fileYX_other <- paste(dataNm, "_1_", paste("x", "y",sep="on"), ".RData", sep="")

  if( substr(dataNm, 1, 4) == "rand"){
    fileXY <- fileXY_rand
    fileYX <- fileYX_rand
  } else{
    fileXY <- fileXY_other
    fileYX <- fileYX_other
  }

  dir(folderLearners)
  load(file=paste(folderLearners, fileXY, sep=""))
  cmfm_learner_xy <- cmemLearnerAux
  load(file=paste(folderLearners, fileYX, sep=""))
  cmfm_learner_yx <- cmemLearnerAux

  sigma.rbf.Y_xy <- cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y$val
  sigma.rbf.Y_yx <- cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.Y$val
  
  
  
  indx_sigma.rbf.Y_xy <- which.min(abs(sigma.rbf.Y_xy-cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y$seq))
  indx_sigma.rbf.Y_yx <- which.min(abs(sigma.rbf.Y_yx-cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.Y$seq))

  sigma.rbf.X_xy <- cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.X$val
  sigma.rbf.X_yx <- cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.X$val
  
  indx_sigma.rbf.X_xy <- which.min(abs(sigma.rbf.X_xy-cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.X$seq))
  indx_sigma.rbf.X_yx <- which.min(abs(sigma.rbf.X_yx-cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.X$seq))

  q.rbf.Y_xy <- sum(sort(as.numeric(dist(cmfm_learner_xy$hyperParams$trainData$y)^2))<1/sigma.rbf.Y_xy)/length(cmfm_learner_xy$hyperParams$trainData$y)
  q.rbf.Y_yx <- sum(sort(as.numeric(dist(cmfm_learner_yx$hyperParams$trainData$y)^2))<1/sigma.rbf.Y_yx)/length(cmfm_learner_yx$hyperParams$trainData$y)
  
  q.rbf.X_xy <- sum(sort(as.numeric(dist(cmfm_learner_xy$hyperParams$trainData$x)^2))<1/sigma.rbf.X_xy)/length(cmfm_learner_xy$hyperParams$trainData$y)
  q.rbf.X_yx <- sum(sort(as.numeric(dist(cmfm_learner_yx$hyperParams$trainData$x)^2))<1/sigma.rbf.X_yx)/length(cmfm_learner_yx$hyperParams$trainData$y)
  
  
  indx_lambda_xy <- which.min(abs(cmfm_learner_xy$hyperParams$data$optimizable$lambda$val-cmfm_learner_xy$hyperParams$data$optimizable$lambda$seq))
  indx_lambda_yx <- which.min(abs(cmfm_learner_yx$hyperParams$data$optimizable$lambda$val-cmfm_learner_yx$hyperParams$data$optimizable$lambda$seq))
  res <- matrix(c(indx_sigma.rbf.Y_xy, indx_sigma.rbf.Y_yx, 
                  indx_sigma.rbf.X_xy, indx_sigma.rbf.X_yx, 
                  q.rbf.Y_xy, q.rbf.Y_yx, 
                  q.rbf.X_xy, q.rbf.X_yx,
                  indx_lambda_xy, indx_lambda_yx), 2, 5)  
  dimnames(res) <- list(dir=c("xy","yx"), par=c("sigma.rbf.Y","sigma.rbf.X", "q.rbf.Y","q.rbf.X","lambda"))
  return(res)
}, simplify="array")
  
  dimnames(indxs)[[3]] <- 1:300
  names(dimnames(indxs))[3] <- "dataset"
  
  table(as.numeric(indxs[,"sigma.rbf.Y",]))
  table(as.numeric(indxs[,"sigma.rbf.X",]))
  table(as.numeric(indxs[,"lambda",]))
  
  return(indxs)
  
}, simplify="array")

names(dimnames(indxs))[4] <- "learner"
dimnames(indxs)
dim(indxs)

table(as.numeric(indxs[,"sigma.rbf.Y",,]))
table(as.numeric(indxs[,"sigma.rbf.X",,]))
table(as.numeric(indxs[,"lambda",,]))

hist(as.numeric(indxs[,"q.rbf.Y",,]))
max(as.numeric(indxs[,"q.rbf.Y",,]))

hist(as.numeric(indxs[,"q.rbf.X",,]))
max(as.numeric(indxs[,"q.rbf.X",,]))

sapply(dimnames(indxs)$learner , function(learner) 
  table(as.numeric(indxs["xy","sigma.rbf.Y",,learner]), as.numeric(indxs["yx","sigma.rbf.Y",,learner])))
  


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

