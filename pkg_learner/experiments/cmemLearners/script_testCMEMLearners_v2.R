# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"


repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)


# simulation functions
source("./pkg_causaLearner/genData/func_getData_v2.R", echo=FALSE)
# learners
source("./pkg_learner/func_learners_v3.R", echo=FALSE)
source("./pkg_learner/ob_cme_learner_v2.R", echo=FALSE)
# for stdrize
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R", echo=FALSE)
library(reshape)
library(ggplot2)

p <- 2
nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
nodes <- rep(list(nodesX),p)
names(nodes) <- c("x","y")
dag <- matrix(c(0,0,1,0),2,2)
colnames(dag) <- rownames(dag) <- c("x","y")

n <- 100

# random function
set.seed(4)
simTest <- simRandSEM(p, n, nodes, sigma=2, sigmaErr=0, dagMat=dag)
plot(getGraph(simTest$dag))




X <- simTest$x
X <- apply(X, 2, stdrize)

apply(X, 2, mean)
apply(X, 2, sd)

x <- X[,"x"]
y <- X[,"y"]
plot(x,y)
plot(y,x)

trainDataXY <- constructData(as.matrix(x), y)
trainDataYX <- constructData(as.matrix(y), x)

# train hyperparameters
cmem_rbf_rbf_L2_none_xy <- setParams(learner=cmem_rbf_rbf_L2_none, trainData=trainDataXY, plot=FALSE)
cmem_rbf_rbf_L2_none_yx <- setParams(learner=cmem_rbf_rbf_L2_none, trainDataYX, plot=FALSE)
# train learn parameters
cmem_rbf_rbf_L2_none_xy <- cmem_rbf_rbf_L2_none$learn(learner=cmem_rbf_rbf_L2_none_xy)
cmem_rbf_rbf_L2_none_yx <- cmem_rbf_rbf_L2_none$learn(cmem_rbf_rbf_L2_none_yx)
# calculate measures
pm <- proc.time()
cmem_rbf_rbf_L2_none_xy$calcMsrs(cmem_rbf_rbf_L2_none_xy)
cmem_rbf_rbf_L2_none_yx$calcMsrs(cmem_rbf_rbf_L2_none_yx)
proc.time() - pm # 0.415 secs


KCDC2_1 <- function(x,y, lambda){
  
  sigma0 <- 1/median(as.numeric(dist(x)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  
  
  n <- length(x)
  L  <- kernlab:::kernelMatrix(kernelX, x) 
  K  <- kernlab:::kernelMatrix(kernelY, y)
  I <- diag(n)
  
  
  Blambda <- base:::solve(as.matrix(L)+n*lambda*I)
  Alambda <- Blambda%*%K%*%Blambda
  
 
  
  LAL <- L%*%Alambda%*%L
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  res <- (c/n) - (b/n)^2
  
  return(res)
}

KCDCrel2_1 <- function(x, y, lambda, numPerms){
  
  sigma0 <- 1/median(as.numeric(dist(x)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(x)
  set.seed(12345)            
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  
  mesr <- KCDC2_1(x, y, lambda)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y,y, lambda)))
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    
    res <- KCDC2_1(x, y[rperm], lambda)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y[rperm],y[rperm], lambda))) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}

n <- length(x)
lambda <- cmem_rbf_rbf_L2_none$hyperParams$data$non_optimizable$lambda$val
numPerms <- cmem_rbf_rbf_L2_none$msrs$KCDCrel$pars$numPerms

pm <- proc.time()
KCDC2_1(x,y, lambda)
KCDC2_1(y,x, lambda)
KCDCrel2_1(x,y,lambda, numPerms)
KCDCrel2_1(y,x,lambda, numPerms)
proc.time() - pm # 0.757 secs



# set lambda parameter according to "Conditional mean embeddings as regressors" Grünewälder et al

pm <- proc.time()
# train hyperparameters
cmem_rbf_rbf_L2_lambda_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainData=trainDataXY, plot=TRUE)
cmem_rbf_rbf_L2_lambda_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataYX, plot=TRUE)
proc.time() - pm # 2.072 secs

getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "sigmaY")

getHyperPar(cmem_rbf_rbf_L2_lambda_yx, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_yx, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_yx, "sigmaY")

# train learn parameters
learner <- cmem_rbf_rbf_L2_lambda_xy
cmem_rbf_rbf_L2_lambda_xy <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_xy)
cmem_rbf_rbf_L2_lambda_yx <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_yx)
# calculate measures
cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
 

vecRegLoss1 <- function(xTr,yTr, xTe=NULL, yTe=NULL, lambda){
  
  sigma0 <- 1/median(as.numeric(dist(xTr)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(yTr)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(xTr)
  
  LTr  <- kernlab:::kernelMatrix(kernelX, xTr) 
  KTr  <- kernlab:::kernelMatrix(kernelY, yTr)
  I <- diag(n)
  
  if(!is.null(xTe)){
    nTe <- length(xTe)
    LTrTe  <- kernlab:::kernelMatrix(kernelX, xTr, xTe) 
    KTrTe  <- kernlab:::kernelMatrix(kernelY, yTr, yTe)
    KTe <- kernlab:::kernelMatrix(kernelY, yTe)
  } else{
    nTe <- n
    LTrTe <- LTr
    KTrTe <- KTr
    KTe <- KTr
  }
  
  
  
  Blambda <- solve(LTr+n*lambda*I)
  Alambda <- Blambda%*%KTr%*%Blambda
  
  LAL <- t(LTrTe)%*%Alambda%*%LTrTe
  LBK <- t(LTrTe)%*%Blambda%*%KTrTe
  
  
  res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}

vecRegLoss1(x,y,lambda=lambda)

CV.parallel2_1 <- function(x, y, numFolds, lambdas, fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  trainData <- cbind(x=x,y=y)

  nParams <- length(lambdas)
  dimnames <- list(as.character(1:numFolds), lambdas)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- lambdas[1]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      lossesTrain <- vecRegLoss1(xTr=curTrain[,"x"], yTr=curTrain[,"y"],  lambda=p)
      lossesTest <- vecRegLoss1(xTr=curTrain[,"x"], yTr=curTrain[,"y"], xTe=curTest[,"x"], yTe=curTest[,"y"], lambda=p)
      res <- c(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, p=lambdas, mc.cores=2, SIMPLIFY="array")
  

  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=lambdas)
  
  #print("exits CV parallel")
  return(losses)
}

lambdas <- cmem_rbf_rbf_L2_lambda$hyperParams$data$optimizable$lambda$seq
numFolds <- cmem_rbf_rbf_L2_lambda$optimizeParams$numFolds

pm <- proc.time()
lossesXY <- CV.parallel2_1(x, y, numFolds, lambdas)
lossesMeanXY <- apply(lossesXY, c("trainTest", "params"), mean)
lossesYX <- CV.parallel2_1(y, x, numFolds, lambdas)
lossesMeanYX <- apply(lossesYX, c("trainTest", "params"), mean)
proc.time() - pm # 4.858 secs

# compare lambdas
getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "lambda")
which.min(lossesMeanXY["test",])
getHyperPar(cmem_rbf_rbf_L2_lambda_yx, "lambda")
which.min(lossesMeanYX["test",])

# compare sigma and gamma
# sigma
getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "sigmaX")
1/median(as.numeric(dist(x)^2))
# gamma
getHyperPar(cmem_rbf_rbf_L2_lambda_xy, "sigmaY")
1/median(as.numeric(dist(y)^2))

# compare KCDC and KCDC rel
# x->y
cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
KCDC2_1(x, y, lambda=as.numeric(names(which.min(lossesMeanXY["test",]))))
KCDCrel2_1(x,y, lambda=as.numeric(names(which.min(lossesMeanXY["test",]))), numPerms)

# x <- y
cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
KCDC2_1(y,x, as.numeric(names(which.min(lossesMeanYX["test",]))))
KCDCrel2_1(y,x, as.numeric(names(which.min(lossesMeanYX["test",]))), numPerms)


# cross validate lambda, sigma 
# train hyperparameters 
pm <- proc.time()
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainData=trainDataXY, plot=FALSE)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainData=trainDataYX, plot=FALSE)
proc.time() - pm # 25 seconds
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaY")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "sigmaY")
# train learn parameters
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
# calculate measures
cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)


vecRegLoss_2 <- function(xTr,yTr, xTe=NULL, yTe=NULL, lambda, sigma){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigma))
  sigma0 <- 1/median(as.numeric(dist(yTr)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(xTr)
  
  LTr  <- kernlab:::kernelMatrix(kernelX, xTr) 
  KTr  <- kernlab:::kernelMatrix(kernelY, yTr)
  I <- diag(n)
  
  if(!is.null(xTe)){
    nTe <- length(xTe)
    LTrTe  <- kernlab:::kernelMatrix(kernelX, xTr, xTe) 
    KTrTe  <- kernlab:::kernelMatrix(kernelY, yTr, yTe)
    KTe <- kernlab:::kernelMatrix(kernelY, yTe)
  } else{
    nTe <- n
    LTrTe <- LTr
    KTrTe <- KTr
    KTe <- KTr
  }
  
  
  
  Blambda <- solve(LTr+n*lambda*I)
  Alambda <- Blambda%*%KTr%*%Blambda
  
  LAL <- t(LTrTe)%*%Alambda%*%LTrTe
  LBK <- t(LTrTe)%*%Blambda%*%KTrTe
  
  
  res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}

vecRegLoss_2(x,y,lambda=lambda, sigma=1)

lambdas <- cmem_rbf_rbf_L2_lambda_kernParsX$hyperParams$data$optimizable$lambda$seq
sigmasXY <-  cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$sigmaX$seq
sigmasYX <-  cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$sigmaX$seq
numFolds <- cmem_rbf_rbf_L2_lambda$optimizeParams$numFolds

log(lambdas,10)

sigmas <- sigmasXY

CV.parallel2_2 <- function(x, y, numFolds, lambdas, sigmas,  fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  trainData <- cbind(x=x,y=y)
  
  params <- constructParams(lambda=lambdas, sigma=sigmas)
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), lambdas, sigmas)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- params[[1]]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      lossesTrain <- vecRegLoss_2(xTr=curTrain[,"x"], yTr=curTrain[,"y"],  lambda=p$lambda, sigma=p$sigma)
      lossesTest <- vecRegLoss_2(xTr=curTrain[,"x"], yTr=curTrain[,"y"], xTe=curTest[,"x"], yTe=curTest[,"y"], lambda=p$lambda, sigma=p$sigma)
      res <- c(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, p=params, mc.cores=2, SIMPLIFY="array")
  
  
  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=names(params))
  
  #print("exits CV parallel")
  return(losses)
}


pm <- proc.time()
lossesXY <- CV.parallel2_2(x, y, numFolds, lambdas, sigmasXY)
lossesMeanXY <- apply(lossesXY, c("trainTest", "params"), mean)
lossesYX <- CV.parallel2_2(y, x, numFolds, lambdas, sigmasYX)
lossesMeanYX <- apply(lossesYX, c("trainTest", "params"), mean)
proc.time() - pm # 25 secs

dim(lossesMeanXY)

gridXY <- cbind(expand.grid(lambda=lambdas, sigma=sigmasXY), loss=lossesMeanXY["test",])
parsXY <- as.numeric(sapply(strsplit(strsplit(names(which.min(lossesMeanXY["test",])), " ")[[1]], "="), function(el) el[2]))
v <- ggplot(gridXY, aes(log(lambda, 10), log(sigma, 10), z = log(loss, 10)))
v <- v + geom_raster(aes(fill = log(loss, 10))) 
v <- v + geom_contour(colour = "white", bins = 10)
v <- v + geom_point(aes(x=log(parsXY[1],10), y=log(parsXY[2],10)), colour="green", size=2)
v

gridYX <- cbind(expand.grid(lambda=lambdas, sigma=sigmasYX), loss=lossesMeanYX["test",])
parsYX <- as.numeric(sapply(strsplit(strsplit(names(which.min(lossesMeanYX["test",])), " ")[[1]], "="), function(el) el[2]))
v <- ggplot(gridYX, aes(log(lambda, 10), log(sigma, 10), z = log(loss, 10)))
v <- v + geom_raster(aes(fill = log(loss, 10))) 
v <- v + geom_contour(colour = "white", bins = 10)
v <- v + geom_point(aes(x=log(parsYX[1],10), y=log(parsYX[2],10)), colour="green", size=2)
v


# compare lambdas
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "lambda"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX"),10)
log(parsXY,10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "lambda"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "sigmaX"),10)
log(parsYX,10) 
# difference due to slighlty different loss calculations due to different kernel implementations
# kernlab vs my implementation

# compare  gamma

# gamma
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaY")
1/median(as.numeric(dist(y)^2))

KCDC2_2 <- function(x,y, lambda, sigma){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigma))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  # sigma0; cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$sigmaY
  
  
  n <- length(x)
  L  <- kernlab:::kernelMatrix(kernelX, x) 
  K  <- kernlab:::kernelMatrix(kernelY, y)
  # plot(L, cmem_rbf_rbf_L2_lambda_kernParsX_yx$learnParams$Lx); abline(a=0,b=1, col="red")
  # plot(K, cmem_rbf_rbf_L2_lambda_kernParsX_yx$learnParams$Ky); abline(a=0,b=1, col="red")
  
  
  I <- diag(n)
  
  cmem_rbf_rbf_L2_lambda_kernParsX_yx$learnParams$Blambda
  Blambda <- base:::solve(as.matrix(L)+n*lambda*I)
  # plot(Blambda, cmem_rbf_rbf_L2_lambda_kernParsX_yx$learnParams$Blambda); abline(a=0,b=1, col="red")
  Alambda <- Blambda%*%K%*%Blambda
  
  
  
  LAL <- L%*%Alambda%*%L
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  res <- (c/n) - (b/n)^2
  
  return(res)
}

KCDCrel2_2 <- function(x, y, lambda, sigma, numPerms){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigma))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(x)
  set.seed(12345)            
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  
  mesr <- KCDC2_2(x, y, lambda, sigma)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y,y, lambda)))
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    
    res <- KCDC2_2(x, y[rperm], lambda, sigma)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y[rperm],y[rperm], lambda))) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}


# compare KCDC and KCDC rel
# x->y
cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
KCDC2_2(x, y, lambda=parsXY[1], sigma=parsXY[2])
KCDCrel2_2(x,y, lambda=parsXY[1], sigma=parsXY[2], numPerms)

# x <- y
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
KCDC2_2(y, x, lambda=parsYX[1], sigma=parsYX[2])
KCDCrel2_2(y, x, lambda=parsYX[1], sigma=parsYX[2], numPerms)
# we see slightly different results because of difference in selected lambda

# cross validate lambda, sigma and gamma
# train hyperparameters 
pm <- proc.time()
cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainData=trainDataXY, plot=FALSE)
cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataYX, plot=FALSE)
proc.time() - pm # 20 secs
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "lambda")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "sigmaY")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "sigmaY")
# train learn parameters
cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
# calculate measures
cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
cmem_rbf_rbf_L2_lambda_kernParsXY_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)

vecRegLoss_3 <- function(xTr,yTr, xTe=NULL, yTe=NULL, lambda, sigmaX, sigmaY){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigmaX))
  
  kernelY <- do.call("rbfdot", list(sigma=sigmaY))
  
  n <- length(xTr)
  
  LTr  <- kernlab:::kernelMatrix(kernelX, xTr) 
  KTr  <- kernlab:::kernelMatrix(kernelY, yTr)
  I <- diag(n)
  
  if(!is.null(xTe)){
    nTe <- length(xTe)
    LTrTe  <- kernlab:::kernelMatrix(kernelX, xTr, xTe) 
    KTrTe  <- kernlab:::kernelMatrix(kernelY, yTr, yTe)
    KTe <- kernlab:::kernelMatrix(kernelY, yTe)
  } else{
    nTe <- n
    LTrTe <- LTr
    KTrTe <- KTr
    KTe <- KTr
  }
  
  
  
  Blambda <- solve(LTr+n*lambda*I)
  Alambda <- Blambda%*%KTr%*%Blambda
  
  LAL <- t(LTrTe)%*%Alambda%*%LTrTe
  LBK <- t(LTrTe)%*%Blambda%*%KTrTe
  
  dKTe <- diag(KTe)
  res <- diag(LAL) + dKTe - 2*diag(LBK)
  res <- sum(res/dKTe)
  
  #res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}

vecRegLoss_3(x,y,lambda=lambda, sigmaX=1, sigmaY=3)

lambdas <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$hyperParams$data$optimizable$lambda$seq
sigmasX_XY <-  cmem_rbf_rbf_L2_lambda_kernParsXY_xy$hyperParams$data$optimizable$sigmaX$seq
sigmasX_YX <-  cmem_rbf_rbf_L2_lambda_kernParsXY_yx$hyperParams$data$optimizable$sigmaX$seq
sigmasY_XY <-  cmem_rbf_rbf_L2_lambda_kernParsXY_xy$hyperParams$data$optimizable$sigmaY$seq
sigmasY_YX <-  cmem_rbf_rbf_L2_lambda_kernParsXY_yx$hyperParams$data$optimizable$sigmaY$seq
numFolds <- cmem_rbf_rbf_L2_lambda$optimizeParams$numFolds



CV.parallel2_3 <- function(x, y, numFolds, lambdas, sigmasX, sigmasY,  fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  trainData <- cbind(x=x,y=y)
  
  params <- constructParams(lambda=lambdas, sigmaX=sigmasX, sigmaY=sigmasY)
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), lambdas, sigmasX, sigmasY)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- params[[1]]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      lossesTrain <- vecRegLoss_3(xTr=curTrain[,"x"], yTr=curTrain[,"y"],  lambda=p$lambda, sigmaX=p$sigmaX, sigmaY=p$sigmaY)
      lossesTest <- vecRegLoss_3(xTr=curTrain[,"x"], yTr=curTrain[,"y"], xTe=curTest[,"x"], yTe=curTest[,"y"], lambda=p$lambda, sigmaX=p$sigmaX, sigmaY=p$sigmaY)
      res <- c(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, p=params, mc.cores=2, SIMPLIFY="array")
  
  
  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=names(params))
  
  #print("exits CV parallel")
  return(losses)
}


pm <- proc.time()
lossesXY <- CV.parallel2_3(x, y, numFolds, lambdas, sigmasX_XY, sigmasY_XY)
lossesMeanXY <- apply(lossesXY, c("trainTest", "params"), mean)
lossesYX <- CV.parallel2_3(y, x, numFolds, lambdas, sigmasX_YX, sigmasY_YX)
lossesMeanYX <- apply(lossesYX, c("trainTest", "params"), mean)
proc.time() - pm # 25 secs

dim(lossesMeanXY)
dimnames(lossesMeanXY)

gridXY <- cbind(expand.grid(lambda=lambdas, sigmaX=sigmasX_XY, sigmaY=sigmasY_XY), loss=lossesMeanXY["test",])
parsXY <- as.numeric(sapply(strsplit(strsplit(names(which.min(lossesMeanXY["test",])), " ")[[1]], "="), function(el) el[2]))
v <- ggplot(gridXY, aes(log(lambda, 10), log(sigmaX, 10), z = log(loss, 10)))
v <- v + geom_raster(aes(fill = log(loss, 10))) 
v <- v + geom_contour(colour = "white", bins = 10)
v <- v + geom_point(aes(x=log(parsXY[1],10), y=log(parsXY[2],10)), colour="green", size=2)
v <- v + facet_wrap(~sigmaY)
v

gridYX <- cbind(expand.grid(lambda=lambdas, sigmaX=sigmasX_YX, sigmaY=sigmasY_YX), loss=lossesMeanYX["test",])
rownames(gridYX) <- seq(nrow(gridYX))
parsYX <- as.numeric(sapply(strsplit(strsplit(names(which.min(lossesMeanYX["test",])), " ")[[1]], "="), function(el) el[2]))
v <- ggplot(gridYX, aes(log(lambda, 10), log(sigmaX, 10), z = log(loss, 10)))
v <- v + geom_raster(aes(fill = log(loss, 10))) 
v <- v + geom_contour(colour = "white", bins = 10)
v <- v + geom_point(aes(x=log(parsYX[1],10), y=log(parsYX[2],10)), colour="green", size=2)
v <- v + facet_wrap(~sigmaY)
v


# compare lambdas, sigmas and gammas
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "lambda"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "sigmaX"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_xy, "sigmaY"),10)
log(parsXY,10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "lambda"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "sigmaX"),10)
log(getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsXY_yx, "sigmaY"),10)
log(parsYX,10) #again prolly numerical instability



KCDC2_3 <- function(x,y, lambda, sigmaX, sigmaY){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigmaX))
  kernelY <- do.call("rbfdot", list(sigma=sigmaY))
  #sigmaY; cmem_rbf_rbf_L2_lambda_kernParsXY_xy$hyperParams$data$optimizable$sigmaY
  #sigmaX; cmem_rbf_rbf_L2_lambda_kernParsXY_xy$hyperParams$data$optimizable$sigmaX
  
  n <- length(x)
  L  <- kernlab:::kernelMatrix(kernelX, x) 
  K  <- kernlab:::kernelMatrix(kernelY, y)
  # plot(L, cmem_rbf_rbf_L2_lambda_kernParsXY_xy$learnParams$Lx); abline(a=0,b=1, col="red")
  # plot(K, cmem_rbf_rbf_L2_lambda_kernParsXY_xy$learnParams$Ky); abline(a=0,b=1, col="red")
  
  
  I <- diag(n)
  
  
  Blambda <- base:::solve(as.matrix(L)+n*lambda*I)
  # plot(Blambda, cmem_rbf_rbf_L2_lambda_kernParsXY_xy$learnParams$Blambda); abline(a=0,b=1, col="red")
  Alambda <- Blambda%*%K%*%Blambda
  
  
  
  LAL <- L%*%Alambda%*%L
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  res <- (c/n) - (b/n)^2
  
  return(res)
}

KCDCrel2_3 <- function(x, y, lambda, sigmaX, sigmaY, numPerms){
  
  
  kernelX <- do.call("rbfdot", list(sigma=sigmaX))
  kernelY <- do.call("rbfdot", list(sigma=sigmaY))
  
  n <- length(x)
  set.seed(12345)            
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  
  mesr <- KCDC2_3(x, y, lambda, sigmaX, sigmaY)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y,y, lambda)))
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    
    res <- KCDC2_3(x, y[rperm], lambda, sigmaX, sigmaY)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y[rperm],y[rperm], lambda))) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}


# compare KCDC and KCDC rel
# x->y
cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
KCDC2_3(x, y, lambda=parsXY[1], sigmaX=parsXY[2], sigmaY=parsXY[3])
KCDCrel2_3(x,y, lambda=parsXY[1], sigmaX=parsXY[2], sigmaY=parsXY[3], numPerms)
# differences due to difference in kernel implementation

# x <- y
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
KCDC2_3(y, x, lambda=parsYX[1], sigmaX=parsYX[2], sigmaY=parsYX[3])
KCDCrel2_3(y, x, lambda=parsYX[1], sigmaX=parsYX[2], sigmaY=parsYX[3], numPerms)
# sligthly different because of slightly different calculations by using one type of kernel matrix
# implementation or the other (kernlab vs mine)



# Lets see for 100 datasets what KCDC and KCDCrel is 

q <- 100
n <- 200
set.seed(4)
p <- 2
(ps <- rep(p, q))
(ns <- rep(n, q))
nodes <- list(dist="runif", pars=list(min=-2, max=2), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- c("x","y")
nodess <- lapply(p, function(p) nodes)
dag <- matrix(c(0,0,1,0),2,2)
colnames(dag) <- rownames(dag) <- c("x","y")

# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
# random functions
set.seed(5)
dataTestList <- simRandSEMs(q, ps, ns, nodess, sigma=5, sigmaErr=0, dagMat=dag)
plotPairsList(dataTestList)

# chosen function
edgL <- vector("list", length=2)
names(edgL) <- c("x","y")
edgL[["x"]] <- list(edges=c("y"))
XtoY <- graphNEL(nodes=c("x","y"), edgeL=edgL, edgemode="directed")
plot(XtoY)
nodesUnif <-  list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                   y=list(dist="runif", pars=list(min=-1, max=1), a=1, b=1))
funcs <- list(fx=function(n) n, fy=function(x, n) sin(x)*n) #sin(x)*n
sem <- list(dag=XtoY, funcs=funcs, simPars=list(n=100, nodes=nodesUnif))
dataTestList <- simSEMs(q, sem)
plotPairsList(dataTestList)


pm <- proc.time()
msrs <- mcmapply(function(el, nm){
  # i <- 9; el <- dataTestList$xs[[i]]; nm <- dataTestList$names[i]
  X <- apply(el, 2, norml)
  
  print(paste("name: ", nm))
  #print("head(X)")
  #print(head(X))
  #print(apply(X, 2, mean))
  #print(apply(X, 2, sd))
  
  x <- X[,1]
  y <- X[,2]
  
  trainDataXY <- constructData(as.matrix(x), y)
  trainDataYX <- constructData(as.matrix(y), x)
  
  # train hyperparameters
  cmem_rbf_rbf_L2_none_xy <- setParams(learner=cmem_quad_quad_L2_none, trainData=trainDataXY, plot=FALSE)
  cmem_rbf_rbf_L2_none_yx <- setParams(learner=cmem_quad_quad_L2_none, trainDataYX, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataYX, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataYX, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataYX, plot=FALSE)
  
  # train learn parameters
  cmem_rbf_rbf_L2_none_xy <- cmem_rbf_rbf_L2_none$learn(cmem_rbf_rbf_L2_none_xy)
  cmem_rbf_rbf_L2_none_yx <- cmem_rbf_rbf_L2_none$learn(cmem_rbf_rbf_L2_none_yx)
  #cmem_rbf_rbf_L2_lambda_xy <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_xy)
  #cmem_rbf_rbf_L2_lambda_yx <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  # calculate measures
  msrsFixXY <- cmem_rbf_rbf_L2_none_xy$calcMsrs(cmem_rbf_rbf_L2_none_xy)
  msrsFixYX <- cmem_rbf_rbf_L2_none_yx$calcMsrs(cmem_rbf_rbf_L2_none_yx)
  #msrsOpt1XY <- cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
  #msrsOpt1YX <- cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
  #msrsOpt2XY <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  #msrsOpt2YX <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #msrsOpt3XY <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #msrsOpt3YX <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  
  #KCDC lambda fix
  KCDCq <- msrsFixYX["KCDC"]/msrsFixXY["KCDC"] 
  #KCDCrel lambda fix
  KCDCrelq <- msrsFixYX["KCDCrel"]/msrsFixXY["KCDCrel"]
  
  #KCDC lambda opt1
  #KCDCq_opt1 <- msrsOpt1YX["KCDC"]/msrsOpt1XY["KCDC"]
  #KCDCrel lambda opt1
  #KCDCrelq_opt1 <- msrsOpt1YX["KCDCrel"]/msrsOpt1XY["KCDCrel"]
  
  #KCDC lambda opt2
  #KCDCq_opt2 <- msrsOpt2YX["KCDC"]/msrsOpt2XY["KCDC"]
  #KCDCrel lambda opt2
  #KCDCrelq_opt2 <- msrsOpt2YX["KCDCrel"]/msrsOpt2XY["KCDCrel"]
  
  #KCDC lambda opt3
  #KCDCq_opt3 <- msrsOpt3YX["KCDC"]/msrsOpt3XY["KCDC"]
  #KCDCrel lambda opt3
  #KCDCrelq_opt3 <- msrsOpt3YX["KCDCrel"]/msrsOpt3XY["KCDCrel"]
  
  res <- c(KCDCq, KCDCrelq)
  names(res) <- c("KCDC","KCDCrel")
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2, KCDCq_opt3, KCDCrelq_opt3)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2","KCDCopt3","KCDCrelopt3")
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2")
  return(res)
}, el=dataTestList$xs, nm=dataTestList$names, 
SIMPLIFY="array", mc.cores=6)
proc.time() - pm 

# 36 mins for none, lambda and lambda_kernParsX
# 90 mins for none, lambda, lambda_kernParsX and lambda_kernParsXY

names(dimnames(msrs)) <- c("measure", "database")

#save("msrs", file="/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_learner/experiments/cmemLearners/testCMEMlearners.RData")

msrsDB <- melt(msrs)

pctRight <- function(x){
  
  res <- sum(x>1)/length(x)*100
  return(res)
}
  

cast(msrsDB, measure~., value="value", fun.aggregate="pctRight")



