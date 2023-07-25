# learner package

library(CVST) # cross validation framework
library(kernlab) # kernels, krr
library(FNN) # knnx.index
library(gptk) #gaussian process toolkit
library(lbfgs) # non-convex optimization for hsic regression
library(numDeriv) #comparing numerical approx of gradient with analytical
library(earth) # mars model for initial residual bandwidth estimation
#library(dHSIC) # dhsic (in getFixedParams)
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_dHSIC/dHSIC.R")
library(parallel) #mcmapply (in CV.parallel)

#############################################################################################################*
# learners

# this method should really be a method of kernel learners but since we are only using a 
# general pseudo-learner class well leave it outside for now
makeKernel <- function(learner, kernelName){
  kernel <- learner$hyperParams$non_data[[kernelName]]
  hyperParams <- learner$getHyperParams(learner)
  hyperParamsNames <- names(hyperParams) 
  kernelParamsNames <- names(formals(learner$hyperParams$non_data[[kernelName]]))
  indx <- match(kernelParamsNames, hyperParamsNames)
  kernelParams <- hyperParams[indx]
  kernelParams <- lapply(kernelParams, function(el){
    if(is.list(el)){
      res <- el$val
    } else{
      res <- el
    }
    return(res)
  })
  kernel <- do.call(kernel, kernelParams)
  return(kernel)
}  


learn.vanilla <- function(learner) {
  return(learner)
}
predict.vanilla <- function(learner, data) {
  
  gy <- as.matrix(data$y)
  x <- data$x
  return(list(x=as.matrix(x), gy=gy, gyh=matrix(0, nrow(gy), ncol(gy))))
}


learn.krr <- function(learner) {
  
  kernelXs <- makeKernel(learner, kernelName="kernelXs")
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  Kxs <- kernelMatrix(kernelXs, x)
  N <- nrow(Kxs)
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  alpha <- solve(Matrix(Kxs + diag(lambda, N))) %*% y
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.krr <- function(learner, data) {
  kernelXs <- makeKernel(learner, "kernelXs")
  
  kxs <- kernelMatrix(kernelXs, data$x, learner$hyperParams$trainData$x)
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}

learn.qhsic <- function (learner) {
  
  kernelXs <- makeKernel(learner, kernelName="kernelXs")
  kernelXb <- makeKernel(learner, kernelName="kernelXb")
  
  trainData <- learner$hyperParams$trainData
  
  mu_y <- mean(trainData$y)
  
  if(nrow(trainData$x) > learner$optimizeParams$maxPoints){
    wishInducing <- makeGrid(trainData$x, learner$optimizeParams$maxPoints)
    indxs <- as.numeric(knnx.index(data=trainData$x, query=wishInducing, k=1))
    trainData <- getSubset(trainData, indxs) 
  } else{
    indxs <- NA
  }
  
  x <- trainData$x
  y <- trainData$y
  Kxs <- kernelMatrix(kernelXs, x)
  Kxb <- kernelMatrix(kernelXb, x)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  alpha <- solve(Matrix(Kxbc%*%Kxs + diag(lambda, N))) %*% Kxbc %*% y
  learner$learnParams$alpha <- alpha
  learner$learnParams$avgy <- mu_y
  learner$learnParams$indxInducing <- indxs
  return(learner)
}
predict.qhsic <- function (learner, data){
  kernelXs <- makeKernel(learner, "kernelXs")
  trainData <- learner$hyperParams$trainData
  
  if(nrow(trainData$x) > learner$optimizeParams$maxPoints){
    trainData <- getSubset(trainData, learner$learnParams$indxInducing) 
  } else{
    indxs <- NA
  }
  
  kxs <- kernelMatrix(kernelXs, data$x, trainData$x)
  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}	

# auxiliary to find data points closest to being linearly spaced
makeGrid <- function(x, num){
  
  n <- nrow(x)
  p <- ncol(x)
  numPerDim <- ceiling(num^(1/p))
  extremes <- apply(x, 2, range)
  markPerDim <- apply(extremes, 2, function(col) seq(col[1], col[2], length.out=numPerDim))
  grid <- expand.grid(as.data.frame(markPerDim))
  grid <- grid[sample(1:nrow(grid),num),]
  return(grid)
}


.hsicRegLoss <- function(alpha, y, Kxs, Kxb, kernelRg, lambda){
  
  rs <- y-Kxs%*%alpha
  Krg <- kernelMatrix(kernelRg, rs)
  N = nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- Kxb%*%H
  Krgc <- Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc)) + lambda*t(alpha)%*%Kxs%*%alpha  
  return(res)
}
.hsicRegLossGrad <- function(alpha, y, Kxs, Kxb, kernelRg, lambda){
  rs <- y-Kxs%*%alpha
  Krg <- kernelMatrix(kernelRg, rs)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  ones <- matrix(1,N,1)
  aux <- Kxbc*Krg*(ones%*%t(rs)-rs%*%t(ones))
  gamma <- kernelRg@kpar$sigma
  
  grad <- apply(diag(N),2, function(ek)t(ones)%*%(aux*(ones%*%t(ek)%*%Kxs-Kxs%*%ek%*%t(ones)))%*%ones)
  
  
  
  grad <- 2*gamma*grad  + 2*lambda*Kxs%*%alpha
  return(grad)
}


learn.hsic <- function (learner) {
  #print("enters hsic learn function")
  kernelXs <- makeKernel(learner, kernelName="kernelXs")
  kernelXb <- makeKernel(learner, kernelName="kernelXb")
  kernelRg <- makeKernel(learner, kernelName="kernelRg")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  Kxs <- kernelMatrix(kernelXs, x)
  Kxb <- kernelMatrix(kernelXb, x)
  N <- nrow(Kxs)
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  max_iterations1 <- learner$optimizeParams$max_iterations1
  max_iterations2 <- learner$optimizeParams$max_iterations2
  num_init <- learner$optimizeParams$num_init
  
  #varias (m) inicializaciones en paralelo y escoger la mejor
  
  alpha_krr <- solve(Matrix(Kxs + diag(lambda, N))) %*% y
  
  ALPHA0 <- matrix(rnorm(num_init*N), N, num_init) 
  # que tengan misma norma 2 que la de alpha_krr las m inicializaciones
  ALPHA0 <- apply(ALPHA0, 2, function(col){
    res <- col/as.numeric(sqrt(t(col)%*%col))*as.numeric(sqrt(t(alpha_krr)%*%alpha_krr))
  })
  ALPHA0 <- as.data.frame(ALPHA0)
  #apply(ALPHA0, 2, function(col) as.numeric(sqrt(t(col)%*%col)))
  #sqrt(t(alpha_krr)%*%alpha_krr)
  
  fxs <- mcmapply(FUN=function(alpha0){
    res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1, max_iterations=max_iterations1)
    return(res$value)
  }, alpha0=as.data.frame(ALPHA0), mc.cores=1)
  
  # escoger la mejor inicialización
  indx <- which.min(fxs)
  alpha0 <- ALPHA0[,indx]
  
  # optimizar con esa inicialización pero mas iteraciones
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1, max_iterations=max_iterations2)
  alpha <- res$par
  
  learner$learnParams$alpha <- alpha
  learner$learnParams$avgy <- mean(y)
  
  #print("exits hsic learn function")
  return(learner)
}
predict.hsic <- function (learner, data) {
  kernelXs <- makeKernel(learner, "kernelXs")
  kxs <- kernelMatrix(kernelXs, data$x, learner$hyperParams$trainData$x)
  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}



learn.qkric <- function (data, params) {
  stopifnot(isRegression(data))
  kernelXs = do.call(params$kernelXs, params$psKernXs)
  kernelXb = do.call(params$kernelXb, params$psKernXb)
  return(.qkric(data$x, kernelXs, kernelXb, data$y, getN(data) * params$lambda, params$mu))
}
predict.qkric <- function (model, newData) {
  stopifnot(isRegression(newData))
  return(as.matrix(.qkric.predict(newData$x, model)))
}	


.hsicRegLossGradGamma <- function(alpha, y, Kxs, Kxb, kernelRg){
  rs <- y-Kxs%*%alpha
  Krg <- kernelMatrix(kernelRg, rs)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  ones <- matrix(1,N,1)
  gamma <- kernelRg@kpar$sigma
  grad <- -gamma*t(ones)%*%(Kxbc*Krg)%*%ones  
  return(grad)
}

.qkric <- function(data, kernelXs, kernelXb, y, lambda, mu) {
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  #Kxsc <- H%*%Kxs%*%H
  I <- diag(N)
  
  alpha = solve(Matrix((mu*H+(1-mu)*Kxbc)%*%Kxs+lambda*I)) %*% (mu*H+(1-mu)*Kxbc) %*% y
  return(list(data=data, kernelXs=kernelXs, kernelXb=kernelXb, alpha=alpha, avg=mean(y)))
}
.qkric.predict <- function(newData, qkric) {  
  kxs = kernelMatrix(qkric$kernelXs, newData, qkric$data)
  #N = nrow(kxs)
  #n <- ncol(kxs)
  #H1 <- diag(N)-matrix(1/N,N,N)
  #H2 <- diag(n)-matrix(1/n,n,n)
  #kxsc <- H1 %*% kxs %*% H2
  pred <- kxs %*% qkric$alpha
  pred <- pred - mean(pred) + qkric$avg
  return(pred)
}


learn.gptk <- function(learner) return(learner)

predict.gptk <- function(learner, data){
  aux <- gpPosteriorMeanVar(learner$hyperParams$data$optimizable$model$val, X=data$x)
  pred <- list(x=data$x, gy=data$y, gyh=aux)
  return(pred)
}



#############################################################################################################*
# Loss functions for cross validation and evaluation of learners

sse <- function(learner, pred){
  #print("enters sse")
  res <- sum((pred$gy - pred$gyh)^2)
  #print("exits sse")
  return(res)
}
mse <- function(learner, pred){
  #print("enters mse")
  res <- mean((pred$gy - pred$gyh)^2)
  #print("exits mse")
  return(res)
}
rmse <- function(learner, pred){
  #print("enters rmse")
  
  res <- sqrt(mean((pred$gy - pred$gyh)^2))
  #print("exits rmse")
  return(res)
}

# only for kernel methods - evaluates regularizer t(alpha)%*%Kxs%*%alpha
regL <- function(learner, pred){
  
  alpha <- learner$learnParams$alpha
  kernelXs <- makeKernel(learner, "kernelXs")
  Kxs <-  kernelMatrix(kernelXs, pred$x)
  return(as.numeric(t(alpha)%*%Kxs%*%alpha))
}

qhsicLoss <- function(learner, pred){
  
  kernelXb <- makeKernel(learner, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, pred$x)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  
  res <- t(pred$gyh)%*%Kxbc%*%pred$gyh + t(pred$gy)%*%Kxbc%*%pred$gy - 2*t(pred$gyh)%*%Kxbc%*%pred$gyh
  return(res)
}
nqhsicLoss <- function(learner, pred){
  
  resids <- pred$gy-pred$gyh
  kernelXb <- makeKernel(learner, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, pred$x)
  N <- nrow(Kxb)
  kernelR <- do.call("vanilladot", list())
  Kr <- kernelMatrix(kernelR, matrix(resids,N,1))
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Krc <- H%*%Kr%*%H
  res <- sum(diag(Kxbc%*%Krc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krc%*%Krc))))
  return(res)
}
hsicLoss  <- function(learner, pred){
  
  #print("enters hsicLoss")
  
  kernelXb <- makeKernel(learner, "kernelXb")
  
  
  Kxb <-  kernelMatrix(kernelXb, pred$x)
  
  resids <- pred$gy-pred$gyh
  kernelRg <- makeKernel(learner, "kernelRg") 
  Krg <- kernelMatrix(kernelRg, resids)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- Kxb%*%H
  Krgc <- Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc)) 
  #print("exits hsicLoss")
  return(res)
}
nhsicLoss  <- function(model, test, learner, param){
  
  #print("enters hsicLoss")
  
  kernelXb <- makeKernel(learner, "kernelXb")
  
  
  Kxb <-  kernelMatrix(kernelXb, pred$x)
  
  resids <- pred$gy-pred$gyh
  kernelRg <- makeKernel(learner, "kernelRg") 
  Krg <- kernelMatrix(kernelRg, resids)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Krgc <- H%*%Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krgc%*%Krgc))))
  #print("exits hsicLoss")
  return(res)
}
qkricLoss <- function(learner, pred){
  
  
  kernelXb <- makeKernel(learner, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, pred$x)
  
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  
  
  Kxbc <- H%*%Kxb%*%H
  
  mseLoss <- sum((pred$gy - test$gyh)^2)
  qhsicLoss <- (t(pred$gyh)%*%Kxbc%*%pred$gyh + t(pred$gy)%*%Kxbc%*%pred$gy - 2*t(pred$gyh)%*%Kxbc%*%pred$gyh)
  hyperParams <- learner$hyperParams$getHyperParams(learner) 
  mu <- hyperParams$mu
  res <- mu*mseLoss + (1-mu)*qhsicLoss
  
  return(res)
}

# Fair Learning / Consisent Regression Regularizer
nhsicReg  <- function(learner, pred){
  
  #print("enters nhsicReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  hyperParams <- learner$hyperParams$getHyperparams(learner)
  indxSens <- hyperParams$indxSens
  
  x <- matrix(pred$x, N, p)
  
  # we want the "real" hsic even though the learner may approximate with linear kernel for Xs and/or Yhs
  # this is a bit hacky but beta and kappa parameters must be in learner even if corresponding kernels are linear
  learnerAux <- learner
  learnerAux$hyperParams$non_data$kernelXb <- "rbfdot"
  learnerAux$hyperParams$non_data$kernelYhk <- "rbfdot"
  
  kernelXb <- makeKernel(learnerAux, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, x[,indxSens])
  
  
  kernelYhk <- makeKernel(learnerAux, "kernelYhk")
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhkc <- H%*%Kyhk%*%H
  res <- sum(diag(Kxbc%*%Kyhkc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Kyhkc%*%Kyhkc))))
  #print("exits nhsicReg")
  return(res)
}
hsicReg  <- function(learner, pred){
  
  #print("enters hsicReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  hyperParams <- learner$hyperParams$getHyperparams(learner)
  
  indxSens <- hyperParams$indxSens
  
  x <- matrix(test$x, N, p)
  
  
  # we want the "real" hsic even though the learner may approximate with linear kernel for Xs and/or Yhs
  # this is a bit hacky but beta and kappa parameters must be in learner even if corresponding kernels are linear
  learnerAux <- learner
  learnerAux$hyperParams$non_data$kernelXb <- "rbfdot"
  learnerAux$hyperParams$non_data$kernelYhk <- "rbfdot"
  
  kernelXb <- makeKernel(learnerAux, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, x[,indxSens])
  
  
  kernelYhk <- makeKernel(learnerAux, "kernelYhk")
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhkc <- H%*%Kyhk%*%H
  res <- sum(diag(Kxbc%*%Kyhkc))
  #print("exits hsicReg")
  return(res)
}
hsicYhReg  <- function(learner, pred){
  
  #print("enters hsicYhReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  
  
  x <- matrix(pred$x, N, p)
  
  learnerAux <- learner
  learnerAux$hyperParams$non_data$kernelYhk <- "rbfdot"
  kernelYhk <- makeKernel(learnerAux, "kernelYhk")
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  
  
  N <- nrow(Kyhk)
  H <- diag(N)-matrix(1/N,N,N)
  Kyhkc <- H%*%Kyhk%*%H
  res <- sqrt(sum(diag(Kyhkc%*%Kyhkc)))
  #print("exits hsicYhReg")
  return(res)
}

# actual regularizer used
nhsicRegA  <- function(learner, pred){
  
  #print("enters nhsicReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  hyperParams <- learner$hyperParams$getHyperparams(learner)
  indxSens <- hyperParams$indxSens
  
  
  x <- matrix(pred$x, N, p)
  
  
  
  kernelXb <- makeKernel(learner, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, as.matrix(x[,indxSens]))
  
  kernelYhk <- makeKernel(learner, "kernelYhk")
  
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhkc <- H%*%Kyhk%*%H
  res <- sum(diag(Kxbc%*%Kyhkc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Kyhkc%*%Kyhkc))))
  #print("exits nhsicReg")
  return(res)
}
hsicRegA  <- function(learner, pred){
  
  #print("enters hsicReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  hyperParams <- learner$hyperParams$getHyperparams(learner)
  indxSens <- hyperParams$indxSens
  
  x <- matrix(test$x, N, p)
  
  
  kernelXb <- makeKernel(learner, "kernelXb")
  Kxb <-  kernelMatrix(kernelXb, as.matrix(x[,indxSens]))
  
  kernelYhk <- makeKernel(learner, "kernelYhk")
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhkc <- H%*%Kyhk%*%H
  res <- sum(diag(Kxbc%*%Kyhkc))
  #print("exits hsicReg")
  return(res)
}
hsicYhRegA  <- function(learner, pred){
  
  #print("enters hsicYhReg")
  if(!is.null(dim(pred$x))){
    p <- dim(pred$x)[2]
    N <- dim(pred$x)[1]
  } else{
    N <- length(pred$x)
    p <- 1
  }
  
  
  
  x <- matrix(pred$x, N, p)
  
  kernelYhk <- makeKernel(learner, "kernelYhk")
  Kyhk <- kernelMatrix(kernelYhk, pred$gyh)
  N <- nrow(Kyhg)
  H <- diag(N)-matrix(1/N,N,N)
  Kyhkc <- H%*%Kyhk%*%H
  res <- sqrt(sum(diag(Kyhkc%*%Kyhkc)))
  #print("exits hsicYhReg")
  return(res)
}


corre  <- function(learner, pred){
  
  #print("enters corre")
  
  res <- cor(pred$gyh, pred$gy)
  #print("exits corre")
  return(res)
}


######################################################################################################*
# Heuristics or initialization of models hyper parameters

getFixedParams <- function(learner, data, indxSens=NULL, indxPred=NULL, plot=FALSE, print=FALSE){
  #print("enters getFixedParams")
  
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    p <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    p <- 1
  }
  if(n > 500){
    data <- getSubset(data, 1:500) 
    n <- 500
  }
  
  if(is.null(indxSens)){
    indxSens <- 1:p
  }
  
  if(is.null(indxPred)){
    indxPred <- 1:p
  }
  
  x <- matrix(data$x, n, p)
  xx <- x[,indxPred, drop=F]
  xs <- x[,indxSens, drop=F]
  
  sigma0 <- 1/median(as.numeric(dist(xx)^2))
  beta0 <- 1/median(as.numeric(dist(xs)^2))
  
  #obtain residuals with a mars model
  model  <- earth(x=xx, y=data$y, nfold=5)
  yh0 <- predict(model, xx)
  kappa0 <- 1/median(as.numeric(dist(yh0)^2))
  res0 <- data$y-yh0
  gamma0 <- 1/median(as.numeric(dist(res0)^2))
  
  
  # fit sigma
  ord <- 10
  sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
  varsHsics.sig <- sapply(sigmas1, function(sd1){
    kernelXs <- do.call("rbfdot", list(sigma=sd1))
    Kxs <- kernelMatrix(kernelXs, xx)
    N <- nrow(Kxs)
    H <- diag(N)-matrix(1/N,N,N)
    Kxsc <- H%*%Kxs%*%H
    distsX <- (Kxsc)[lower.tri(Kxsc)]
    res <- var(distsX)
    #plot(distsX2, distsX, main=res, ylim=c(-0.6,0.6))
    #hist(distsX,100, main=res)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.sig)
  sigmaVar <- sigmas1[indxMaxVar]
  
  # obtain dynamic range of sigma
  spl <- spline(log(sigmas1,10), varsHsics.sig)
  splf <- splinefun(log(sigmas1,10), varsHsics.sig)
  dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
  tol <- 1e-2
  DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
  DR.sig <- range(log(DR.sig,10))
  
  # obtain sigma saturation point
  sigmas2 <- 10^seq(DR.sig[1],DR.sig[2],1)
  dhsics.sig <- sapply(sigmas2, function(sd1){
    kernelXs <- do.call("rbfdot", list(sigma=sd1))
    kernelY <- do.call("vanilladot", list())
    Kxs <- kernelMatrix(kernelXs, xx)
    Ky <- kernelMatrix(kernelY, matrix(data$y,length(data$y),1))
    
    Ks <- vector("list", 2)
    Ks[[1]] <- Kxs
    Ks[[2]] <- Ky
    Kxx <- vector("list", 2)
    Kxx[[1]] <- Kxs
    Kxx[[2]] <- Kxs
    Kyy <- vector("list", 2)
    Kyy[[1]] <- Ky
    Kyy[[2]] <- Ky
    
    dhsicXY <- dhsic(K=Ks)$dHSIC
    dhsicXX <- dhsic(K=Kxx)$dHSIC
    dhsicYY <- dhsic(K=Kyy)$dHSIC
    
    #print("*************************")
    #print(paste("sigma: ", sd1 ,sep=""))
    #print(paste("dhsicXY: ", dhsicXY, sep=""))
    #print(paste("dhsicXX: ", dhsicXX, sep=""))
    #print(paste("dhsicYY: ", dhsicYY, sep=""))
    
    res <- (dhsicXY)/(sqrt(dhsicXX)*sqrt(dhsicYY))
    return(res)
  })
  dhsics.sig_fac <- dhsics.sig/max(dhsics.sig[which(dhsics.sig<Inf)])
  
  #indxSat <- which(dhsics.sig_fac > 0.998)
  #aux <- which.max(sigmas2[indxSat])
  #indxSat <- indxSat[aux]
  #sigmaSat <- sigmas2[indxSat]
  sigmaSat <- sigmas2[which.max(dhsics.sig)]
  
  sigmas <- c(med=sigma0, var=sigmaVar, sat=sigmaSat)
  
  # fit beta	
  betas1 <- (10^seq(-ord,ord,1)) # *beta0
  varsHsics.bet <- sapply(betas1, function(sd1){
    kernelXb <- do.call("rbfdot", list(sigma=sd1))
    Kxb <- kernelMatrix(kernelXb, xs)
    N <- nrow(Kxb)
    H <- diag(N)-matrix(1/N,N,N)
    Kxbc <- H%*%Kxb%*%H
    distsX <- (Kxbc)[lower.tri(Kxbc)]
    res <- var(distsX)
    #plot(distsX2, distsX, main=res, ylim=c(-0.6,0.6))
    #hist(distsX,100, main=res)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.bet)
  betaVar <- betas1[indxMaxVar]
  
  betas <- c(med=beta0, var=betaVar)
  
  #fit kappa
  
  kappas1 <- (10^seq(-ord,ord,1)) # *gamma0
  varsHsics.kap <- sapply(kappas1, function(sd1){
    kernelYhk <- do.call("rbfdot", list(sigma=sd1))
    Kyhk <- kernelMatrix(kernelYhk, yh0)
    N <- nrow(Kyhk)
    H <- diag(N)-matrix(1/N,N,N)
    Kyhkc <- H%*%Kyhk%*%H
    distsYh <- (Kyhkc)[lower.tri(Kyhkc)]
    res <- var(distsYh)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.kap)
  kappaVar <- kappas1[indxMaxVar]
  
  kappas <- c(med=kappa0, var=kappaVar)
  
  #fit gamma
  
  gammas1 <- (10^seq(-ord,ord,1)) # *gammar
  varsHsics.gam <- sapply(gammas1, function(sd1){
    kernelRg <- do.call("rbfdot", list(sigma=sd1))
    Krg <- kernelMatrix(kernelRg, res0)
    N <- nrow(Krg)
    H <- diag(N)-matrix(1/N,N,N)
    Krgc <- H%*%Krg%*%H
    distsR <- (Krgc)[lower.tri(Krgc)]
    res <- var(distsR)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.gam)
  gammaVar <- gammas1[indxMaxVar]
  
  gammas <- c(med=gamma0, var=gammaVar)
  
  if(print){
    print(paste("sigmaVar: ", sigmaVar, sep=""))
    print(paste("sigmaSat: ", sigmaSat, sep=""))
    print(paste("simgaMed: ", sigma0, sep=""))
    print(paste("betaVar: ", betaVar, sep=""))
    print(paste("betaMed: ", beta0, sep=""))
    print(paste("kappaVar: ", kappaVar, sep=""))
    print(paste("kappaMed: ", kappa0, sep=""))
    print(paste("gammaVar: ", gammaVar, sep=""))
    print(paste("gammaMed: ", gamma0, sep=""))
  }
  
  if(plot){
    
    # sigma
    
    # beta
    dhsics.bet <- sapply(betas1, function(sd1){
      kernelXb <- do.call("rbfdot", list(sigma=sd1))
      kernelYh <- do.call("vanilladot", list())
      Kxb <- kernelMatrix(kernelXb, xs)
      Kyh <- kernelMatrix(kernelYh, matrix(yh0,length(data$y),1))
      Ks <- vector("list", 2)
      Ks[[1]] <- Kxb 
      Ks[[2]] <- Kyh
      Kxx <- vector("list", 2)
      Kxx[[1]] <- Kxb 
      Kxx[[2]] <- Kxb 
      Kyhh <- vector("list", 2)
      Kyhh[[1]] <- Kyh
      Kyhh[[2]] <- Kyh
      
      dhsicXYh <- dhsic(K=Ks)$dHSIC
      dhsicXX <- dhsic(K=Kxx)$dHSIC
      dhsicYhYh <- dhsic(K=Kyhh)$dHSIC
      
      #print("*************************")
      #print(paste("sigma: ", sd1 ,sep=""))
      #print(paste("dhsicXYh: ", dhsicXYh, sep=""))
      #print(paste("dhsicXX: ", dhsicXX, sep=""))
      #print(paste("dhsicYhYh: ", dhsicYhYh, sep=""))
      
      res <- (dhsicXYh)/(sqrt(dhsicXX)*sqrt(dhsicYhYh))
      return(res)
    })
    dhsics.bet_fac <- dhsics.bet/max(dhsics.bet[which(dhsics.bet<Inf)])
    indxSat <- which(dhsics.bet_fac > 0.95)
    aux <- which.min(dhsics.bet_fac[indxSat])
    indxSat <- indxSat[aux]
    betaSat <- betas1[indxSat]
    
    # kappa
    dhsics.kap <- sapply(kappas1, function(sd1){
      kernelX <- do.call("vanilladot", list())
      kernelYh <- do.call("rbfdot", list(sigma=sd1))
      Kx <- kernelMatrix(kernelX, matrix(xs, length(data$y),1))
      Kyhk <- kernelMatrix(kernelYh, matrix(yh0,length(data$y),1))
      Ks <- vector("list", 2)
      Ks[[1]] <- Kx 
      Ks[[2]] <- Kyhk
      Kxx <- vector("list", 2)
      Kxx[[1]] <- Kx 
      Kxx[[2]] <- Kx 
      Kyhh <- vector("list", 2)
      Kyhh[[1]] <- Kyhk
      Kyhh[[2]] <- Kyhk
      
      dhsicXYh <- dhsic(K=Ks)$dHSIC
      dhsicXX <- dhsic(K=Kxx)$dHSIC
      dhsicYhYh <- dhsic(K=Kyhh)$dHSIC
      
      #print("*************************")
      #print(paste("sigma: ", sd1 ,sep=""))
      #print(paste("dhsicXYh: ", dhsicXYh, sep=""))
      #print(paste("dhsicXX: ", dhsicXX, sep=""))
      #print(paste("dhsicYhYh: ", dhsicYhYh, sep=""))
      
      res <- (dhsicXYh)/(sqrt(dhsicXX)*sqrt(dhsicYhYh))
      return(res)
    })
    dhsics.kap_fac <- dhsics.kap/max(dhsics.kap[which(dhsics.kap<Inf)])
    indxSat <- which(dhsics.kap_fac > 0.95)
    aux <- which.min(dhsics.kap_fac[indxSat])
    indxSat <- indxSat[aux]
    kappaSat <- kappas1[indxSat]
    
    # gamma
    dhsics.gam <- sapply(gammas1, function(sd1){
      kernelX <- do.call("vanilladot", list())
      kernelRg <- do.call("rbfdot", list(sigma=sd1))
      Kx <- kernelMatrix(kernelX, matrix(xx, length(data$y),1))
      Krg <- kernelMatrix(kernelRg, matrix(res0,length(data$y),1))
      Ks <- vector("list", 2)
      Ks[[1]] <- Kx 
      Ks[[2]] <- Krg
      Kxx <- vector("list", 2)
      Kxx[[1]] <- Kx 
      Kxx[[2]] <- Kx 
      Krr <- vector("list", 2)
      Krr[[1]] <- Krg
      Krr[[2]] <- Krg
      
      dhsicXR <- dhsic(K=Ks)$dHSIC
      dhsicXX <- dhsic(K=Kxx)$dHSIC
      dhsicRR <- dhsic(K=Krr)$dHSIC
      
      #print("*************************")
      #print(paste("sigma: ", sd1 ,sep=""))
      #print(paste("dhsicXR: ", dhsicXR, sep=""))
      #print(paste("dhsicXX: ", dhsicXX, sep=""))
      #print(paste("dhsicRR: ", dhsicRR, sep=""))
      
      res <- (dhsicXR)/(sqrt(dhsicXX)*sqrt(dhsicRR))
      return(res)
    })
    dhsics.gam_fac <- dhsics.gam/max(dhsics.gam[which(dhsics.gam<Inf)])
    indxSat <- which(dhsics.gam_fac > 0.95)
    aux <- which.min(dhsics.gam_fac[indxSat])
    indxSat <- indxSat[aux]
    gammaSat <- gammas1[indxSat]
    
    
    
    print(paste("betaSat: ", betaSat, sep=""))
    print(paste("kappaSat: ", kappaSat, sep=""))
    print(paste("gammaSat: ", gammaSat, sep=""))
    
    par(mfrow=c(2,2))
    # sigma
    plot(log(sigmas2,10), dhsics.sig, main="fit sigma", xlab="log(sigma)", ylab="hsic/var", xlim=range(log(sigmas1,10)),ylim=range(dhsics.sig[which(dhsics.sig<Inf)], varsHsics.sig))
    lines(log(sigmas1,10), varsHsics.sig, col="red", type="p")
    abline(v=DR.sig, col="purple")
    abline(v=log(sigma0,10), col="red")
    abline(v=log(sigmaSat,10), col="green")
    abline(v=log(sigmaVar,10), col="blue")
    # beta
    plot(log(betas1,10), dhsics.bet, main="fit beta", xlab="log(betas)", ylab="hsic/var", ylim=range(dhsics.bet[which(dhsics.bet<Inf & dhsics.bet>-Inf)], varsHsics.bet))
    lines(log(betas1,10), varsHsics.bet, col="red", type="p")
    abline(v=log(beta0,10), col="red")
    abline(v=log(betaSat,10), col="green")
    abline(v=log(betaVar,10), col="blue")
    # gamma
    plot(log(gammas1,10), dhsics.gam, type="p", main="gammas fit", xlab="log(gamma)",ylab="hsic/var", ylim=range(dhsics.gam[which(dhsics.gam<Inf & dhsics.gam>-Inf)], varsHsics.gam))
    lines(log(gammas1,10), varsHsics.gam, col="red", type="p")
    abline(v=log(gamma0,10), col="red")
    abline(v=log(gammaSat,10), col="green")
    abline(v=log(gammaVar,10), col="blue")
    
    # beta-gamma
    par(mfrow=c(1,1))
    dhsics.gamBet <- sapply(betas1, function(sd1) sapply(gammas1, function(sd2){
      kernelXb <- do.call("rbfdot", list(sigma=sd1))
      kernelYhg <- do.call("rbfdot", list(sigma=sd2))
      Kxb <- kernelMatrix(kernelXb, xs)
      Kyhg <- kernelMatrix(kernelYhg, matrix(yh0, length(data$y),1))
      
      Ks <- vector("list", 2)
      Ks[[1]] <- Kxb 
      Ks[[2]] <- Kyhg
      Kxx <- vector("list", 2)
      Kxx[[1]] <- Kxb 
      Kxx[[2]] <- Kxb 
      Kyhh <- vector("list", 2)
      Kyhh[[1]] <- Kyhg
      Kyhh[[2]] <- Kyhg
      
      dhsicXYh <- dhsic(K=Ks)$dHSIC
      dhsicXX <- dhsic(K=Kxx)$dHSIC
      dhsicYhYh <- dhsic(K=Kyhh)$dHSIC
      
      #print("*************************")
      #print(paste("sigma: ", sd1 ,sep=""))
      #print(paste("dhsicXYh: ", dhsicXYh, sep=""))
      #print(paste("dhsicXX: ", dhsicXX, sep=""))
      #print(paste("dhsicYhYh: ", dhsicYhYh, sep=""))
      
      res <- (dhsicXYh)/(sqrt(dhsicXX)*sqrt(dhsicYhYh))
      return(res)
    }))
    dimnames(dhsics.gamBet) <- list(gamma=gammas1, beta=betas1)
    
    df <- melt(dhsics.gamBet)
    df <- df[which(df$value<Inf),]
    df <- df[which(df$value>-Inf),]
    v <- ggplot(df, aes(log(beta), log(gamma), z = value))
    v <- v + geom_raster(aes(fill = value)) + geom_contour(colour = "white", bins = 10)
    v <- v + geom_point(aes(x=log(betaVar), y=log(gammaVar)), colour="blue", size=2)
    v <- v + geom_point(aes(x=log(betaSat), y=log(gammaSat)), colour="green", size=2)
    print(v)
    
  }
  
  sigmaType  <- learner$hyperParams$data$non_optimizable$sigma$type
  betaType   <- learner$hyperParams$data$non_optimizable$beta$type
  gammaType  <- learner$hyperParams$data$non_optimizable$gamma$type
  kappaType  <- learner$hyperParams$data$non_optimizable$kappa$type
  
  
  #print("exits getFixedParams")
  return(list(sigma=sigmas[sigmaType], beta=betas[betaType], kappa=kappas[kappaType], gamma=gammas[gammaType], indxPred=list(indxPred), indxSens=list(indxSens)))
}

myGpCreate <- function(learner, trainData){
  options <- gpOptions(approx=learner$optimizeParams$approx)
  options$kern$comp <- learner$hyperParams$non_data$kernel  
  options$numActive <- min(learner$optimizeParams$numActive, nrow(trainData$x))
  
  if(learner$optimizeParams$fixInducing){
    options$fixInducing <- TRUE
    wishInducing <- makeGrid(trainData$x, options$numActive)
    
    indxs <- knnx.index(data=trainData$x, query=wishInducing, k=1)
    
    options$fixIndices <- as.numeric(indxs)
  }
  
  model <- gpCreate(q=dim(trainData$x)[2], d=1, X=trainData$x, y=as.matrix(trainData$y), options)
  return(list(modelInit=model))
}


######################################################################################################*
# optimization functions for learners' hyperparameters

CV.parallel <- function(learner, params, fac=1, verbose=TRUE) {
  stopifnot(class(learner) == "emley.learner" && class(params) == "CVST.params")
  
  #print("enters CV parallel")
  
  numFolds <- learner$optimizeParams$numFolds
  trainData <- learner$hyperParams$trainData
  
  
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), names(params))
  
  n <- getN(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- params[[6]]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- getSubset(trainData, -validationIndex)
      curTest <- getSubset(trainData, validationIndex)
      # either mean squared error or mean classification error
      
      
      learnerAux <- learner
      learnerAux$hyperParams$trainData <- curTrain
      
      nmsPars <- names(learnerAux$hyperParams$data$optimizable)
      for(pr in nmsPars){
        learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
      }
      
      #learner <- learnerAux
      learnerAux <- try(learnerAux$learn(learnerAux))
      if(class(learnerAux)=="try-error"){
        print(paste("singular matrix for parameters = ", paste(paste(names(p), p, collapse="-") ,collapse=", "), "for fold = ", f))
        res <- rep(NA, length(learner$optimizeParams$losses))
        names(res) <- names(learner$optimizeParams$losses)
        res <- cbind(res, res)
        return(res)
      }
      
      
      predTrain <- learnerAux$predict(learner=learnerAux, data=curTrain)
      predTest <- learnerAux$predict(learnerAux, data=curTest)
      
      
      
      lossesTrain <- sapply(names(learnerAux$optimizeParams$losses), function(func){
        do.call(func, list(learner=learnerAux, pred=predTrain))
      })
      
      lossesTest <- sapply(names(learnerAux$optimizeParams$losses), function(func){
        do.call(func, list(learner=learnerAux, pred=predTest))
      })
      
      
      
      res <- cbind(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, params, mc.cores=2, SIMPLIFY="array")
  
  
  
  losses <- aperm(losses, c(2,3,4,1))
  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=names(params), var=names(learner$optimizeParams$losses))
  
  #print("exits CV parallel")
  return(losses)
}


optHP.CV <- function(learner, plot=FALSE, fac=1){
  #print("enters optHP.CV")
  
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  numParams <- length(optimHyperParams)
  
  numEach <- lapply(optimHyperParams, function(el) length(el$seq))
  
  paramsList <- lapply(optimHyperParams, function(el){ 
    if(el$log){
      res <- 10^el$seq
    } else{
      res <- el$seq
    }
  })
  names(paramsList) <- names(optimHyperParams)
  params <- c(paramsList, lapply(learner$hyperParams$data$non_optimizable, function(el) el$val), learner$hyperParams$non_data)
  
  
  
  params <- do.call(constructParams, params)
  
  #params <- params2
  grid <- CV.parallel(learner, params, fac=fac)
  
  
  grid <- apply(grid, c(1,3,4), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  
  # reshape back into one dimension per hyperparameter
  
  
  
  dimnms <- dimnames(grid)
  
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], numEach , dim(grid)[3])
  
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  
  dimnames(grid) <- dimnms
  
  # test whether the dimensions are right (only when lambda and sigma are passed for CV)
  #iLambda <- 25
  #jSigma  <- 17
  #dimnames(grid2)$lambda[iLambda]
  #dimnames(grid2)$sigma[jSigma]
  #k <- (jSigma-1)*numEach$lambda + iLambda
  #dimnames(grid)$params[k]
  #grid[2,k,2]
  #grid2[2,iLambda,jSigma,2]
  
  # obtain best hyperparameters
  
  testGrid <- as.array(apply(grid, names(paramsList), function(el) el["test",1])) #grid["test",,,1]
  dimnames(testGrid) <- dimnms[which(names(dimnms) %in% names(paramsList))]
  
  minTestGrid <- min(testGrid, na.rm=T)
  #we add a tiny bit of noise to get exactly one minimum
  testGrid <- testGrid+ rnorm(length(testGrid),mean=0, sd=max(minTestGrid,1e-10)*1e-10)
  minTestGrid <- min(testGrid, na.rm=T)
  optMinBool <- testGrid==minTestGrid
  
  opt <- list()
  
  #print("optimal values")
  for(nm in names(optimHyperParams)){
    opt[[nm]] <- paramsList[[nm]][which(apply(optMinBool, nm, any, na.rm=T))]
    #print(opt[[nm]])
  }
  
  
  
  
  # check (only when lanbda and sigma are passed for CV)
  #grid["test",which.min(abs(as.numeric(dimnames(grid)$lambda)-opt$lambda)),which.min(abs(as.numeric(dimnames(grid)$sigma)-opt$sigma)),1]
  #min(grid["test",,,1])
  
  
  res <-  list(opt=opt, grid=grid)
  
  if(plot){
    if(numParams<=2){
      plot.optHP.CV(res)
    } else{
      print("cannot plot for more than 2 cross validated hyperparams at the moment. ")
    }
    
  } 
  
  #print("exits optHP.CV")
  return(res)
}

plot.optHP.CV <- function(opt){
  
  params <- names(opt$opt)
  numParams <- length(params)
  
  if(numParams==1){
    
    numLosses <- dim(opt$grid)[3]
    
    numRows <- ceiling(sqrt(numLosses))
    numCols <- ceiling(numLosses/numRows)
    
    valsPar <- as.numeric(dimnames(opt$grid)[[2]])
    
    par(mfrow=c(numRows, numCols))
    for(loss in dimnames(opt$grid)$var){
      
      plot(log(valsPar, 10), opt$grid["train",,loss], ylim=range(opt$grid[,,loss],na.rm=T), type="b", xlab=names(opt$opt), ylab="", main=loss)
      lines(log(valsPar,10), opt$grid["test",,loss], type="b", col="red")
      lines(rep(log(opt$opt[[1]],10), 2), opt$grid[,which.min(abs(valsPar - opt$opt[[1]])),loss], col="red" , cex=1.5, type="p")
    }
    par(mfrow=c(1,1))
  }
  
  if(numParams==2){
    testGrid <- opt$grid["test",,,1]
    df <- melt(testGrid)
    
    indxMax <- which.max(df$value)
    df$value[which(df$value > quantile(df$value, 0.95))] <- quantile(df$value, 0.95) 
    
    v <- ggplot()
    v <- v + geom_raster(aes_string(paste("log(", params[1],",10)"), paste("log(",params[2],",10)"), fill = "value"), data=df) 
    v <- v + geom_contour(aes_string(paste("log(",params[1],", 10)"), paste("log(",params[2],", 10)"), z="value"),colour = "white", bins = 5, data=df)
    v <- v + geom_point(aes(x=log(opt$opt[[1]], 10), y=log(opt$opt[[2]], 10)), colour="green", size=2)
    v <- v + geom_point(aes_string(x=paste("log(",params[1],", 10)"), y=paste("log(",params[2],", 10)")), colour="red", size=2, data=df[indxMax,])
    print(v)
    
  }
}


myGpOptimise <- function(learner, plot){
  
  model <- gpOptimise(model=learner$hyperParams$data$non_optimizable$modelInit$val, display=F, iters=learner$optimizeParams$iters, gradcheck= FALSE)
  
  #if(plot) gpPlot(model=model, Xstar=learner$hyperParams$trainData$x)
  
  return(list(opt=list(model=model), grid=NULL)) 
}

#needed to change something from one or more of following gptk package functions to get gpOptimise to work
gpCreate <- function (q, d, X, y, options){
  if (dim(X)[2] != q) 
    stop("Input matrix X does not have dimension ", q)
  if (dim(y)[2] != d) 
    stop("Input matrix y does not have dimension ", d)
  if (any(is.nan(y)) && !options$isMissingData) 
    stop("NaN values in y, but no missing data declared.")
  if (options$isMissingData && options$isSpherical) 
    stop("If there is missing data, spherical flag cannot be set.")
  y = as.matrix(y)
  X = as.matrix(X)
  model <- list(type = "gp", y = y, X = X, approx = options$approx, 
                beta = options$beta, learnScales = options$learnScales, 
                scaleTransform = optimiDefaultConstraint("positive"), 
                optimiseBeta = options$optimiseBeta, betaTransform = optimiDefaultConstraint("positive"), 
                q = dim(X)[2], d = dim(y)[2], N = dim(y)[1])
  model$optimiser = options$optimiser
  model$isSpherical = options$isSpherical
  model$isMissingData = options$isMissingData
  model$scale = matrix(1, 1, model$d)
  if (!model$isMissingData) {
    model$bias = colMeans(y)
  }
  else {
    for (i in 1:model$d) {
      model$indexPresent[[i]] = which(!is.nan(y[, i]))
      if (length(model$indexPresent[[i]]) == 0) {
        model$bias[i] = 0
      }
      else {
        model$bias[i] = mean(model$y[model$indexPresent[[i]], 
                                     i])
      }
    }
  }
  if (("scale2var1" %in% names(options)) && (options$scale2var1)) {
    model$scale = sd(model$y)
    model$scale[which(model$scale == 0)] = 1
    if (model$learnScales) 
      warning("Both learn scales and scale2var1 set for GP")
    if ("scaleVal" %in% names(options)) 
      warning("Both scale2var1 and scaleVal set for GP")
  }
  if ("scaleVal" %in% names(options)) 
    model$scale = kronecker(matrix(1, 1, model$d), options$scaleVal)
  model$m = gpComputeM(model)
  model$computeS = FALSE
  if (options$computeS) {
    model$computeS = TRUE
    model$S = model$m %*% t(model$m)
    if (model$approx != "ftc") 
      stop("If compute S is set, approximation type must be 'ftc'.")
  }
  if (is.list(options$kern) && ("nParams" %in% options$kern)) 
    model$kern = options$kern
  else model$kern = kernCreate(model$X, options$kern)
  if (options$approx == "ftc") {
    model$k = 0
    model$X_u = list()
    if (model$optimiseBeta && length(options$beta) == 0) 
      stop("options.beta cannot be empty if it is being optimised.")
  }
  ##############################################################################################################################
  # Approximating Gps such as FITC    
  ##############################################################################################################################    
  else if (options$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    
    model$k = options$numActive
    model$fixInducing = options$fixInducing
    ###############################################
    # In case you want to specify inducing points
    ###############################################
    if (options$fixInducing) {
      if (length(options$fixIndices) != options$numActive) {
        stop(paste("Length of indices for fixed inducing variables must ", 
                   "match number of inducing variables"))
      }
      model$X_u = model$X[options$fixIndices, ,drop=FALSE] #added drop=FALSE
      model$inducingIndices = options$fixIndices
    }
    else {
      ind = sample(1:model$N, size = model$N)
      ind = ind[1:model$k]
      model$X_u = model$X[ind, , drop = FALSE]
    }
  } # end of if options$approx %in% c("dtc", "dtcvar", "fitc", "pitc")
  if (model$k > model$N) 
    stop("Number of active points cannot be greater than number of data.")
  if (model$approx == "pitc") {
    numBlocks = ceiling(model$N/model$k)
    numPerBlock = ceiling(model$N/numBlocks)
    startVal = 1
    endVal = model$k
    model$blockEnd = matrix(0, 1, numBlocks)
    for (i in 1:numBlocks) {
      model$blockEnd[i] = endVal
      endVal = numPerBlock + endVal
      if (endVal > model$N) 
        endVal = model$N
    }
  }
  initParams = gpExtractParam(model)
  model = gpExpandParam(model, initParams)
  return(model)
}

gpOptimise <- function(model, display = TRUE, iters = 2000, gradcheck = FALSE){
  params = gpExtractParam(model)
  options = optimiDefaultOptions()
  options$display = FALSE
  if (display) {
    options$display = TRUE
    if ((length(params) <= 100) && gradcheck) 
      options$gradcheck = TRUE
  }
  options$maxit = iters
  if ("optimiser" %in% names(model)){ 
    optim = get(paste(model$optimiser, "optim", sep = ""), mode = "function")
  } else {
    optim = get("CGoptim", mode = "function")
  }
  fn = get("gpObjective", mode = "function")
  grad = get("gpGradient", mode = "function")
  #print("a")
  params = optim(params, fn, grad, options, model)
  #print("b")
  model = gpExpandParam(model, params)
  return(model)
}

gpObjective <- function (params, model) {
  #print("entra gpObjective")
  model = gpExpandParam(model, params)
  f = -gpLogLikelihood(model)
  #print("sale gpObjective")
  return(f)
}

gpGradient <- function (params, model){
  #print("entra gpGradient")
  model = gpExpandParam(model, params)
  g = -gpLogLikeGradients(model)
  #print("sale gpGradient")
  return(g)
}

gpLogLikelihood <- function (model){
  if (model$approx == "ftc") {
    if ("S" %in% names(model)) {
      ll = -0.5 * (model$d * model$logDetK_uu + sum(model$invK_uu * 
                                                      model$S))
      return(ll)
    }
    ll = 0
    for (i in 1:dim(model$m)[2]) {
      if ((!"isSpherical" %in% names(model)) || model$isSpherical) 
        ll = ll - 0.5 * model$logDetK_uu - 0.5 * t(model$m[, 
                                                           i, drop = FALSE]) %*% model$invK_uu %*% model$m[, 
                                                                                                           i, drop = FALSE]
      else {
        if (model$isMissingData) 
          m = model$m[model$indexPresent[[i]], i]
        else m = model$m[, i]
        ll = ll - 0.5 * model$logDetK_uu[i] - 0.5 * t(m) %*% 
          model$invK_uu[[i]] %*% m
      }
    }
  }
  else if (model$approx %in% c("dtc", "dtcvar")) {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      E = model$K_uf %*% model$m
      EET = E %*% t(E)
      if (length(model$beta) == 1) {
        ll = -0.5 * (model$d * (-(model$N - model$k) * 
                                  log(model$beta) - model$logDetK_uu + model$logDetA) - 
                       (sum(model$Ainv * EET) - sum(model$m * model$m)) * 
                       model$beta)
        if (model$approx == "dtcvar") 
          ll = ll - model$d * 0.5 * sum(model$diagD)
      }
      else stop("Not implemented variable length beta yet.")
    }
    else {
      ll = 0
      for (i in 1:model$d) {
        ind = gpDataIndices(model, i)
        e = model$K_uf[, ind, drop = FALSE] %*% model$m[ind, 
                                                        i, drop = FALSE]
        if (length(model$beta) == 1) {
          ll = ll - 0.5 * ((-(model$N - model$k) * log(model$beta) - 
                              model$logDetK_uu + model$logDetA[i]) - (t(e) %*% 
                                                                        model$Ainv[[i]] %*% e - t(model$m[ind, i, 
                                                                                                          drop = FALSE]) %*% model$m[ind, i, drop = FALSE]) * 
                             model$beta)
          if (is.nan(ll)) 
            stop("Log likelihood is NaN")
          if (model$approx == "dtcvar") 
            stop("Not implemented dtcvar for non-spherical yet.")
        }
        else stop("Not implemented variable length beta yet.")
      }
    }
  }
  else if (model$approx == "fitc") {
    #print("entra aqui")
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      if (length(model$beta) == 1) {
        if (FALSE) {
          Dinvm = model$Dinv %*% model$m
          K_ufDinvm = model$K_uf %*% Dinvm
          ll = -0.5 * (model$d * (sum(log(model$diagD)) - 
                                    (model$N - model$k) * log(model$beta) + model$detDiff) + 
                         (sum(Dinvm * model$m) - sum((model$Ainv %*% 
                                                        K_ufDinvm) * K_ufDinvm)) * model$beta)
          ll = ll - 0.5 * model$N * model$d * log(2 * 
                                                    pi)
        }
        else {
          ll = -model$d * (sum(log(diag(model$Lm))) + 
                             0.5 * (-(model$N - model$k) * log(model$beta) + 
                                      (model$N * log(2 * pi) + sum(log(model$diagD)))))
          for (i in 1:model$d) ll = ll - 0.5 * model$beta * 
              (t(model$scaledM[, i, drop = FALSE]) %*% 
                 model$scaledM[, i, drop = FALSE] - t(model$bet[, 
                                                                i, drop = FALSE]) %*% model$bet[, i, drop = FALSE])
        }
      }
      else stop("Variable length Beta not implemented yet.")
    }
    else {
      if (length(model$beta) == 1) {
        if (FALSE) {
          ll = 0
          for (i in 1:model$d) {
            ind = gpDataIndices(model, i)
            Dinvm = model$Dinv[[i]] %*% model$m[ind, 
                                                i, drop = FALSE]
            K_ufDinvm = model$K_uf[, ind, drop = FALSE] %*% 
              Dinvm
            ll = ll - 0.5 * (sum(log(model$diagD[[i]])) - 
                               (length(ind) - model$k) * log(model$beta) + 
                               model$detDiff[i] + (sum(Dinvm * model$m[ind, 
                                                                       i, drop = FALSE]) - sum((model$Ainv[[i]] %*% 
                                                                                                  K_ufDinvm) * K_ufDinvm)) * model$beta + 
                               length(ind) * log(2 * pi))
          }
        }
        else {
          ll = 0
          for (i in 1:model$d) {
            ind = gpDataIndices(model, i)
            ll = ll - (sum(log(diag(model$Lm[[i]]))) + 
                         0.5 * (-(length(ind) - model$k) * log(model$beta) + 
                                  (length(ind) * log(2 * pi) + sum(log(model$diagD[[i]])))))
            ll = ll - 0.5 * model$beta * (t(model$scaledM[[i]]) %*% 
                                            model$scaledM[[i]] - t(model$bet[[i]]) %*% 
                                            model$bet[[i]])
          }
        }
      }
      else stop("Variable length Beta not implemented yet.")
    }
  }
  else if (model$approx == "pitc") {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      if (length(model$beta) == 1) {
        ll = model$d * (model$logDetA - model$logDetK_uu + 
                          model$k * log(model$beta))
        K_ufDinvm = matrix(0, model$k, model$d)
        Dinvm = list()
        for (i in 1:length(model$blockEnd)) {
          ind = gpBlockIndices(model, i)
          Dinvm[[i]] = model$Dinv[[i]] %*% model$m[ind, 
                                                   , drop = FALSE]
          K_ufDinvm = K_ufDinvm + model$K_uf[, ind, drop = FALSE] %*% 
            Dinvm[[i]]
        }
        ll = ll - model$beta * sum((model$Ainv %*% K_ufDinvm) * 
                                     K_ufDinvm)
        for (i in 1:length(model$blockEnd)) {
          ind = gpBlockIndices(model, i)
          ll = ll + model$d * (model$logDetD[i] - length(ind) * 
                                 log(model$beta))
          +model$beta * sum(Dinvm[[i]] * model$m[ind, 
                                                 , drop = FALSE])
        }
        ll = -0.5 * ll
        ll = ll - 0.5 * model$N * model$d * log(2 * pi)
      }
      else stop("Variable Length Beta not implemented yet.")
    }
    else {
      if (length(model$beta) == 1) {
        ll = 0
        Dinvm = matrix(0, model$blockEnd, model$d)
        Dinvm = lapply(split(Dinvm, row(Dinvm)), split, 
                       1:model$d)
        for (j in 1:model$d) {
          ll = ll + model$logDetA[j] - model$logDetK_uu + 
            model$k * log(model$beta)
          K_ufDinvm = matrix(0, model$k, 1)
          for (i in 1:length(model$blockEnd)) {
            ind = gpDataIndices(model, j, i)
            Dinvm[[i]][[j]] = model$Dinv[[i]][[j]] %*% 
              model$m[ind, j, drop = FALSE]
            K_ufDinvm = K_ufDinvm + model$K_uf[, ind] %*% 
              Dinvm[[i]][[j]]
          }
          ll = ll - model$beta * sum((model$Ainv[[i]] %*% 
                                        K_ufDinvm) * K_ufDinvm)
          for (i in 1:length(model$blockEnd)) {
            ind = gpDataIndices(model, j, i)
            ll = ll + model$logDetD[i, j] - length(ind) * 
              log(model$beta)
            +model$beta * sum(Dinvm[[i]][[j]] * model$m[ind, 
                                                        j, drop = FALSE])
            ll = ll + length(ind) * log(2 * pi)
          }
        }
        ll = -0.5 * ll
      }
      else stop("Variable Length Beta not implemented yet.")
    }
  }
  if (model$learnScales) 
    ll = ll - sum(log(model$scale))
  ll = ll - model$d * model$N/2 * log(2 * pi)
  return(ll)
}

gpLogLikeGradients <- function (model, X = model$X, M, X_u, gX_u.return = FALSE, gX.return = FALSE, g_beta.return = FALSE){
  if (missing(X_u)) {
    X_u = list()
    if ("X_u" %in% names(model)) 
      X_u = model$X_u
    if (missing(M) && (!"S" %in% names(model))) 
      M = model$m
  }
  gX_u = list()
  gX = list()
  g_scaleBias = gpScaleBiasGradient(model)
  g_meanFunc = list()
  if ("meanFunction" %in% names(model) && length(model$meanFunction) > 0) 
    g_meanFunc = gpMeanFunctionGradient(model)
  if (model$approx == "ftc") {
    if (gX_u.return && gX.return) {
      gKX = kernGradX(model$kern, X, X)
      gKX = gKX * 2
      dgKX = kernDiagGradX(model$kern, X)
      for (i in 1:model$N) gKX[i, , i] = dgKX[i, ]
      gX = matrix(0, model$N, model$q)
    }
    g_param = matrix(0, 1, model$kern$nParams)
    g_beta = list()
    if ("beta" %in% names(model)) 
      g_beta = 0
    if ("S" %in% names(model)) {
      gK = localSCovarianceGradients(model)
      if (gX_u.return && gX.return) {
        counter = 0
        for (i in 1:model$N) {
          counter = counter + 1
          for (j in 1:model$q) gX[i, j] = gX[i, j] + 
              t(gKX[, j, i, drop = FALSE]) %*% gK[, counter, 
                                                  drop = FALSE]
        }
      }
      g_param = g_param + kernGradient(model$kern, X, gK)
    }
    else {
      for (k in 1:model$d) {
        gK = localCovarianceGradients(model, M[, k], 
                                      k)
        if (gX_u.return && gX.return) {
          ind = gpDataIndices(model, k)
          counter = 0
          for (i in ind) {
            counter = counter + 1
            for (j in 1:model$q) gX[i, j] = gX[i, j] + 
                gKX[ind, j, i, drop = FALSE] %*% gK[, counter, 
                                                    drop = FALSE]
          }
        }
        if (model$isMissingData) {
          g_param = g_param
          +kernGradient(model$kern, X[model$indexPresent[[k]], 
                                      ], gK)
        }
        else g_param = g_param + kernGradient(model$kern, 
                                              X, gK)
      }
      if ("beta" %in% names(model) && model$optimiseBeta) {
        model$beta = as.matrix(model$beta)
        if (dim(model$beta)[1] == 1) 
          g_beta = g_beta + sum(diag(gK))
        else if (dim(model$beta)[2] == 1 && dim(model$beta)[1] == 
                 model$N) 
          g_beta = g_beta + diag(gK)
        else if (dim(model$beta)[2] == model$d && dim(model$beta)[1] == 
                 model$N) 
          g_beta[, k] = diag(gK)
        else stop("Unusual dimensions for model$beta.")
      }
    }
  }
  else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    gK = gpCovGrads(model, M)
    gK_uu = gK$gK_uu
    gK_uf = gK$gK_uf
    gK_star = gK$g_Lambda
    g_beta = gK$gBeta
    gParam_u = kernGradient(model$kern, X_u, gK_uu)
    gParam_uf = kernGradient(model$kern, X_u, X, gK_uf)
    g_param = gParam_u + gParam_uf
    gKX = kernGradX(model$kern, X_u, X_u)
    gKX = gKX * 2
    dgKX = kernDiagGradX(model$kern, X_u)
    for (i in 1:model$k) gKX[i, , i] = dgKX[i, ]
    if (!model$fixInducing || gX_u.return || gX.return || g_beta.return) {
      #print("entra a aqui 1")
      gX_u = matrix(0, model$k, model$q)
      #print("I")
      for (i in 1:model$k) {
        for (j in 1:model$q) gX_u[i, j] = t(gKX[, j, i]) %*% gK_uu[, i, drop = FALSE]
      }
      #print("II")
      gKX_uf = kernGradX(model$kern, X_u, X)
      #print("III")
      for (i in 1:model$k) {
        for (j in 1:model$q){
          #print(paste("(i,j)= ", i, j))
          #print(paste("dim(gX_u)", paste(dim(gX_u), collapse="-")))
          #print(paste("class(gX_u[i,j])", paste(class(gX_u[i,j]), collapse="-")))
          #print(paste("gX_u[i,j]", gX_u[i,j]))
          #print(paste("dim(gKX_uf)", paste(dim(gKX_uf), collapse="-")))
          #print(paste("dim(gK_uf)", paste(dim(gK_uf), collapse="-")))
          #print(paste("dim(gKX_uf[, j, i])", paste(dim(gKX_uf[,j,i]), collapse="-")))
          #print(paste("dim(gK_uf[i, , drop = FALSE])", paste(dim(gK_uf[i,,drop=FALSE]), collapse="-")))
          
          #print(paste("t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE])", t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE])))
          #print(paste("class(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE]))", class(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE]))))
          
          #print(gX_u[i, j])
          #print(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE]))
          #print(gX_u[i, j] + as.matrix(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE])))
          #print(gX_u[i, j] + as.numeric(as.matrix(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE]))))
          gX_u[i, j] = gX_u[i, j] + as.numeric(as.matrix(t(gKX_uf[, j, i]) %*% t(gK_uf[i, , drop = FALSE])))
        }
      }
    }
    #print("IV")
    if (gX_u.return && gX.return) {
      print("entra aqui 2")
      gX = matrix(0, model$N, model$q)
      gKX_uf = kernGradX(model$kern, X, X_u)
      for (i in 1:model$N) {
        for (j in 1:model$q) gX[i, j] = t(gKX_uf[, j, i, drop = FALSE]) %*% gK_uf[, i, drop = FALSE]
      }
    }
  }
  else stop("Unknown model approximation.")
  if (model$approx == "ftc") {
  }
  else if (model$approx == "dtc") {
  }
  else if (model$approx %in% c("fitc", "dtcvar")) {
    if (gX_u.return && gX.return) {
      gKXdiag = kernDiagGradX(model$kern, X)
      for (i in 1:model$N) gX[i, ] = gX[i, ] + gKXdiag[i, 
                                                       ] %*% gK_star[i]
    }
    g_param = g_param + kernDiagGradient(model$kern, X, gK_star)
  }
  else if (model$approx == "pitc") {
    if (gX_u.return && gX.return) {
      startVal = 1
      for (i in 1:length(model$blockEnd)) {
        endVal = model$blockEnd[i]
        ind = startVal:endVal
        gKXblock = kernGradX(model$kern, X[ind, , drop = FALSE], 
                             X[ind, , drop = FALSE])
        gKXblock = gKXblock * 2
        dgKXblock = kernDiagGradX(model$kern, X[ind, 
                                                , drop = FALSE])
        for (j in 1:length(ind)) gKXblock[j, , j] = dgKXblock[j, 
                                                              ]
        for (j in ind) {
          for (k in 1:model$q) {
            subInd = j - startVal + 1
            gX[j, k] = gX[j, k] + t(gKXblock[, k, subInd, 
                                             drop = FALSE]) %*% gK_star[[i]][, subInd, 
                                                                             drop = FALSE]
          }
        }
        startVal = endVal + 1
      }
    }
    for (i in 1:length(model$blockEnd)) {
      ind = gpBlockIndices(model, i)
      g_param = g_param + kernGradient(model$kern, X[ind, 
                                                     , drop = FALSE], gK_star[[i]])
    }
  }
  else stop("Unrecognised model approximation")
  if (!(gX_u.return && gX.return && g_beta.return)) {
    if ((!"optimiseBeta" %in% names(model) && model$approx != 
         "ftc") || model$optimiseBeta) 
      gParam = unlist(c(g_param, g_meanFunc, g_scaleBias, 
                        g_beta))
    else gParam = unlist(c(g_param, g_meanFunc, g_scaleBias))
  }
  else gParam = unlist(c(g_param, g_meanFunc, g_scaleBias))
  if (!(gX_u.return || gX.return || g_beta.return)) 
    gParam = c(gX_u, gParam)
  return(as.numeric(gParam))
}


######################################################################################################*

optLambdas4Mus <- function(trainData, learner, numFolds, mus, parallel=FALSE, plot=FALSE){
  
  fixedParams <- do.call(learner$getFixedParams, list(data=trainData))
  
  print("enters optLambdas4Mus")
  # for each mu get optimal lambdas
  opts <- lapply(mus, function(mu){
    fixedParams$mu <- mu
    opt  <- optLambda(trainData, learner,   numFolds=numFolds, parallel=parallel, plot=plot, fixedParams=fixedParams)
    return(opt$opt[[1]])
  })
  print("exits optLambdas4Mus")
  return(opts)
}


findRelevantMus <- function(trainData, learner, numFolds, numMus){
  # we are interested in finding the mus for which qhsic(mu) and mse(mu) are more or less equal
  # more or less meaning that the quotient is between 0.5 and 2 and we distribute points between these
  # two points equally in the  log-scale
  qn <- seq(log(0.5,10), log(2,10), length.out=numMus)
  
  
  mu_posts <- seq(log(0.001/0.999,10),log(0.999/0.001,10), length.out=11)
  mu_posts <- 10^(mu_posts)/(1+10^(mu_posts))
  
  mu_posts <- c(0,mu_posts,1)
  numPosts <- length(mu_posts)
  
  opt_posts <- optLambdas4Mus(trainData, learner, numFolds, mu_posts, parallel=TRUE, plot=FALSE)
  mse_posts <- sapply(opt_posts, function(opt) mse(learner$learn(trainData, opt), trainData, learner, opt))
  qhsic_posts <- sapply(opt_posts, function(opt) qhsicLoss(learner$learn(trainData, opt), trainData, learner, opt))  
  logDiff_posts <- log((1-mu_posts)*qhsic_posts,10)-log((mu_posts)*mse_posts,10)
  plot(mu_posts, logDiff_posts)
  abline(h=0, col="red")
  
  tol <- 0.01
  
  dists <- sapply(qn, function(q) sapply(logDiff_posts, function(post) abs(q-post)))
  # get best mus_post for each qn
  indx <- apply(dists, 2, which.min)
  mus <- mu_posts[indx]
  
  # get worst approx
  minDists <- apply(dists,2,min)
  maxMinDist <- max(minDists)
  
  # get posts with unacceptable distances
  indx_improve <- which(minDists>tol)
  indx_improve <- unique(indx[indx_improve])
  indx_improve <- sort(union(setdiff(indx_improve, length(mu_posts)), setdiff(indx_improve-1,0)))
  
  count <- 0
  
  while(maxMinDist>tol & count<= 15){
    count <- count +1
    print("********")
    print(paste("count: ", count))
    print(paste("maxMinDist: ", maxMinDist))
    
    # calculate new mus_posts and log diffs_posts
    mu_posts_new <- (mu_posts[2:numPosts]+mu_posts[1:(numPosts-1)])/2
    
    # use only mid points where we need closer approx, otherwise takes long ting
    mu_posts_new <- mu_posts_new[indx_improve]
    print(paste("length(mu_posts_new): ", length(mu_posts_new)))
    
    opt_posts_new <- optLambdas4Mus(trainData, learner, numFolds, mu_posts_new, parallel=TRUE, plot=TRUE)
    mse_posts_new <- sapply(opt_posts_new, function(opt) mse(learner$learn(trainData, opt), trainData, learner, opt))
    qhsic_posts_new <- sapply(opt_posts_new, function(opt) qhsicLoss(learner$learn(trainData, opt), trainData, learner, opt))  
    logDiff_posts_new <- log((1-mu_posts_new)*qhsic_posts_new,10)-log(mu_posts_new*mse_posts_new,10)
    mu_posts <- c(mu_posts, mu_posts_new)
    logDiff_posts <- c(logDiff_posts, logDiff_posts_new)
    indxPos <- order(mu_posts)
    mu_posts <- mu_posts[indxPos]
    logDiff_posts <- logDiff_posts[indxPos]
    numPosts <- length(mu_posts)
    
    plot(mu_posts, logDiff_posts)
    abline(h=0, col="red")
    
    
    # 1 update, dists
    dists <- sapply(qn, function(q) sapply(logDiff_posts, function(post) abs(q-post)))
    # a) update bets mus
    indx <- apply(dists, 2, which.min)
    mus <- mu_posts[indx]
    # b) update worst approx
    minDists <- apply(dists,2,min)
    maxMinDist <- max(minDists)
    # c) update posts to add
    indx_improve <- which(minDists>tol)
    indx_improve <- unique(indx[indx_improve])
    indx_improve <- sort(union(setdiff(indx_improve, length(mu_posts)), setdiff(indx_improve-1,0)))
    print(paste("length(indx_improve): ", length(indx_improve)))
    
  }
  
  return(sort(mus))
}

plot.optLambas4Mus <- function(opts, trainData, learner, numFolds){      
  
  mus <- sapply(opts, function(el) el$mu)
  lambdas <- sapply(opts, function(el) el$lambda)
  #plot(mus, lambdas)
  
  mses <- sapply(opts, function(opt) mse(learner$learn(trainData, opt), trainData, learner, opt))
  qhsics <- sapply(opts, function(opt) qhsicLoss(learner$learn(trainData, opt), trainData, learner, opt))
  
  facs <- (mus*mses)/((1-mus)*qhsics)
  spl <- spline(facs, mus)
  #plot(facs, mus)
  #lines(spl, col="blue")
  splf <- splinefun(facs, mus)
  facs.x <- seq(min(facs),max(facs),length.out=100)
  mus.x <- splf(facs.x, deriv=0)
  dMuDfac.x <- splf(facs.x, deriv=1)
  # plot(facs.x, dMuDfac.x)
  
  ########################################*
  # || alfa(mu) ||
  ########################################*
  
  alphas <- sapply(opts, function(opt) as.numeric(learner$learn(trainData, opt)$alpha))
  normAlphas <- apply(alphas, 2, function(col) sqrt(sum(col^2)))  
  
  ########################################*
  # || d alfa / d fac ||
  ########################################*
  {
    # obtain lambda(mu) and lambda'(mu) approximations
    spl <- spline(mus, lambdas)
    # plot(mus,lambdas)
    # lines(spl, col="blue")
    
    splf <- splinefun(mus, lambdas)
    mus.x <- seq(0,1,length.out=100)
    lambda.x <- splf(mus.x, deriv=0)
    dLambdaDmu.x <- splf(mus.x, deriv=1)
    # plot(mus.x, dLambdaDmu.x)
    
    # calculate kernel matrices
    fixedParams <- learner$getFixedParams(trainData)
    Kxb <- kernelMatrix(do.call(fixedParams$kernelXb, fixedParams$psKernXb), trainData$x)
    Kxs <- kernelMatrix(do.call(fixedParams$kernelXs, fixedParams$psKernXs), trainData$x)
    n <- nrow(Kxs)
    H <- diag(n)-matrix(1/n,n,n)
    I <- diag(n)
    Kxbc <- H%*%Kxb%*%H
    
    # calculate || d alfa / d fac ||
    dAlphaDfacNorm <- mapply(FUN=function(mu, lambda, dLambdaDmu, dMuDfac){
      # alpha = Amu^-1 %*% B %*% y
      # where 
      #        Amu = (mu*H + (1-mu)*Kxbc) %*% Kxs + lambda*I
      #          B = mu*H + (1-mu)*Kxbc
      B <- mu*H + (1-mu)*Kxbc
      Amu <- (mu*H+(1-mu)*Kxbc)%*%Kxs + lambda*I
      
      # d alfa / d mu = (d(Amu^-1) / d mu) %*% B %*% y 
      #   where d(Amu^-1) / d mu = -Amu^-1 %*% d Amu / d mu %*% Amu^-1,
      #     where d Amu / d mu = (dB / dmu) %*% Kxs + lambda'(mu)*I
      #       where dB/dmu = H - Kxbc
      dBDmu <- H-Kxbc
      dAmuDmu <- dBDmu%*%Kxs + dLambdaDmu*I
      AmuInv <- Matrix(solve(Amu))
      dAmuInvDmu <- -AmuInv %*% dAmuDmu %*% AmuInv
      dAlphaDmu <- dAmuInvDmu%*%B%*%  trainData$y 
      
      # d alfa / d fac = d alfa / d mu * d mu / d fac
      dAlphaDfac <- dAlphaDmu * dMuDfac
      dAlphaDfacNorm <- sqrt(as.numeric(t(dAlphaDfac)%*%dAlphaDfac))
      return(dAlphaDfacNorm)
    }, 
    mu=mus.x, lambda=lambda.x, dLambdaDmu=dLambdaDmu.x, dMuDfac=dMuDfac.x)
    # plot(facs.x, dAlphaDfacNorm)
    dAlphaDfacNorm_rel <- dAlphaDfacNorm/normAlphas
    # plot(facs.x, log(dAlphaDfacNorm_rel))
  }
  
  
  
  #####################################################*
  # || alfa(mu) - alpha(0)||, alpha(0) = alpha_qhsic
  # || alfa(mu) - alpha(1)||, alpha(1) = alpha_krr
  #####################################################*
  {
    # alpha_krr
    fixedParams1 <- fixedParams
    fixedParams1$mu <- 1
    
    opt  <- optLambda(trainData, learner, numFolds=numFolds, parallel=TRUE, plot=FALSE, fixedParams=fixedParams1)
    alpha_krr <- as.numeric(learner$learn(trainData, opt$opt[[1]])$alpha)
    normAlpha_krr <- sqrt(sum(alpha_krr^2))
    # alpha_qhsic
    fixedParams1$mu <- 0
    opt  <- optLambda(trainData, learner, numFolds=numFolds, parallel=TRUE, plot=FALSE, fixedParams=fixedParams1)
    alpha_qhsic <- as.numeric(learner$learn(trainData, opt$opt[[1]])$alpha)
    normAlpha_qhsic <- sqrt(sum(alpha_qhsic^2))
    
    # || alfa(mu) - alpha(1)||
    difNormKr <- alphas-alpha_krr
    difNormKr <- apply(difNormKr,2, function(col) sqrt(sum(col^2)))
    difNormKr_rel <- difNormKr / normAlpha_krr  
    # || alfa(mu) - alpha(0)||
    difNormQh <- alphas-alpha_qhsic
    difNormQh <- apply(difNormQh,2, function(col) sqrt(sum(col^2)))
    difNormQh_rel <- difNormQh / normAlpha_qhsic
    
    # plot(facs, difNormKr_rel)
    # plot(facs, difNormQh_rel)
  }
  
  #####################################################*
  # angle(alfa(mu), alpha(0)), alpha(0) = alpha_qhsic
  # angle(alfa(mu) - alpha(1)), alpha(1) = alpha_krr
  #####################################################*
  {
    # <alfa(mu), alpha(1)> = t(alfa(mu)) %*% alpha(1)
    prodIntKr <- alphas*alpha_krr
    prodIntKr <- apply(prodIntKr,2, sum)
    # <alfa(mu), alpha(1)> = t(alfa(mu)) %*% alpha(0)
    prodIntQh <- alphas*alpha_qhsic
    prodIntQh <- apply(prodIntQh,2, sum)
    # angle(alfa(mu), alpha(1))= arcosine(<alfa(mu), alpha(1)>/(||alpha(mu)||*||alpha(1)||))
    angleKr <- acos(prodIntKr/(normAlphas*normAlpha_krr))/pi*180
    # angle(alfa(mu), alpha(0))= arcosine(<alfa(mu), alpha(0)>/(||alpha(mu)||*||alpha(0)||))
    angleQh <- acos(prodIntQh/(normAlphas*normAlpha_qhsic))/pi*180
    
    # plot(facs, angleKr)
    # plot(facs, angleQh)
  }
  
  
  #####################################################*
  # Make plots 
  #####################################################*	                
  
  #plot(mus, lambdas)
  #lines(spl, col="blue")
  par(mfrow=c(2,2))
  plot(facs, normAlphas/((normAlpha_krr + normAlpha_qhsic) / 2))
  plot(facs.x, dAlphaDfacNorm_rel)
  plot(facs, difNormKr_rel)
  #plot(facs, difNormKr_rel, ylim=range(difNormKr_rel, difNormQh_rel))
  #lines(facs, difNormQh_rel, type="p",col="red")
  plot(facs, angleKr )
  #plot(facs, angleKr , ylim=range(angleKr, angleQh))
  #lines(facs, angleQh, type="p",col="red")
  
  
}


constructLearner <- function (hyperParams, learnParams, optimizeParams, hueristicSet, optimizeSet, learn, predict){
  # hyperParams should be a list with 2 named lists: data-dependent and non data-dependent
  stopifnot(is.list(hyperParams) && all(c("data", "non_data") %in% names(hyperParams)))
  # hyperParams$data should have 2 named lists: optimizable and non-optimizable parameters
  stopifnot(is.list(hyperParams$data) && all(c("optimizable", "non_optimizable") %in% names(hyperParams$data)))
  # optimizable, non_optimizible and non_data should all be lists
  stopifnot(is.list(hyperParams$data$optimizable) && is.list(hyperParams$data$non_optimizable) && is.list(hyperParams$non_data))            
  # learnParams should be a list, possibly to remain always empty (eg. bayesian learner with MCMC opt)
  stopifnot(is.list(learnParams))
  # optimizeParams is a list, it should at least include a vector of valid loss functions to be evaluated, the
  # first is considered the optimization loss. Could also include optimization method parameters (burn in, iterations,
  # convergence, etc)
  stopifnot(is.list(optimizeParams) && "losses" %in% names(optimizeParams))
  validLossFuncs <- c("sse", "mse", "rmse", "qhsicLoss", "hsicLoss","nqhsicLoss", "nhsicLoss", "qkricLoss", "nhsicReg", "hsicReg", 
                      "hsicYhReg", "nhsicRegA", "hsicRegA", "hsicYhRegA", "corre")
  stopifnot( all(sapply(optimizeParams$losses, function(el) is.function(el))) && all(names(optimizeParams$losses) %in% validLossFuncs))
  # heuristicSet, optimizeSet, learn and predict should all be functions
  stopifnot(is.function(learn) && is.function(predict) && is.function(heuristicSet) && is.function(optimizeSet))
  
  # getyperParams simply appends all hyperparameters- data_optimizable, data_non_optimizable and non_data in a list
  getHyperParams <- function(learner){
    return(c(learner$hyperParams$data$optimizable, learner$hyperParams$data$non_optimizable, learner$hyperParams$non_data))
  }
  
  learner <- list(hyperParams = hyperParams, learnParams=learnParams, optimizeParams=optimizeParams, 
                  getHyperParams= getHyperParams, 
                  heuristicSet=heuristicSet, optimizeSet=optimizeSet, learn=learn, predict=predict)
  class(learner) = "emley.learner"
  
  return(learner)
}


rbfSigma <- function(sigma) do.call("rbfdot", list(sigma=sigma))
rbfBeta <-  function(beta) do.call("rbfdot", list(sigma=beta))
rbfKappa <- function(kappa) do.call("rbfdot", list(sigma=kappa))
rbfGamma <- function(gamma) do.call("rbfdot", list(sigma=gamma))


setParams <- function(learner, trainData, plot=FALSE){
  #learner should be an object of class "emley.learner"
  stopifnot(class(learner)=="emley.learner")
  stopifnot(class(trainData)=="CVST.data")
  
  
  
  
  dataNonOptimParams <- learner$heuristicSet(learner, trainData)
  
  
  
  paramsNonOpt <- names(learner$hyperParams$data$non_optimizable)
  
  
  
  for(nm in paramsNonOpt){
    # nm <- paramsNonOpt[1]
    learner$hyperParams$data$non_optimizable[[nm]]$val <- dataNonOptimParams[[nm]]
  }
  
  
  learner$hyperParams$trainData <- trainData
  
  
  
  dataOptimParams <- learner$optimizeSet(learner, plot=plot)
  
  
  paramsOpt <- names(learner$hyperParams$data$optimizable) 
  for(nm in paramsOpt){
    # nm <- paramsOpt[1]
    learner$hyperParams$data$optimizable[[nm]]$val <- dataOptimParams$opt[[nm]]
  }
  
  
  
  return(learner)
  
} 

plot.emeley.1D <- function(predList){
  if(dim(predList[[1]]$x)[2]==1){
    
    if(length(predList)>1){
      aux <- sapply(predList, function(el){
        cbind(el$x, el$gy)
      }, simplify="array")
      
      aux <- aux[,,2:length(predList)]-aux[,,1:(length(predList)-1)]
    } else{
      aux <- 0
    }
    
    if(all(aux==0)){
      plot(predList[[1]]$x, predList[[1]]$gy, xlab="x",ylab="y")
      for(i in 1:length(predList)){
        xx <- predList[[i]]$x
        yy <- predList[[i]]$gyh
        o <- order(xx)
        xx <- xx[o]
        yy <- yy[o]
        lines(xx,yy, col="red", lwd=2)
      }
      legend("bottomleft", legend=names(predList), lwd=2, col="red")
      
    } else{
      print("data for different predictors not the same")
    }
    
    
    
  } else{
    print("this function is only for 1d regression")
  }
}


