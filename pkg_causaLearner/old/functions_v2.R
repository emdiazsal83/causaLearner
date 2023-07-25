library(CVST) # cross validation framework
library(kernlab) # kernels, krr
library(lbfgs) # non-convex optimization for hsic regression
library(numDeriv) #comparing numerical approx of gradient with analytical

library(graph) # graphNEL
library(pcalg) # unifDAG
library(kpcalg)
library(bnlearn) # count.graphs
library(gRbase) # topoSort
library(gtools) # permutations
library(dHSIC) # dhsic.test
library(earth) # mars model for initial residual bandwidth estimation
library(abind) #bind arrays
library(ggplot2)
library(reshape)
library(parallel) #mcmapply, mlapply
library(R.matlab)
library(gptk) #gaussian process toolkit

library(FNN) # knn.dist
library(fields) # rdist
library(pracma) # Toeplitz(a, b)
library(pROC) # roc, auc

library(WeightedROC) # WeightedROC, WeightedAUC
library(boot) # boot


#############################################################################################################
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
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
	Kxs <- kernelMatrix(kernelXs, x)
  Kxb <- kernelMatrix(kernelXb, x)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  alpha <- solve(Matrix(Kxbc%*%Kxs + diag(lambda, N))) %*% Kxbc %*% y
  learner$learnParams$alpha <- alpha
  learner$learnParams$avgy <- mean(y)

	return(learner)
}
predict.qhsic <- function (learner, data){
  kernelXs <- makeKernel(learner, "kernelXs")
  kxs <- kernelMatrix(kernelXs, data$x, learner$hyperParams$trainData$x)
  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
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



#############################################################################################################
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


######################################################################################################
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
  }
  
  if(is.null(indxSens)){
    indxSens <- 1:p
  }
  
  if(is.null(indxPred)){
    indxPred <- 1:p
  }
  
  x <- matrix(data$x, n, p)
  xx <- x[,indxPred]
  xs <- x[,indxSens]
  
  sigma0 <- 1/median(as.numeric(dist(xx)^2))
  beta0 <- 1/median(as.numeric(dist(xs)^2))
  
  #obtain residuals with a mars model
  model  <- earth(x=xx, y=data$y, nfold=5)
  yh0 <- as.numeric(predict(model, xx))
  kappa0 <- 1/median(as.numeric(dist(yh0)^2))
  
  res0 <- data$y-as.numeric(predict(model, xx))
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
  model <- gpCreate(q=dim(trainData$x)[2], d=1, X=trainData$x, y=as.matrix(trainData$y), options)
  return(list(modelInit=model))
}


######################################################################################################
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
      model$X_u = model$X[options$fixIndices, ]
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
  if ("optimiser" %in% names(model)) 
    optim = get(paste(model$optimiser, "optim", sep = ""), mode = "function")
  else optim = get("CGoptim", mode = "function")
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


######################################################################################################

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
      
  ########################################
  # || alfa(mu) ||
  ########################################
  
  alphas <- sapply(opts, function(opt) as.numeric(learner$learn(trainData, opt)$alpha))
  normAlphas <- apply(alphas, 2, function(col) sqrt(sum(col^2)))  
	               
  ########################################
  # || d alfa / d fac ||
  ########################################
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
    
  
  
  #####################################################
  # || alfa(mu) - alpha(0)||, alpha(0) = alpha_qhsic
  # || alfa(mu) - alpha(1)||, alpha(1) = alpha_krr
  #####################################################
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
  
  #####################################################
  # angle(alfa(mu), alpha(0)), alpha(0) = alpha_qhsic
  # angle(alfa(mu) - alpha(1)), alpha(1) = alpha_krr
  #####################################################
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
  
  
	#####################################################
	# Make plots 
	#####################################################	                
						
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
        lines(xx,yy, col=i)
      }
      legend("bottomleft", legend=names(predList), lwd=2, col=1:length(predList))
      
    } else{
      print("data for different predictors not the same")
    }
    
    
    
  } else{
    print("this function is only for 1d regression")
  }
}


########################################################################################
# Functions for applying hsic regression to causality
########################################################################################


# simulate from a given sem defined as a dag, set of functions and noise distributions
# defined as n = a*dist(n_iid, pars)^b, dag could be obtained from functions but we
# allow for this redundancy

simSEM <- function(sem){
  if(all(names(sem$simPars$nodes)!=nodes(sem$dag))) stop("names (sem,ns) dont match")
  if(all(names(sem$funcs)!=paste("f",nodes(sem$dag),sep=""))) stop("names (nodes, funcs) dont match")
  
  n <- sem$simPars$n
  p <- length(sem$simPars$nodes)
  
  
  
  ns <- mapply(FUN=function(dist, pars, a, b){
    
    a*(do.call(dist, c(list(n=n), pars)))^b
    
    
  }, dist=lapply(sem$simPars$nodes, function(el) el$dist), 
  pars=lapply(sem$simPars$nodes, function(el) el$pars),
  a=lapply(sem$simPars$nodes, function(el) el$a),
  b=lapply(sem$simPars$nodes, function(el) el$a))
  
  colnames(ns) <- nodes(sem$dag)
  
  dagAmat <- amat(as.bn(sem$dag))
  Nodes <- nodes(sem$dag)
  topoNodes <- topoSort(sem$dag)
  
  x <- matrix(NA, n, p)
  colnames(x) <- Nodes
  indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
  x[,Nodes[indxNoParents]] <- ns[,Nodes[indxNoParents]]
  
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    
    indxParents <- which(dagAmat[,nd]==1)
    Parents <- Nodes[indxParents]
    argsNd <- cbind(x[,Parents], ns[,nd])
    argsNd <- as.list(as.data.frame(argsNd))
    names(argsNd) <- c(Parents,"n")
    
    x[,nd] <- do.call(sem$funcs[[paste("f",nd,sep="")]], argsNd)
  }
  
  return(list(dag=dagAmat,x=x,n=ns))
}

simRandAddSEM <- function(p, n, nodes){
  # p <- 4; n <- 100 
  # nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
  # nodes <- rep(list(nodes),4)
  dag <- unifDAG(p)
  
  
  ns <- mapply(FUN=function(dist, pars, a, b){
    
    a*(do.call(dist, c(list(n=n), pars)))^b
    
    
  }, dist=lapply(nodes, function(el) el$dist), 
  pars=lapply(nodes, function(el) el$pars),
  a=lapply(nodes, function(el) el$a),
  b=lapply(nodes, function(el) el$a))
  
  colnames(ns) <- nodes(dag)
  
  dagAmat <- amat(as.bn(dag))
  Nodes <- nodes(dag)
  topoNodes <- topoSort(dag)
  
  x <- matrix(NA, n, p)
  colnames(x) <- Nodes
  indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
  x[,Nodes[indxNoParents]] <- ns[,Nodes[indxNoParents]]
  
  kernelX <- do.call("rbfdot", list(sigma=1))
  
  
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # (nd <- setdiff(topoNodes, Nodes[indxNoParents])[1])
    indxParents <- which(dagAmat[,nd]==1)
    Parents <- Nodes[indxParents]
    
    K <- kernelMatrix(kernelX, x[,Parents], matrix(rnorm(n*length(Parents), 0, 1), n, length(Parents)))
    alpha <- rnorm(n, 0 , 1)    
    fx <- as.numeric(K %*% alpha)
    
    # plot(x[,Parents[2]],fx)
    
    x[,nd] <- fx + ns[,nd]
    
    
  }
  
  return(list(dag=dagAmat,x=x,n=ns))
}


simRandAddSEMs <- function(q, p, n, nodes, nms=1:q){
  
  
  res <- mapply(FUN=function(p, n, nodes){
    simRandAddSEM(p, n, nodes)
  }, p=p, n=n, nodes=nodes, SIMPLIFY=F)
  
  dags <- lapply(res, function(el) el$dag)
  xs <- lapply(res, function(el) el$x)
  ns <- lapply(res, function(el) el$n)
  
  names(dags) <- nms
  names(xs) <- nms
  names(ns) <- nms
  
  return(list(dags=dags, xs=xs, ns=ns, names=nms))
  
}


# given train data and a graph, estimate corresponding non-linear additive SEM returning a  named list of models
fitSEMgivenDAG <- function(G, trainData, learner){
  
  
  if(class(G)=="graphNEL"){
    G <- amat(as.bn(G))
    Nodes <- colnames(G)
  } else if(class(G)=="matrix"){
    Nodes <- colnames(G) 
  } else{
    stop("G must be a graph or a matrix")
  }
  
  
  if(all(colnames(trainData)!=Nodes)) stop("names (G, trainData) dont match")
  
  
  
  topoNodes <- topoSort(G)
  p <- length(Nodes)
  nTrain <- nrow(trainData)
  
  learnerList <- list()
  learnerList[Nodes] <- list(NULL)
  
  
  indxNoParents <- which(apply(G,2, function(col) all(col==0)))
  
  
  print(paste("no parent nodes: ", Nodes[indxNoParents], sep=""))
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # nd <- "x"
    print("*********************************")
    print(paste("node: ", nd,sep=""))
    indxParents <- which(G[,nd]==1)
    Parents <- Nodes[indxParents]
    
    trainDataO <- constructData(as.matrix(trainData[,Parents]), trainData[,nd])
    
    
    learnerAux <- setParams(learner, trainDataO)
    learnerAux <- learnerAux$learn(learnerAux)
    
    learnerList[[nd]] <- learnerAux
    
  }
  
  return(learnerList)
  
}

# given data a graph and a list of models for each node in the graph, predict node variables
# obtain residuals and return predictions and residuals

predictSEM <- function(G, data, learnerList, plot=TRUE){
  
  
  if(class(G)=="graphNEL"){
    G <- amat(as.bn(G))
    Nodes <- colnames(G)
  } else if(class(G)=="matrix"){
    Nodes <- colnames(G) 
  } else{
    stop("G must be a graph or a matrix")
  }
  
  
  if(all(colnames(data)!=Nodes)) stop("names (G, data) dont match")
  if(all(names(learnerList)!=Nodes)) stop("names (learnerList, G) dont match")
  
  topoNodes <- topoSort(G)
  p <- length(Nodes)
  n <- nrow(data)
  
  
  
  indxNoParents <- which(apply(G,2, function(col) all(col==0)))
  
  
  preds <- matrix(NA, nrow(data), ncol(data))
  resids <- matrix(NA, nrow(data), ncol(data))
  colnames(preds) <- colnames(data)
  colnames(resids) <- colnames(data)
  
  
  preds[,indxNoParents] <- 0
  resids[,indxNoParents] <- data[,indxNoParents]
  
  print(paste("no parent nodes: ", Nodes[indxNoParents], sep=""))
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # nd <- "x"
    print("*********************************")
    print(paste("node: ", nd,sep=""))
    indxParents <- which(G[,nd]==1)
    Parents <- Nodes[indxParents]
    
    learner <- learnerList[[nd]]
    
    dataOb <- constructData(as.matrix(data[,Parents]), data[,nd])
    
    pred <- learner$predict(learner, data=dataOb)
    
    if(plot){
      predList <- list(pred)
      names(predList) <- nd  
      plot.emeley.1D(predList=predList)
    } 
    
    preds[,nd] <- pred$gyh
    resids[,nd] <- pred$gyh - pred$gy
    
  }
  

  return(list(preds=preds, resids=resids))
  
}


# obtains all possible unique regressions for each node (ie for all possible m-node dags)
getUniqueRegs <- function(allDags){
  # slice <- "w"
	uniqueRegs <- sapply(dimnames(allDags)[[2]], function(slice) unique(allDags[,slice,], MARGIN=2), simplify="array")
	dimnames(uniqueRegs) <- list(nodeFrom=dimnames(uniqueRegs)[[1]], numReg=1:(dim(uniqueRegs)[2]), nodeTo=dimnames(uniqueRegs)[[3]])
	return(uniqueRegs)
}

# obtains the unique regressions for a certain set of candidate dags
getUniqueRegsList <- function(dags){
  # slice <- "w"
  uniqueRegsList <- lapply(dimnames(dags)[[2]], function(slice){ 
    # slice <- "1"
    res <- unique( adrop(dags[,slice,,drop=FALSE], drop=2), MARGIN=2)
    dimnames(res) <- list(nodeFrom=dimnames(res)[[1]], numReg=1:(dim(res)[2]))
    return(res)
    })
  names(uniqueRegsList) <-c(dimnames(dags)[[2]])
  return(uniqueRegsList)
}


# obtains the models for all possible unique regressions implied in all m-node dags (output of getUniqueRegs)
fitSEMClassGivenDAGClass <- function(uniqueRegs, trainData, learner){
	

	nodes <- dimnames(uniqueRegs)[[1]]
	numRegs <- dimnames(uniqueRegs)[[2]]
	numRegsTot <- (length(numRegs)-1)*length(nodes)
	
	pm0 <- proc.time()
	learnerList <- lapply(nodes,  function(nodeTo){
		# nodeTo <- nodes[1]
	  print("*********************")
		print(paste(nodeTo, " regressions:", sep=""))
		 res <- lapply(numRegs, function(numReg){
			 # numReg <- numRegs[2]
		   print("******")
			 print(paste(nodeTo,", reg # ", numReg, sep=""))
			 indxPreds <- which(uniqueRegs[,numReg,nodeTo]==1)
			 	if(length(indxPreds)>0){
					
			 	  trainDataO <- constructData(as.matrix(trainData[,nodes[indxPreds]]), trainData[,nodeTo])
			 	  learnerAux <- setParams(learner, trainDataO)
			 	  learnerAux <- learnerAux$learn(learnerAux)
			 	  
					print("estimated time to completion:")
					numRegsLeft <- ((length(nodes)-match(nodeTo, nodes))*(length(numRegs)-1)+(length(numRegs)-match(numReg,numRegs)))
					numRegsDone <- numRegsTot -numRegsLeft
					print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
				} else{
					learnerAux <- list(NULL)
				}
				return(learnerAux)
			})
		  names(res) <- numRegs
		  return(res)
		})
	print(proc.time()-pm0) # 5.2 mins for krr1, 250 train points 
	
	names(learnerList) <- nodes
	
	
	return(learnerList)
}

# uses models to obtain predictions and residuals for all possible unique regression implied in all m-node dags
predictSEMClass <- function(uniqueRegs, data, learnerList, plot=TRUE){
  
  dims <- dimnames(uniqueRegs)[-1]
  n <- dim(data)[1]
  dims$sim <- 1:n
  dims$predResid <- c("pred","resid")
  dims <- dims[c(3,4,1,2)]
  
  nodes <- dimnames(uniqueRegs)[[1]]
  numRegs <- dimnames(uniqueRegs)[[2]]
  numRegsTot <- (length(numRegs)-1)*length(nodes)
  
  pm0 <- proc.time()
  predResid <- sapply(nodes,  function(nodeTo){
    print("*********************")
    print(paste(nodeTo, " regressions:", sep=""))
    sapply(numRegs, function(numReg){
      # nodeTo <- "w"; numReg <- "2"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegs[,numReg,nodeTo]==1)
      if(length(indxPreds)>0){
        
        dataOb <- constructData(as.matrix(data[,nodes[indxPreds]]), data[,nodeTo])
        learnerAux <- learnerList[[nodeTo]][[numReg]]
        pred <- learnerAux$predict(learnerAux, data=dataOb)
        if(plot){
          predList <- list(pred)
          names(predList) <- nodeTo
          plot.emeley.1D(predList=predList)
        } 
        preds  <- pred$gyh
        resids <- pred$gyh - pred$gy
        res <- cbind(preds, resids)
        
        
        print("estimated time to completion:")
        numRegsLeft <- ((length(nodes)-match(nodeTo, nodes))*(length(numRegs)-1)+(length(numRegs)-match(numReg,numRegs)))
        numRegsDone <- numRegsTot -numRegsLeft
        print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        res <- cbind(rep(0,n), data[,nodeTo])
      }
      return(res)
    }, simplify="array")
  }, simplify="array")
  print(proc.time()-pm0) # 30 mins for krr, 
  
  dimnames(predResid) <- dims 
  
  return(predResid)
}

# obtains the models for unique regressions implied in a set of dags (output of getUniqueRegsList)
fitSEMSetGivenDAGSet <- function(uniqueRegsList, trainData, learner){
  
  
  n <- dim(trainData)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  #count <- 0
  #pm0 <- proc.time()
  learnerList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "4"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "8"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>0){
        #count <- count + 1
        
        trainDataO <- constructData(as.matrix(trainData[,nodes[indxPreds]]), trainData[,nodeTo])
        learnerAux <- setParams(learner, trainDataO)
        learnerAux <- learnerAux$learn(learnerAux)
        
        
        
        #print("estimated time to completion:")
        #numRegsLeft <- numRegsTot - count
        #numRegsDone <- numRegsTot -numRegsLeft
        #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        learnerAux <- list(NULL)
      }
      return(learnerAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
    })
  #print(proc.time()-pm0) #  
  
  names(learnerList) <- nodes
  
  return(learnerList)
}

# uses models to obtain predictions and residuals for set of candidate dags
predictSEMSet <- function(uniqueRegsList, data, learnerList, plot=TRUE){
  
  
  n <- dim(data)[1]
  dims <- lapply(uniqueRegsList, function(el){
    # el <- uniqueRegsList[[1]]
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  #count <- 0
  #pm0 <- proc.time()
  predResids <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    
    res <- sapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "2"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>0){
        #count <- count + 1
        
        dataO <- constructData(as.matrix(data[,nodes[indxPreds]]), data[,nodeTo])
        learnerAux <- learnerList[[nodeTo]][[numReg]]
        pred <- learnerAux$predict(learnerAux, data=dataO)
        if(plot){
          predList <- list(pred)
          names(predList) <- nodeTo
          plot.emeley.1D(predList=predList)
        } 
        preds  <- pred$gyh
        resids <- pred$gyh - pred$gy
        res <- cbind(pred=preds, resid=resids)
        
        #print("estimated time to completion:")
        #numRegsLeft <- numRegsTot - count
        #numRegsDone <- numRegsTot -numRegsLeft
        #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        res <- cbind(pred=rep(0,n), resid=data[,nodeTo])
      }
      return(res)
    }, simplify="array")
    
    
    dimnames(res) <- list(numObs=1:n, predResid=c("pred","resid"), numReg=numRegs[[nodeTo]])
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(predResids) <- nodes
  
  return(predResids)
}


# obtains the p-values of fitted residuals in all m-mode dags
# NOTA!: CREO QUE VALE LA PENA PONER CMPLX FUNCS, Y PARS EN UNA LISTA Y SI SON VARIAS EN UNA LISTA DE LISTAS

complexityScore <- function(allDags, uniqueRegs, vars, cmplxScores){
	
	#indxDag <- 543
	#allDags[,,indxDag]
	
  matlabFunctions <- c("Shannon_KDP",  "Shannon_vME")
  # Functions which are in matlab but which I have now programmed in R
  # "Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", 
  # "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", "Shannon_spacing_LL", "Shannon_PSD_SzegoT", 
  #"Shannon_Edgeworth", "Shannon_MaxEnt1", "Shannon_MaxEnt2", "Shannon_expF",
  
  matlab <- rep(FALSE, length(cmplxScores))
  indxEntropy <- which(sapply(cmplxScores, function(el) el$func) == "score_sumMarginalEntropies")
  
  type <- matlab
  
  entropyEstimators <- sapply( cmplxScores[indxEntropy], function(el) el$pars$type)
  
  if(length(indxEntropy)>0) type[indxEntropy] <- entropyEstimators
  
  indxMatlab <- which(entropyEstimators %in% matlabFunctions)
  matlab[indxEntropy][indxMatlab] <- TRUE
  
  data.frame(sapply(cmplxScores, function(el) el$func), matlab, type)
  
  if(any(matlab)){
    Matlab$startServer(port=9998)
    ## Connect to Matlab session
    matlabSession <- Matlab(port=9998)
    setVerbose(matlabSession, threshold=200000)
    open(matlabSession)
    # Load ITE library into Matlab session
    evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/matlab'))")
    #evaluate(matlab, "help HShannon_kNN_k_initialization")
    
  } else{
    matlabSession <- NULL
  }
  
  cmplxFuncs <- sapply(cmplxScores, function(el) el$func)
  pars <- lapply(cmplxScores, function(el) el$pars)
	
	nodes <- dimnames(allDags)[[1]]
	numDags <- dim(allDags)[3]
	count <- 0
	pm0 <- proc.time()
	scores <- apply(allDags, "dag", function(dag){
		#get residuals for the four regressions implied in each column
		# dag <- allDags[,,1]
	  count <<- count + 1
		print("*******************")
		print(paste("dag # ",count, sep=""))
		
		# obtain residuals/vars to evaluate for each dag
		varsDag <- sapply(nodes, function(nodeTo){
			children <- dag[,nodeTo]
			indxReg <- which(apply(uniqueRegs[,,nodeTo],2, function(col) all(col==children)))
			rsNode <- vars[,indxReg, nodeTo]
			return(rsNode)
		}, simplify="array")
		
	
		#score <- dhsic.test(residsDag, matrix.input=TRUE, pairwise=FALSE)[[type]]
		
		#PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
		# ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
		# ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
		# COMBINAR COMO CORRESPONDE
		
		
		
		score <- mapply(FUN=function(func, p, m, matlabSession){
		    # i<-4; func <- cmplxFuncs[i]; p <- pars[[i]]; m <- matlab[i]
		    print(paste("function: ", func))
		    print("parameters:")
		    print(p)
		    if(m){
		      ps <- c(p, list(vars=varsDag), matlabSession=list(matlabSession))
		    } else{
		      ps <- c(p, list(vars=varsDag))
		    }
		    do.call(func, ps)
		  }, func=cmplxFuncs, p=pars, m=matlab, MoreArgs=list(matlabSession=matlabSession)) 
		
	
		if(length(indxEntropy)>0) names(score)[indxEntropy] <- entropyEstimators
		
		print("estimated time until completion")
		print((proc.time()-pm0)*(numDags-count)/count)
		return(score)
	})
	print(proc.time()-pm0) #
	
	if(any(matlab)){
	  close(matlabSession)
	}
	

	dims <- dimnames(scores)
	dims <- list(cmplxFunc=dims[[1]], dag=dims$dag)
	dimnames(scores) <- dims
	scores <- aperm(scores, c(2,1))
	
	return(scores)
}

complexityScoreList <- function(dags, uniqueRegsList, vars, cmplxScores, prnt=FALSE){
  
  
  matlabFunctions <- c("Shannon_KDP",  "Shannon_vME")
  # Functions which are in matlab but which I have now programmed in R
  # "Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", 
  # "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", "Shannon_spacing_LL", "Shannon_PSD_SzegoT", 
  #"Shannon_Edgeworth", "Shannon_MaxEnt1", "Shannon_MaxEnt2", "Shannon_expF",
  
  matlab <- rep(FALSE, length(cmplxScores))
  indxEntropy <- which(sapply(cmplxScores, function(el) el$func) == "score_sumMarginalEntropies")
  type <- matlab
  entropyEstimators <- sapply( cmplxScores[indxEntropy], function(el) el$pars$type)
  if(length(indxEntropy)>0) type[indxEntropy] <- entropyEstimators
  indxMatlab <- which(entropyEstimators %in% matlabFunctions)
  matlab[indxEntropy][indxMatlab] <- TRUE
  data.frame(sapply(cmplxScores, function(el) el$func), matlab, type)
  if(any(matlab)){
    Matlab$startServer(port=9998)
    ## Connect to Matlab session
    matlabSession <- Matlab(port=9998)
    setVerbose(matlabSession, threshold=200000)
    open(matlabSession)
    # Load ITE library into Matlab session
    evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/matlab'))")
    #evaluate(matlab, "help HShannon_kNN_k_initialization")
    
  } else{
    matlabSession <- NULL
  }
  
  cmplxFuncs <- sapply(cmplxScores, function(el) el$func)
  pars <- lapply(cmplxScores, function(el) el$pars)
  
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,31]
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
      
    # obtain residuals/vars to evaluate for each dag
    varsDag <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      children <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==children)))
      rsNode <- vars[[nodeTo]][,indxReg]
      return(rsNode)
    }, simplify="array")
    
    #PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
    # ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
    # ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
    # COMBINAR COMO CORRESPONDE
    
    score <- mapply(FUN=function(func, p, m, matlabSession){
      # i<-7; func <- cmplxFuncs[i]; p <- pars[[i]]; m <- matlab[i]
      if(prnt){
        print(paste("function: ", func))
        print("parameters:")
        print(p)
      }
      if(m){
        ps <- c(p, list(vars=varsDag), matlabSession=list(matlabSession))
      } else{
        ps <- c(p, list(vars=varsDag))
      }
      return(do.call(func, ps))
    }, func=cmplxFuncs, p=pars, m=matlab, MoreArgs=list(matlabSession=matlabSession)) 
    
    if(length(indxEntropy)>0) names(score)[indxEntropy] <- entropyEstimators
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  })
  if(prnt) print(proc.time()-pm0) #
  
  if(any(matlab)){
    close(matlabSession)
  }
  
  
  dims <- dimnames(scores)
  dims <- list(cmplxFunc=dims[[1]], dag=dims$dag)
  dimnames(scores) <- dims
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

# Score functions for residuals (model based) or variables (non-parametric)

score_pvalHSIC <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- score$p.value
  return(score)
}

score_pHSIC <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- -log(score$p.value)
  return(score)
}

score_HSIC <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- score$statistic
  return(score)
}

score_HSIC_fix <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method, kernel="gaussian.fixed", bandwidth=1)
  score <- score$statistic
  return(score)
}

score_pvalUnifPart <- function(vars, method, numParts){
  
  n <- nrow(vars)
  
  
  
  # n <- 150; numParts <- 8
  numPerPart <- ceil(n/numParts)
  smpl <- rep(1:numParts, numPerPart)[1:n]
  smpl <- smpl[sample(1:n)]
  
  # we sample which residuals to use for each hsic pvalue but we need at least two in each batch
  #smpl <- c(rep(1:numParts, 2) , sample(1:numParts, size=n-numParts*2, replace=TRUE))
  

  pvalsPart <- sapply(1:numParts, function(part) {
    indxPart <- which(smpl==part)
    
    
    pval <- dhsic.test(vars[indxPart,], matrix.input=TRUE, pairwise=FALSE, method=method)
    pval <- pval$p.value
    return(pval)
  })
  
  
  
  pvalsPart <- unique(na.omit(pvalsPart))
  if(length(pvalsPart)>2){
    ks_test <- ks.test(pvalsPart, punif)  
    score <- ks_test$p.value
  } else{
    hsic_test <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
    score <- hsic_test$p.value
  } 
  
  
  
  return(score)
  
}

score_pvalUnifBoot <- function(vars, method, numSmpls){
  
  n <- nrow(vars)
  pvalsBoot <- sapply(1:numSmpls, function(i) {
    indxPart <- sample(1:n, replace=T)
    pval <- dhsic.test(vars[indxPart,], matrix.input=TRUE, pairwise=FALSE, method=method)
    pval <- pval$p.value
    return(pval)
  })
  ks_test <- ks.test(pvalsBoot, punif)
  score <- ks_test$statistic
  
  return(score)
  
}


score_sumMarginalEntropies <- function(vars, type=c("Shannon_1sp", "Shannon_Gauss","Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", 
                                                    "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", 
                                                    "Shannon_spacing_LL", "Shannon_KDP", "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1", 
                                                    "Shannon_MaxEnt2", "Shannon_expF", "Shannon_vME"), ...){
  
  matlabFunctions <- c(  "Shannon_KDP",  "Shannon_vME")
  
  # Las siguientes son funciones q están en matlab pero q ya programamos en R también
  # "Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", "Shannon_spacing_Vpconst",
  # "Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", "Shannon_spacing_LL",
  # "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1", "Shannon_MaxEnt2", "Shannon_expF",
  
  pars <- list(...)
  
  func <- type 
  if(type %in% matlabFunctions){
    func <- "genericMatlabEntropy"
    pars <- c(pars, type=type)
  }
  
  score <- apply(vars, 2, function(col) do.call(func, c(list(x=col), pars)))
  score <- sum(score)
  
  return(score)
  
}


# meta score functions

score_pUnifHSIC <- function(scoreMat){
  # scoreMat <- scores[[1]][[1]]
  newVar <- "score_pUnifHSIC"; expr <- "score_pvalUnifPart * score_pHSIC"
  newScore <- as.matrix(within(as.data.frame(scoreMat), eval(parse(text=paste(newVar, expr, sep=" <- ")))))
  
  dimnms <- dimnames(scoreMat)
  dimnms$cmplxFunc <- c(dimnms$cmplxFunc, newVar)
  
  dimnames(newScore) <- dimnms
  
  return(newScore)
}

score_bestpUnif_pHSIC <- function(scoreArr){
  if(dim(scoreArr)[1]>1){
    diffPvals <- apply(scoreArr[,"score_pvalUnifPart",], "recipe", function(x) min(x[-which.min(x)])-min(x))
  } else{
    diffPvals <- scoreArr[1,"score_pvalUnifPart",]
  }
  
  res <- scoreArr[,"score_pHSIC",which.max(diffPvals), drop=FALSE]    
  dimnames(res)$cmplxFunc <- "score_bestpUnif_pHSIC"
  res <- adrop(res, drop=3)
  return(res)
}


# Ranking and probability transformation function  score -> rank, score -> prob/conf

# function to put each score on gaussian 0,1 scale
gaussianize <- function(x,...){
  
  Fn <- ecdf(x)
  y <- Fn(x)
  nug <- 0.001
  y[which(y==1)] <- 1- nug
  z <- qnorm(y)
  #hist(z)
  
  return(z)
}

# function to express each hypothesis as a probability
scoreToProb <- function(x,...){
  
  Fn <- ecdf(x)
  y <- 1-Fn(x)
  
  #plot(x, y)
  
  nug <- 0.001
  y[which(y==0)] <- nug
  z <- y / sum(y)
  # plot(x, z)
  
  return(z)
}

# function to correct a given score for fact many hypothesis are "similar" (share many edges and non edges)
correctScoreToAdd <- function(x, hyps){
  
  p <- dim(hyps)[2]
  offDiag <- !diag(p)
  
  # take mean score for each edge and non edge
  
  unwrappedEdges <- t(apply(hyps, 3, function(col) col[offDiag]))
  unwrappedNonEdges <- t(apply((!hyps)*1, 3, function(col) col[offDiag]))
  
  
  meanScorePerEdge <- apply(unwrappedEdges*x,2,mean)
  meanScorePerNonEdge <- apply(unwrappedNonEdges*x,2,mean)
  
  
  # sum score of each edge  and non-edge in hypothesis 
  # to get corrected hypothesis score
  
  corrX <- as.numeric(unwrappedEdges %*% meanScorePerEdge + unwrappedNonEdges %*% meanScorePerNonEdge)
  
  return(corrX)
}




# Entropy functions

# univariate
Shannon_1sp <- function(x){
  # From Mooij et al. 2016 pg 26
  res <- sort(x)
  res <- unique(res)
  n <- length(x)
  res <- sum(log(abs(res[2:n]-res[1:(n-1)])))
  #res <- res/(n-1)
  #res <- res + digamma(n) - digamma(1)
  return(res)
}

# this is the univariate version
Shannon_Gauss <- function(x){
  res <- log(var(x))
  #res <- 0.5(log(2*pi*exp(1))+res)
  return(res)
}

# multivariate
Shannon_kNN_k_matlab <- function(x, matlabSession, k){
  if(class(x) != "matrix"){
    x <- matrix(x, length(x), 1)
  }
  
  
  # pass sample to matlab
  setVariable(matlabSession, x=x)
  #evaluate(matlab, "size(x)")
  
  # calculate entropy in matlab - this code varies by type of entropy
  evaluate(matlabSession, "mult=1;")
  parmInitString <- paste("co = HShannon_kNN_k_initialization(mult, {'k',",k,"});")
  evaluate(matlabSession, parmInitString)
  evaluate(matlabSession, "H = HShannon_kNN_k_estimation(x, co);")
  
  # Bring back to R
  res <- getVariable(matlabSession, c("H"))$H
  
  
  
  return(as.numeric(res))
}

# R version
# multivariate
Shannon_kNN_k <- function(x, k){
  if(class(x) != "matrix"){
    x <- matrix(x, length(x), 1)
  }
  d <- ncol(x)
  n <- nrow(x)
  squared_distances <- FNN:::knn.dist(x, k)^2
  V <- pi^(d/2) / gamma(d/2+1)
  H <-  log(n-1) - digamma(k) + log(V) + d / n * sum(log(sqrt(squared_distances[,k]))) 
  return(H)
  
}

# univariate
Shannon_spacing_V <- function(x){
  
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  x_sorted <- c(rep(x_sorted[1], m), x_sorted, rep(x_sorted[n], m))
  diffs <- x_sorted[(2*m+1):(n+2*m)] - x_sorted[1:n]
  H = mean(log (n / (2*m) * diffs))
  return(H)
}

# univariate
Shannon_spacing_Vb <- function(x){
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  
  A <- x_sorted[(1+m):n]
  B <- x_sorted[1:(n-m)]
  
  diffs <- A - B;
  
  b <- sum(1/(m:n)) + log(m/(n+1)) # bias correction
  
  H <- mean(log((n+1)/m*diffs)) + b
  return(H)
}

# univariate
Shannon_spacing_Vpconst <- function(x){
  
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  x_sorted <- c(rep(x_sorted[1], m), x_sorted, rep(x_sorted[n], m))
  differs <- x_sorted[(2*m+1):(n+2*m)] - x_sorted[1:n]
  c <- c(rep(1,m), rep(2, n-2*m), rep(1,m)) # piecewise constant correction
  H <- mean(log (n / m * differs/c))
  return(H)
}

# univariate
Shannon_spacing_Vplin <- function(x){
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  x_sorted <- c(rep(x_sorted[1], m), x_sorted, rep(x_sorted[n],m))
  diffs <- x_sorted[(2*m+1):(n+2*m)] - x_sorted[1:n]
  C <- c(1+((1:m)-1)/m, rep(2,n-2*m), 1+(n-((n-m+1):n))/m) # piecewise linear correction
  H <- mean(log (n/ m * diffs/C))
  return(H)
}

# univariate
Shannon_spacing_Vplin2 <- function(x){
  n <- length(x)
  m <- floor(sqrt(n)) 
  x_sorted <- sort(x);
  x_sorted <- c(rep(x_sorted[1],m), x_sorted, rep(x_sorted[n],m)) #%with the smallest (left) and largest (right) element
  diffs <- x_sorted[(2*m+1):(n+2*m)] - x_sorted[1:n]
  c1 <- 1 + ((1:m)+1)/m - (1:m)/(m^2)
  c2 <- rep(2,n-2*m-1)
  c3 <- 1 + (n-((n-m):n))/(m+1)
  C <- c(c1,c2,c3)
  H <- mean(log (n / m * diffs/C))
  return(H)
}

# univariate
Shannon_spacing_VKDE <- function(x){
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  
  stdX <- sd(x) # sample standard deviation
  h <- 1.06 * stdX * n^(-1/5)
  A <- x_sorted[c(1:m, (n-m+1):n)]/h 
  B <- x_sorted/h
  sD <- fields:::rdist(A,B)^2 #squared distances between each real i.e. distance corresponds  absolute value
  s13 <- apply(exp(-sD/2),2, mean) / (sqrt(2*pi)*h) 
  s2 <- (2*m/n) / (x_sorted[(2*m+1):n] - x_sorted[1:(n-2*m)])
    
  H <- -mean(log(c(s13, s2)))
  return(H)
}

# univariate
Shannon_spacing_LL <- function(x){
  n <- length(x)
  m <- floor(sqrt(n))
  x_sorted <- sort(x)
  x_sorted <- c(rep(x_sorted[1],m), x_sorted, rep(x_sorted[n],m))
  
  
  # we need rolling indices which are padded when we are near the edges of a vector and 
  # there are not enough values to the left or to the right
  
  # even numbered vectors shd be centered on the left of the two central indices ex: if winLen = 6
  # then for vector x = c(1, 2 ,3, 4, 5, 6, 7, 8, 9, 10) the rolling windows would be
  
  # c(0, 0, 1,  2,  3,  4 ); center <- 1
  # c(0, 1, 2,  3,  4,  5 ); center <- 2
  # c(1, 2, 3,  4,  5,  6 ); center <- 3
  # c(2, 3, 4,  5,  6,  7 ); center <- 4
  # c(3, 4, 5,  6,  7,  8 ); center <- 5
  # c(4, 5, 6,  7,  8,  9 ); center <- 6
  # c(5, 6, 7,  8,  9,  10); center <- 7
  # c(6, 7, 8,  9,  10, 0 ); center <- 8
  # c(7, 8, 9,  10, 0,  0 ); center <- 9
  # c(8, 9, 10, 0,  0,  0 ); center <- 10
  
  
  len <- length(x_sorted)
  winLen <- 2*m+1
  
  b <- sapply(1:length(x_sorted), function(i){
    indx <- (i-m):(i+m)
    indx[which(indx<=0 | indx> len)] <- NA
    res <- x_sorted[indx]
    res[which(is.na(res))] <- 0
    res <- res - mean(res)
    
    res <- as.numeric(t(-m:m) %*% res) / sum(res^2)
    
    return(res)
  })
  
  
  # b <- rep(0, len)
  # 
  # B = matrix(0, winLen, len)
  # 
  # for(i in 1:len){
  #   print(i)
  #   mid <- i + ((winLen+1)%%2)/2 
  #   CEIL <- ceiling(mid-floor(winLen/2))
  #   padLeft <- sum((CEIL:1)<=0)
  #   ini <- max(1, CEIL)
  #   FLOOR <- floor(mid + floor(winLen/2))
  #   padRight <- sum(((len):FLOOR)>len)
  #   fin <- min(len, FLOOR)
  #   z <- x_sorted[ini:fin]
  #   z <- c(rep(0, padLeft), z, rep(0, padRight)) 
  #   B[,i] <- z
  # }
  #     
  # 
  # for(i in 1:len){
  #   aux <- B[,i] - mean(B[,i])
  #   w <- (t(-m:m) %*% aux) / sum(aux^2)
  #   b[i] <- w
  # }
  
  
  
  b <- b[(m+1):(len-m)]
  H <- -mean(log(b/n))
  return(H)
  
}

# univariate
Shannon_PSD_SzegoT <- function(x){
  p_max <- 10 
  K <- 100 
  #x <- rnorm(10000)
  
  a <- 1/(2*max(abs(x)))

  
  x <- x*a
  
  ###################################################
  #        Estimation of the AR parameters          #
  ###################################################
  
  #all_a_AR <- matrix(0, p_max-1,p_max)
  #all_e <- rep(0,p_max-1)
  
  all_a_AR <- matrix(0, p_max, p_max +1)
  all_e <- rep(0,p_max)
  
  # i is complex number so phi_x is too
  phi_x <- apply(exp(1i*2*pi*(0:p_max)%*%t(x)),1, mean)

 for (order in 1:p_max){       
    # LONG AR
    # complex number matrix
    # order <- p_max
    #print(order)
    
    if(order== 1){
      R <- matrix(phi_x[1], 1, 1)
    } else{
      R <- Toeplitz(phi_x[1:order], phi_x[1:order])
    }
    r <- phi_x[2:(order+1)]
    Delta <- diag((1:order)^4)
    aux <- solve(-(R + (1e-5)*Delta), r)
    a_AR <- c(1 , aux)    
    all_a_AR[order, 1:(order+1)] <- a_AR
    all_e[order] <- as.numeric(Re(t(a_AR) %*% Conj(phi_x[1:(order+1)])))
  }
  # all_aAR is complex
  # all_e is real

  MDL <- length(x)*log(all_e) + (1:length(all_e))*log(length(x))

  p <- which.min(MDL)
  min_value <- MDL[p]

  a_AR <- all_a_AR[p,1:(p+1)]
  
  ###################################################
  #     Extrapolation of the autocorrelation        #
  ###################################################


  L <- K*(p+1)

  phi_x_tilde <- rep(0, L)
  phi_x_tilde[1:(p+1)] <- phi_x[1:(p+1)]

  for(kk in (p+2):L){
      # kk <- 20
      phi_x_tilde[kk] <- - as.complex(t(a_AR[2:length(a_AR)]) %*% phi_x_tilde[(kk-1):(kk-p)])
  }
  
  R <- Toeplitz(phi_x_tilde, phi_x_tilde)

  ######################################################
  #  Estimation of the Entropy using Szeg�'s Theorem   #
  ######################################################

  lambdas <- pracma:::eig(R)
    
  H <- -mean(Re(log2(lambdas)*lambdas))
  
  #H2 <- H - log2(a)
  #H2 <- H2 / log2(exp(1))
  #H3 <- H/log2(exp(1)) - log(a)
  #H4 <- log(2^H) - log(a)
  
  H <- log(2^H / a)
  
  #H; H2; H3; H4
  
  return(H)
  
}

# multivariate
Shannon_Edgeworth <- function(x){
  #x <- matrix(rnorm(10000*3), 10000, 10)
  #x <- matrix(rnorm(10000*1), 10000, 1)
  if(class(x) != "matrix"){
    x <- matrix(x, length(x), 1)
  }
  n <- nrow(x)
  d <- ncol(x)
  
  
  # normalize Y to have zero mean and unit std:
  x <- apply(x, 2, function(col) col-mean(col)) # E=0, this step does not change the Shannon entropy of the variable
    
  s <- sqrt(apply(x^2, 2, sum)/(n-1))
  x <- t(t(x) / s)
  
  
  
  H_whiten <- log(prod(s)) #we will take this scaling into account via the entropy transformation rule [ H(Wz) = H(z)+log(|det(W)|) ] at the end
      
  dete <- det(cov(x))
  #print(dete)
  dete <- max(dete, 1e-300)
  
  H_normal <- log(dete)/2 + d/2 * log(2*pi) + d/2 #Shannon entropy of a normal variable with cov(Y.') covariance.
                         
  #t1b:
    # t1b <- 0
    # for(i in 1:d){ #d terms
    #   kappa_iii <- mean(x[,i]^3) 
    #   t1b <- t1b + kappa_iii^2
    # }
  
  t1 <- sum(apply(x^3, 2, mean)^2)
  #t1-t1b              
           
  #t2:
   # pm <- proc.time()
   # t2b <- 0
   # for(i in 1:d){
   #   for(j in setdiff(1:d,i)){ #j\ne i; 2*nchoosek(d,2) terms
   #     kappa_iij <- mean(x[,i]^2 * x[,j])
   #     t2b <- t2b + kappa_iij^2
   #   }
   # }
   # proc.time()-pm
  
  #pm <- proc.time()
  t2 <- sum(((t(x^2) %*% x/n)[diag(d)==0])^2) 
  #proc.time()-pm
  #t2-t2b
  
  t2 <- 3 * t2
                         
  #t3:
  # pm <- proc.time()
  # t3b <- 0
  # for(i in 1:(d-2)){ #i<j<k; nchoosek(d,3) terms
  #   for(j in (i+1):(d-1)){
  #     for(k in (j+1):d){
  #       kappa_ijk <- mean(x[,i]*x[,j]*x[,k])
  #       t3b <- t3b + kappa_ijk^2
  #     }
  #   }
  # }
  # proc.time()-pm
  
  
  #pm <- proc.time()
  indx <- expand.grid(i=1:max(1,d-2), j=1:max(1, d-1), k=1:d)
  indx <- indx[which(indx$j>indx$i & indx$k>indx$j),]
  t3 <- sum(apply(indx, 1, function(row) mean(x[,row[1]]*x[,row[2]]*x[,row[3]])^2))
  #proc.time() - pm
  #t3-t3b
  
  
  
  
  t3 <- t3 / 6
  
  H <- (H_normal - (t1+t2+t3) / 12) + H_whiten
  return(H)
}

# univariate
Shannon_MaxEnt1 <- function(x){
  n <- length(x)
  
  # normalize Y to have zero mean and unit std:
  x <- x - mean(x) # E=0, this step does not change the Shannon entropy of the variable
  
  s <- sqrt(sum(x^2)/(n-1))
  x <- x/s
  
  H_whiten <- log(s) #we will take this scaling into account via the entropy transformation rule [ H(wz) = H(z)+log(|w|) ] at the end
  
  # H1,H2 -> H:
  H1 <- ( 1+log(2*pi) ) / 2 # %=H[N(0,1)]
  
  #H2:
  k1 <- 36 / ( 8*sqrt(3) - 9 )
  k2a <- 1 / ( 2 - 6/pi )
  H2 <- k1 * mean(x * exp(-x^2/2))^2 + k2a * (mean(abs(x)) - sqrt(2/pi))^2
  H <- H1 - H2
  
  #take into account the 'std=1' pre-processing:
  H <- H + H_whiten
  return(H)
}

# univariate
Shannon_MaxEnt2 <- function(x){
  n <- length(x)
  
  # normalize Y to have zero mean and unit std:
  x <- x - mean(x) # E=0, this step does not change the Shannon entropy of the variable
  
  s <- sqrt(sum(x^2)/(n-1))
  x <- x/s
  
  H_whiten <- log(s) #we will take this scaling into account via the entropy transformation rule [ H(wz) = H(z)+log(|w|) ] at the end
  
  # H1,H2 -> H:
  H1 <- ( 1+log(2*pi) ) / 2 # %=H[N(0,1)]
  
  #H2:
  k1 <- 36 / ( 8*sqrt(3) - 9 )
  k2b <- 24 / (16*sqrt(3) - 27)
  H2 <- k1 * mean(x * exp(-x^2/2))^2 + k2b * (mean(exp(-x^2/2)) - sqrt(1/2))^2
  H <- H1 - H2
  
  #take into account the 'std=1' pre-processing:
  H <- H + H_whiten
  return(H)
}

# multivariate
# for now multivariate gaussian entropy
Shannon_expF <- function(x){
  # n <- 10000; d<- 10; x <- matrix(rnorm(n*d), n, d)
  if(class(x) != "matrix"){
    x <- matrix(x, length(x), 1)
  }
  n <- nrow(x)
  d <- ncol(x)
  m <- apply(x, 2, mean)
  C <- cov(x)
  invC <- solve(C)
  t1 <- as.numeric(invC %*% m)
  t2 <- invC / 2
  In <- solve(t2)
  term1 <- sum(diag(In %*% t1 %*% t(t1))) /4 - log(det(t2)) / 2 + d * log(pi) / 2
  
  
  s <- In %*% t1
  gradF.t1 <- s / 2
  gradF.t2 <- -In/2 - (s  %*% t(s)) / 4
 

  term2 <- sum(t1 * gradF.t1) + sum(t2 * gradF.t2)

  H <-  term1 - term2 #assumption: k, the carrier measure is zero
  return(H)
  
}

# multivariate
Shannon_vME <- function(x){
  # n <- 10000; d<- 10; x <- matrix(rnorm(n*d), n, d)
  if(class(x) != "matrix"){
    x <- matrix(x, length(x), 1)
  }
  n <- nrow(x)
  d <- ncol(x)
  stop("currently not implemented")
  
  return(H) 
}

# Faltan 2: "Shannon_KDP",  "Shannon_vME"

# Falta agregar referencia ITE-package y referencia original a cada función

# Faltan comparar version matlab y R para asegurar q las programe bien


genericMatlabEntropy <- function(x, matlabSession, type=c("Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", 
                                                          "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", 
                                                          "Shannon_spacing_LL", "Shannon_KDP", "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1", 
                                                          "Shannon_MaxEnt2", "Shannon_expF", "Shannon_vME") , ...){
  
  pars <- list(...)
  x <- matrix(x, 1, length(x))
  
  # pass sample to matlab
  setVariable(matlabSession, x=x)
  #evaluate(matlab, "size(x)")
  
  # calculate entropy in matlab - this code varies by type of entropy
  evaluate(matlabSession, "mult=1;")
  initFuncString <- paste("H", type, "_initialization", sep="")
  
  if(length(pars)>0){
    parmString <- paste("'", names(pars) ,"',", sapply(pars, function(el) el), sep="")
    parmString <- paste(",{", parmString, "}", sep="")  
    parmInitString <- paste("co = ", initFuncString,"(mult", parmString,");", sep="")
  } else{
    parmInitString <- paste("co = ", initFuncString,"(mult);", sep="")
  }
  evaluate(matlabSession, parmInitString)
  estFuncString <- paste("H", type, "_estimation", sep="")
  evaluate(matlabSession, "H = ", estFuncString,"(x, co);", sep="")
  
  # Bring back to R
  res <- getVariable(matlabSession, c("H"))$H
  
  return(as.numeric(res))
}



#Matlab$startServer(port=9997)
#matlabSession <- Matlab(port=9997)
#setVerbose(matlabSession, threshold=200000)
#open(matlabSession)
# Load ITE library into Matlab session
#evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/matlab'))")
#x <- rnorm(100)
#Shannon_kNN_k(x, matlabSession, k=5)
#genericMatlabEntropy(x, matlabSession, type="Shannon_kNN_k", k=5)
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_V")
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vb") 
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vpconst")
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vplin")
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vplin2")
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_VKDE")
#genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_LL")
#genericMatlabEntropy(x, matlabSession, type="Shannon_KDP")
#genericMatlabEntropy(x, matlabSession, type="Shannon_PSD_SzegoT")
#genericMatlabEntropy(x, matlabSession, type="Shannon_Edgeworth")
#genericMatlabEntropy(x, matlabSession, type="Shannon_MaxEnt1")
#genericMatlabEntropy(x, matlabSession, type="Shannon_MaxEnt2")
#genericMatlabEntropy(x, matlabSession, type="Shannon_expF")
#genericMatlabEntropy(x, matlabSession, type="Shannon_vME")
#close(matlabSession)

udag2dags <- function(uDAG, checkValid=TRUE){
  numEdges <- sum(uDAG)/2
  numDags <- 2^numEdges
  numNodes <- nrow(uDAG)
  V <- colnames(uDAG)
  
  aux <- uDAG  
  aux[upper.tri(aux)] <- 0
  aux <- melt(aux, stringsAsFactors=F)
  aux <- aux[which(aux$value>0),1:2]
  aux$X1 <- as.character(aux$X1)
  aux$X2 <- as.character(aux$X2)
  
  
  lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
  upperOnes <- (lowerOnes==0)*1
  
  lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
  upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
  
  matOnes <- mapply(FUN=function(a,b){
    aux1 <- a
    colnames(aux1) <- c("x1","x2")
    aux2 <- b
    colnames(aux2) <- c("x2","x1")
    base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
  
  res <- matrix(0, numNodes, numNodes)
  rownames(res) <- V
  colnames(res) <- V
  
  dags <- sapply(matOnes, function(mat){
    # mat <- matOnes[[1]]
    mat2 <- mat[which(mat[,3]==1),]
    indxRow <- match(mat2[,1], V)
    indxCol <- match(mat2[,2], V)
    indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
    
    res2 <- res
    res2[indxMat] <-   1
    return(res2)
    
  }, simplify="array")
  
  # check for and delete graphs with cycles
  
  if(checkValid){
    indxValid <- which(apply(dags, c(3), function(mat){ 
      isValidGraph(mat, "dag")
    }))
    dags <- dags[,,indxValid]
  }
  
  return(dags)
  
}

pdag2dags <- function(pDAG){
  
  # undirected edges adjacency matrix
  uDAG <- (pDAG==1 & t(pDAG)==1)*1
  
  # number of undirected edges
  numUEdges <- sum(uDAG)/2
  
  # number of dags we can get from cpdag
  numDags <- 2^numUEdges
  numNodes <- nrow(pDAG)
  V <- colnames(pDAG)
  
  # if there is at least one undirected edge 
  # we permute those edges
  
  if(numDags>1){
  
    aux <- uDAG  
    aux[upper.tri(aux)] <- 0
    aux <- melt(aux, stringsAsFactors=F)
    aux <- aux[which(aux$value>0),1:2]
    aux$X1 <- as.character(aux$X1)
    aux$X2 <- as.character(aux$X2)
  
  
    lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
    upperOnes <- (lowerOnes==0)*1
  
    lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
    upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
  
    matOnes <- mapply(FUN=function(a,b){
      aux1 <- a
      colnames(aux1) <- c("x1","x2")
      aux2 <- b
      colnames(aux2) <- c("x2","x1")
      base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
  
    res <- pDAG - uDAG
  
  
    dags <- sapply(matOnes, function(mat){
      # mat <- matOnes[[1]]
      mat2 <- mat[which(mat[,3]==1),]
      indxRow <- match(mat2[,1], V)
      indxCol <- match(mat2[,2], V)
      indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
    
      res2 <- res
      res2[indxMat] <-   1
      return(res2)
    
    }, simplify="array")
  } else{
    dags <- pDAG
    dim(dags) <- c(dim(dags),1)
  }
  
  
  
  return(dags)
  
}

mindag2dags <- function(minDAG, checkValid=TRUE){
  
  noLinks <- (minDAG==0 & t(minDAG)==0)*1
  diag(noLinks) <- 0
  
  numEdges <- sum(noLinks)/2
  numDags <- 3^numEdges
  numNodes <- nrow(minDAG)
  V <- colnames(minDAG)
  
  # if there is at least one pair of nodes with no edges
  # we permute those
  if(numDags>1){
  
    aux <- noLinks  
    aux[upper.tri(aux)] <- 0
    aux <- melt(aux, stringsAsFactors=F)
    aux <- aux[which(aux$value>0),1:2]
    aux$X1 <- as.character(aux$X1)
    aux$X2 <- as.character(aux$X2)
  
  
    lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
    lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
  
  
    res <-  minDAG*0
  
    noLinks2 <- sapply(lowerOnes, function(mat){
    # mat <- matOnes[[1]]
    mat2 <- mat[which(mat[,3]==1),]
    indxRow <- match(mat2[,1], V)
    indxCol <- match(mat2[,2], V)
    indxMat1 <- matrix(c(indxRow, indxCol), length(indxRow), 2)
    indxMat2 <- matrix(c(indxCol, indxRow), length(indxRow), 2)
    
    res2 <- res
    res2[indxMat1] <-   1
    res2[indxMat2] <-   1
    return(res2)
    
  }, simplify="array")
  
    dagss <- sapply(1:(dim(noLinks2)[3]), function(i){
    noLinks <- noLinks2[,,i]
    if(sum(noLinks)==0){
      dags <- minDAG
      dim(dags) <- c(dim(dags),1)
    } else{
      aux <- noLinks  
      aux[upper.tri(aux)] <- 0
      aux <- melt(aux, stringsAsFactors=F)
      aux <- aux[which(aux$value>0),1:2]
      aux$X1 <- as.character(aux$X1)
      aux$X2 <- as.character(aux$X2)
      
      
      lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
      upperOnes <- (lowerOnes==0)*1
      
      lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
      upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
      
      matOnes <- mapply(FUN=function(a,b){
        aux1 <- a
        colnames(aux1) <- c("x1","x2")
        aux2 <- b
        colnames(aux2) <- c("x2","x1")
        base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
      
      res <- minDAG
      
      dags <- sapply(matOnes, function(mat){
        # mat <- matOnes[[1]]
        mat2 <- mat[which(mat[,3]==1),]
        indxRow <- match(mat2[,1], V)
        indxCol <- match(mat2[,2], V)
        indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
        
        res2 <- res
        res2[indxMat] <-   1
        return(res2)
        
      }, simplify="array")
    }
    return(dags)
    
  }, simplify="array")
  
    dagss <- do.call(what=abind, args=c(dagss, along=3) )
  
    if(checkValid){
      indxValid <- which(apply(dagss, c(3), function(mat){ 
        isValidGraph(mat, "dag")
      }))
      dagss <- dagss[,,indxValid]
    }
  } else{
    dagss <- minDAG
    dim(dagss) <- c(dim(dagss),1)
  }
  
  return(dagss)
  
}


intToBitss <- function(x, maxInt) paste(rev(strsplit(substr(paste(as.integer(intToBits(x-1)), collapse=""),1, ceiling(log(maxInt,2))), split="")[[1]]), collapse="")

BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+1

dagIDtoMatrix <- function(dagID, p){
  dagBin <- intToBitss(dagID, 2^(p^2))
  dag <- matrix(as.numeric(sapply(1:nchar(dagBin), function(i) substr(dagBin, i, i))), p , p)
  return(dag)
}


BinToDec("000")
BinToDec("001")
BinToDec("010")
BinToDec("011")
BinToDec("100")
BinToDec("101")
BinToDec("110")
BinToDec("111")


intToBitss(1, 8)
intToBitss(2, 8)
intToBitss(3, 8)
intToBitss(4, 8)
intToBitss(5, 8)
intToBitss(6, 8)
intToBitss(7, 8)
intToBitss(8, 8)



# obtain binary id string for set of hypotheses
getHypID <- function(allHyps){
  
  hypIDs <- data.frame(bin=apply(allHyps, 3, function(mat) paste(as.character(as.numeric(mat)), collapse="")))
  hypIDs$id <- sapply(hypIDs$bin, BinToDec)
  hypIDs <- hypIDs[,c("id","bin")]
  
  return(hypIDs)
}  


# generate all the DAGS with m-nodes
genAllDAGS <- function(m){
  # possible topological orderings
  topOrd <- permutations(m,m)
  numTopOrds <- nrow(topOrd)
  
  # obtain all permutations of DAG matrix for one topological ordering
  maxEdges <- m*(m-1)/2
  numEdgePerms = 2^maxEdges
  
  # create binary strings indicating whether a certain directed edge is switched on or not for all possible configs
  edgePerms <- sapply(0:(numEdgePerms-1), function(i) as.numeric(strsplit(substr(paste(as.integer(intToBits(i)), collapse=""),1,maxEdges),"")[[1]]))
  
  # now we put in adjacency matrix form by placing in 
  dags1topOrd <- array(0, dim=c(m,m,numEdgePerms))
  indMat <- matrix(1:m^2,m,m)
  dags1topOrd[array(upper.tri(indMat), dim=c(m,m,numEdgePerms))] <- as.numeric(edgePerms)
  
  # we now create for all topological orders by permuting
  dagsAlltopOrd <- sapply(1:numTopOrds, function(i){
    perm <-  topOrd[i,]
    aux <- dags1topOrd[,perm,]
    return(aux[perm,,])
  }, simplify="array")
  # stack adjacency matrices in one direction
  
  dim(dagsAlltopOrd) <- c(m,m,numEdgePerms*numTopOrds) 
  
  # remove duplicates
  dagsAlltopOrd <- unique(dagsAlltopOrd, MARGIN=3)
  return(dagsAlltopOrd)
  
  
}


causalOrderingMooij <- function(Xtr, Xte, learner,  method, alpha=0){
  print("enters causalOrderingMooij")
  d <- ncol(Xtr)
  S <- colnames(Xtr) # 1:d
  ord <- as.numeric()
  
  totRegs <- d*(d+1)/2 -1
  #count <- 0
  # we'll choose one node at a time to be next in back to front causal ordering
  #pm0 <- proc.time()
  for (iter in d:2){
    # iter <- 3
    #print("**********************")
    #print(paste("iter: ",iter))
    
    # for each  nodes not chosen that are left we use the other nodes not chosen 
    # and see if we get independent residuals. We choose "most" independent
    # regression and take the dependent variable to be next in back to front
    # causal ordering
    pvals <- sapply(S, function(effectNodeTry){
      # effectNodeTry <- S[1]
      #print(paste("iter: ", iter,"effectNodeTry: ", effectNodeTry))
      #count <<- count  + 1
      causeNodesTry <- setdiff(S, effectNodeTry)
      
      trainDataO <- constructData(as.matrix(Xtr[,causeNodesTry]), Xtr[,effectNodeTry])
      
      learnerAux <- setParams(learner, trainData=trainDataO)
      learnerAux <- learnerAux$learn(learnerAux)
      testDataO <- constructData(as.matrix(Xte[,causeNodesTry]), Xte[,effectNodeTry])
      pred <- learnerAux$predict(learnerAux, data=testDataO)
      rs <- pred$gyh - pred$gy
      
      
      
      pval <- dhsic.test(Xte[, causeNodesTry], rs, method=method)$p.value
      #print(paste("p-value: ", pval))
      
      #timePast <- proc.time()-pm0
      #avgTime <- timePast/count
      #regsLeft <- totRegs - count
      #timeLeft <- regsLeft*avgTime
      
      #print(paste("Estimated time left causal ordering: ", round(timeLeft[3]/60,2), " mins."))
      
      return(pval)
    })
    
    mostEffectNode <- S[which.max(pvals)]
    if(pvals[mostEffectNode] < alpha){
      print("no consitent DAGs")
      return(NULL)
    }
    ord <- c(mostEffectNode, ord)
    S <- setdiff(S, mostEffectNode)
    #print(paste("causal order so far: ", paste(ord, collapse="-> ")))
  }
  ord <- c(S, ord)
  
  print(paste("final causal order: ", paste(ord, collapse="-> ")))
  
  print("exits causalOrderingMooij")
  return(ord)
}

minimalDAG <- function(Xtr, Xte, learner, method, alpha=0.2){
  print("enters minimalDAG")
  d <- ncol(Xtr)
  print("calculating causal ordering")
  
  ord <- causalOrderingMooij(Xtr, Xte, learner, method=method, alpha=0)
  
  if(is.null(ord)){
    print("no consitent DAGs")
    return(NULL)
  }  
  
  totRegs <- d*(d-1)/2
  #count <- 0
  
  print("calculating parents")
  # for each node starts with all the parents and removes certain parents if we keep
  # independence of residuals with reduced inputs
  #pm0 <- proc.time()
  parents <- lapply(1:d, function(j){
    node <- ord[j]
    #print(paste("node: ", node))
    if(j  == 1){
      prnts_node <- as.numeric()
      #print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    } else{
      
      prnts_node <- ord[1:(j-1)]
      # we will try and take off as many inputs as we can maintaining independence
      for(k in 1:(j-1)){
        #print(paste("for node ", node, "testing taking off input: ", ord[k]))
        #count <<- count + 1
        notParent_node_try <- ord[k]
        prnts_node_try <- setdiff(prnts_node, notParent_node_try)
        
        if(length(prnts_node_try)==0){
          # if theres only one potential parent then to test if we should eliminate it we
          # check if the response variable "node" is independent of the explanatory variable
          # "prnts_node"
          pval <- dhsic.test(Xte[, prnts_node], Xte[,node])$p.value
        } else{
          
          trainDataO <- constructData(as.matrix(Xtr[,prnts_node_try]), Xtr[,node])
          learnerAux <- setParams(learner, trainDataO)
          learnerAux <- learnerAux$learn(learnerAux)
          testDataO <- constructData(as.matrix(Xte[,prnts_node_try]), Xte[,node])
          pred <- learnerAux$predict(learnerAux, data=testDataO)
          rs <- pred$gyh - pred$gy
          
          pval <- dhsic.test(Xte[,prnts_node_try], rs, method=method)$p.value
          
          
          #timePast <- proc.time()-pm0
          #avgTime <- timePast/count
          #regsLeft <- totRegs - count
          #timeLeft <- regsLeft*avgTime
          
          #print(paste("Estimated time left parents: ", round(timeLeft[3]/60,2), " mins."))
        }
        
        # if we maintain independence of residuals and inputs (p-value big) while not
        # using input k we take it off the parents list
        #print(paste("pval: ", pval))
        if(pval >= alpha){
          #print(paste("for node ", node,"taking off input: ", notParent_node_try))
          prnts_node <- setdiff(prnts_node, notParent_node_try)
        }
      }
      
    }
    print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    return(prnts_node)
  })
  
  # make into a matrix and plot
  V <- ord
  parents2 <- sapply(parents, function(el) (V %in% el)*1)
  
  rownames(parents2) <- V
  colnames(parents2) <- V
  
  parents3 <- getGraph(parents2)
  #plot(parents3)
  
  print("exits minimalDAG")
  return(parents3)
}

dataRegime <- function(data, type=c("recycle", "holdout")){
  dataReg <- switch(type, 
                    recycle={
                      xTrain <- data
                      xTest <- data
                      return(list(train=xTrain, test=xTest))
                    },
                    holdout={
                      midpoint <- ceiling(nrow(data)/2)
                      xTrain <- data[1:midpoint,]
                      xTest <- data[(midpoint+1):nrow(data),]
                      return(list(train=xTrain, test=xTest))
                    })
  return(dataReg)
}


# approx dagset methods
exhaustive <- function(data, pars){
  p <- ncol(data)
  DAGset <- genAllDAGS(p)
  return(DAGset)
}

forceEdge <- function(data, pars){
  p <- ncol(data)
  DAGset <- genAllDAGS(p)
  DAGset <- DAGset[,,-1]
  return(DAGset)
}

pcSkeleton <- function(data, pars){
  V <- colnames(data)
  uDAG <- pcalg:::skeleton(suffStat=list(data=data, ic.method=pars$ic.method), indepTest=eval(parse(text=pars$indepTest)), alpha=pars$alpha, labels=V)
  uDAG <- amat(as.bn(uDAG))
  DAGset <- udag2dags(uDAG)
  return(DAGset)
}

pcMarkovEquiv <- function(data, pars){
  V <- colnames(data)
  pDAG <- pc(suffStat=list(data=data, ic.method=pars$ic.method), indepTest=eval(parse(text=pars$indepTest)), alpha=pars$alpha, labels=V)
  pDAG <- amat(as.bn(pDAG))
  DAGset <- pdag2dags(pDAG)
  return(DAGset)
}

minDAG <- function(data, pars){
  dataReg <- dataRegime(data, type=pars$dataRegime)
  xTrain <- dataReg$train
  xTest <- dataReg$test
  minDAG <- minimalDAG(Xtr=xTrain, Xte=xTest, learner=eval(parse(text=pars$learner)), method=pars$method)
  # obtain generating set of "super-DAGs" (those obtained by adding edges but not taking away)
  # 1. identify non-present edges
  # 2. obtain all super-DAGs
  minDAG <- amat(as.bn(minDAG))
  DAGset <- mindag2dags(minDAG)
  return(DAGset)
}



########################################################################################
# Functions for CauseEffectPairs testing
########################################################################################
genLearnRecipes <- function(learnerList, approxDagSetMethodList, dataSettingList){
  # obtain N different recipes
  
  # divide model and non-model based learner-complexityFunc pairs
  indxNonModel <- which(learnerList=="vanilla")
  indxModel <- which(learnerList!="vanilla")  
  nonModelBased <- expand.grid(learner=learnerList[indxNonModel], approxDagSetMethod=approxDagSetMethodList, dataSetting=dataSettingList ,stringsAsFactors=F)
  modelBased <- expand.grid(learner=learnerList[indxModel], approxDagSetMethod=approxDagSetMethodList, dataSetting=dataSettingList,  stringsAsFactors=F)
  
  
  recipe <- rbind(modelBased, nonModelBased)
  
  # recipe$approxDagSetMethod <- mapply(FUN=function(el, lrnr, dtaReg){
  #   # el <- recipe$approxDagSetMethod[[1]]
  #   res <- el
  #   if(res$name=="minDAG"){
  #     res$pars$learner <- lrnr
  #     res$pars$dataRegime <- dtaReg
  #   }
  #   return(res)
  # }, recipe$approxDagSetMethod, recipe$learner, recipe$dataSetting, SIMPLIFY=FALSE)
  
  #all(sapply(recipe$approxDagSetMethod, function(el) el$name) == sapply(aux, function(el) el$name))
  
  # funcAux <- function(el){
  #   if(is.null(names(el$pars))){
  #     res <- "bla"
  #   } else{
  #     res <- names(el$pars)
  #   }
  #   res <- paste(res, collapse=".")
  #   return(res)
  # }
  
  #all(sapply(recipe$approxDagSetMethod, funcAux)==sapply(aux, funcAux))
  
  return(recipe)
  
}

setdiffRecipes <- function(recipe, recipeSubract){
  srce <- c(rep(1, nrow(recipe)), rep(2, nrow(recipeSubtract)))
  res <- rbind(recipe, recipeSubtract)
  d1 <- duplicated(res) 
  d2 <- rev(duplicated(res[nrow(res):1,]))
  del <- which(d1 | d2 ) #| srce==2
  res <- res[-del,]
  return(res)
}

# version with aprox dag set for each learner
applyLearnRecipesDeprecated <- function(recipe, cmplxScoreList, dataList){
  
  cmplxScorePack <- eval(parse(text=cmplxScoreList))
  cmplxScoresList <- cmplxScorePack$score
  indxMtch <- match(sapply(cmplxScorePack$score, function(el) el$func), cmplxScorePack$tab$cmplxFunc)
  argType <- cmplxScorePack$tab[indxMtch,"argType"]
  
  
  # for each learner x dataSetting get union of all hypotheses to be considered (save time by not running 
  # same hypothesis for different approxDagSets)
  # apply learner under all data settings, and get all complexity scores for all q data sets
  # form N x q recipe scores
  
  learnerDataSetting <- with(recipe, paste(learner, dataSetting, sep="."))
  approxDagMethod <- recipe$approxDagSetMethod
  #approxDagMethod <- with(recipe, sapply(approxDagSetMethod, function(el) el$name))
  apprxDagMthsPerDataSet <- sapply(unique(learnerDataSetting), function(el) recipe$approxDagSetMethod[which(learnerDataSetting==el)])
  
  # I CAN IMPROVE BY RUNNING APPROX DAG METHOD ONLY ONCE PER DATA SET INSTEAD OF FOR EACH LEARNER
  # AND DATA SETTING GIVEN THAT IM NOT ALLOWING THE PARAMETERS OF THE APPROX METHOD TO DEPEND ON THE GENERAL
  #LEARNER AND DATA SETTING ANYMORE (BEFORE I WAS ALLOWING THE MIN DAG METHOD TO TAKE ITS LEARNER AND DATA
  # SETTING FROM THE GENERAL LEARNER AND DATA SETTING)
  
  firstDag <- 1
  lastDag <- length(dataList$dags)
  
  # for each data set ...
  scores <- mapply(FUN=function(dag, x, nm) {
    # i <- 1; dag <- dataList$dags[[i]]; x <- dataList$xs[[i]]; nm <- dataList$names[i]  ;plot(getGraph(dag))
    print("*****************************************")
    print(paste("data set:  ", nm))
    print(head(x))
    print("true dag")
    print(dag)
    
    # get the learner and data regime for each element of the following mapply loop
    aux <- strsplit(names(apprxDagMthsPerDataSet), split="\\.")
    lrnrs <- sapply(aux, function(el) el[[1]])
    dataRegs <- sapply(aux, function(el) el[[2]])
    
    # ... get union of hypothesis for each learner-dataSetting
    scoreList <- mapply(function(apprxDagMths, lrnr, dataReg){
      # i <- 3; apprxDagMths <- apprxDagMthsPerDataSet[[i]]; apprxDagMths; lrnr <- lrnrs[i]; dataReg <- dataRegs[i]
      print("*****************************************")
      print(paste("learner: ", lrnr, " data-regime: ", dataReg))
      
      print("get union of hypotheses:")
      hypsList <- lapply(apprxDagMths, function(apprxDagMth){
        # (apprxDagMth <- apprxDagMths[2])
        print("************************************")
        print(paste("apprx dag set method: ", apprxDagMth))
        #DAGset <- getApproxDAGset(data=x, method=apprxDagMth)
        apprxPack <- eval(parse(text=apprxDagMth))
        
        # for minDAG function, if "learner" and "dataRegime" parameters are NULL, take from recipe
        if(apprxPack$func=="minDAG"){
          if(is.null(apprxPack$pars$learner)) apprxPack$pars$learner <- lrnr
          if(is.null(apprxPack$pars$dataRegime)) apprxPack$pars$dataRegime <- dataReg
        }
        
        pars <- list(data=x, pars=apprxPack$pars)
        DAGset <- do.call(apprxPack$func, pars)
        
        print(paste("number of dags in set: ", dim(DAGset)[3]))
        return(DAGset)
      })
      hypsListID <- lapply(hypsList, getHypID)
      hypArray <- do.call("abind", hypsList)
      # apply union to eliminate duplicates
      hypArray <- unique(hypArray, MARGIN=3)
      
      print(paste("number of hypotheses to be scored: ", dim(hypArray)[3]))
      
      # apply learner-data setting to each hypothesis to get residuals
      dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
      uniqueRegsList <- getUniqueRegsList(dags=hypArray)
      # get data
      dataRg <- dataRegime(x, type=dataReg)
      xTrain <- dataRg$train
      xTest <- dataRg$test 
      # get learner
      learner <- eval(parse(text=lrnr))
      
      # learn
      print(paste("fitting regressions needed to evaluate all hypotheses"))
      modsClass <- fitSEMSetGivenDAGSet(uniqueRegsList, trainData=xTrain, learner=learner)
      # predict
      print(paste("obtaining residuals needed to evaluate all hypotheses"))
      predsClass <- predictSEMSet(uniqueRegsList, data=xTest, learnerList=modsClass, plot=FALSE)
      
      
      
      
      if(lrnr=="vanilla"){
        cmplxScores <- cmplxScoresList[which(argType=="variables")]
      } else{
        cmplxScores <- cmplxScoresList[which(argType!="variables")]
      }
      print(paste("scoring all hypotheses"))
      lapply(predsClass, dimnames)
      predsSet <- lapply(predsClass, function(el) adrop(el[,"resid", ,drop=FALSE], drop=2))
      scores <- complexityScoreList(dags=hypArray, uniqueRegsList=uniqueRegsList, vars=predsSet, cmplxScores=cmplxScores, prnt=FALSE)
      
      # arrange scores back according to the approx dag set method, put max (worse) score on dags outside the approx set
      # also add rank and prob
      
      indxPerApprxDagMth <- lapply(hypsListID, function(ids) match(ids$id, as.numeric(dimnames(hypArray)$dag)))
      scoreList <- lapply(indxPerApprxDagMth, function(indx) scores[indx,,drop=FALSE])
      names(scoreList) <- paste(nm, lrnr, dataReg, apprxDagMths, sep=".")
      return(scoreList)
    }, 
    apprxDagMths=apprxDagMthsPerDataSet, lrnr=lrnrs, dataReg=dataRegs)
    
    
    
    # unfold lists so that we have one element per row in the recipe
    scoreList <- unlist(scoreList, recursive=FALSE) 
    
    
    names(scoreList) <- paste(rep(lrnrs, sapply(apprxDagMthsPerDataSet, length)), unlist(apprxDagMthsPerDataSet), rep(dataRegs, sapply(apprxDagMthsPerDataSet, length)), sep=".")
    
    
    return(scoreList)
    
  }, 
  dag=dataList$dags[firstDag:lastDag], x=dataList$xs[firstDag:lastDag], nm=dataList$names[firstDag:lastDag], SIMPLIFY=FALSE)
  
  names(scores) <- dataList$names[firstDag:lastDag]
  
  return(scores)
  
}


applyLearnRecipes <- function(recipe, cmplxScoreList, dataList){
  
  cmplxScorePack <- eval(parse(text=cmplxScoreList))
  cmplxScoresList <- cmplxScorePack$score
  indxMtch <- match(sapply(cmplxScorePack$score, function(el) el$func), cmplxScorePack$tab$cmplxFunc)
  argType <- cmplxScorePack$tab[indxMtch,"argType"]
  


  approxDagMethod <- recipe$approxDagSetMethod
  
  dataSetting <- recipe$dataSetting
  
  #  we will run each approx dag set method only once
  # (save time by only running once and not once for each learner x data setting since approx dag set method independent of learner and data setting-
  # it has its own data setting parameter)
  
  apprxDagMthsUni <-  unique(recipe$approxDagSetMethod)

  
    
  # then we need to join the hypothesis sets according to learner data setting methods so as to not score the same hypothesis twice
  learnerDataSetting <- with(recipe, paste(learner, dataSetting, sep="."))
  # here we dont take "unique" coz that wd mean that there are identical rows in recipe... we shd check that there isnt when 
  # constructing the recipe
  apprxDagMthsPerDataSet <- sapply(unique(learnerDataSetting), function(el) recipe$approxDagSetMethod[which(learnerDataSetting==el)])
  
  # get the learner and data regime for each element of the following mapply loop
  aux <- strsplit(names(apprxDagMthsPerDataSet), split="\\.")
  lrnrs <- sapply(aux, function(el) el[[1]])
  dataRegs <- sapply(aux, function(el) el[[2]])
  
  
  # from which indices i in apprxDagMthsUni[i] can we find each element in apprxDagMthsPerDataSet
  
  indxsI <- sapply(apprxDagMthsPerDataSet, function(apprxDagMth) match(apprxDagMth, apprxDagMthsUni))
  
  
  firstDag <- 1
  lastDag <- length(dataList$dags)
  
  # for each data set ...
  scores <- mapply(FUN=function(dag, x, nm) {
    # i <- 7; dag <- dataList$dags[[i]]; x <- dataList$xs[[i]]; nm <- dataList$names[i]  ;plot(getGraph(dag))
    print("*****************************************")
    print(paste("data set:  ", nm))
    print(head(x))
    print("true dag")
    print(dag)
    
    
    # get the approximate hypothesis set for each approx method x data setting : the approx method, with the specific data used
    # must be independent of the learner used. This also prevents approx methods with certain randomnmess (ex permutation based hsic used)
    # being run more than once (for diff learners) and resulting in different hypothesis sets despite having exact same parameters.
    
    print("run approximate dag set methods for each datasetting")
    
    hypsListPerApprxDagMth <- lapply(apprxDagMthsUni, function(apprxDagMth){
      # (apprxDagMth <- apprxDagMthsUni[3])
      print("************************************")
      print(paste("apprx dag set method: ", apprxDagMth))
      
      apprxPack <- eval(parse(text=apprxDagMth))
      
      pars <- list(data=x, pars=apprxPack$pars)
      
      DAGset <- do.call(apprxPack$func, pars)
      print(paste("number of dags in set: ", dim(DAGset)[3]))
      return(DAGset)
    })
    
    # ... get union of hypothesis for each learner-dataSetting
    print("get union of hypotheses:")
    hypListList <- mapply(FUN=function(indxI, AMnm){
      # i <- 1; indxI <- indxsI[[i]]
      hypsList <- hypsListPerApprxDagMth[indxI]
   
      hypsListID <- lapply(hypsList, getHypID)
      hypArray <- do.call("abind", hypsList)
      # apply union to eliminate duplicates
      hypArray <- unique(hypArray, MARGIN=3)
      dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
      print(paste("number of hypotheses to be scored for learner-dataSetting ", AMnm, " is ",dim(hypArray)[3]))
      return(list(hypsListID=hypsListID, hypArray=hypArray))
    }, indxI=indxsI, AMnm=names(indxsI), SIMPLIFY=FALSE)
    
  
    # get scores for each hypothesis list
    scoreList <- mapply(function(hypList, apprxDagMths, lrnr, dataReg){
      # i <- 1; hypList <- hypListList[[i]]; lrnr <- lrnrs[i]; dataReg <- dataRegs[i]
      print("*****************************************")
      print(paste("learner: ", lrnr, " data-regime: ", dataReg))
      

      hypArray <- hypList$hypArray
      hypsListID <- hypList$hypsListID
      
      # apply learner-data setting to each hypothesis to get residuals
      uniqueRegsList <- getUniqueRegsList(dags=hypArray)
      
      # get data
      dataRg <- dataRegime(x, type=dataReg)
      xTrain <- dataRg$train
      xTest <- dataRg$test 
      # get learner
      learner <- eval(parse(text=lrnr))
      
      # learn
      print(paste("fitting regressions needed to evaluate all hypotheses"))
      modsClass <- fitSEMSetGivenDAGSet(uniqueRegsList, trainData=xTrain, learner=learner)
      # predict
      print(paste("obtaining residuals needed to evaluate all hypotheses"))
      predsClass <- predictSEMSet(uniqueRegsList, data=xTest, learnerList=modsClass, plot=FALSE)
      
      # ESTA PATE EVENTUALMENTE TENDRA QUE CAMBIAR PARA CAUSAL LEARNERS NON-MODEL BASED COMO KERNEL DEVIANCE
      cmplxScores <- cmplxScoresList[which(argType=="residuals")]
      
      print(paste("scoring all hypotheses"))
      lapply(predsClass, dimnames)
      predsSet <- lapply(predsClass, function(el) adrop(el[,"resid", ,drop=FALSE], drop=2))
      scores <- complexityScoreList(dags=hypArray, uniqueRegsList=uniqueRegsList, vars=predsSet, cmplxScores=cmplxScores, prnt=FALSE)
      
      # arrange scores back according to the approx dag set method
      # also add rank and prob
      
      indxPerApprxDagMth <- lapply(hypsListID, function(ids) match(ids$id, as.numeric(dimnames(hypArray)$dag)))
      scoreList <- lapply(indxPerApprxDagMth, function(indx) scores[indx,,drop=FALSE])
      names(scoreList) <- paste(nm, lrnr, dataReg, apprxDagMths, sep=".")
      return(scoreList)
    }, 
    hypList=hypListList, apprxDagMths=apprxDagMthsPerDataSet, lrnr=lrnrs, dataReg=dataRegs, SIMPLIFY=FALSE)
    
    
    
    # unfold lists so that we have one element per row in the recipe
    scoreList <- unlist(scoreList, recursive=FALSE) 
    
    
    names(scoreList) <- paste(rep(lrnrs, sapply(apprxDagMthsPerDataSet, length)), unlist(apprxDagMthsPerDataSet), rep(dataRegs, sapply(apprxDagMthsPerDataSet, length)), sep=".")
    
    # scoreList[["krr1.minDAG_gptk.recycle"]]; scoreList[["gptk.minDAG_gptk.recycle"]] 
    
    return(scoreList)
    
  }, 
  dag=dataList$dags[firstDag:lastDag], x=dataList$xs[firstDag:lastDag], nm=dataList$names[firstDag:lastDag], SIMPLIFY=FALSE)
  
  names(scores) <- dataList$names[firstDag:lastDag]
  
  return(scores)
  
}

# to apply a function at the data-recipe level you need, the scores, the function to apply,
# you might also need:
# the true dag of each dag, the complexity pack used in order to transform scores, 
# a function at the dataset level to get info about the data or dag,
# a function ath the dataset_recipe level to get info about the scores (which hypotheses, which scores, etc)
deliverFuncToScore <- function(scores, func, dags=rep(NA, length(scores)), cmplxPack=NULL, 
                               getDagPars=function(dag) return(list()), getDataLevelPars=function(scoreMat, dag) return(list()), ...){
  
  pars <- list(...)
  # pars <- list(groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc=aggFunc)
  # pars <- list(rank2Func="probRD")
  #print("length(pars)")
  #print(length(pars))
  #print("names(pars)")
  #print(names(pars))
  
  pars <- c(pars, cmplxPack=cmplxPack)
  # pars <- list()
  
  res <- mapply(FUN=function(data_score, dag){
    # i <- 2; data_score <- scores[[i]]; dag <- dags[[i]]
    #print("***********")
    #print("true Dag")
    #print(dag)
    parsDag <- getDagPars(dag) 
    res <- lapply(data_score, function(data_score_recipe){
      # data_score_recipe <- data_score[[1]]
      
      parsDataScoreRecipe <- getDataLevelPars(data_score_recipe, dag)
      pars2 <- c(pars, parsDag, parsDataScoreRecipe)
      pars2$scoreMat <- data_score_recipe
      
      names(pars2)
      
      # cmplxPack <- pars2$cmplxPack; scoreMat <- pars2$scoreMat  
      # rank2Func <- pars2$rank2Func; scoreMat <- pars2$scoreMat  
      # groupFuncs <- pars2$groupFuncs; groupFuncPars <- pars2$groupFuncPars; aggFunc <- pars2$aggFunc; trueDag <- pars2$trueDag; dagHyps <- pars2$dagHyps; scoreMat <- pars2$scoreMat
      res <- do.call(func, pars2)
      
      return(res)
    })
    return(res)
  }, data_score=scores, dag=dags, SIMPLIFY=FALSE)
  return(res)
}

# usually the getDataLevel Pars will be "getHypots" which obtains hypotheses in id form which have been scored 
# transforms to matrix form
getHypots <- function(scoreMat, dag){
  p <- nrow(dag)
  ids <- as.numeric(rownames(scoreMat))
  dagHyps <- sapply(ids, function(id) dagIDtoMatrix(id,p), simplify="array")
  return(list(dagHyps=dagHyps))
}


# functions corresponding to "func" argument of deliverFuncToScore function
#a print infinite, not-a-numbers or nas
printErrors <- function(scoreMat){
  if(any(is.infinite(scoreMat) | is.nan(scoreMat) | is.na(scoreMat))){
    print(paste("THERE ARE ", sum(is.infinite(scoreMat)) + sum(is.nan(scoreMat)) + sum(is.na(scoreMat))," NON VALID SCORES:"))
    print(scoreMat)
  } 
  return(NULL)
}
#b correct infinite, not-a-numbers or nas
correctErrors <- function(scoreMat, funcReplace){
  
  scoreMat2 <- apply(scoreMat, 2, function(col){
    indxError <- which(is.infinite(col) | is.nan(col) | is.na(col))
    if(length(indxError)>0){
      print("indxError")
      print(indxError)
      col[indxError] <- do.call(funcReplace, list(col[-indxError]))
    }
    return(col)
  })
  dim(scoreMat2) <- dim(scoreMat)
  dimnames(scoreMat2) <- dimnames(scoreMat)
  #rownames(scoreMat2) <- rownames(scoreMat)
  #colnames(scoreMat2) <- colnames(scoreMat)
  return(scoreMat2)
}
#c normalize scores to prob
correctToProb <- function(scoreMat, dagHyps, trueDag){
  scoreT <- apply(scoreMat, 2, correctScoreToAdd, hyps=dagHyps)
  dim(scoreT) <- dim(scoreMat)
  scoreT <- apply(scoreT, 2, scoreToProb)
  dim(scoreT) <- dim(scoreMat)
  dimnames(scoreT) <- dimnames(scoreMat)
  return(scoreT)
}
# se also aggregateScores and rankedDecisions below

# unwrap scores list structure (data list, then recipe list, each with score matrix) into a database long format
unwrapScoreFunc <- function(scoresFunc, dataList, ws){
  innerDims <- names(dimnames(scoresFunc[[1]][[1]]))
  scoresFuncs <- mapply(FUN=function(scrFuncs, dta, w, nms){
    # i <- 1; scrFuncs <- scoresFunc[[i]]; nms <- names(scoresFunc)[i]; dta <- dataList$xs[[i]]
    res <- mapply(FUN=function(scrFunc,  nm){
      # j <- 1; scrFunc <- scrFuncs[[j]]; nm <- names(scrFunc)[j]
      res <- melt(scrFunc)
      res$recipe <- nm
      return(res)
    }, scrFunc=scrFuncs, nm=names(scrFuncs), SIMPLIFY=FALSE)
    
    res <- do.call(rbind, res)
    res$dataset <- nms
    res$p <- ncol(dta)
    res$n <- nrow(dta)
    res$w <- w
    res$weighted_value <- res$w*res$value
    
    colnames(res) <- c(innerDims,"value","recipe","dataset","p","n", "w","weighted_value")
    res <- res[,c("dataset", "p", "n", "w", "recipe", rev(innerDims),"value", "weighted_value")]
    return(res)
  }, scrFuncs=scoresFunc, dta=dataList$xs, w=ws, nms=names(scoresFunc), SIMPLIFY=FALSE)
  
  scoresFuncs <- do.call(rbind, scoresFuncs)
  rownames(scoresFuncs) <- 1:nrow(scoresFuncs)
  
  scoresFuncs$recipeFull <- apply(scoresFuncs[,c("recipe",innerDims[2:length(innerDims)])], 1, paste, collapse=".")
  aux <- strsplit(scoresFuncs$recipe, "\\.")
  scoresFuncs$learner <- sapply(aux, function(el) el[1])
  scoresFuncs$approxMethod <- sapply(aux, function(el) el[2])
  scoresFuncs$dataReg <- sapply(aux, function(el) el[3])
  scoresFuncs <- scoresFuncs[,c("dataset", "p", "n", "w","recipeFull", "recipe", "learner", "approxMethod", "dataReg", rev(innerDims),"value", "weighted_value")]
  return(scoresFuncs)
}

# uwrap scores with lis structure and with same number of hypothesis into an array format
scoreDBToArray <- function(scoreDB, names_datasets, value=c("value", "weighted_value")){
  
  value <- match.arg(value)
  
  scoreArr <- sapply(unique(scoreDB$recipeFull), function(reci){
    # reci <- unique(scoreDB$recipeFull)[212]
    indxReci <- which(scoreDB$recipeFull==reci)
    mat <- scoreDB[indxReci,]
    mat <- cast(mat, dataset~hypothesis, value="value")
    
    if(nrow(mat) < length(names_datasets)){
      print(paste("for recipe ", reci, " there are only scores for", nrow(mat), " datasets" ))
    }
    
    hyps <- colnames(mat)[c(2,3)]
    # meta learners are not available for all datasets
    mat <- mat[match(as.character(names_datasets), mat$dataset),]
    mat <- as.matrix(mat[, c("hyp0pos","hyp1pos")])
    dimnames(mat) <- list(dataset=names_datasets, hypothesis=hyps)
    return(mat)
  }, simplify="array")
  
  return(scoreArr)
}


# This combines recipes of any kind as long as they have some hypotheses in common, so you can combine recipes
# which use different approximation method, learner or data setting which might be quite weird. Also because
# any recipe can be combined as long as they have a hypothesis set scored which has some intersection 
# (actually one must be a subset of the other as is currently implemented) then you end up combining different
# recipes for different data sets which is a hassle as you have to deal with NAs for some datasets 
addMetaLearner2 <- function(scores, func, ...){
  pars <- list(...)
  # pars <- list()
  
  newScores <- mapply(FUN=function(data_score, nm){
    # i <- 24; data_score <- scores[[i]]; nm <- names(scores)[i]
    
    # order the list of recipe-scores in a data-recipes-scores according to the number of hypotheses
    data_score_aux <- data_score[order(sapply(data_score, function(el) length(dimnames(el)$dag)))]
    
    print("*****************************************")
    print(paste("data set:  ", nm))
    
    # which recipes to merge?
    # we can only merge recipes that have common hypotheses.
    # we try and form combinations of 2, 3,..., n number of recipes scored 
    recipeMerge <- lapply(2:length(data_score_aux), function(m){
      # m <- 2
      #print("**************")
      #print(paste("m: ", m))
      
      # obtener los indices (cuales columnas de la matriz de scores) de las m-combinaciones posibles
      combos <- combn(1:length(data_score_aux), m) 
      
      # cuales de esas combinaciones de recipes si tienen hipotesis en común
      # es decir, cuales si podemos combinar?
      subsetBool <- apply(combos,2, function(indx){
        # indx <- combos[,1]
        
        # obtenemos una lista con listas de las hipotesis (dags) q corresponden a cada una de las recipes
        # que estamos considerando combinar
        dagList <- lapply(indx, function(i) dimnames(data_score_aux[[i]])$dag)  
        dagList <- dagList[order(sapply(dagList, length))]
        
        # we will only merge two recipes if one is a subset of the other. 
        # since we ordered the recipes by the number of hypothesis from low to high, and
        # the combo indices are from low to high, we need only check that the hypotheses of a recipe
        # given by a given index is a subset of the hypotheses of the recipe given by a subsequent recipe in
        # the m-combination. 
        
        # To compare "subsequent" sets of hypothesis we take the dagList which is of size m=2,3,..., length(recipe)
        # and put elements 1,...,m-1 and 2,...,m so that by comparing element by element we are comparing sets of 
        # hypotheses which are subsequent accoridng to indx
        dagList1 <- dagList[1:(length(dagList)-1)]
        dagList2 <- dagList[2:(length(dagList))]
        
        # now we check that all earlier elements (elements in dagList1) are in later elements (elements in dagList2)
        res <- mapply(FUN=function(dag1, dag2){
          # i <- 1; dag1 <- dagList1[[i]]; dag2 <- dagList2[[i]]
          all(dag1 %in% dag2)
        }, dag1=dagList1, dag2=dagList2)    
        res <- all(res)
        return(res)
      })
      
      # we take the indices of only those combinations which we found that all smaller size sets are subsets of larger sets
      mergableCombos <- combos[,which(subsetBool), drop=FALSE]
      # get the recipe names to merge
      mergable <- apply(mergableCombos,2, function(indx) names(data_score_aux)[indx])
      mergable <- as.list(as.data.frame(mergable))
      
      # for each combination of recipes which was found to be "mergeable" (smaller size hypothesis sets are a subset of larger size ones)
      # we take the smaller (or equal) set of each mergeable Combo
      dagList <- lapply(mergableCombos[1,], function(i) dimnames(data_score_aux[[i]])$dag)
      
      # if there are 2 or more m-combinatins with the same "common" hypothesis we want to merge all the corresponding 
      # recipes so we take the union of the hypothesis set lists to merge
      dagUniList <- unique(dagList)
      # ok now we know what sets of hypothesis will result from merging certain recipes, but which ones should we merge?
      # for each unique hypothesis set we look for the corresponding recipes to merge, we first identify for each mergeable
      # of the m-combinations which is the unique hypothesis set which will result from merging. Now for all that
      # have the same we form a union of recipes to merge
      indxUnion <- match(dagList, dagUniList)
      recipeMerge <- lapply(unique(indxUnion), function(i){
        # i <- 1
        indx <- which(indxUnion==i)
        unique(unlist(mergable[indx]))
      })
      
      # return the recipes that will be merged and the resulting set of hypotheses which will result 
      return(list(recipe=recipeMerge, dags=dagUniList))
    })
    
    # recipeMerge is a list of lists. For each m indicating the size of the combinations there is a list
    # indicating which recipes can be merged and the corresponding resulting hypothesis set, but 
    # accross different m, there may be repeated pairs of (recipes to merge x resulting hypothesis set)
    # so we need to elimnate them
    
    # we first put the recipes and resulting hypothesis sets of different m's on the same level 
    recipeList <- lapply(recipeMerge, function(el) el$recipe)
    recipeList <- unlist(recipeList, recursive=FALSE)
    dagList <- lapply(recipeMerge, function(el) el$dags)   
    dagList <- unlist(dagList, recursive = FALSE)
    
    
    
    #indxOrder <- order(sapply(dagList, length))
    #dagList <- dagList[indxOrder]
    #recipeList <- recipeList[indxOrder]
    
    # we take the union of resulting hypothesis sets and corresponding union of recipes to merge
    dagUniList <- unique(dagList)
    indxUnion <- match(dagList, dagUniList)
    recipeMerge <- lapply(unique(indxUnion), function(i){
      # i <- 1
      indx <- which(indxUnion==i)
      unique(unlist(recipeList[indx]))
    })
    
    # we then form a list of extra-scores, one for each list of recipes to merge and its corresponding resulting
    # hypothesis set
    extra_score <- mapply(FUN=function(dagList, mixRecipe){
      # i <- 2; dagList <- dagUniList[[i]]; mixRecipe <- recipeMerge[[i]]
      print("*****************************************")
      print(paste("mixing recipes:  ", paste(mixRecipe, collapse=" & ")))
      print(paste("these apply to hypotheses-dags:", paste(dagList, collapse=", ")))
      
      # we get the indices of the columns to mix
      indxMix <- which(names(data_score_aux) %in% mixRecipe)
      # we form an array with the recipes to merge and the corresponding hypotheses
      scoreArr <- sapply(data_score_aux[indxMix], function(el) el[dagList,,drop=FALSE], simplify="array")
      names(dimnames(scoreArr))[3] <- "recipe"
      
     
      
      pars2 <- pars
      pars2$scoreArr <- scoreArr
      
      res <- do.call(func, pars2)
      return(res)
    }, dagList=dagUniList, mixRecipe=recipeMerge, SIMPLIFY=FALSE)
    
    # obtain names for the merged recipes
    splitzies <- lapply(recipeMerge, strsplit, split="\\.")
    nmsEx <- sapply(splitzies, function(el){
      # el <- splitzies[[1]]
      lrnrs <- sapply(el, function(el2) el2[1])
      aprxMet <- sapply(el, function(el2) el2[2])
      datReg <- sapply(el, function(el2) el2[3])
      paste(paste(unique(lrnrs), collapse="_"), paste(unique(aprxMet), collapse="_"), paste(unique(datReg), collapse="_") ,sep=".")
    })
    
    names(extra_score) <- nmsEx
    
    newScore <- c(data_score, extra_score)
    return(newScore)
  }, data_score=scores, nm=names(scores))
  
  names(newScores) <- names(scores)
  return(newScores)
  
}

#This one only combines recipes where the approx method (and parameters!) and data setting are the same 
# and only the learner changes: guaranteeing that the same hypotheses will be scored by all learners

addMetaLearner <- function(scores, func, ...){
  pars <- list(...)
  # pars <- list()
  
  
  data_score <- scores[[1]]
  
  recipeNms <- names(data_score)
  aux <- strsplit(recipeNms, split="\\.")
  lrnrs <- sapply(aux, function(el) el[1])
  approxMets <- sapply(aux, function(el) el[2])
  dataSetts <-  sapply(aux, function(el) el[3]) 
  approxMet_dataSetts <- paste(approxMets, dataSetts, sep=".")
  approxMet_dataSetts_tab <- table(approxMet_dataSetts)
  approxMet_dataSetts_tab <- approxMet_dataSetts_tab[which(approxMet_dataSetts_tab>1)]
  indxMergeList <- sapply(names(approxMet_dataSetts_tab), function(el) which(el == approxMet_dataSetts))
  
  print("*****************************************")
  print(paste("mixing ", length(indxMergeList), " recipes per data:"))
  
  aux <- lapply(1:length(indxMergeList), function(i){ 
    print(paste("for approx dag method-datasetting ",names(indxMergeList)[i], " mixing following recipes (", length(indxMergeList[[i]]), ", total)"))
    print(paste(names(data_score)[indxMergeList[[i]]], collapse=" & "))
    print("")
    return(NULL)
    })
  
  
  if(length(approxMet_dataSetts_tab)>0){
    newScores <- mapply(FUN=function(data_score, nm){
    # i <- 7; data_score <- scores[[i]]; nm <- names(scores)[i]
    
    #print("*****************************************")
    #print(paste("data set:  ", nm))
    
    # we then form a list of extra-scores, one for each list of recipes to merge and its corresponding resulting
    # hypothesis set
    extra_score <- lapply(indxMergeList, FUN=function(indxMerge){
      # i <- 1; indxMerge <- indxMergeList[[i]] 
      recipeMerge <- recipeNms[indxMerge]
      
      dagsScored <- sapply(recipeMerge, function(el) rownames(data_score[[el]]))
      
      if(is.list(dagsScored)) stop("error, trying to merge recipes with different size hypothesis sets scored")
      
      dim(dagsScored) <- c(nrow(data_score[[ indxMerge[1] ]]), length(indxMerge))
      
      chkSame <- apply(dagsScored, 1, function(row) all(row==row[1]))
      if(any(!chkSame)) stop("error, trying to merge recipes with different hypothesis set scored")
      #print("*****************************************")
      #print(paste("mixing recipes:  ", paste(recipeMerge, collapse=" & ")))
      #print(paste("these apply to hypotheses-dags:", paste(dagsScored[,1], collapse=", ")))
      
      # we get the indices of the columns to mix
      indxMix <- which(names(data_score) %in% recipeMerge)
      # we form an array with the recipes to merge and the corresponding hypotheses
      scoreArr <- sapply(data_score[indxMix], function(el) el, simplify="array")
      names(dimnames(scoreArr))[3] <- "recipe"
      
      
      
      pars2 <- pars
      pars2$scoreArr <- scoreArr
      
      res <- do.call(func, pars2)
      return(res)
    })
    names(extra_score) <- sapply(indxMergeList, function(el) paste( paste(lrnrs[el], collapse="_"), unique(approxMets[el]), unique(dataSetts[el]), sep="."))
    
    newScore <- c(data_score, extra_score)
    return(newScore)
  
    }, data_score=scores, nm=names(scores), SIMPLIFY=FALSE)
    names(newScores) <- names(scores)
  } else{
    newScores <- scores
  }
  
  return(newScores)
  
}



numDataPairs <- function(pairs, folder){
  num <- sapply(pairs, function(pair){ 
    
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    dat <- read.csv(fileFolder, sep="", header=F)
    n <- nrow(dat)
    return(n)
  })
  return(num)
}

summaryDataPairs <- function(meta, pairs, folder){
  sapply(pairs, function(pair){ 
    indxPair <- which(meta$pairNumber==pair)
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    dat <- read.csv(fileFolder, sep="", header=F)
    print("********************")
    print(paste("pair:", pair))
    print("summary cause:")
    print(summary(dat[,meta[indxPair, "firstCauseCol"]]))
    print("summary effect:")
    print(summary(dat[,meta[indxPair,"firstEffectCol"]]))
    return(NULL)
  })
  return(NULL)
}

createTCEPList <- function(pairs, folder){
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- pairs
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  dags <- lapply(fileFolders, function(el) dag) 
  xs <- fileFolders
  ns <- lapply(fileFolders, function(el) NA)
  res <- list(dags=dags, xs=xs, ns=ns, names=pairs)
  return(res)
}


# Old recipe-style functions


scorePairs <- function(meta, pairs, learner, numFolds, folder, fac=1.25){
  # pairs <- c(1,2,3)
  count <- 0
  numPairs <- length(pairs)
  
  pm0 <- proc.time()
  pvalues <- sapply(pairs, function(pair){ 
    count <<- count + 1
    print(paste("count: ", count))
    # read in a cause-effect pair
    # pair <-  69, 98 ,70, 17
    
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    dat <- read.csv(fileFolder, sep="", header=F)
    n <- nrow(dat)
    # we mix incase they are ordered in some way so as to obtain similar train/test samples
    maxData <- 1000
    set.seed(1)
    dat <- dat[sample(n, size=min(maxData,n)),]
    
    # standardize
    dat <- apply(as.matrix(dat), 2, function(col){
        res <- col-mean(col)
        #print(mean(res))
        res <- res/sd(col)
        #print(sd(res))
        return(res)
    }) 
    
    print("col means and sds")
    print(apply(dat, 2, mean))
    print(apply(dat, 2, sd))
    
    
    n <- nrow(dat)
    print("**************************")
    print(paste("pair #: ", pair))
    print(paste("number of observations: ", n))
    print(head(dat))
    indx <- which(meta$pairNumber==pair)
    smpl <- 1:(min(n,100))
    plot(dat[smpl, meta[indx,"firstCauseCol"]],dat[smpl, meta[indx,"firstEffectCol"]])
    
    midpoint <- ceiling(n/2)
    xTrain <- dat[1:midpoint,]
    xTest <- dat[(midpoint+1):n,]
    
    
    
    # estimate V1 -> V2 model and obtain p-value using gamma approx
    pvals <- sapply(as.numeric(meta[which(meta$pairNumber==pair), c("firstCauseCol", "firstEffectCol")]), function(indx){
      # indX <- 1; indY <- 2
      # indX <- 2; indY <- 1
      indX <- indx
      indY <- setdiff(c(1,2), indx)
      rs <- Residuals(Xtr=xTrain[,indX], Ytr=xTrain[,indY], Xte=xTest[,indX], Yte=xTest[,indY], learner, numFolds, fac=fac)
      
      # under Null hypothesis of independence (i.e for true causal direction) p-value of independence test
      # should be distributed uniformly. Lets: 
      #   a) divide test residuals in 10
      #   b) apply hsic independence test to each subset
      #   c) apply a k.s. test for uniformity to those p-values
      numTest <- length(rs$rsTe)
      
      numParts <- 10
      smpl <- sample(1:numParts, size=numTest, replace=TRUE)
      pvalsPart <- sapply(1:numParts, function(part) {
        indxPart <- which(smpl==part)
        res <- dhsic.test(xTest[indxPart,indX], rs$rsTe[indxPart], method="gamma")$p.value
      })
      
      ks_test <- ks.test(pvalsPart, punif)
      
      statUnif <- ks_test$statistic
      pvalUnif <- ks_test$p.value
      
      hsic_test <- dhsic.test(xTest[,indX], rs$rsTe, method="gamma")
      
      statIndep <- hsic_test$statistic
      pvalIndep <- hsic_test$p.value
      
      #bootstrap
      pvalsBoot <- sapply(1:100, function(i) {
        indxPart <- sample(1:length(rs$rsTe), replace=T)
        res <- dhsic.test(xTest[indxPart,indX], rs$rsTe[indxPart], method="gamma")$p.value
      })
      # hist(pvalsBoot)
      
      bootUnifTest <- ks.test(pvalsBoot, punif)
      statUnifBoot <- bootUnifTest$statistic
      pvalUnifBoot <- bootUnifTest$p.value
      
      pval <- dhsic.test(xTest[,indX], rs$rsTe, method="gamma")$p.value
      
      res <- c(indepStat=statIndep, unifStat=statUnif, indepPval=pvalIndep, unifPval=pvalUnif, unifBootPval=pvalUnifBoot, unifBootStat=statUnifBoot)
      return(res)
    })
    
    timePerPair <- ((proc.time()-pm0)[3])/count
    timeTillComp <- timePerPair*(numPairs-count)
    print(paste("estimated time till completion: ", round(timeTillComp/60,2), " mins"))
    
    return(pvals)
  }, simplify="array")
  
  
  dimnames(pvalues) <- list(type=c("indepStat","unifStat","indepPval","unifPval","unifBootPval","unifBootStat"), direction=c("cause", "effect"), pair=pairs)
  
  return(pvalues)
}

scoreSimPairs <- function(X, learner, numFolds, fac=1.25){
  count <- 0
  numPairs <- length(dimnames(X)$ps)*length(dimnames(X)$funcs)
  
  pm0 <- proc.time()
  # for all pairs associated to each function f and each p 
  pvalues <- sapply(dimnames(X)$ps, function(p) sapply(dimnames(X)$funcs, function(f){ 
    # p <- "1"; f <- "1"
    count <<- count + 1
    print(paste("count: ", count))
    print(paste("(p,f): ", paste(p,f,sep="-")))
    # read in a cause-effect pair
    
    dat <- rbind(X[,f,p,"train",], X[,f,p,"test",])
    
    # standardize
    dat <- apply(as.matrix(dat), 2, function(col){
      res <- col-mean(col)
      #print(mean(res))
      res <- res/sd(col)
      #print(sd(res))
      return(res)
    }) 
    
    print("col means and sds")
    print(apply(dat, 2, mean))
    print(apply(dat, 2, sd))
    
    
    n <- nrow(dat)
    print("**************************")
    print(paste("pair #: ", count))
    print(paste("number of observations: ", n))
    print(head(dat))
    
    smpl <- 1:(min(n,100))
    plot(dat[smpl, "x"],dat[smpl, "y"])
    
    midpoint <- ceiling(n/2)
    xTrain <- dat[1:midpoint,]
    xTest <- dat[(midpoint+1):n,]
    
    
    # estimate V1 -> V2 model and obtain p-value using gamma approx
    pvals <- sapply(1:2, function(indx){
      # indX <- 1; indY <- 2
      # indX <- 2; indY <- 1
      # indx <- 2
      indX <- indx
      indY <- setdiff(c(1,2), indx)
      rs <- Residuals(Xtr=xTrain[,indX], Ytr=xTrain[,indY], Xte=xTest[,indX], Yte=xTest[,indY], learner, numFolds, fac=fac)
      
      # plot(xTest[,indX], rs$Te)
      
      # under Null hypothesis of independence (i.e for true causal direction) p-value of independence test
      # should be distributed uniformly. Lets: 
      #   a) divide test residuals in 10
      #   b) apply hsic independence test to each subset
      #   c) apply a k.s. test for uniformity to those p-values
      numTest <- length(rs$rsTe)
      numParts <- 10
      smpl <- sample(1:numParts, size=numTest, replace=TRUE)
      pvalsPart <- sapply(1:numParts, function(part) {
        indxPart <- which(smpl==part)
        res <- dhsic.test(xTest[indxPart,indX], rs$rsTe[indxPart], method="gamma")$p.value
      })
      #hist(pvalsPart)
      
      ks_test <- ks.test(pvalsPart, punif)
      
      statUnif <- ks_test$statistic
      pvalUnif <- ks_test$p.value
      
      #bootstrap
      pvalsBoot <- sapply(1:100, function(i) {
        indxPart <- sample(1:length(rs$rsTe), replace=T)
        res <- dhsic.test(xTest[indxPart,indX], rs$rsTe[indxPart], method="gamma")$p.value
      })
      # hist(pvalsBoot)
      
      bootUnifTest <- ks.test(pvalsBoot, punif)
      statUnifBoot <- bootUnifTest$statistic
      pvalUnifBoot <- bootUnifTest$p.value
      
      pval <- dhsic.test(xTest[,indX], rs$rsTe, method="gamma")$p.value
      
      res <- c(indepStat=statIndep, unifStat=statUnif, indepPval=pvalIndep, unifPval=pvalUnif, unifBootPval=pvalUnifBoot, unifBootStat=statUnifBoot)
      return(res)
    })
    
    timePerPair <- ((proc.time()-pm0)[3])/count
    timeTillComp <- timePerPair*(numPairs-count)
    print(paste("estimated time till completion: ", round(timeTillComp/60,2), " mins"))
    
    return(pvals)
  }, simplify="array"), simplify="array")
  
  
  
  dimnames(pvalues) <- list(type=c("indepStat","unifStat","indepPval","unifPval","unifBootPval","unifBootStat"), direction=c("cause", "effect"), funcs=dimnames(X)$funcs, ps=dimnames(X)$ps)
  
  return(pvalues)
}


########################################################################################
# Performance functions for causal learning methods
########################################################################################
# overview

# Performance evaluation

# combine accross learners (but asame approx method and data setting to keep same hypothesis space) or within learners (combine different scores
# for same learner-aprroxDagset-recycle)

measuresDagDistEdge <- function(scoreMat, dagHyps, trueDag, cmplxPack){
  
  p <- dim(trueDag)[1]
  offDiag <- !diag(p)
  
  # Transform scores to weights according to complexity pack chosen
  tab <- eval(parse(text=cmplxPack))$tab
  entropyFuncs <- eval(parse(text=cmplxPack))$entropyFuncs
  cmplxFunc <- colnames(scoreMat)
  entropIndx <- which(cmplxFunc %in% entropyFuncs)
  cmplxFunc[entropIndx] <- "score_sumMarginalEntropies"
  indx <- match(cmplxFunc, tab$cmplxFunc)
  rankFunc <- tab$rankFunc[indx]
  probFunc <- tab$probFunc[indx]
  ws <- sapply(1:ncol(scoreMat), function(i){
    res <- do.call(rankFunc[i], list(scoreMat[,i], dagHyps))
    res <- do.call(probFunc[i], list(res, dagHyps))
  }, simplify="matrix")
  dim(ws) <- dim(scoreMat)
  
  
  # Best dag per score
  indxWinner <- apply(scoreMat, 2, which.min)
  chosenDag <- dagHyps[,,indxWinner, drop=FALSE]
  
  
  # edge distances best score
  edgeDists <- apply(chosenDag, 3, edgeDist, dagTrue=trueDag)/choose(p,2)
  nonEdgeDists <- apply(chosenDag, 3, nonEdgeDist, dagTrue=trueDag)/choose(p,2)
  totalEdgeDists <- (edgeDists + nonEdgeDists)/2
  #apply(chosenDag, 3, totalEdgeDist, dagTrue=trueDag)/(2*choose(p,2)) 
  
  #weighted edgeDistances
  edgeDistsW <- apply(dagHyps, 3, edgeDist, dagTrue=trueDag)/choose(p,2)
  nonEdgeDistsW <- apply(dagHyps, 3, nonEdgeDist, dagTrue=trueDag)/choose(p,2)
  #totalEdgeDistsW <- apply(dagHyps, 3, totalEdgeDist, dagTrue=trueDag)/(2*choose(p,2))
  edgeDistsW <- apply(ws, 2, function(col) sum(col*edgeDistsW))
  nonEdgeDistsW <- apply(ws, 2, function(col) sum(col*nonEdgeDistsW))
  #totalEdgeDistsW <- apply(ws, 2, function(col) sum(col*totalEdgeDistsW))
  totalEdgeDistsW <- (edgeDistsW + nonEdgeDistsW)/2
  
  resEdgeDists <- rbind(edgeDists, nonEdgeDists, totalEdgeDists, edgeDistsW, nonEdgeDistsW, totalEdgeDistsW)
  rownames(resEdgeDists) <- c("edgeD", "nonEdgeD", "totEdgeD","WedgeD", "WnonEdgeD", "WtotEdgeD")
  colnames(resEdgeDists) <- colnames(scoreMat)
  
  # edgewise contingency table based measures
  conTabs <- apply(chosenDag, 3, function(pred) contingencyTable(pred[offDiag], trueDag[offDiag]))
  ccrs  <- apply(conTabs, 2, correctCR)
  msrs  <- apply(conTabs, 2, misCR)
  ppps  <- apply(conTabs, 2, posPP) 
  npps  <- apply(conTabs, 2, negPP)
  senss <- apply(conTabs, 2, sensitivity)
  specs <- apply(conTabs, 2, specificity) 
  fprs  <- apply(conTabs, 2, fpr)
  fnrs  <- apply(conTabs, 2, fnr)
  tsss  <- apply(conTabs, 2, tss) 
  
  # weighted edgewise contingency table based measures
  conTabs <- apply(dagHyps, 3, function(pred) contingencyTable(pred[offDiag], trueDag[offDiag]))
  ccrsW  <- apply(conTabs, 2, correctCR)
  msrsW  <- apply(conTabs, 2, misCR)
  pppsW  <- apply(conTabs, 2, posPP) 
  nppsW  <- apply(conTabs, 2, negPP)
  senssW <- apply(conTabs, 2, sensitivity)
  specsW <- apply(conTabs, 2, specificity) 
  fprsW  <- apply(conTabs, 2, fpr)
  fnrsW  <- apply(conTabs, 2, fnr)
  tsssW  <- apply(conTabs, 2, tss) 
  ccrsW  <- apply(ws, 2, function(col) sum(col*ccrsW))
  msrsW  <- apply(ws, 2, function(col) sum(col*msrsW))
  pppsW  <- apply(ws, 2, function(col) sum(col*pppsW)) 
  nppsW  <- apply(ws, 2, function(col) sum(col*nppsW))
  senssW <- apply(ws, 2, function(col) sum(col*senssW))
  specsW <- apply(ws, 2, function(col) sum(col*specsW)) 
  fprsW  <- apply(ws, 2, function(col) sum(col*fprsW))
  fnrsW  <- apply(ws, 2, function(col) sum(col*fnrsW))
  tsssW  <- apply(ws, 2, function(col) sum(col*tsssW))
  
  #AUC
  preds <- apply(dagHyps, 3, function(mat) as.numeric(mat[offDiag])) %*% ws
  
  obs <- trueDag[offDiag]
  numTrials <- length(obs)
  numCasesGoal <- numTrials/2
  numCasesNow <- sum(obs)
  numCasesAdd <- numCasesGoal - numCasesNow
  indx0 <- which(obs==0)
  indx1 <- which(obs==1)
  smpl <- as.numeric()
  
  #print("numCasesAdd")
  #print(numCasesAdd)
  
  if(numCasesAdd>0) smpl <- sample(indx0, size=numCasesAdd)
  if(numCasesAdd<0) smpl <- sample(indx1, size=-numCasesAdd)
  
  obsNew <- obs
  obsNew[smpl] <- (!obs[smpl])*1
  predsNew <- preds
  predsNew[smpl,] <- 1-preds[smpl,] 
  
  
  aucs <- apply(predsNew, 2, function(pred) as.numeric(roc(obsNew, pred)$auc))
  
  
  resEdgeCont <- rbind(ccrs,  msrs,  ppps,  npps,  senss,  specs,  fprs,  fnrs, tsss, 
                       ccrsW, msrsW, pppsW, nppsW, senssW, specsW, fprsW, fnrsW, tsssW, aucs)
  rownames(resEdgeCont) <- c("ccr",  "msr",  "ppp", "npp",  "sens",  "spec", "fpr", "fnr", "tss",
                             "Wccr", "Wmsr", "Wppp","Wnpp", "Wsens", "Wspec","Wfpr","Wfnr","Wtss", "edgeAUC")
  colnames(resEdgeCont) <- colnames(scoreMat)
  
  res <- rbind(resEdgeDists, resEdgeCont)
  
  names(dimnames(res)) <- c("measure","score")
  
  return(res)
}


# I. Dag/Edge wise evaluation: Distance from true hyp to chosen hyp 

edgeDist <- function(dagChosen, dagTrue){
  difs <- dagChosen - dagTrue
  difs[which(difs<0)] <- 0
  res <- sum(difs)
  return(res)
}    
nonEdgeDist <- function(dagChosen, dagTrue){
  difs <- dagTrue - dagChosen
  difs[which(difs<0)] <- 0
  res <- sum(difs)
  return(res)
}
totalEdgeDist <- function(dagChosen, dagTrue) sum(abs(dagTrue -  dagChosen))


#   1) number of edges, non-edges and total diff between true and chosen hyp



#   2) get distancest to all hyp then weight by prob/conf

# II. Edge wise evaluation

contingencyTable <- function(pred, obs, ws=rep(1, length(pred))){
  
  
  ws <- ws/sum(ws)*length(pred)
  tab <- table(factor(pred, levels=c(0,1)), factor(obs, levels=c(0,1)))
  rsums <- apply(tab, 2, sum)
  tab <- rbind(tab, rsums)
  csums <- apply(tab, 1, sum)
  tab <- cbind(tab, csums)
  res <- as.numeric(tab)
  names(res) <- c("tn","fp","f","fn","tp","t","n","p","N")
  
  return(res)
}

correctCR <- function(conTab) (conTab["tn"] + conTab["tp"])/conTab["N"]
misCR <- function(conTab) (conTab["fn"] + conTab["fp"])/conTab["N"]
posPP <- function(conTab){
  if(conTab["p"]==0){
    res <- 1
  } else{
    res <- conTab["tp"]/conTab["p"]
  }
  return(res)
}
negPP <- function(conTab){
  if(conTab["n"]==0){
    res <- 1
  } else{
    res <- conTab["tn"]/conTab["n"]
  }
  return(res)
}
sensitivity <- function(conTab){
  if(conTab["t"]==0){
    # if t==0, tp==0, so you got 100% of trues 
    res <- 1
  } else{
    res <- conTab["tp"]/conTab["t"]
  }
  return(res)
}
specificity <- function(conTab){
  if(conTab["f"]==0){
    # if f==0, tn==0 so you got 100% the falses
    res <- 1
  } else{
    res <- conTab["tn"]/conTab["f"]
  }
  return(res)
}
fpr <- function(conTab){
  if(conTab["f"]==0){
    res <- 0
  } else{
    res <- conTab["fp"]/conTab["f"]
  }
  return(res)
}
fnr <- function(conTab){
  if(conTab["t"]==0){
    res <- 0
  } else{
    res <- conTab["fn"]/conTab["t"]
  }
  return(res)
}
tss <- function(conTab) sensitivity(conTab) + specificity(conTab) - 1


#  1) compare true hyp to a single chosen hyp
#    a) accuracy, ROC/AUC, TSS, recall, precison, etc 
#  2) compare true hyp to all hyps and then weight measures by prob
#    b) accuracy, ROC/AUC, TSS, recall, precison, etc



# III. Dag-set wise evaluation

# Ranking of two hypotheses as in Mooij et al, pg 32-33
quot <- function(score_hA, score_hB){
  rnk <- 1/min(score_hA, score_hB)*c(-1,1)[(score_hA < score_hB)*1+1]
  return(rnk)
}

compQuot <- function(score_hA, score_hB){
  rnk <- 1/(1-min(score_hA, score_hB))*c(-1,1)[(score_hA < score_hB)*1+1]
  return(rnk)
}

differ <- function(score_hA, score_hB){
  return(score_hA - score_hB)
}

probRD <- function(score_hA, score_hB){
  rnk <- max(score_hA, score_hB)*c(-1,1)[(score_hA < score_hB)*1+1]
  return(rnk)
}

# 1. divide hypothesis space in two - part with true hyp and other

aggregateScores <- function(scoreMat, dagHyps, trueDag, groupFuncs, groupFuncPars, aggFunc){
  
  # construct names for the different bipartition methods groupFuncs and their parameters groupFuncPars
  nms <- mapply(function(func, pars){
    # i <- 2; func <- groupFuncs[i]; pars <- groupFuncPars[[i]]
    aux <- unlist(pars)
    res <- func
    if(!is.null(aux)){
      bit1 <- c("ED","NED","TED")[match(aux["distFun"], c("edgeDist","nonEdgeDist","totalEdgeDist"))]
      bit2 <- paste("k", aux["k"], sep="")
      res <- paste(res, bit1, bit2, sep="")
    }
    return(res)
  }, func=groupFuncs, pars=groupFuncPars, USE.NAMES = F)
  
  
  # Aggregate score matrix
  #print("aggregating score matrix")
  aggScoreMat <- mapply(function(grpFunc, grpFuncPars){
    # i <- 1; grpFunc <- groupFuncs[i]; grpFuncPars <- groupFuncPars[[i]]
    #print("*************")
    #print(paste("group function: ", grpFunc))
    pars <- grpFuncPars
    pars$dagHyps <- dagHyps
    pars$trueDag <- trueDag
    grp <- do.call(grpFunc, pars)
    indx0 <- which(grp==0)
    indx1 <- which(grp==1)
    res <- apply(scoreMat, 2, aggFunc, indx1=indx0, indx2=indx1)
    return(res)
  }, grpFunc=groupFuncs, grpFuncPars=groupFuncPars, SIMPLIFY="array")
  
  names(dimnames(aggScoreMat))[c(1,3)] <- c("hypothesis", "bipartCrit")
  dimnames(aggScoreMat)$bipartCrit <- nms
  
  return(aggScoreMat)
  
}


# a ) true hyp vs else

agg_truevsRest <- function(dagHyps, trueDag){
  grp <- apply(dagHyps, 3, function(dag) any(duplicated(abind(dag, trueDag, along=3), MARGIN=3)))*1
  return(grp)
}

# b )  +- k edges 

agg_distKvsRest <- function(dagHyps, trueDag, distFun, k){
  # is this function the same as agg_truevsRest if I pass k=0??
  distss <- apply(dagHyps, 3, distFun, dagTrue=trueDag)
  grp <- (distss <= k)*1  
  return(grp)
}

#     i) true hyp +- k edges vs else
#     ii) true hyp + k edges vs else (less independencies)
#     iii) true hyp - k edges vs else (more independencies)
# c)   Independencies
#    i) Markov equivalence class vs else

agg_MECvsRest <- function(dagHyps, trueDag){
  
  p <- nrow(trueDag)
  MEdagTrue <- pdag2allDags(dag2cpdag(trueDag))$dags
  MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
  grp <- (duplicated(abind(MEdagTrue, dagHyps ,along=3), MARGIN=3)[(dim(MEdagTrue)[3]+1):(dim(dagHyps)[3]+dim(MEdagTrue)[3])])*1
  return(grp)
}

#    ii) Markov equivalence class +- k edges

agg_distKMECvsRest <- function(dagHyps, trueDag, distFun, k){
  # is this the same as agg_distKvsRest if I pass k=0??
  p <- nrow(trueDag)
  MEdagTrue <- pdag2allDags(dag2cpdag(trueDag))$dags
  MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
  distss <- apply(MEdagTrue, 3, function(dag) apply(dagHyps, 3, distFun, dagTrue=dag))
  dim(distss) <- c(dim(dagHyps)[3], dim(MEdagTrue)[3])
  grp <- apply(distss,1, function(row) any(row <= k)*1)
  return(grp)
}

#     A) true hyp +- k edges vs else
#     B) true hyp + k edges vs else (less independencies)
#     C) true hyp - k edges vs else (more independencies)        

# 2: aggregate scores by hypothesis group (max?, mean?, dircetly on score and on prob)

aggMin <- function(x, indx1, indx2){
  
  if(length(indx1)==0){
    res1 <- max(x)+100 #Inf
  } else{
    res1 <- min(x[indx1])
  }
  
  if(length(indx2)==0){
    res2 <- max(x)+100 # Inf
  } else{
    res2 <- min(x[indx2])
  }
  
  res <- c(res1, res2)
  names(res) <- c("0","1")
  return(res)
  
}

aggSum <- function(x, indx1, indx2){
  
 
  if(length(indx1)==0){
    res1 <- 0 
  } else{
    res1 <- sum(x[indx1])
  }
  
  if(length(indx2)==0){
    res2 <- 0
  } else{
    res2 <- sum(x[indx2])
  }
  
  res <- c(res1, res2)
  names(res) <- c("0","1")
  return(res)
  
}


# 3. Add metascore -> combine scores from same recipe to obtain new score

# 4. Add metalearner -> combine scores from different learner-dataSetting 
# (but same approxDagSet method to keep hypothesis space the same) 
# eg. use  uniformity p-value score accross learners to decide which to use

# 5.  Rank transformation on score (This is only for 2 hypotheses so must be done after dividing hypothesis space)

rankedDecisions <- function(scoreMat, cmplxPack, rank2Func=NULL){
  
  # For Ranked decisions
  
  # Transform scores  according to complexity pack chosen
  #print("obtaining rank decision functions")
  if(! is.null(rank2Func)){
    rank2Func <- rep(rank2Func, ncol(scoreMat))
  } else{
    tab <- eval(parse(text=cmplxPack))$tab
    entropyFuncs <- eval(parse(text=cmplxPack))$entropyFuncs
    cmplxFunc <- colnames(scoreMat)
    entropIndx <- which(cmplxFunc %in% entropyFuncs)
    cmplxFunc[entropIndx] <- "score_sumMarginalEntropies"
    indx <- match(cmplxFunc, tab$cmplxFunc)
    rank2Func <- tab$rank2Func[indx]
    
  }
  
  # we add both ways because later we will use one or the other to balance "positive" cases
  #print("applying rank decision functions")
  rnkDec <- sapply(1:dim(scoreMat)[2], function(i){
    res <- apply(scoreMat[,i,], 2, function(hyps){
      res1 <- do.call(rank2Func[i], list(hyps[1], hyps[2]))
      res2 <- do.call(rank2Func[i], list(hyps[2], hyps[1]))
      res <- c(res1, res2)
      names(res) <- c("hyp1pos","hyp0pos")
      return(res)
    })
  }, simplify="array")
  
  rnkDec <- aperm(rnkDec, c(1,3,2))
  dimnames(rnkDec)[[2]] <- dimnames(scoreMat)[[2]]
  names(dimnames(rnkDec))[2] <- names(dimnames(scoreMat))[2]
  
  return(rnkDec)
}


# 6: Calculate 2-class performance measures: accuracy, ROC/AUC, TSS, recall, precison, etc

statsArrayToDB <- function(statsArr){
  statss <- melt(statsArr)
  
  aux <- strsplit(as.character(statss$fullRecipe), "\\.")
  statss$learner <- sapply(aux, function(el) el[1])
  statss$approxMethod <- sapply(aux, function(el) el[2])
  statss$dataReg <- sapply(aux, function(el) el[3])
  statss$cmplxFunc <- sapply(aux, function(el) el[4])
  statss$bipartCrit <- sapply(aux, function(el) el[5])
  statss <- statss[,c("fullRecipe", "learner", "approxMethod", "dataReg", "cmplxFunc","bipartCrit","measure","statistic","value")]
  return(statss)  
  
}

measures2DagClass <- function(rankDecArray, posNeg, ws, n.boot){
  
  indxMat <- cbind(1:n, posNeg+1)
  pm <- proc.time()
  statss <- apply(rankDecArray, 3, function(mat){
    # mat <- rnkDecsArr[,,1]
    predHard <- (apply(mat,1, which.max)-1)
    predHard[which(!posNeg)] <- (!predHard[which(!posNeg)])*1 
    predSoft <- mat[indxMat]
    contBased <- getStats(posNeg, predHard, ws, n.boot, stat="cntgncyStats")
    rankBased <- getStats(posNeg, predSoft, ws, n.boot, stat="aucStat")
    stats <- cbind(contBased, rankBased)
    return(stats)
  })
  proc.time()-pm # 33 mins for 2000 bootstrap samples for 50 data sets
  
  nmsFullRecipe <- dimnames(statss)[2]
  
  # inner stats matrix which we can't access from outer environment has 7 x 9:
  # 7 stats  5 quantiles, the mean and the original 
  # 10 measures: ccr, mcr, ppp, npp, fpr, fnr, sens, spec, tss, auc
  
  mat <- rnkDecsArr[,,1]
  predHard <- (apply(mat,1, which.max)-1)
  predHard[which(!posNeg)] <- (!predHard[which(!posNeg)])*1 
  predSoft <- mat[indxMat]
  contBased <- getStats(posNeg, predHard, ws, n.boot, stat="cntgncyStats")
  rankBased <- getStats(posNeg, predSoft, ws, n.boot, stat="aucStat")
  stats <- cbind(contBased, rankBased)
  
  dim(statss) <- c(dim(stats), dim(statss)[2])
  dimnames(statss) <- c(dimnames(stats), nmsFullRecipe)
  names(dimnames(statss)) <- c("statistic", "measure", "fullRecipe")
  
  return(statss)
}

getStats <- function(obs, pred, ws, n.boot, stat){
  data <- cbind(obs=obs, pred=pred, ws=ws)
  boots <- boot(data, statistic=eval(parse(text=stat)), R=n.boot, sim="ordinary", stype="i")
  qs <- apply(boots$t, 2 , quantile, probs=c(0.025, 0.18, 0.5, 0.84, 0.975))
  mu <- apply(boots$t,2, mean)
  dim(mu) <- c(1, ncol(qs))
  stats <- rbind(qs[1:2,,drop=F], mu, qs[3,,drop=F], boots$t0, qs[4:5,,drop=F])
  rownames(stats) <- c("q025", "q18", "mu", "q50", "orig", "q84", "q975")
  return(stats)
}

cntgncyStats <- function(data, indx=1:nrow(data)){
  
  pred <- data[indx,"pred"]
  obs <- data[indx, "obs"]
  ws <- data[indx, "ws"]
  
  
  # don't name the same as something outside otherwise it screws up (in this case)
  tab <- contingencyTable(pred, obs, ws)
  
  ccr <- correctCR(tab)
  mcr <- misCR(tab) 
  ppp <- posPP(tab)
  npp <- negPP(tab)
  sens <- sensitivity(tab)
  spec <- specificity(tab)
  fpr <- fpr(tab)
  fnr <- fnr(tab)
  tss <- tss(tab)
  res <- c(ccr, mcr, ppp, npp, sens, spec, fpr, fnr, tss)
  names(res) <- c("ccr", "mcr", "ppp", "npp", "sens", "spec", "fpr", "fnr", "tss")
  return(res)
} 

aucStat <- function(data, indx=1:nrow(data)){
  res <- WeightedAUC(WeightedROC(data[indx,"pred"], data[indx,"obs"], data[indx,"ws"]))
  names(res) <- "auc"
  return(res)
} 

###############################################################################################################
# modifications to functions from dHSIC package be able to use hsic test by passing it pre-calculated kernels
dhsic <- function (X, Y, K, kernel = "gaussian", bandwidth = 1, matrix.input = FALSE) {
  if (missing(K)) {
    if (!missing(Y)) {
      X <- list(X, Y)
    }
    if (matrix.input) {
      if (!is.matrix(X)) {
        stop("X needs to be a matrix if matrix.input=TRUE")
      }
      else {
        X <- split(X, rep(1:ncol(X), each = nrow(X)))
      }
    }
    d <- length(X)
    for (j in 1:d) {
      if (is.matrix(X[[j]]) == FALSE) {
        X[[j]] <- as.matrix(X[[j]])
      }
    }
    len <- nrow(X[[1]])
    if (len < 2 * d) {
      warning("Sample size is smaller than twice the number of variables. dHSIC is trivial.")
      result = list(dHSIC = 0, time = c(GramMat = 0, HSIC = 0))
      return(result)
    }
    if (length(kernel) < d) {
      kernel <- rep(kernel[1], d)
    }
    if (length(bandwidth) < d) {
      bandwidth <- rep(bandwidth[1], d)
    }
    median_bandwidth <- function(x) {
      bandwidth <- dHSIC:::median_bandwidth_rcpp(x[sample(1:len), , drop = FALSE], len, ncol(x))
      if (bandwidth == 0) {
        bandwidth <- 0.001
      }
      return(bandwidth)
    }
    custom_grammat <- function(x, fun) {
      KX <- matrix(nrow = len, ncol = len)
      for (i in 1:len) {
        for (j in i:len) {
          KX[i, j] <- match.fun(fun)(x[i, ], x[j, ])
          KX[j, i] <- KX[i, j]
        }
      }
      return(KX)
    }
    K <- vector("list", d)
    ptm <- proc.time()
    for (j in 1:d) {
      if (kernel[j] == "gaussian") {
        bandwidth[j] <- median_bandwidth(X[[j]])
        K[[j]] <- dHSIC:::gaussian_grammat_rcpp(X[[j]], bandwidth[j], len, ncol(X[[j]]))
      }
      else if (kernel[j] == "gaussian.fixed") {
        K[[j]] <- dHSIC:::gaussian_grammat_rcpp(X[[j]], bandwidth[j], len, ncol(X[[j]]))
      }
      else if (kernel[j] == "discrete") {
        bandwidth[j] <- NA
        K[[j]] <- discrete_grammat_rcpp(X[[j]], len, ncol(X[[j]]))
      }
      else {
        bandwidth[j] <- NA
        K[[j]] <- custom_grammat(X[[j]], kernel[j])
      }
    }
    timeGramMat <- as.numeric((proc.time() - ptm)[1])
  }
  else {
    if (is.list(K)) {
      d <- length(K)      
      #len <- nrow(K) ### AQUI CAMBIE ALGO!!!!!!!!!!!!!
      len <- nrow(K[[1]]) # originally len <- nrow(K)
    }
    else {
      stop("K needs to be a list of matrices")
    }
  }
  ptm <- proc.time()
  term1 <- 1
  term2 <- 1
  term3 <- 2/len
  for (j in 1:d) {
    term1 <- term1 * K[[j]]
    term2 <- 1/len^2 * term2 * sum(K[[j]])
    term3 <- 1/len * term3 * colSums(K[[j]])
  }
  term1 <- sum(term1)
  term3 <- sum(term3)
  dHSIC = 1/len^2 * term1 + term2 - term3
  timeHSIC <- as.numeric((proc.time() - ptm)[1])   ### AQUI CAMBIE ALGO!!!!!!!!!!!!!
  #result = list(dHSIC = dHSIC, time = c(GramMat = timeGramMat,HSIC = timeHSIC), bandwidth = bandwidth)
  result = list(dHSIC = dHSIC, bandwidth = bandwidth)
  return(result)
}

dhsic.test <- dHSIC:::dhsic.test

dhsic.testK <- function (X, Y, K, alpha = 0.05, method = "permutation", kernel = "gaussian", 
                         B = 1000, pairwise = FALSE, bandwidth = 1, matrix.input = FALSE){
  if (pairwise & missing(K)) {
    if (matrix.input) {
      if (!is.matrix(X)) {
        stop("X needs to be a matrix if matrix.input=TRUE")
      }
      else {
        X <- split(X, rep(1:ncol(X), each = nrow(X)))
      }
    }
    if (!missing(Y) || !is.list(X) || length(X) <= 2) {
      stop("pairwise only makes sense if number of variables is greater than 2")
    }
    d <- length(X)
    for (j in 1:d) {
      if (is.matrix(X[[j]]) == FALSE) {
        X[[j]] <- as.matrix(X[[j]])
      }
    }
    len <- nrow(X[[1]])
    if (len < 2 * d) {
      warning("Sample size is smaller than twice the number of variables. Test is trivial.")
      test = list(statistic = 0, crit.value = Inf, p.value = 1, 
                  time = c(GramMat = 0, dHSIC = 0, CritVal = 0))
      return(test)
    }
    pValVec <- rep(0, d - 1)
    statVec <- rep(0, d - 1)
    critVec <- rep(0, d - 1)
    ptm <- proc.time()
    for (j in (d:2)) {
      resTmp <- dhsic.test(do.call("cbind", X[1:(j - 1)]), 
                           X[[j]], alpha = alpha/(d - 1), method = method, 
                           kernel = kernel, B = B, pairwise = FALSE, bandwidth = bandwidth)
      pValVec[d - j + 1] <- resTmp$p.value
      statVec[d - j + 1] <- resTmp$statistic
      critVec[d - j + 1] <- resTmp$crit.value
    }
    ind <- which.min(pValVec)
    stat <- statVec[ind]
    critical_value <- critVec[ind]
    p_value = min(pValVec[ind] * (d - 1), 1)
    timeCritVal <- as.numeric((proc.time() - ptm)[1])
    timeGramMat <- as.numeric(resTmp$time[1])
    timeHSIC <- as.numeric(resTmp$time[2])
  }
  else {
    if (missing(K)) {
      if (matrix.input) {
        if (!is.matrix(X)) {
          stop("X needs to be a matrix if matrix.input=TRUE")
        }
        else {
          X <- split(X, rep(1:ncol(X), each = nrow(X)))
        }
      }
      if (!missing(Y)) {
        X <- list(X, Y)
      }
      d <- length(X)
      for (j in 1:d) {
        if (is.matrix(X[[j]]) == FALSE) {
          X[[j]] <- as.matrix(X[[j]])
        }
      }
      len <- nrow(X[[1]])
      if (len < 2 * d) {
        warning("Sample size is smaller than twice the number of variables. Test is trivial.")
        test = list(statistic = 0, crit.value = Inf, 
                    p.value = 1, time = c(GramMat = 0, dHSIC = 0, 
                                          CritVal = 0), bandwidth = NA)
        return(test)
      }
      if (length(kernel) < d) {
        kernel <- rep(kernel[1], d)
      }
      if (length(bandwidth) < d) {
        bandwidth <- rep(bandwidth[1], d)
      }
      median_bandwidth <- function(x) {
        bandwidth <- median_bandwidth_rcpp(x[sample(len, 
                                                    replace = FALSE), , drop = FALSE], len, ncol(x))
        if (bandwidth == 0) {
          bandwidth <- 0.001
        }
        return(bandwidth)
      }
      custom_grammat <- function(x, fun) {
        KX <- matrix(nrow = len, ncol = len)
        for (i in 1:len) {
          for (j in i:len) {
            KX[i, j] <- match.fun(fun)(x[i, ], x[j, ])
            KX[j, i] <- KX[i, j]
          }
        }
        return(KX)
      }
      K <- vector("list", d)
      ptm <- proc.time()
      for (j in 1:d) {
        if (kernel[j] == "gaussian") {
          bandwidth[j] <- median_bandwidth(X[[j]])
          K[[j]] <- gaussian_grammat_rcpp(X[[j]], bandwidth[j], 
                                          len, ncol(X[[j]]))
        }
        else if (kernel[j] == "gaussian.fixed") {
          K[[j]] <- gaussian_grammat_rcpp(X[[j]], bandwidth[j], 
                                          len, ncol(X[[j]]))
        }
        else if (kernel[j] == "discrete") {
          bandwidth[j] <- NA
          K[[j]] <- discrete_grammat_rcpp(X[[j]], len, 
                                          ncol(X[[j]]))
        }
        else {
          bandwidth[j] <- NA
          K[[j]] <- custom_grammat(X[[j]], kernel[j])
        }
      }
      timeGramMat <- as.numeric((proc.time() - ptm)[1])
    }
    else {
      if (is.list(K)) {
        d <- length(K)
        #len <- nrow(K)    # ### AQUI CAMBIE ALGO!!!!!!!!!!!!!
        len <- nrow(K[[1]]) # originally len <- nrow(K)
      }
      else {
        stop("K needs to be a list of matrices")
      }
    }
    ptm <- proc.time()
    term1 <- 1
    term2 <- 1
    term3 <- 2/len
    for (j in 1:d) {
      term1 <- term1 * K[[j]]
      term2 <- 1/len^2 * term2 * sum(K[[j]])
      term3 <- 1/len * term3 * colSums(K[[j]])
    }
    term1 <- sum(term1)
    term3 <- sum(term3)
    dHSIC = 1/len^2 * term1 + term2 - term3
    timeHSIC <- as.numeric((proc.time() - ptm)[1])
    ptm <- proc.time()
    if (method == "permutation") {
      dhsic_perm_fun <- function(i) {
        term1 <- K[[1]]
        term2 <- 1/len^2 * sum(K[[1]])
        term3 <- 2/len^2 * colSums(K[[1]])
        for (j in 2:d) {
          perm <- sample(0:(len - 1), replace = FALSE)
          #K_perm <- shuffle_grammat_rcpp(K[[j]], perm, len)   ### AQUI CAMBIE ALGO!!!!
          K_perm <- dHSIC:::shuffle_grammat_rcpp(K[[j]], perm, len)
          term1 <- term1 * K_perm
          term2 <- 1/len^2 * term2 * sum(K_perm)
          term3 <- 1/len * term3 * colSums(K_perm)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2 * term1 + term2 - term3)
      }
      dHSIC_perm <- sapply(1:B, dhsic_perm_fun)
      sortdHSIC <- sort(len * dHSIC_perm)
      Bind <- sum(len * dHSIC == sortdHSIC) + ceiling((1 - 
                                                         alpha) * (B + 1))
      if (Bind <= B) {
        critical_value <- sortdHSIC[Bind]
      }
      else {
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_perm >= dHSIC) + 1)/(B + 1)
    }
    if (method == "bootstrap") {
      dhsic_boot_fun <- function(i) {
        term1 <- K[[1]]
        term2 <- 1/len^2 * sum(K[[1]])
        term3 <- 2/len^2 * colSums(K[[1]])
        for (j in 2:d) {
          boot <- sample(0:(len - 1), replace = TRUE)
          # K_boot <- shuffle_grammat_rcpp(K[[j]], boot, len)  ### AQUI CAMBIE ALGO!!!
          K_boot <- dHSIC:::shuffle_grammat_rcpp(K[[j]], boot, len)
          term1 <- term1 * K_boot
          term2 <- 1/len^2 * term2 * sum(K_boot)
          term3 <- 1/len * term3 * colSums(K_boot)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2 * term1 + term2 - term3)
      }
      dHSIC_boot <- sapply(1:B, dhsic_boot_fun)
      sortdHSIC <- sort(len * dHSIC_boot)
      Bind <- sum(len * dHSIC == sortdHSIC) + ceiling((1 - 
                                                         alpha) * (B + 1))
      if (Bind <= B) {
        critical_value <- sortdHSIC[Bind]
      }
      else {
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_boot >= dHSIC) + 1)/(B + 1)
    }
    else if (method == "gamma") {
      est.a <- vector("numeric", d)
      est.b <- vector("numeric", d)
      est.c <- vector("numeric", d)
      for (j in 1:d) {
        est.a[j] <- 1/(len^2) * sum(K[[j]])
        est.b[j] <- 1/(len^2) * sum(K[[j]]^2)
        est.c[j] <- 1/len^3 * sum(rowSums(K[[j]])^2)
      }
      prod.a <- prod(est.a)
      prod.b <- prod(est.b)
      prod.c <- prod(est.c)
      oneoutprod.a <- vector("numeric", d)
      oneoutprod.b <- vector("numeric", d)
      oneoutprod.c <- vector("numeric", d)
      for (j in 1:d) {
        oneoutprod.a[j] <- prod.a/est.a[j]
        oneoutprod.b[j] <- prod.b/est.b[j]
        oneoutprod.c[j] <- prod.c/est.c[j]
      }
      est.d <- est.a^2
      prod.d <- prod.a^2
      oneoutprod.d <- oneoutprod.a^2
      exp.est <- (1 - sum(oneoutprod.a) + (d - 1) * prod.a)/len
      term1 <- prod.b
      term2 <- (d - 1)^2 * prod.d
      term3 <- 2 * (d - 1) * prod.c
      term4 <- 0
      term5 <- 0
      term6 <- 0
      term7 <- 0
      for (r in 1:(d - 1)) {
        term4 <- term4 + est.b[r] * oneoutprod.d[r]
        term5 <- term5 + est.b[r] * oneoutprod.c[r]
        term6 <- term6 + est.c[r] * oneoutprod.d[r]
        for (s in (r + 1):d) {
          term7 <- term7 + 2 * est.c[r] * est.c[s] * 
            oneoutprod.d[r]/est.d[s]
        }
      }
      term4 <- term4 + est.b[d] * oneoutprod.d[d]
      term5 <- -2 * (term5 + est.b[d] * oneoutprod.c[d])
      term6 <- -2 * (d - 1) * (term6 + est.c[d] * oneoutprod.d[d])
      factor1 <- len - 2 * d
      factor2 <- len * (len - 1) * (len - 2)
      for (j in 1:(2 * d - 3)) {
        factor1 <- factor1 * (len - 2 * d - j)
        factor2 <- factor2 * (len - 2 - j)
      }
      var.est = 2 * factor1/factor2 * (term1 + term2 + 
                                         term3 + term4 + term5 + term6 + term7)
      a <- (exp.est^2)/var.est
      b <- len * var.est/exp.est
      critical_value <- qgamma(1 - alpha, shape = a, scale = b)
      p_value <- pgamma(len * dHSIC, shape = a, scale = b, 
                        lower.tail = FALSE)
    }
    else if (method == "eigenvalue") {
      est1 <- vector("numeric", d)
      est2 <- matrix(NA, len, d)
      for (j in 1:d) {
        est1[j] <- 1/(len * (len - 1)) * sum(K[[j]] - 
                                               diag(diag(K[[j]])))
        est2[, j] <- 1/len * rowSums(K[[j]])
      }
      est1_prod <- prod(est1)
      est2_prod <- apply(est2, 1, prod)
      a1 <- matrix(1, len, len)
      a2 <- (d - 1)^2 * est1_prod
      a3 <- (d - 1) * est2_prod
      a5 <- matrix(0, len, len)
      a6 <- matrix(0, len, len)
      a8 <- matrix(0, len, len)
      a9 <- matrix(0, len, 1)
      for (j in 1:d) {
        a1 <- a1 * K[[j]]
        a5 <- a5 + K[[j]] * est1_prod/est1[j]
        a6 <- a6 + K[[j]] * matrix(rep(est2_prod/est2[, 
                                                      j], len), len, len)
        a9 <- a9 + est2[, j, drop = FALSE] * est1_prod/est1[j]
        i <- j + 1
        while (i <= d) {
          a8 <- a8 + est2[, j, drop = FALSE] %*% t(est2[, 
                                                        i, drop = FALSE]) * est1_prod/(est1[j] * 
                                                                                         est1[i]) + est2[, i, drop = FALSE] %*% t(est2[, 
                                                                                                                                       j, drop = FALSE]) * est1_prod/(est1[j] * 
                                                                                                                                                                        est1[i])
          i <- i + 1
        }
      }
      a3 <- matrix(rep(a3, len), len, len)
      a4 <- t(a3)
      a6 <- -a6
      a7 <- t(a6)
      a8 <- a8
      a9 <- matrix(rep((1 - d) * a9, len), len, len)
      a10 <- t(a9)
      H2 <- 1/(d * (2 * d - 1)) * (a1 + a2 + a3 + a4 + 
                                     a5 + a6 + a7 + a8 + a9 + a10)
      eigenvalues <- eigen(H2, only.values = TRUE, symmetric = TRUE)$values/len
      M <- 5000
      Z <- rnorm(M * len)^2
      chi_dist <- d * (2 * d - 1) * colSums(matrix(Z * 
                                                     eigenvalues, len, M))
      critical_value <- as.numeric(quantile(chi_dist, 1 - 
                                              alpha))
      p_value <- 1 - ecdf(chi_dist)(dHSIC * len)
    }
    timeCritVal <- as.numeric((proc.time() - ptm)[1])
    stat <- dHSIC * len
  }
  #test = list(statistic = stat, crit.value = critical_value, 
  #            p.value = p_value, time = c(GramMat = timeGramMat, dHSIC = timeHSIC, 
  #                                        CritVal = timeCritVal), bandwidth = bandwidth)
  
  test = list(statistic = stat, crit.value = critical_value, 
              p.value = p_value, bandwidth = bandwidth)
  
  
  return(test)
}

#################################################################################################################
# Miscellaneous

# for pairs plots, to plot smoother and calculate 1-pval_hsic as a dependence strength
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- 1-dhsic.test(list(x,y))$p.value
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) 
    cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r) 
}

