library(CVST) # cross validation framework
library(kernlab) # kernels, krr
library(lbfgs) # non-convex optimization for hsic regression
library(numDeriv) #comparing numerical approx of gradient with analytical

library(pcalg) # unifDAG, graphNELs
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
########################################################################################
# functions HSICreg
########################################################################################

constructLearner <- function (learn, predict, loss, getFixedParams){
    stopifnot(is.function(learn) && is.function(predict) && is.function(loss) && is.function(getFixedParams))
    learner = list(learn = learn, predict = predict, loss=loss, getFixedParams=getFixedParams)
    class(learner) = "CVST.learner"
    return(learner)
}

learn.krr <- function(data, params) {
	stopifnot(isRegression(data))
    kernelXs = do.call(params$kernelXs, params$psKernXs)
    return(.krr(data$x, kernelXs, data$y, getN(data) * params$lambda))
}
predict.krr <- function(model, newData) {
	stopifnot(isRegression(newData))
    return(as.matrix(.krr.predict(newData$x, model)))
}
    
learn.qhsic <- function (data, params) {
    stopifnot(isRegression(data))
    kernelXs = do.call(params$kernelXs, params$psKernXs)
	kernelXb = do.call(params$kernelXb, params$psKernXb)
	return(.qhsic(data$x, kernelXs, kernelXb, data$y, getN(data) * params$lambda))
}
predict.qhsic <- function (model, newData) {
    stopifnot(isRegression(newData))
    return(as.matrix(.qhsic.predict(newData$x, model)))
}	


# start with alpha = 0
learn.hsic <- function (data, params) {
    stopifnot(isRegression(data))  
	#print("enters hsic learn function")
    kernelXs = do.call(params$kernelXs, params$psKernXs)
	kernelXb = do.call(params$kernelXb, params$psKernXb)
	kernelRg = do.call(params$kernelRg, params$psKernRg)
    #print("exits hsic learn function")
	return(.hsic(data$x, kernelXs, kernelXb, kernelRg, data$y, getN(data) * params$lambda))
}
predict.hsic <- function (model, newData) {
    stopifnot(isRegression(newData))
    return(as.matrix(.hsic.predict(newData$x, model)))
}

# start with alpha_qhsic
learn.hsic2 <- function (data, params) {
  stopifnot(isRegression(data))  
  #print("enters hsic learn function")
  kernelXs = do.call(params$kernelXs, params$psKernXs)
  kernelXb = do.call(params$kernelXb, params$psKernXb)
  kernelRg = do.call(params$kernelRg, params$psKernRg)
  #print("exits hsic learn function")
  return(.hsic2(data$x, kernelXs, kernelXb, kernelRg, data$y, getN(data) * params$lambda))
}

# start with gamma close to zero and alpha_qhsic
learn.hsic3 <- function (data, params) {
  stopifnot(isRegression(data))  
  #print("enters hsic learn function")
  kernelXs = do.call(params$kernelXs, params$psKernXs)
  kernelXb = do.call(params$kernelXb, params$psKernXb)
  kernelRg = do.call(params$kernelRg, params$psKernRg)
  #print("exits hsic learn function")
  return(.hsic3(data$x, kernelXs, kernelXb, kernelRg, data$y, getN(data) * params$lambda))
}

# use estimate of alpha from previous lambda -> look at graduated convexity to solve with this approach
learn.hsic4 <- function (data, params) {
  stopifnot(isRegression(data))  
  #print("enters hsic learn function")
  kernelXs = do.call(params$kernelXs, params$psKernXs)
  kernelXb = do.call(params$kernelXb, params$psKernXb)
  kernelRg = do.call(params$kernelRg, params$psKernRg)
  #print("exits hsic learn function")
  return(.hsic4(data$x, kernelXs, kernelXb, kernelRg, data$y, getN(data) * params$lambda, params$ini))
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



.krr <- function (data, kernelXs, y, lambda){
    Kxs = kernelMatrix(kernelXs, data)
    N = nrow(Kxs)
    alpha = solve(Matrix(Kxs + diag(lambda, N))) %*% y
    return(list(data = data, kernelXs = kernelXs, alpha = alpha))
}
.krr.predict <- function (newData, krr) {
    kxs = kernelMatrix(krr$kernelXs, newData, krr$data)
    return(kxs %*% krr$alpha)
}

.qhsic <- function(data, kernelXs, kernelXb, y, lambda) {
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  alpha <- solve(Matrix(Kxbc%*%Kxs + diag(lambda, N))) %*% Kxbc %*% y
  return(list(data=data, kernelXs=kernelXs, kernelXb=kernelXb, alpha=alpha, avg=mean(y)))
}
.qhsic.predict <- function(newData, qhsic) {  
  kxs = kernelMatrix(qhsic$kernelXs, newData, qhsic$data)
  pred <- kxs %*% qhsic$alpha
  pred <- pred - mean(pred) + qhsic$avg
  return(pred)
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


.hsic <- function(data, kernelXs, kernelXb, kernelRg, y, lambda) {
  #print("enters .hsic")
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N  <- nrow(Kxs)
  
  #alpha0= rep(0, N)
  
  
  qhsic <- constructLearner(learn.qhsic, predict.qhsic, qhsicLoss, getFixedParams)
  psAux <- list(kernelXs="rbfdot", kernelXb="rbfdot", psKernXs=list(sigma=kernelXs@kpar$sigma),  psKernXb=list(sigma=kernelXb@kpar$sigma),lambda=lambda)
  
  
  dataAux <- constructData(x=data,y=y)
  modAux    <- qhsic$learn(dataAux, psAux)
  alpha0 <- as.numeric(modAux$alpha)
  
  
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1, max_iterations=100)
  
  alpha <- res$par
  #print("exits .hsic")
  return(list(data=data, kernelXs=kernelXs, kernelRg=kernelRg, alpha=alpha, avg=mean(y)))
}
.hsic.predict <- function(newData, hsic) {
  #require(kernlab)
  kxs = kernelMatrix(hsic$kernelXs, newData, hsic$data)
  pred <- kxs %*% hsic$alpha
  pred <- pred - mean(pred) + hsic$avg
  return(pred)
}
.hsic2 <- function(data, kernelXs, kernelXb, kernelRg, y, lambda) {
  #print("enters .hsic")
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N  <- nrow(Kxs)
  
  #alpha0= rep(0, N)
  
  
  qhsic <- constructLearner(learn.qhsic, predict.qhsic, qhsicLoss, getFixedParams)
  psAux <- list(kernelXs="rbfdot", kernelXb="rbfdot", psKernXs=list(sigma=kernelXs@kpar$sigma),  psKernXb=list(sigma=kernelXb@kpar$sigma),lambda=lambda/(2*1e-10))
  
  
  dataAux <- constructData(x=data,y=y)
  modAux    <- qhsic$learn(dataAux, psAux)
  alpha0 <- as.numeric(modAux$alpha)
  
  
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1)
  
  alpha <- res$par
  #plot(alpha0, alpha, main="before after")
  #print("exits .hsic")
  return(list(data=data, kernelXs=kernelXs, kernelRg=kernelRg, alpha=alpha, avg=mean(y)))
}
.hsic3 <- function(data, kernelXs, kernelXb, kernelRg, y, lambda) {
  #print("enters .hsic")
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N  <- nrow(Kxs)
  
  #alpha0= rep(0, N)
  
  #ord <- 10
  #gammas <- 10^seq(-ord,ord,1)
  
  gammaChosen <- kernelRg@kpar$sigma
  
  gammas <- gammaChosen*seq(0.01, 1, length.out=10)
  
  #indx <- which(gammas>gammaChosen)[1]
  #gammas <- gammas[1:indx]
  
  gamma0 <- gammas[1]
  

  
  qhsic <- constructLearner(learn.qhsic, predict.qhsic, qhsicLoss, getFixedParams)
  psAux <- list(kernelXs="rbfdot", kernelXb="rbfdot", psKernXs=list(sigma=kernelXs@kpar$sigma),  psKernXb=list(sigma=kernelXb@kpar$sigma),lambda=lambda/(2*gamma0))
  
  
  dataAux <- constructData(x=data,y=y)
  modAux    <- qhsic$learn(dataAux, psAux)
  alphaFirst <- as.numeric(modAux$alpha)
  
  kernelRg_t <- rbfdot(sigma=gamma0)
  
  print("first")
  pm <- proc.time()
  alphat <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alphaFirst, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg_t, lambda=lambda, invisible=1)$par
  print(proc.time()-pm)
  
  #plot(alphaFirst,alphat,main=paste("first, lambda= ", lambda))
  #abline(a=0,b=1, col="red")
  
  
  maxIter <- 10
  count <- 0
  stepConst <- (gammaChosen-gamma0)/maxIter
  gammat <- gamma0
  
  varsDists <- as.numeric() 
  times <- as.numeric()
  grads <- as.numeric()
  
  for(gamma in gammas){
    count <- count +1
    print(paste("count: ", count))
    gradAlpha <- .hsicRegLossGrad(alphat, y, Kxs, Kxb, kernelRg_t, lambda)
    gradGamma <- .hsicRegLossGradGamma(alphat, y, Kxs, Kxb, kernelRg_t)
    grad <- c(gradAlpha, gradGamma)
    normGrad <- sqrt(t(grad)%*%grad)
    #stepSize <- stepConst/normGrad
    #stepSize <- max(exp(-normGrad/gammat)*stepConst,gamma0)
    stepSize <- 10^(-(maxIter-count))
	gammat <- min(gammaChosen, gammat+stepSize)
    
    kernelRg_t <- rbfdot(sigma=gamma)
    #print(paste("gammat: ", gammat))
    #print(paste("step size %:", (stepSize/gammaChosen)*100))
    alphaLast <- alphat
    
	pm <- proc.time()
    res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alphaFirst, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg_t, lambda=lambda, invisible=1)
	times <- c(times, (proc.time()-pm)[3])
	alphat <- res$par
	
	grads <- c(grads, sqrt(sum((alphat-alphaLast)^2)))
	
	resids <- y - Kxs%*%alphat
	Krg <- as.matrix(kernelMatrix(kernelRg_t, resids))
	varsDists <- c(varsDists, var(Krg[lower.tri(Krg)]))
	
	#print(proc.time()-pm)
    
	#plot(alphaLast,alphat,main=paste("count: ", count, ", gamma: ", gammat, ", lambda: ", lambda))
    #abline(a=0,b=1, col="red")
  }
  
  maxVar <- gammas[which.max(varsDists)]
  maxGrad <- gammas[which.max(grads)]
  
  
  #print("last")
  pm <- proc.time()
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alphaFirst, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1)
  totTime <- (proc.time()-pm)[3]
  
  alpha <- res$par
  
  par(mfrow=c(2,2))
  plot(gammas, varsDists)
  abline(v=c(maxVar, maxGrad), col=c("red","green"))
  
  plot(gammas, grads)
  abline(v=c(maxVar, maxGrad), col=c("red","green"))
  
  plot(gammas, times, ylim=c(0,max(c(times, totTime/length(gammas)))))
  abline(v=c(maxVar, maxGrad), col=c("red","green"))
  abline(h=c(mean(times),totTime), col=c("blue","purple"))
  
  plot(alphat, alpha, main=.hsicRegLoss(alpha, y, Kxs, Kxb, kernelRg, lambda)-.hsicRegLoss(alphat, y, Kxs, Kxb, kernelRg, lambda))
  abline(a=0, b=1, col="red")
  
  par(mfrow=c(1,1))
  
  
  
  
  #plot(alphat,alpha, main=paste("last, gamma: ", gammat, ", lambda: ", lambda))
  #abline(a=0,b=1, col="red")
  #print("exits .hsic")
  return(list(data=data, kernelXs=kernelXs, kernelRg=kernelRg, alpha=alpha, avg=mean(y)))
}
.hsic4 <- function(data, kernelXs, kernelXb, kernelRg, y, lambda, alpha0) {
  #print("enters .hsic")
  Kxs <- kernelMatrix(kernelXs, data)
  Kxb <- kernelMatrix(kernelXb, data)
  N  <- nrow(Kxs)
  
  #alpha0= rep(0, N)
  
  
  #qhsic <- constructLearner(learn.qhsic, predict.qhsic, qhsicLoss, getFixedParams)
  #psAux <- list(kernelXs="rbfdot", kernelXb="rbfdot", psKernXs=list(sigma=kernelXs@kpar$sigma),  psKernXb=list(sigma=kernelXb@kpar$sigma),lambda=lambda)
  
  
  #dataAux <- constructData(x=data,y=y)
  #modAux    <- qhsic$learn(dataAux, psAux)
  #alpha0 <- as.numeric(modAux$alpha)
  
  
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1)
  
  alpha <- res$par
  #print("exits .hsic")
  return(list(data=data, kernelXs=kernelXs, kernelRg=kernelRg, alpha=alpha, avg=mean(y)))
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


CV <- function(data, learner, params, fold=5, fac=1, verbose=TRUE) {
  stopifnot(class(learner) == "CVST.learner" &&
            class(data) == "CVST.data" &&
            class(params) == "CVST.params")
  
  print("enters CV")
 
  nParams <- length(params)
  dimnames <- list(as.character(1:fold), names(params))

  lossMatTest  <- matrix(0, fold, nParams, dimnames=dimnames)
  lossMatTrain  <- matrix(0, fold, nParams, dimnames=dimnames)
  mseL <- matrix(0, fold, nParams, dimnames=dimnames)
  qhsicL <- matrix(0, fold, nParams, dimnames=dimnames)
  hsicL <- matrix(0, fold, nParams, dimnames=dimnames)
  regL <- matrix(0, fold, nParams, dimnames=dimnames)
  
  n <- getN(data)
  size =  ceiling(n / fold)

  #alphat <- matrix(0,size*(fold-1),fold)
  
  for (ind in 1:nParams) {
	  #print("*************************************")
	  #print(ind)
    p <- params[[ind]]
	
    for (f in 1:fold) {
      validationIndex <- seq((f-1)*size + 1, min(f*size, n))
      curTrain <- getSubset(data, -validationIndex)
      curTest <- getSubset(data, validationIndex)
      # either mean squared error or mean classification error
		
	  #p$ini <- alphat[,f]
	  model <- learner$learn(curTrain, p)
	  alpha <- model$alpha
  	  
	  #plot(alphat[,f], alpha)
	  #abline(a=0, b=1, col="red")
	  
	  #if(f < fold){
	  #	  alphat[,f+1] <- alpha
	  #} else{
	  #	  alphat[,1] <- alpha
	  #}
	  
  	  kernelXs <- do.call(p$kernelXs, p$psKernXs)
	  Kxs <-  kernelMatrix(kernelXs, curTrain$x)
      lossMatTest[f, ind] <- do.call(learner$loss, list(model=model, test=curTest, learner=learner, param=p)) 
	  lossMatTrain[f, ind] <- do.call(learner$loss, list(model=model, test=curTrain, learner=learner, param=p)) 
	  
	  mseL[f, ind] <- do.call("mse", list(model=model, test=curTest, learner=learner, param=p)) 
	  qhsicL[f, ind] <- do.call("qhsicLoss", list(model=model, test=curTest, learner=learner, param=p)) 
	  hsicL[f, ind] <- do.call("hsicLoss", list(model=model, test=curTest, learner=learner, param=p)) 
	  regL[f, ind] <- as.numeric(t(alpha)%*%Kxs%*%alpha)

	  
    }
    if (FALSE) {
      cat(names(params)[ind], "(", mean(lossMatTest[, ind]), ")\n")
    }
  }
  
  
  
  means <- apply(lossMatTest, 2, mean)
  lambdas <- sapply(params, function(x) x$lambda)
  meansFac <- means/min(means)
  indx <- which(meansFac <= fac)
  maxLambda <- max(lambdas[indx])
  winner <- which(lambdas==maxLambda)
  
  
  gridL <- abind(lossMatTest, lossMatTrain, mseL, qhsicL, hsicL, regL, along=3)
  dimnames(gridL) <- list(fold=1:fold, params=names(params), var=c("lossTest","lossTrain","mse","qhsicL","hsicL","regL"))
  
  print("exits CV")
  return(list(opt=params[winner], grid=gridL))
}

CV.parallel <- function(data, learner, params, fold=5, fac=1, verbose=TRUE) {
  stopifnot(class(learner) == "CVST.learner" &&
              class(data) == "CVST.data" &&
              class(params) == "CVST.params")
  
  print("enters CV parallel")
  
  nParams <- length(params)
  dimnames <- list(as.character(1:fold), names(params))
  
  n <- getN(data)
  size <- ceiling(n / fold)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- params[[6]]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- getSubset(data, -validationIndex)
      curTest <- getSubset(data, validationIndex)
      # either mean squared error or mean classification error
      
      
      model <- try(learner$learn(data=curTrain, params=p))
      if(class(model)=="try-error"){
        print(paste("singular matrix for lambda = ", p$lambda, " fold = ", f))
        res <- c(lossTest=NA, lossTrain=NA, rmseL=NA, qhsicL=NA, hsicL=NA, regL=NA, nhsicRegL=NA, hsicRegL=NA, hsicYhRegL=NA,
                 nhsicRegAL=NA, hsicRegAL=NA, hsicYhRegAL=NA, corr=NA)
        return(res)
      }
        
      alpha <- model$alpha
      
      kernelXs <- do.call(p$kernelXs, p$psKernXs)
      Kxs <-  kernelMatrix(kernelXs, curTrain$x)
      
      lossTest <- do.call(learner$loss, list(model=model, test=curTest, learner=learner, param=p)) 
      lossTrain <- do.call(learner$loss, list(model=model, test=curTrain, learner=learner, param=p)) 
      rmseL <- do.call("rmse", list(model=model, test=curTest, learner=learner, param=p)) 
      qhsicL <- do.call("qhsicLoss", list(model=model, test=curTest, learner=learner, param=p)) 
      hsicL <- do.call("hsicLoss", list(model=model, test=curTest, learner=learner, param=p)) 
      regL <- as.numeric(t(alpha)%*%Kxs%*%alpha)
      
      nhsicRegL <- do.call("nhsicReg", list(model=model, test=curTest, learner=learner, param=p))
      hsicRegL <- do.call("hsicReg", list(model=model, test=curTest, learner=learner, param=p))
      hsicYhRegL <- do.call("hsicYhReg", list(model=model, test=curTest, learner=learner, param=p))
      
      nhsicRegAL <- do.call("nhsicRegA", list(model=model, test=curTest, learner=learner, param=p))
      hsicRegAL <- do.call("hsicRegA", list(model=model, test=curTest, learner=learner, param=p))
      hsicYhRegAL <- do.call("hsicYhRegA", list(model=model, test=curTest, learner=learner, param=p))
      corr <- do.call("corre", list(model=model, test=curTest, learner=learner, param=p))
      
      res <- c(lossTest=lossTest, lossTrain=lossTrain, rmseL=rmseL, qhsicL=qhsicL, hsicL=hsicL, regL=regL, 
               nhsicRegL=nhsicRegL, hsicRegL=hsicRegL, hsicYhRegL=hsicYhRegL, 
               nhsicRegAL=nhsicRegAL, hsicRegAL=hsicRegAL, hsicYhRegAL=hsicYhRegAL, corr=corr)
      return(res)
    }, f=1:fold, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, params, mc.cores=2, SIMPLIFY="array")
  
  
  means <- apply(losses["lossTest",,], 2, mean) #, na.rm=T (add if its ok for some folds to have uninvertible matrices)
  lambdas <- sapply(params, function(x) x$lambda)
  meansFac <- means/min(means, na.rm=T)
  indx <- which(meansFac <= fac)
  maxLambda <- max(lambdas[indx])
  winner <- which(lambdas==maxLambda)
  
  losses <- aperm(losses, c(2,3,1))
  dimnames(losses) <- list(fold=1:fold, params=names(params), var=c("lossTest","lossTrain","rmse","qhsicL","hsicL","regL",
                                                                    "nhsicRegL", "hsicRegL", "hsicYhRegL",
                                                                    "nhsicRegAL", "hsicRegAL", "hsicYhRegAL","corr"))
  
  print("exits CV parallel")
  return(list(opt=params[winner], grid=losses))
}



sse <- function(model, test, learner, param){
	#print("enters sse")
	pred <- learner$predict(model, test)
	
	res <- sum((pred - test$y)^2)
	#print("exits sse")
	return(res)
}
mse <- function(model, test, learner, param){
  #print("enters mse")
  pred <- learner$predict(model, test)
  
  res <- mean((pred - test$y)^2)
  #print("exits mse")
  return(res)
}

rmse <- function(model, test, learner, param){
  #print("enters rmse")
  pred <- learner$predict(model, test)
  
  res <- sqrt(mean((pred - test$y)^2))
  #print("exits rmse")
  return(res)
}


qhsicLoss <- function(model, test, learner, param){

	pred <- learner$predict(model, test)
	#qhsic(x,r) = a'KxKnKxa + yn'Knyn -2a'Kxyn
	
	
	kernelXs <- do.call(param$kernelXs, param$psKernXs)
	Kxs <-  kernelMatrix(kernelXs, test$x)
	kernelXb <- do.call(param$kernelXb, param$psKernXb)
	Kxb <-  kernelMatrix(kernelXb, test$x)
	
	N <- nrow(Kxs)
	H <- diag(N)-matrix(1/N,N,N)
	
	#Kxsc <- H%*%Kxs%*%H
	Kxbc <- H%*%Kxb%*%H
	
	res <- t(pred)%*%Kxbc%*%pred + t(test$y)%*%Kxbc%*%test$y - 2*t(pred)%*%Kxbc%*%test$y
	return(res)
}
nqhsicLoss <- function(model, test, learner, param){

	pred <- learner$predict(model, test)
	#qhsic(x,r) = a'KxKnKxa + yn'Knyn -2a'Kxyn
	resids <- test$y-pred
	
	
	kernelXb <- do.call(param$kernelXb, param$psKernXb)
	Kxb <-  kernelMatrix(kernelXb, test$x)
	N <- nrow(Kxb)
	kernelR <- do.call("vanilladot", list())
	Kr <- kernelMatrix(kernelR, matrix(resids,N,1))
	

	H <- diag(N)-matrix(1/N,N,N)
	
	Kxbc <- H%*%Kxb%*%H
	Krc <- H%*%Kr%*%H
	
	res <- sum(diag(Kxbc%*%Krc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krc%*%Krc))))
	
	
	
	return(res)
}
hsicLoss  <- function(model, test, learner, param){
	
	#print("enters hsicLoss")
	kernelXs <- do.call(param$kernelXs, param$psKernXs)
	kernelXb <- do.call(param$kernelXb, param$psKernXb)
	
	Kxs <-  kernelMatrix(kernelXs, test$x)
	Kxb <-  kernelMatrix(kernelXb, test$x)
		
	pred <- learner$predict(model, test)
	resids <- test$y-pred
	kernelRg <- do.call(param$kernelRg, param$psKernRg)
	Krg <- kernelMatrix(kernelRg, resids)
	N <- nrow(Kxs)
	H <- diag(N)-matrix(1/N,N,N)
	Kxbc <- Kxb%*%H
	Krgc <- Krg%*%H
	res <- sum(diag(Kxbc%*%Krgc)) 
	#print("exits hsicLoss")
	return(res)
}
nhsicLoss  <- function(model, test, learner, param){
  
  #print("enters hsicLoss")
  kernelXs <- do.call(param$kernelXs, param$psKernXs)
  kernelXb <- do.call(param$kernelXb, param$psKernXb)
  
  Kxs <-  kernelMatrix(kernelXs, test$x)
  Kxb <-  kernelMatrix(kernelXb, test$x)
  
  pred <- learner$predict(model, test)
  resids <- test$y-pred
  kernelRg <- do.call(param$kernelRg, param$psKernRg)
  Krg <- kernelMatrix(kernelRg, resids)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Krgc <- H%*%Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krgc%*%Krgc))))
  #print("exits hsicLoss")
  return(res)
}
qkricLoss <- function(model, test, learner, param){
  
  pred <- learner$predict(model, test)
  
  kernelXs <- do.call(param$kernelXs, param$psKernXs)
  Kxs <-  kernelMatrix(kernelXs, test$x)
  kernelXb <- do.call(param$kernelXb, param$psKernXb)
  Kxb <-  kernelMatrix(kernelXb, test$x)
  
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  
  
  Kxbc <- H%*%Kxb%*%H
  
  mseLoss <- sum((pred - test$y)^2)
  qhsicLoss <- (t(pred)%*%Kxbc%*%pred + t(test$y)%*%Kxbc%*%test$y - 2*t(pred)%*%Kxbc%*%test$y)
  mu <- param$mu
  res <- mu*mseLoss + (1-mu)*qhsicLoss
  
  return(res)
}

# Fair Learning / Consisent Regression Regularizer
nhsicReg  <- function(model, test, learner, param){
  
  #print("enters nhsicReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
  if(!is.null(model$indxSens)){
    indxSens <- model$indxSens  
  } else{
    indxSens <- 1:p
  }
  
  
  
  x <- matrix(test$x, N, p)
  
  
  kernelXb <- do.call("rbfdot", param$psKernXb)
  Kxb <-  kernelMatrix(kernelXb, x[,indxSens])
  
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("rbfdot", param$psKernYhg)
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhgc <- H%*%Kyhg%*%H
  res <- sum(diag(Kxbc%*%Kyhgc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Kyhgc%*%Kyhgc))))
  #print("exits nhsicReg")
  return(res)
}
hsicReg  <- function(model, test, learner, param){
  
  #print("enters hsicReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
  if(!is.null(model$indxSens)){
    indxSens <- model$indxSens  
  } else{
    indxSens <- 1:p
  }
  
  x <- matrix(test$x, N, p)
  
  #kernelXs <- do.call(param$kernelXs, param$psKernXs)
  kernelXb <- do.call("rbfdot", param$psKernXb)
  
  #Kxs <-  kernelMatrix(kernelXs, x)
  Kxb <-  kernelMatrix(kernelXb, x[,indxSens])
  
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("rbfdot", param$psKernYhg)
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhgc <- H%*%Kyhg%*%H
  res <- sum(diag(Kxbc%*%Kyhgc))
  #print("exits hsicReg")
  return(res)
}
hsicYhReg  <- function(model, test, learner, param){
  
  #print("enters hsicYhReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
  #if(!is.null(model$indxSens)){
  #  indxSens <- model$indxSens  
  #} else{
  #  indxSens <- 1:p
  #}
  
  x <- matrix(test$x, N, p)
  
  #kernelXs <- do.call(param$kernelXs, param$psKernXs)
  #kernelXb <- do.call("rbfdot", param$psKernXb)
  
  #Kxs <-  kernelMatrix(kernelXs, x)
  #Kxb <-  kernelMatrix(kernelXb, x[,indxSens])
  
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("rbfdot", param$psKernYhg)
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kyhg)
  H <- diag(N)-matrix(1/N,N,N)
  Kyhgc <- H%*%Kyhg%*%H
  res <- sqrt(sum(diag(Kyhgc%*%Kyhgc)))
  #print("exits hsicYhReg")
  return(res)
}

# actual regularizer used
nhsicRegA  <- function(model, test, learner, param){
  
  #print("enters nhsicReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
  if(!is.null(model$indxSens)){
    indxSens <- model$indxSens  
  } else{
    indxSens <- 1:p
  }
  
  x <- matrix(test$x, N, p)
  
  if(class(model$kernelXb)=="vanillakernel"){
    kernelXb <- do.call(vanilladot, list())
    } else{
    kernelXb <- model$kernelXb
    }
  
  
  Kxb <-  kernelMatrix(kernelXb, as.matrix(x[,indxSens]))
  
  
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("vanilladot", list())
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhgc <- H%*%Kyhg%*%H
  res <- sum(diag(Kxbc%*%Kyhgc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Kyhgc%*%Kyhgc))))
  #print("exits nhsicReg")
  return(res)
}
hsicRegA  <- function(model, test, learner, param){
  
  #print("enters hsicReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
  if(!is.null(model$indxSens)){
    indxSens <- model$indxSens  
  } else{
    indxSens <- 1:p
  }
  
  x <- matrix(test$x, N, p)
  
  
  if(class(model$kernelXb)=="vanillakernel"){
    kernelXb <- do.call(vanilladot, list())
  } else{
    kernelXb <- model$kernelXb  
  }
  
  
  Kxb <-  kernelMatrix(kernelXb, as.matrix(x[,indxSens]))
  
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("vanilladot", list())
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Kyhgc <- H%*%Kyhg%*%H
  res <- sum(diag(Kxbc%*%Kyhgc))
  #print("exits hsicReg")
  return(res)
}
hsicYhRegA  <- function(model, test, learner, param){
  
  #print("enters hsicYhReg")
  if(!is.null(dim(test$x))){
    p <- dim(test$x)[2]
    N <- dim(test$x)[1]
  } else{
    N <- length(test$x)
    p <- 1
  }
  
 
  
  x <- matrix(test$x, N, p)
  
 
  pred <- learner$predict(model, test)
  yh <- pred
  kernelYhg <- do.call("vanilladot", list())
  Kyhg <- kernelMatrix(kernelYhg, yh)
  N <- nrow(Kyhg)
  H <- diag(N)-matrix(1/N,N,N)
  Kyhgc <- H%*%Kyhg%*%H
  res <- sqrt(sum(diag(Kyhgc%*%Kyhgc)))
  #print("exits hsicYhReg")
  return(res)
}



corre  <- function(model, test, learner, param){
  
  #print("enters corre")
  
  pred <- learner$predict(model, test)
 
  res <- cor(pred, test$y)
  #print("exits corre")
  return(res)
}


getFixedParams <- function(data, plot=FALSE){
	#print("enters getSigmaR")
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
  } else{
    n <- length(data$x)
  }
  if(n > 500){
    data <- getSubset(data, 1:500) 
  }
  
	sigma0 <- 1/median(as.numeric(dist(data$x)^2))
	beta0 <- sigma0
	
	#obtain residuals with a mars model
	model  <- earth(x=data$x, y=data$y, nfold=5)
	res0 <- data$y-as.numeric(predict(model, data$x))
	gamma0 <- 1/median(as.numeric(dist(res0)^2))
	
	# fit sigma
	ord <- 10
	sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
	varsHsics.sig <- sapply(sigmas1, function(sd1){
	  kernelXb <- do.call("rbfdot", list(sigma=sd1))
	  Kxb <- kernelMatrix(kernelXb, data$x)
	  N <- nrow(Kxb)
	  H <- diag(N)-matrix(1/N,N,N)
	  Kxbc <- H%*%Kxb%*%H
	  distsX <- (Kxbc)[lower.tri(Kxbc)]
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
	  Kxs <- kernelMatrix(kernelXs, data$x)
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
	
	
	# fit beta	
	betas1 <- (10^seq(-ord,ord,1)) # *beta0
	varsHsics.bet <- sapply(betas1, function(sd1){
		kernelXb <- do.call("rbfdot", list(sigma=sd1))
		Kxb <- kernelMatrix(kernelXb, data$x)
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
	
	#fit gamma
	
	gammas1 <- (10^seq(-ord,ord,1)) # *gamma0
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

	print(paste("sigmaVar: ", sigmaVar, sep=""))
	print(paste("sigmaSat: ", sigmaSat, sep=""))
	print(paste("simgaMed: ", sigma0, sep=""))
	print(paste("betaVar: ", betaVar, sep=""))
	print(paste("betaMed: ", beta0, sep=""))
	print(paste("gammaVar: ", gammaVar, sep=""))
	print(paste("gammaMed: ", gamma0, sep=""))
	
	
	if(plot){
	  
	  # sigma
	  
	  # beta
	  dhsics.bet <- sapply(betas1, function(sd1){
	    kernelXb <- do.call("rbfdot", list(sigma=sd1))
	    kernelR <- do.call("vanilladot", list())
	    Kxb <- kernelMatrix(kernelXb, data$x)
	    Kr <- kernelMatrix(kernelR, matrix(res0,length(data$y),1))
	    Ks <- vector("list", 2)
	    Ks[[1]] <- Kxb 
	    Ks[[2]] <- Kr
	    Kxx <- vector("list", 2)
	    Kxx[[1]] <- Kxb 
	    Kxx[[2]] <- Kxb 
	    Krr <- vector("list", 2)
	    Krr[[1]] <- Kr
	    Krr[[2]] <- Kr
	    
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
	  dhsics.bet_fac <- dhsics.bet/max(dhsics.bet[which(dhsics.bet<Inf)])
	  indxSat <- which(dhsics.bet_fac > 0.95)
	  aux <- which.min(dhsics.bet_fac[indxSat])
	  indxSat <- indxSat[aux]
	  betaSat <- betas1[indxSat]
	  
	  # gamma
	  dhsics.gam <- sapply(gammas1, function(sd1){
	    kernelX <- do.call("vanilladot", list())
	    kernelRg <- do.call("rbfdot", list(sigma=sd1))
	    Kx <- kernelMatrix(kernelX, matrix(data$x, length(data$y),1))
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
	    kernelRg <- do.call("rbfdot", list(sigma=sd2))
	    Kxb <- kernelMatrix(kernelXb, data$x)
	    Krg <- kernelMatrix(kernelRg, matrix(res0, length(data$y),1))
	    
	    Ks <- vector("list", 2)
	    Ks[[1]] <- Kxb 
	    Ks[[2]] <- Krg
	    Kxx <- vector("list", 2)
	    Kxx[[1]] <- Kxb 
	    Kxx[[2]] <- Kxb 
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
	
		
	#print("exits getSigmaR")
	return(list(kernelXs="rbfdot", kernelXb="rbfdot", kernelRg="rbfdot", psKernXs=list(sigma=sigmaSat), psKernXb=list(sigma=betaVar), psKernRg=list(sigma=gammaVar)))
}
getFixedParams2 <- function(data){
	#print("enters getSigmaR")
	sigma0 <- 1/median(as.numeric(dist(data$x)^2))
	beta0 <- sigma0
	
	
	# fit sigma
	ord <- 10
	sigmas1 <- sigma0*(10^seq(-ord,ord,1))
	dhsics <- sapply(sigmas1, function(sd1){
		kernelXs <- do.call("rbfdot", list(sigma=sd1))
		kernelY <- do.call("vanilladot", list())
		Kxs <- kernelMatrix(kernelXs, data$x)
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
	dhsics_fac <- dhsics/max(dhsics[which(dhsics<Inf)])
	indxSat <- which(dhsics_fac > 0.998)
	aux <- which.max(sigmas1[indxSat])
	indxSat <- indxSat[aux]
	sigmaSat <- sigmas1[indx]
	varsHsics <- sapply(sigmas1, function(sd1){
		kernelXb <- do.call("rbfdot", list(sigma=sd1))
		Kxb <- kernelMatrix(kernelXb, data$x)
		N <- nrow(Kxb)
		H <- diag(N)-matrix(1/N,N,N)
		Kxbc <- H%*%Kxb%*%H
		distsX <- (Kxbc)[lower.tri(Kxbc)]
		res <- var(distsX)
		#plot(distsX2, distsX, main=res, ylim=c(-0.6,0.6))
		#hist(distsX,100, main=res)
		return(res)
	})
	indxMaxVar <- which.max(varsHsics)
	sigmaVar <- sigmas1[indxMaxVar]
	
	par(mfrow=c(2,2))
	plot(log(sigmas1,10), dhsics, main="fit sigma", xlab="log(sigma)", ylab="hsic/var", ylim=range(dhsics, varsHsics))
	lines(log(sigmas1,10), varsHsics, col="red", type="p")
	abline(v=log(sigma0,10), col="red")
	abline(v=log(sigmaSat,10), col="green")
	abline(v=log(sigmaVar,10), col="blue")
	print(paste("sigmaSat: ", sigmaSat, sep=""))
	print(paste("sigmaVar: ", sigmaVar, sep=""))
	print(paste("median precision: ", sigma0, sep=""))
	
	
	# fit beta	
	betas1 <- beta0*(10^seq(-ord,ord,1))
	dhsics <- sapply(betas1, function(sd1){
		kernelXb <- do.call("rbfdot", list(sigma=sd1))
		kernelR <- do.call("vanilladot", list())
		Kxb <- kernelMatrix(kernelXb, data$x)
		Kr <- kernelMatrix(kernelR, matrix(res1,length(data$y),1))
		Ks <- vector("list", 2)
		Ks[[1]] <- Kxb 
		Ks[[2]] <- Kr
		Kxx <- vector("list", 2)
		Kxx[[1]] <- Kxb 
		Kxx[[2]] <- Kxb 
		Krr <- vector("list", 2)
		Krr[[1]] <- Kr
		Krr[[2]] <- Kr
		
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
	dhsics_fac <- dhsics/max(dhsics[which(dhsics<Inf)])
	indxSat <- which(dhsics_fac > 0.95)
	aux <- which.min(dhsics_fac[indxSat])
	indxSat <- indxSat[aux]
	betaSat <- betas1[indxSat]
	varsHsics <- sapply(betas1, function(sd1){
		kernelXb <- do.call("rbfdot", list(sigma=sd1))
		Kxb <- kernelMatrix(kernelXb, data$x)
		N <- nrow(Kxb)
		H <- diag(N)-matrix(1/N,N,N)
		Kxbc <- H%*%Kxb%*%H
		distsX <- (Kxbc)[lower.tri(Kxbc)]
		res <- var(distsX)
		#plot(distsX2, distsX, main=res, ylim=c(-0.6,0.6))
		#hist(distsX,100, main=res)
		return(res)
	})
	indxMaxVar <- which.max(varsHsics)
	betaVar <- betas1[indxMaxVar]
	plot(log(betas1,10), dhsics, main="fit beta", xlab="log(betas)", ylab="hsic/var", ylim=range(dhsics, varHsics))
	lines(log(betas1,10), varsHsics, col="red", type="p")
	abline(v=log(beta0,10), col="red")
	abline(v=log(betaSat,10), col="green")
	abline(v=log(betaVar,10), col="blue")
	print(paste("betaSat: ", betaSat, sep=""))
	print(paste("betaVar: ", betaVar, sep=""))
	print(paste("median precision: ", beta0, sep=""))
	
	
	par(mfrow=c(1,1))
		
	#print("exits getSigmaR")
	return(list(kernelXs="rbfdot", kernelXb="rbfdot", kernelRg="vanilladot", psKernXs=list(sigma=sigmaSat), psKernXb=list(sigma=betaVar), psKernRg=list()))
}
getFixedParamsCR <- function(data, indxSens=NULL, indxPred=NULL, plot=FALSE){
  #print("enters getSigmaR")
  
  
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
  gamma0 <- 1/median(as.numeric(dist(yh0)^2))
  
  res0 <- data$y-as.numeric(predict(model, xx))
  gamma0r <- 1/median(as.numeric(dist(res0)^2))
  
  
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
  
  #fit gamma
  
  gammas1 <- (10^seq(-ord,ord,1)) # *gamma0
  varsHsics.gam <- sapply(gammas1, function(sd1){
    kernelYhg <- do.call("rbfdot", list(sigma=sd1))
    Kyhg <- kernelMatrix(kernelYhg, yh0)
    N <- nrow(Kyhg)
    H <- diag(N)-matrix(1/N,N,N)
    Kyhgc <- H%*%Kyhg%*%H
    distsYh <- (Kyhgc)[lower.tri(Kyhgc)]
    res <- var(distsYh)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.gam)
  gammaVar <- gammas1[indxMaxVar]
  
  #fit gammaR
  
  gammas1 <- (10^seq(-ord,ord,1)) # *gamma0r
  varsHsics.gamR <- sapply(gammas1, function(sd1){
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
  gammaRVar <- gammas1[indxMaxVar]
  
  print(paste("sigmaVar: ", sigmaVar, sep=""))
  print(paste("sigmaSat: ", sigmaSat, sep=""))
  print(paste("simgaMed: ", sigma0, sep=""))
  print(paste("betaVar: ", betaVar, sep=""))
  print(paste("betaMed: ", beta0, sep=""))
  print(paste("gammaVar: ", gammaVar, sep=""))
  print(paste("gammaMed: ", gamma0, sep=""))
  print(paste("gammaRVar: ", gammaRVar, sep=""))
  print(paste("gammaRMed: ", gamma0r, sep=""))
  
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
    
    # gamma
    dhsics.gam <- sapply(gammas1, function(sd1){
      kernelX <- do.call("vanilladot", list())
      kernelRg <- do.call("rbfdot", list(sigma=sd1))
      Kx <- kernelMatrix(kernelX, matrix(xs, length(data$y),1))
      Kyhg <- kernelMatrix(kernelYh, matrix(yh0,length(data$y),1))
      Ks <- vector("list", 2)
      Ks[[1]] <- Kx 
      Ks[[2]] <- Kyhg
      Kxx <- vector("list", 2)
      Kxx[[1]] <- Kx 
      Kxx[[2]] <- Kx 
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
    })
    dhsics.gam_fac <- dhsics.gam/max(dhsics.gam[which(dhsics.gam<Inf)])
    indxSat <- which(dhsics.gam_fac > 0.95)
    aux <- which.min(dhsics.gam_fac[indxSat])
    indxSat <- indxSat[aux]
    gammaSat <- gammas1[indxSat]
    
    # gammaR
    dhsics.gamR <- sapply(gammas1, function(sd1){
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
    dhsics.gamR_fac <- dhsics.gamR/max(dhsics.gamR[which(dhsics.gamR<Inf)])
    indxSat <- which(dhsics.gamR_fac > 0.95)
    aux <- which.min(dhsics.gamR_fac[indxSat])
    indxSat <- indxSat[aux]
    gammaRSat <- gammas1[indxSat]
    
    
    
    print(paste("betaSat: ", betaSat, sep=""))
    print(paste("gammaSat: ", gammaSat, sep=""))
    print(paste("gammaRSat: ", gammaRSat, sep=""))
    
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
  
  
  #print("exits getSigmaR")
  return(list(kernelXs="rbfdot", kernelXb="rbfdot", kernelRg="rbfdot", kernelYhg="rbfdot", psKernXs=list(sigma=sigma0), 
              psKernXb=list(sigma=sigma0), psKernYhg=list(sigma=gammaVar), psKernRg=list(sigma=gammaVar), indxPred=list(indxPred), indxSens=list(indxSens)))
}


krr <- constructLearner(learn.krr, predict.krr, sse, getFixedParamsCR)
qhsic <- constructLearner(learn.qhsic, predict.qhsic, qhsicLoss, getFixedParamsCR)
hsic <- constructLearner(learn.hsic, predict.hsic, hsicLoss, getFixedParamsCR)
qkric <- constructLearner(learn.qkric, predict.qkric, qkricLoss, getFixedParamsCR)
hsic2 <- constructLearner(learn.hsic2, predict.hsic, hsicLoss, getFixedParamsCR)
hsic3 <- constructLearner(learn.hsic3, predict.hsic, hsicLoss, getFixedParamsCR)
hsic4 <- constructLearner(learn.hsic4, predict.hsic, hsicLoss, getFixedParamsCR)



plot.optLambda <- function(opt){
	
	indx <- which.min(opt$grid[,"lossTest"])
	par(mfrow=c(3,2))
	plot(opt$lambda, opt$grid[,"lossTest"], ylim=range(opt$grid[,c("lossTest","lossTrain")],na.rm=T), type="b", xlab="lambda", ylab="", main="opt. loss")
	lines(opt$lambda, opt$grid[,"lossTrain"], type="b", col="red")
	lines(opt$lambda[indx], opt$grid[indx, "lossTest"], col="red" , cex=1.5, type="p")
	plot(opt$lambda, opt$grid[,"rmse"], type="b", main="root mean square error", xlab="lambda", ylab="")
	lines(opt$lambda[indx], opt$grid[indx, "rmse"], col="red" , cex=1.5, type="p")
	plot(opt$lambda, opt$grid[,"qhsicL"], type="b", main="quasi hsic", xlab="lambda", ylab="")
	lines(opt$lambda[indx], opt$grid[indx, "qhsicL"], col="red" , cex=1.5, type="p")
	plot(opt$lambda, opt$grid[,"hsicL"], type="b", main="hsic", xlab="lambda", ylab="")
	lines(opt$lambda[indx], opt$grid[indx, "hsicL"], col="red" , cex=1.5, type="p")
	plot(opt$lambda, opt$grid[,"regL"], type="b", main="reg", xlab="lambda", ylab="")
	lines(opt$lambda[indx], opt$grid[indx, "regL"], col="red" , cex=1.5, type="p")
	par(mfrow=c(1,1))
}

plot.optLambdaLog <- function(opt){
  
  indx <- which.min(opt$grid[,"lossTest"])
  par(mfrow=c(3,2))
  plot(log(opt$lambda, 10), opt$grid[,"lossTest"], ylim=range(opt$grid[,c("lossTest","lossTrain")],na.rm=T), type="b", xlab="lambda", ylab="", main="opt. loss")
  lines(log(opt$lambda,10), opt$grid[,"lossTrain"], type="b", col="red")
  lines(log(opt$lambda[indx],10), opt$grid[indx, "lossTest"], col="red" , cex=1.5, type="p")
  plot(log(opt$lambda,10), opt$grid[,"rmse"], type="b", main="root mean square error", xlab="lambda", ylab="")
  lines(log(opt$lambda[indx],10), opt$grid[indx, "rmse"], col="red" , cex=1.5, type="p")
  plot(log(opt$lambda,10), opt$grid[,"qhsicL"], type="b", main="quasi hsic", xlab="lambda", ylab="")
  lines(log(opt$lambda[indx],10), opt$grid[indx, "qhsicL"], col="red" , cex=1.5, type="p")
  plot(log(opt$lambda,10), opt$grid[,"hsicL"], type="b", main="hsic", xlab="lambda", ylab="")
  lines(log(opt$lambda[indx],10), opt$grid[indx, "hsicL"], col="red" , cex=1.5, type="p")
  plot(log(opt$lambda,10), opt$grid[,"regL"], type="b", main="reg", xlab="lambda", ylab="")
  lines(log(opt$lambda[indx],10), opt$grid[indx, "regL"], col="red" , cex=1.5, type="p")
  par(mfrow=c(1,1))
}


optLambda <- function(trainData, learner, numFolds, parallel=FALSE, plot=FALSE, fixedParams=NULL, fac=1.25){
	print("enters optLambda")
	
	if(is.null(fixedParams)) fixedParams <- do.call(learner$getFixedParams, list(data=trainData))
	
	ord <- 10
	lambdas1 <- 10^seq(-ord,ord,1)
	
	params <- fixedParams
	params$lambda <- lambdas1
	params <- do.call(constructParams, params)
	
	if(parallel){
	  opt <- CV.parallel(data=trainData, learner, params, fold=numFolds, fac=fac)
	} else{
	  opt <- CV(data=trainData, learner, params, fold=numFolds, fac=fac)
	}
	
	#par(mfrow=c(2,1))
	plot(log(lambdas1,10), opt$grid[1,,"lossTest"], ylim=range(opt$grid[,,"lossTest"], na.rm=T), main="first CV", xlab="log(lambda)", ylab="test loss")
	for(i in 2:(dim(opt$grid)[1])) lines(log(lambdas1,10), opt$grid[i,,"lossTest"], type="p")
	lines(log(lambdas1,10), apply(opt$grid[,,"lossTest"],2,mean, na.rm=T), col="red", type="b")
	abline(v=log(opt$opt[[1]]$lambda,10), col="red")
	
	
	lambdas1 <- c(10^(-ord-1), lambdas1, 10^(ord+1))
	indx <- match(opt$opt[[1]]$lambda, lambdas1)
	lambdas2 <- seq(lambdas1[indx-1],lambdas1[indx+1], length.out=50)
	params <- fixedParams
	
	params$lambda <- lambdas2
	params <- do.call(constructParams, params)
	
	if(parallel){
	  opt <- CV.parallel(trainData, learner, params, fold=numFolds, fac=1)
	} else{
	  opt <- CV(trainData, learner, params, fold=numFolds, fac=1)
	}
	plot(lambdas2, opt$grid[1,,"lossTest"], ylim=range(opt$grid[,,"lossTest"], na.rm=T), main="second CV", xlab="lambda", ylab="test loss")
	for(i in 2:(dim(opt$grid)[1])) lines(lambdas2, opt$grid[i,,"lossTest"], type="p")
	lines(lambdas2, apply(opt$grid[,,"lossTest"],2,mean, na.rm=T), col="red", type="b")
	abline(v=opt$opt[[1]]$lambda, col="red")
	par(mfrow=c(1,1))
	
	
	
	grid <- opt$grid
	dimnames(grid) <- list(fold=1:numFolds, lambda=lambdas2, var=c("lossTest","lossTrain","rmse","qhsicL","hsicL","regL", 
	                                                               "nhsicRegL", "hsicRegL", "hsicYhRegL","corre"))
	grid <- apply(grid, c(2,3), mean, na.rm=T)
	opt$grid <- grid
	
	opt$lambdas <- lambdas2
	
	if(plot) plot.optLambda(opt)
	
	print("exits optLambda")
	return(opt)
}

optLambdaLog <- function(trainData, learner, numFolds, parallel=FALSE, plot=FALSE, fixedParams=NULL, fac=1.25){
  print("enters optLambdaLog")
  
  if(is.null(fixedParams)) fixedParams <- do.call(learner$getFixedParams, list(data=trainData))
  
  ord <- 10
  #lambdas1 <- 10^seq(-ord,ord,1)
  lambdas1 <- c(0, 10^seq(-8,8, length.out=50))
  
  params <- fixedParams
  params$lambda <- lambdas1
  params <- do.call(constructParams, params)
  
  if(parallel){
    opt <- CV.parallel(data=trainData, learner, params, fold=numFolds, fac=fac)
  } else{
    opt <- CV(data=trainData, learner, params, fold=numFolds, fac=fac)
  }
  
  #par(mfrow=c(2,1))
  plot(log(lambdas1,10), opt$grid[1,,"lossTest"], ylim=range(opt$grid[,,"lossTest"], na.rm=T), main="first CV", xlab="log(lambda)", ylab="test loss")
  for(i in 2:(dim(opt$grid)[1])) lines(log(lambdas1,10), opt$grid[i,,"lossTest"], type="p")
  lines(log(lambdas1,10), apply(opt$grid[,,"lossTest"],2,mean, na.rm=T), col="red", type="b")
  abline(v=log(opt$opt[[1]]$lambda,10), col="red")

  
  grid <- opt$grid
  dimnames(grid) <- list(fold=1:numFolds, lambda=lambdas1, var=c("lossTest","lossTrain","rmse","qhsicL","hsicL","regL", 
                                                                 "nhsicRegL", "hsicRegL", "hsicYhRegL",
                                                                 "nhsicRegAL", "hsicRegAL", "hsicYhRegAL","corre"))
  grid <- apply(grid, c(2,3), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  opt$grid <- grid
  
  opt$lambdas <- lambdas1
  
  if(plot) plot.optLambdaLog(opt)
  
  print("exits optLambdaLog")
  return(opt)
}



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


########################################################################################
# Functions for applying hsic regression to causality
########################################################################################

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

# simulate from a given sem and noise
simSEM <- function(sem, ns){
	if(all(colnames(ns)!=nodes(sem$dag))) stop("names (sem,ns) dont match")
	if(all(names(sem$funcs)!=paste("f",nodes(sem$dag),sep=""))) stop("names (nodes, funcs) dont match")
	
	dagAmat <- amat(as.bn(sem$dag))
	n <- nrow(ns)
	Nodes <- nodes(sem$dag)
	topoNodes <- topoSort(sem$dag)
	p <- length(Nodes)
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
	
	return(x)
}

Residuals <- function(Xtr, Ytr, Xte, Yte, learner, numFolds, fac=1.25){
  nsTr  <- constructData(x=Xtr,y=Ytr)
  nsTe  <- constructData(x=Xte, y=Yte)
  
  #optL  <- optLambda(trainData=nsTr, learner,   numFolds, parallel=TRUE, fac=fac)
  optL  <- optLambdaLog(trainData=nsTr, learner,   numFolds, parallel=TRUE, fac=fac)
  
  model <- try(learner$learn(nsTr, optL[[1]][[1]]))
  while(class(model)=="try-error"){
    print("cannot invert full data matrix for given lambda: scaling it up by factor of 10")
    optL[[1]][[1]]$lambda <- optL[[1]][[1]]$lambda*10
    model <- try(learner$learn(nsTr, optL[[1]][[1]]))
  } 
  
  rsTr  <- nsTr$y-learner$predict(model, nsTr)
  rsTe  <- nsTe$y-learner$predict(model, nsTe)
  
  return(list(rsTr=rsTr, rsTe=rsTe))
}


# given data and a graph, estimate corresponding non-linear additive SEM returning residuals
getResiduals <- function(G, xTrain, xTest, learner, numFolds){
	if(all(colnames(xTrain)!=nodes(G))) stop("names (G,xTrain) dont match")
	if(all(colnames(xTrain)!=colnames(xTest))) stop("names (xTrain, xTest) dont match")
	
	dagAmat <- amat(as.bn(G))
	Nodes <- nodes(G)
	topoNodes <- topoSort(G)
	p <- length(Nodes)
	nTrain <- nrow(xTrain)
	nTest <- nrow(xTest)
	rsTrain <- matrix(NA, nTrain, p)
	rsTest <- matrix(NA, nTest, p)
	colnames(rsTrain) <- Nodes
	colnames(rsTest) <- Nodes
	indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
	rsTrain[,Nodes[indxNoParents]] <- xTrain[,Nodes[indxNoParents]]
	rsTest[,Nodes[indxNoParents]] <- xTest[,Nodes[indxNoParents]]
	
	print(paste("no parent nodes: ", Nodes[indxNoParents],sep=""))
	for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
		
		print("*********************************")
		print(paste("node: ", nd,sep=""))
		indxParents <- which(dagAmat[,nd]==1)
		Parents <- Nodes[indxParents]
		
		
		rs <- Residuals(Xtr=xTrain[,Parents], Ytr=xTrain[,nd], Xte=xTest[,Parents] , Yte=xTest[,nd], 
		                learner=learner, numFolds=numFolds)
		
		rsTrain[,nd] <- rs$rsTr
		rsTest[,nd]  <- rs$rsTe

	}
	
	return(list(train=rsTrain, test=rsTest))
	
}

getResidualsMat <- function(dagAmat, xTrain, xTest, learner, numFolds){
	if(all(colnames(xTrain)!=colnames(dagAmat))) stop("names (dagAmat,xTrain) dont match")
	if(all(colnames(xTrain)!=colnames(xTest))) stop("names (xTrain, xTest) dont match")
	
	
	Nodes <- colnames(dagAmat)
	topoNodes <- topoSort(dagAmat)
	p <- length(Nodes)
	nTrain <- nrow(xTrain)
	nTest <- nrow(xTest)
	rsTrain <- matrix(NA, nTrain, p)
	rsTest <- matrix(NA, nTest, p)
	colnames(rsTrain) <- Nodes
	colnames(rsTest) <- Nodes
	indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
	rsTrain[,Nodes[indxNoParents]] <- xTrain[,Nodes[indxNoParents]]
	rsTest[,Nodes[indxNoParents]] <- xTest[,Nodes[indxNoParents]]
	
	print(paste("no parent nodes: ", Nodes[indxNoParents],sep=""))
	for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
		
		print("*********************************")
		print(paste("node: ", nd,sep=""))
		indxParents <- which(dagAmat[,nd]==1)
		Parents <- Nodes[indxParents]
		
		
		rs <- Residuals(Xtr=xTrain[,Parents], Ytr=xTrain[,nd], Xte=xTest[,Parents] , Yte=xTest[,nd], 
		                learner=learner, numFolds=numFolds)
		
		rsTrain[,nd] <- rs$rsTr
		rsTest[,nd]  <- rs$rsTe
	
	}
	
	return(list(train=rsTrain, test=rsTest))
	
}


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
            bandwidth <- median_bandwidth_rcpp(x[sample(1:len), , drop = FALSE], len, ncol(x))
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
                K[[j]] <- gaussian_grammat_rcpp(X[[j]], bandwidth[j], len, ncol(X[[j]]))
            }
            else if (kernel[j] == "gaussian.fixed") {
                K[[j]] <- gaussian_grammat_rcpp(X[[j]], bandwidth[j], len, ncol(X[[j]]))
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
          K_perm <- shuffle_grammat_rcpp(K[[j]], perm, 
                                         len)
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
          K_boot <- shuffle_grammat_rcpp(K[[j]], boot, 
                                         len)
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
    res <- unique(dags[,slice,], MARGIN=2)
    dimnames(res) <- list(nodeFrom=dimnames(res)[[1]], numReg=1:(dim(res)[2]))
    return(res)
    })
  names(uniqueRegsList) <-c(dimnames(dags)[[2]])
  return(uniqueRegsList)
}


# obtains the residuals for all possible unique regressions implied in all m-node dags (output of getUniqueRegs)
getUniqueResids <- function(uniqueRegs, xTrain, xTest, learner, numFolds){
	dims <- dimnames(uniqueRegs)[-1]
	n <- dim(xTrain)[1]
	dims$sim <- 1:n
	dims$set <- c("train","test")
	
	nodes <- dimnames(uniqueRegs)[[1]]
	numRegs <- dimnames(uniqueRegs)[[2]]
	numRegsTot <- (length(numRegs)-1)*length(nodes)
	pm0 <- proc.time()
	resids <- sapply(nodes,  function(nodeTo){
		print("*********************")
		print(paste(nodeTo, " regressions:", sep=""))
		 sapply(numRegs, function(numReg){
			 print("******")
			 print(paste(nodeTo,", reg # ", numReg, sep=""))
			 indxPreds <- which(uniqueRegs[,numReg,nodeTo]==1)
			 	if(length(indxPreds)>0){
					
			 	  
					rs <- Residuals(Xtr=xTrain[,nodes[indxPreds]], Ytr=xTrain[,nodeTo], 
					                Xte=xTest[,nodes[indxPreds]] , Yte=xTest[,nodeTo], 
					                learner=learner, numFolds=numFolds)
					
					rs.train <- rs$rsTr
					rs.test  <- rs$rsTe
					
					
					print("estimated time to completion:")
					numRegsLeft <- ((length(nodes)-match(nodeTo, nodes))*(length(numRegs)-1)+(length(numRegs)-match(numReg,numRegs)))
					numRegsDone <- numRegsTot -numRegsLeft
					print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
				} else{
					rs.train     <- xTrain[,nodeTo]
					rs.test      <- xTest[,nodeTo]
				}
				return(cbind(rs.train, rs.test))
			}, simplify="array")
		}, simplify="array")
	print(proc.time()-pm0) # 30 mins for krr, 
	dimnames(resids) <- list(sim=1:n, set=c("train","test"), numReg=dimnames(uniqueRegs)$numReg, nodes)
	resids <- aperm(resids, c(1,3,4,2))
	return(resids)
}

# obtains the residuals for unique regressions implied in candidate dagas (output of getUniqueRegsList)
getUniqueResidsList <- function(uniqueRegsList, xTrain, xTest, learner, numFolds){
  
  
  n <- dim(xTrain)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims$set <- c("train","test")
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  
  
  
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  count <- 0
  pm0 <- proc.time()
  resids <- lapply(nodes,  function(nodeTo){
    print("*********************")
    print(paste(nodeTo, " regressions:", sep=""))
    #nodeTo <- "x"
    res <- sapply(numRegs[[nodeTo]], function(numReg){
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>0){
        count <- count + 1
        
        rs <- Residuals(Xtr=xTrain[,nodes[indxPreds]], Ytr=xTrain[,nodeTo], 
                        Xte=xTest[,nodes[indxPreds]] , Yte=xTest[,nodeTo], 
                        learner=learner, numFolds=numFolds)
        
        rs.train <- rs$rsTr
        rs.test  <- rs$rsTe
        
        print("estimated time to completion:")
        numRegsLeft <- numRegsTot - count
        numRegsDone <- numRegsTot -numRegsLeft
        print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        rs.train     <- xTrain[,nodeTo]
        rs.test      <- xTest[,nodeTo]
      }
      return(cbind(rs.train, rs.test))
    }, simplify="array")
    res <- aperm(res, c(1,3,2))
    dimnames(res) <- dims[[nodeTo]]
    return(res)
    })
  print(proc.time()-pm0) #  
  
  
  names(resids) <- nodes
  
  return(resids)
}

# obtains the p-values of fitted residuals in all m-mode dags
getPvals <- function(allDags, uniqueRegs, uniqueResids){
	
	#indxDag <- 543
	#allDags[,,indxDag]
	

	nodes <- dimnames(allDags)[[1]]
	numDags <- dim(allDags)[3]
	count <- 0
	pm0 <- proc.time()
	pvals <- apply(allDags, "dag", function(dag){
		#get residuals for the four regressions implied in each column
		count <<- count + 1
		print("*******************")
		print(paste("dag # ",count, sep=""))
		residsDag <- sapply(nodes, function(nodeTo){
			children <- dag[,nodeTo]
			indxReg <- which(apply(uniqueRegs[,,nodeTo],2, function(col) all(col==children)))
			rsNode <- uniqueResids[,indxReg, nodeTo,]
			return(rsNode)
		}, simplify="array")
		
		residsDag <- aperm(residsDag, c(1,3,2))
		dims <- dimnames(residsDag)
		dims <- list(sim=dims$sim, node=dims[[2]], set=dims$set)
		dimnames(residsDag) <- dims
		pval <- apply(residsDag,"set", function(resids) dhsic.test(resids, matrix.input=TRUE, pairwise=FALSE)$p.value)
		print("estimated time until completion")
		print((proc.time()-pm0)*(numDags-count)/count)
		return(pval)
	})
	print(proc.time()-pm0) #
	

	dims <- dimnames(pvals)
	dims <- list(set=dims[[1]], dag=dims$dag)
	dimnames(pvals) <- dims
	pvals <- aperm(pvals, c(2,1))
	
	return(pvals)
}


# obtains the p-values of fitted residuals in candidate dags
getPvalsList <- function(dags, uniqueRegsList, uniqueResidsList){
  
  #indxDag <- 543
  #allDags[,,indxDag]
  
  
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  pvals <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,2]
    count <<- count + 1
    print("*******************")
    print(paste("dag # ",count, sep=""))
    residsDag <- sapply(nodes, function(nodeTo){
      # nodeTo <- "x"
      children <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegsList[[nodeTo]],2, function(col) all(col==children)))
      rsNode <- uniqueResidsList[[nodeTo]][,indxReg,]
      return(rsNode)
    }, simplify="array")
    
    residsDag <- aperm(residsDag, c(1,3,2))
    dims <- dimnames(residsDag)
    dims <- list(sim=dims$sim, node=dims[[2]], set=dims$set)
    dimnames(residsDag) <- dims
    pval <- apply(residsDag,"set", function(resids) dhsic.test(resids, matrix.input=TRUE, pairwise=FALSE)$p.value)
    print("estimated time until completion")
    print((proc.time()-pm0)*(numDags-count)/count)
    return(pval)
  })
  print(proc.time()-pm0) #
  
  
  dims <- dimnames(pvals)
  dims <- list(set=dims[[1]], dag=dims$dag)
  dimnames(pvals) <- dims
  pvals <- aperm(pvals, c(2,1))
  
  return(pvals)
}

causalOrderingMooij <- function(Xtr, Xte, learner, numFolds, alpha){
  print("enters causalOrderingMooij")
  d <- ncol(Xtr)
  S <- colnames(Xtr) # 1:d
  order <- as.numeric()
  
  totRegs <- d*(d+1)/2 -1
  count <- 0
  # we'll choose one node at a time to be next in back to front causal ordering
  pm0 <- proc.time()
  for (iter in d:2){
    print("**********************")
    print(paste("iter: ",iter))
    
    # for all the nodes not chosen left we use the other nodes not chosen 
    # and see if we get independent residuals. We choose "most" independent
    # regression and take the dependent variable to be next in back to front
    # causal ordering
    pvals <- sapply(S, function(effectNodeTry){
      print(paste("iter: ", iter,"effectNodeTry: ", effectNodeTry))
      count <<- count  + 1
      causeNodesTry <- setdiff(S, effectNodeTry)
      rs <- Residuals(Xtr=Xtr[,causeNodesTry], Ytr=Xtr[,effectNodeTry], 
                      Xte=Xte[,causeNodesTry], Yte=Xte[,effectNodeTry], learner, numFolds)
      pval <- dhsic.test(Xte[, causeNodesTry], rs$rsTe)$p.value
      print(paste("p-value: ", pval))
      
      timePast <- proc.time()-pm0
      avgTime <- timePast/count
      regsLeft <- totRegs - count
      timeLeft <- regsLeft*avgTime
      
      print(paste("Estimated time left causal ordering: ", round(timeLeft[3]/60,2), " mins."))
      
      return(pval)
    })
    
    mostEffectNode <- S[which.max(pvals)]
    if(pvals[mostEffectNode] < alpha){
      print("no consitent DAGs")
      return(NULL)
    }
    order <- c(mostEffectNode, order)
    S <- setdiff(S, mostEffectNode)
    print(paste("causal order so far: ", paste(order, collapse="-> ")))
  }
  order <- c(S, order)
  
  print(paste("final causal order: ", paste(order, collapse="-> ")))
  
  print("exits causalOrderingMooij")
  return(order)
}

mooijDAG <- function(Xtr, Xte, learner, numFolds, alpha){
  print("enters mooijDAG")
  d <- ncol(Xtr)
  print("calculating causal ordering")
  
  order <- causalOrderingMooij(Xtr, Xte, learner, numFolds, alpha)
  
  if(is.null(order)){
    print("no consitent DAGs")
    return(NULL)
  }  
  
  totRegs <- d*(d-1)/2
  count <- 0
  
  print("calculating parents")
  # for each node starts with all the parents and removes certain parents if we keep
  # independence of residuals with reduced inputs
  pm0 <- proc.time()
  parents <- lapply(1:d, function(j){
    node <- order[j]
    print(paste("node: ", node))
    if(j  == 1){
      prnts_node <- as.numeric()
    } else{
      
      prnts_node <- order[1:(j-1)]
      # we will try and take off as many inputs as we can maintaining independence
      for(k in 1:(j-1)){
        print(paste("for node ", node, "testing taking off input: ", order[k]))
        count <<- count + 1
        notParent_node_try <- order[k]
        prnts_node_try <- setdiff(prnts_node, notParent_node_try)
        if(length(prnts_node_try)==0){
          pval <- dhsic.test(Xte[, prnts_node], Xte[,node])$p.value
        } else{
          rs <- Residuals(Xtr=Xtr[,prnts_node_try], Ytr=Xtr[,node],
                          Xte=Xte[,prnts_node_try], Yte=Xte[,node], learner, numFolds)
          
          pval <- dhsic.test(Xte[,prnts_node_try], rs$rsTe)$p.value
          
          
          timePast <- proc.time()-pm0
          avgTime <- timePast/count
          regsLeft <- totRegs - count
          timeLeft <- regsLeft*avgTime
          
          print(paste("Estimated time left parents: ", round(timeLeft[3]/60,2), " mins."))
        }
        
        # if we maintain independence of residuals and inputs (p-value small) while not
        # using input k we take it off the parents list
        print(paste("pval: ", pval))
        if(pval >= alpha){
          print(paste("for node ", node,"taking off input: ", notParent_node_try))
          prnts_node <- setdiff(prnts_node, notParent_node_try)
        }
      }
      
    }
    return(prnts_node)
  })
  
  # make into a matrix and plot
  V <- colnames(Xtr)
  parents <- sapply(parents, function(el) (V %in% el)*1)
  
  rownames(parents) <- V
  colnames(parents) <- V
  
  parents <- getGraph(parents)
  
  print("exits mooijDAG")
  return(parents)
}

########################################################################################
# Functions for CauseEffectPairs testing
########################################################################################



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

