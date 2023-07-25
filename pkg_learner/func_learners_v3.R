# learner package

print("func_learners_v3 pkg")

library(CVST) # cross validation framework

library(gptk) #gaussian process toolkit
library(lbfgs) # non-convex optimization for hsic regression (in learn.hsic)
library(earth) # mars model for initial residual bandwidth estimation (in getFixedParams)
#library(dHSIC) # dhsic (in getFixedParams)
source("./pkg_dHSIC/dHSIC.R")
library(kernlab) # kqr
source("./pkg_kernels/func_kernel_pkg.R")
library(parallel) #mcmapply (in CV.parallel)
library(opera) # loss(x,y, loss.type=list(name="pinball", tau=0.95))
library(MLmetrics) # LogLoss in negCE
library(glmnet) # learn.logReg
library(neuralnet) #learn.nnc

#############################################################################################################*
# learners - general learner functions 
#############################################################################################################*

constructData <- function (x, y){
  stopifnot(is.list(x) || is.vector(x) || is.matrix(x))
  stopifnot(is.list(y) || is.vector(y) || is.matrix(y))
  data = list(x = x, y = y)
  class(data) = "CVST.data"
  return(data)
}

constructLearner <- function (hyperParams, learnParams, optimizeParams, hueristicSet, optimizeSet, learn, predict, ...){
  addParams <- list(...)
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
                      "hsicYhReg", "nhsicRegA", "hsicRegA", "hsicYhRegA", "corre", "cmem_L2",
                      "negLogLik", "hsicLoss2", "pinball", "negCE","KCDC","KCRDC","KCMC")
  stopifnot( all(sapply(optimizeParams$losses, function(el) is.function(el))) && all(names(optimizeParams$losses) %in% validLossFuncs))
  # heuristicSet, optimizeSet, learn and predict should all be functions
  stopifnot(is.function(learn) && is.function(predict)) #&& is.function(heuristicSet) && is.function(optimizeSet)
  
  # getyperParams simply appends all hyperparameters- data_optimizable, data_non_optimizable and non_data in a list
  getHyperParams <- function(learner){
    return(c(learner$hyperParams$data$optimizable, learner$hyperParams$data$non_optimizable, learner$hyperParams$non_data))
  }
  
  learner <- list(hyperParams = hyperParams, learnParams=learnParams, optimizeParams=optimizeParams, 
                  getHyperParams= getHyperParams, 
                  heuristicSet=heuristicSet, optimizeSet=optimizeSet, learn=learn, predict=predict)
  learner <- c(learner, addParams)
  class(learner) = "emley.learner"
  
  return(learner)
}

# set heuristic based and optimization (CV, max-likelihood etc) parameters
setParams <- function(learner, trainData, plot=FALSE){
  #learner should be an object of class "emley.learner"
  stopifnot(class(learner)=="emley.learner")
  stopifnot(class(trainData)=="CVST.data")
  
  #dataNonOptimParams <- learner$heuristicSet(learner, trainData)
  dataNonOptimParams <- lapply(learner$heuristicSet, function(el){
    # el <- learner$heuristicSet[1]
    #print(el)
    do.call(el, list(learner=learner, data=trainData))
    })
  
  dataNonOptimParams <- list(val=unlist(lapply(dataNonOptimParams, function(el) el$val), recursive=FALSE),
                             seq=unlist(lapply(dataNonOptimParams, function(el) el$seq), recursive=FALSE))
  
  #dataNonOptimParams <- unlist(dataNonOptimParams, recursive=FALSE)
  
  paramsNonOpt <- names(learner$hyperParams$data$non_optimizable)
  
  if(length(paramsNonOpt)>0){
    indxNullVal <- which(sapply(lapply(learner$hyperParams$data$non_optimizable, function(el) el$val), is.null))
    paramsNonOpt <- paramsNonOpt[indxNullVal]
  }
  
  for(nm in paramsNonOpt){
    # nm <- paramsNonOpt[1]
    learner$hyperParams$data$non_optimizable[[nm]]$val <- dataNonOptimParams$val[[nm]]
  }
  
  
  learner$hyperParams$trainData <- trainData
  
  
  paramsOpt <- names(learner$hyperParams$data$optimizable) 
  
  # for those optimizable parameters, if their seq value is NULL we first have to get the seq from
  # dataNonOptimParams
  if(length(paramsOpt)>0){
    # from parameters to optimize (paramsOpt), get index of those which have a null seq value

    indxNullSeq <- which(sapply(lapply(learner$hyperParams$data$optimizable, function(el) el$seq), is.null))
    for(nm in paramsOpt[indxNullSeq]){
      # nm <- paramsOpt[indxNullSeq][1]
      learner$hyperParams$data$optimizable[[nm]]$seq <- dataNonOptimParams$seq[[nm]]
    }
  }
  
  #dataOptimParams <- learner$optimizeSet(learner, plot=plot)
  dataOptimParams <- do.call(learner$optimizeSet, list(learner=learner, plot=plot))
  
  #dataOptimParams$opt
  #dataOptimParams$grid["test",,,"negLogLik"]
  
  for(nm in paramsOpt){
    # nm <- paramsOpt[1]
    learner$hyperParams$data$optimizable[[nm]]$val <- dataOptimParams$opt[[nm]]
    
  }
  
  learner$hyperParams$data$optimizable$grid <- dataOptimParams$grid
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
        # i <- 1
        xx <- predList[[i]]$x
        o <- order(xx)
        xx <- xx[o]
        yy <- predList[[i]]$gyh
        yy <- yy[o,,drop=FALSE]
        
        
        for(j in 1:ncol(yy)){
          # j <- 1
          lines(xx,yy[,j], col=j, lwd=2)
        }
      }
      legend("bottomleft", legend=names(predList), lwd=2, col="red")
      
    } else{
      print("data for different predictors not the same")
    }
    
    
    
  } else{
    print("this function is only for 1d regression")
  }
}


# this method should really be a method of kernel learners but since we are only using a 
# general pseudo-learner class we'll leave it outside for now
getKernelPars <- function(learner, kernelName){
  # get kernel name
  
  hyperParams <- learner$getHyperParams(learner)
  hyperParamsNames <- names(hyperParams) 
  kernelParamsNames <- learner$hyperParams$non_data[[kernelName]]$pars
  kernelParamsFuncArgs <- names(kernelParamsNames)
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
  names(kernelParams) <- kernelParamsFuncArgs
  
  return(kernelParams)
}  

getHyperPar <- function(learner, parName){
  hyperParams <- learner$getHyperParams(learner)
  hyperParamsNames <- names(hyperParams)
  indx <- which(hyperParamsNames==parName)
  hyperParam <- hyperParams[[indx]]
  if(is.list(hyperParam)){
    hyperParam <- unlist(hyperParam$val)
  } 
  
  names(hyperParam) <- parName
  
  return(hyperParam)
}

######################################################################################################*
# Heuristics or initialization of models hyper parameters

# data <- trainData
getFixedParams_rbf <- function(learner, data, indxSens=NULL, indxPred=NULL, plot=FALSE, print=FALSE){
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
  xx <- x[,indxPred, drop=FALSE]
  xs <- x[,indxSens, drop=FALSE]
  
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
    # i <- 1; sd1 <- sigmas1[i]
    Kxs <- kernelMatrix("kern_rbf", x=xx, pars=list(sigma=sd1))
    N <- nrow(Kxs)
    H <- diag(N)-matrix(1/N,N,N)
    Kxsc <- H%*%Kxs%*%H
    distsX <- (Kxsc)[lower.tri(Kxsc)]
    res <- var(distsX)
    return(res)
  })
  indxMaxVar <- which.max(varsHsics.sig)
  sigmaVar <- sigmas1[indxMaxVar]
  
  # obtain dynamic range of sigma
  spl <- spline(log(sigmas1,10), varsHsics.sig)
  splf <- splinefun(log(sigmas1,10), varsHsics.sig)
  dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
  tol <- 1e-3
  DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
  DR.sig <- range(log(DR.sig,10))
  
  # obtain sigma saturation point
  sigmas2 <- 10^seq(DR.sig[1],DR.sig[2],1)
  dhsics.sig <- sapply(sigmas2, function(sd1){
    
    Kxs <- kernelMatrix("kern_rbf", x=xx, pars=list(sigma=sd1))
    Ky <- kernelMatrix("kern_lin", x=matrix(data$y,length(data$y),1), pars=list(offset=0))
    
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
    
    Kxb <- kernelMatrix("kern_rbf", x=xs, pars=list(sigma=sd1))
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
    
    Kyhk <- kernelMatrix("kern_rbf", x=yh0, pars=list(sigma=sd1))
    
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
    
    Krg <- kernelMatrix("kern_rbf", x=res0, pars=list(sigma=sd1))
    
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
      
      Kxb <- kernelMatrix("kern_rbf", x=xs, pars=list(sigma=sd1))
      Kyh <- kernelMatrix("kern_lin", x=yh0, pars=list(offset=0))
      
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
      
      Kx <- kernelMatrix("kern_lin", x=xs, pars=list(offset=0))
      Kyhk <- kernelMatrix("kern_rbf", x=yh0, pars=list(sigma=sd1))
      
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
      
      Kx <- kernelMatrix("kern_lin", x=xx, pars=list(offset=0))
      Krg <- kernelMatrix("kern_rbf", x=res0, pars=list(sigma=sd1))
      
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
      
      
      Kxb <- kernelMatrix("kern_rbf", x=xs, pars=list(sigma=sd1))
      Kyhg <- kernelMatrix("kern_rbf", x=yh0, pars=list(sigma=sd2))
      
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
  return(list(val=list(sigma=sigmas[sigmaType], beta=betas[betaType], kappa=kappas[kappaType], 
              gamma=gammas[gammaType], indxPred=list(indxPred), indxSens=list(indxSens)),
              seq=list()))
}

getFixedParams_rbf2 <- function(learner, data){
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:1000) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  sigmaX <- 1/median(as.numeric(dist(unique(x))^2))
  sigmaY <- 1/median(as.numeric(dist(unique(y))^2))
  
  numX <- learner$hyperParams$data$optimizable$sigma.rbf.X$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$sigma.rbf.Y$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangeRbf(x, length.out=numX)
  rngY <- getRangeRbf(y, length.out=numY)
  
  #rng <- range(rngX, rngY)
  #rngX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numX)
  #rngY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numY)
  
  return(list(val=list(sigma.rbf.X=sigmaX, sigma.rbf.Y=sigmaY), seq=list(sigma.rbf.X=rngX, sigma.rbf.Y=rngY)))
}
getFixedParams_laplace1 <- function(learner, data){
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:1000) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  scaleX <- median(as.numeric(dist(unique(x))^2))
  scaleY <- median(as.numeric(dist(unique(y))^2))
  
  numX <- learner$hyperParams$data$optimizable$scaleX$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$scaleY$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangeLaplace1(x, length.out=numX)
  rngY <- getRangeLaplace1(y, length.out=numY)
  
  #rng <- range(rngX, rngY)
  #rngX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numX)
  #rngY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numY)
  
  return(list(val=list(scale.laplace1.X=scaleX, scale.laplace1.Y=scaleY), seq=list(scale.laplace1.X=rngX, scale.laplace1.Y=rngY)))
}
getFixedParams_laplace2 <- function(learner, data){
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:1000) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  scaleX <- median(as.numeric(dist(unique(x))^2))
  scaleY <- median(as.numeric(dist(unique(y))^2))
  
  numX <- learner$hyperParams$data$optimizable$scaleX$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$scaleY$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangeLaplace2(x, length.out=numX)
  rngY <- getRangeLaplace2(y, length.out=numY)
  
  #rng <- range(rngX, rngY)
  #rngX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numX)
  #rngY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numY)
  
  return(list(val=list(scale.laplace2.X=scaleX, scale.laplace2.Y=scaleY), seq=list(scale.laplace2.X=rngX, scale.laplace2.Y=rngY)))
}

getFixedParams_quad <- function(learner, data){
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:500) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  offsetX <- 1
  offsetY <- 1
  
  numX <- learner$hyperParams$data$optimizable$offsetX$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$offsetY$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangeQuad(x, length.out=numX)
  rngY <- getRangeQuad(y, length.out=numY)
  
  return(list(val=list(offset.quad.X=offsetX, offset.quad.Y=offsetY), seq=list(offset.quad.X=rngX, offset.quad.Y=rngY)))
}
getFixedParams_log <- function(learner, data){
  
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:500) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  degreeX <- 1
  degreeY <- 1
  
  numX <- learner$hyperParams$data$optimizable$degreeX$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$degreeY$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangeLog(x, length.out=numX)
  rngY <- getRangeLog(y, length.out=numY)
  
  return(list(val=list(degree.log.X=degreeX, degree.log.Y=degreeY), seq=list(degree.log.X=rngX, degree.log.Y=rngY)))
}
getFixedParams_poly <- function(learner, data){
  # degree, scale, offset
  if(!is.null(dim(data$x))){
    n <- dim(data$x)[1]
    px <- dim(data$x)[2]
  } else{
    n <- length(data$x)
    px <- 1
  }
  
  if(!is.null(dim(data$y))){
    n <- dim(data$y)[1]
    py <- dim(data$y)[2]
  } else{
    n <- length(data$y)
    py <- 1
  }
  
  if(n > 1000){
    data <- getSubset(data, 1:500) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  degreeX <- 2
  degreeY <- 2
  scaleX <- 1
  scaleY <- 1
  offsetX <- 0
  offsetY <- 0
  
  numX <- learner$hyperParams$data$optimizable$degreeX$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$degreeY$length.out
  if(is.null(numY)) numY <- 10
  
  rngX <- getRangePoly(x, length.out=numX)
  rngY <- getRangePoly(y, length.out=numY)
  
  return(list(val=list(offset.poly.X=offsetX,     offset.poly.Y=offsetY,     scale.poly.X=scaleX,     scale.poly.Y=scaleY,     degree.poly.X=degreeX, degree.poly.Y=degreeY), 
              seq=list(offset.poly.X=rngX$offset, offset.poly.Y=rngY$offset, scale.poly.X=rngX$scale, scale.poly.Y=rngY$scale, degree.poly.X=rngX$degree, degree.poly.Y=rngY$degree)))
}


getRangeRbf_dep <- function(x, length.out=10){
  sigma1 <- 1/quantile(as.numeric(dist(x)^2), 0.4)  
  sigmas <- 10^seq(-16,16,by=1)
  minPct <- 0.001
  Ks <- sapply(sigmas, function(sd) kern_rbf(unique(x), sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #sigmas.x[c(indxStart, indxMax.x, indxFinish)]
  #log(sigmas.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  #vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  sigmas <- 10^seq(log(sigmas.x[indxStart],10), log(sigmas.x[indxFinish],10), length.out=10)
  #sigmas <- 10^seq(log(sigmas.x[indxStart],10), log(sigma1,10), length.out=10)
  minPct <- 0.001
  Ks <- sapply(sigmas, function(sd) kern_rbf(unique(x), sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  
  ini <- log(sigmas.x, 10)[indxStart]
  fin <- log(sigmas.x, 10)[indxFinish] #log(sigma1, 10)
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}
getRangeLog_dep <- function(x, length.out=10){
  degrees <- c(seq(0.01,0.99, length.out=10), round(seq(1,50, length.out=10)))
  Ks <- sapply(degrees, function(dg) kern_log(x, degree=dg), simplify="array")
  num0s <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    sum(mat2==0, na.rm=T)/2
  })
  deg1st0 <- degrees[which(num0s>0)[1]]
  #plot(degrees, num0s, ylab="num0s") 
  #abline(h=0, v=deg1st0, col="red")
  seq <- seq(1, deg1st0, length.out=length.out)
  return(seq)
}
getRangeQuad_dep <- function(x, length.out=10){
  offsets <- c(10^seq(-9,0, length.out=10), seq(1,100, length.out=10))
  Ks <- sapply(offsets, function(of) kern_quad(x, offset=of), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  minPct <- 0.1
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(offsets[indx], vars[indx])
  splf <- splinefun(offsets[indx], vars[indx])
  offsets.x <- unique(sort(c(offsets[indxMax],10^seq(min(log(offsets[indx], 10)), max(log(offsets[indx], 10)), length.out=100))))
  vars.x <- splf(offsets.x, deriv=0)
  indxMax.x <- which(offsets.x == offsets[indxMax]); offsets.x[indxMax.x]
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  #plot(offsets, vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(offsets.x, vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=offsets.x[c(indxMax.x, indxFinish)], col="orange")
  
  ini <- vars.x[indxMax.x]
  fin <- vars.x[indxFinish]
  seq <- seq(ini, fin, length.out=length.out)
  return(seq)
}

getRangeRbf <- function(x, length.out=10){
  
  sigma0 <- 1/quantile(as.numeric(dist(x)^2), (length(x)-1)/length(x))
  sigma1 <- 1/quantile(as.numeric(dist(x)^2), 0.6)  
  
  sigmas <- 10^seq(-16,16,by=1)
  minPct <- 0.01
  Ks <- sapply(sigmas, function(sd) kern_rbf(x, sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    #hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- min(100,indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct)))
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #abline(v=log(c(sigma0, sigma1),10), col="green")
  #sigmas.x[c(indxStart, indxMax.x, indxFinish)]
  #log(sigmas.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  #vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  sigmas <- 10^seq(log(sigmas.x[indxStart],10), log(sigmas.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(sigmas, function(sd) kern_rbf(x, sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- min(100, indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct)))
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #abline(v=log(c(sigma0, sigma1),10), col="green")
  
  ini <- log(sigmas.x, 10)[indxStart]
  fin <- log(sigmas.x, 10)[indxFinish]
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}

getRangeRbf2 <- function(x, length.out=10){
  dists2 <- as.numeric(dist(x)^2)
  dists2 <- dists2[which(dists2>0)]
  
  sigma1 <- 1/quantile(dists2, 0.9)  
  sigma2 <- 1/quantile(dists2, 0.1)
  
  #res <- 10^seq(log(sigma1,10), log(sigma2,10), length.out=length.out)
  res <- c(sigma1, sigma2)
  
  return(res)
}
getRangeLaplace1 <- function(x, length.out=10){
  
  # scale0 <- quantile(as.numeric(dist(x)^2), (length(x)-1)/length(x))
  # scale1 <- quantile(as.numeric(dist(x)^2), 0.6)  
  
  scales <- 10^seq(-10,10,by=1)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace1(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    #hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  # plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  # lines(spl, col="green")
  # lines(log(scales.x, 10), vars.x, col="blue", type="l")
  # abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  # abline(v=log(c(scale0, scale1),10), col="green")
  # scales.x[c(indxStart, indxMax.x, indxFinish)]
  # log(scales.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  # vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  scales <- 10^seq(log(scales.x[indxStart],10), log(scales.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace1(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  # plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  # lines(spl, col="green")
  # lines(log(scales.x, 10), vars.x, col="blue", type="l")
  # abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  # abline(v=log(c(scale0, scale1),10), col="green")
  
  ini <- log(scales.x, 10)[indxStart]
  fin <- log(scales.x, 10)[indxFinish]
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}
getRangeLaplace2 <- function(x, length.out=10){
  
  # scale0 <- quantile(as.numeric(dist(x)^2), (length(x)-1)/length(x))
  # scale1 <- quantile(as.numeric(dist(x)^2), 0.6)  
  
  scales <- 10^seq(-10,10,by=1)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace2(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    #hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  # plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  # lines(spl, col="green")
  # lines(log(scales.x, 10), vars.x, col="blue", type="l")
  # abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  # abline(v=log(c(scale0, scale1),10), col="green")
  # scales.x[c(indxStart, indxMax.x, indxFinish)]
  # log(scales.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  # vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  scales <- 10^seq(log(scales.x[indxStart],10), log(scales.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace2(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  # plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  # lines(spl, col="green")
  # lines(log(scales.x, 10), vars.x, col="blue", type="l")
  # abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  # abline(v=log(c(scale0, scale1),10), col="green")
  
  ini <- log(scales.x, 10)[indxStart]
  fin <- log(scales.x, 10)[indxFinish]
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}
getRangeLog <- function(x, length.out=10){
  
  
  degrees1 <- 10^seq(-9,-1, length.out=round(length.out/2))
  degrees2 <- round(seq(2,50, length.out=round(length.out/2)))
  
  
  Ks1 <- sapply(degrees1, function(dg) kern_log(unique(x), degree=dg), simplify="array")
  Ks2 <- sapply(degrees2, function(dg) kern_log(unique(x), degree=dg), simplify="array")
  num0s <- apply(Ks2, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    sum(mat2==0, na.rm=T)/2
  })
  vars <- apply(Ks1, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    #hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  indxNonZeros <- which(num0s>0)
  if(length(indxNonZeros)>0){
    deg1st0 <- indxNonZeros[1]
  } else{
    deg1st0 <- length(degrees2)
  }
  
  indxNonZeroVar <- which(vars>1e-10) 
  if(length(indxNonZeroVar)>0){
    degG0 <- indxNonZeroVar[1]
  } else{
    degG0 <- 1
  }
  
  
  #plot(degrees, num0s, ylab="num0s") 
  #abline(h=0, v=degrees[c(deg1st0, degG0)], col="red")
  #plot(log(degrees, 10), log(vars,10), ylab="num0s") 
  #abline(h=0, v=log(degrees[c(deg1st0, degG0)],10), col="red")
  
  degrees1 <- 10^seq(log(degrees1, 10)[degG0],-1, length.out=round(length.out/2)-1)
  degrees2 <- round(seq(2,degrees2[deg1st0], length.out=round(length.out/2)-1))
  
  
  degrees <- unique(sort(c(degrees1, degrees2, 0.5, 0.98)))
  
  return(degrees)
}
getRangeQuad <- function(x, length.out=10){
  
  
  offsets <- 10^seq(-16,16,by=1)
  minPct <- 0.01
  Ks <- sapply(offsets, function(off) kern_quad(x, offset=off), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    #hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(offsets[indx],10), vars[indx])
  splf <- splinefun(log(offsets[indx],10), vars[indx])
  offsets.x <- unique(sort(c(offsets[indxMax],10^seq(min(log(offsets[indx], 10)), max(log(offsets[indx], 10)), length.out=100))))
  vars.x <- splf(log(offsets.x, 10), deriv=0)
  indxMax.x <- which(offsets.x == offsets[indxMax]); offsets.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- min(100, indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct)))
  
  #plot(log(offsets,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(offsets.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(offsets.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #offsets.x[c(indxStart, indxMax.x, indxFinish)]
  #log(offsets.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  #vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  offsets <- 10^seq(log(offsets.x[indxStart],10), log(offsets.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(offsets, function(off) kern_quad(x, offset=off), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    #hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  #spl <- spline(log(offsets[indx],10), vars[indx])
  splf <- splinefun(log(offsets[indx],10), vars[indx])
  offsets.x <- unique(sort(c(offsets[indxMax],10^seq(min(log(offsets[indx], 10)), max(log(offsets[indx], 10)), length.out=100))))
  vars.x <- splf(log(offsets.x, 10), deriv=0)
  indxMax.x <- which(offsets.x == offsets[indxMax]); offsets.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- min(100, indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct)))
  
  #plot(log(offsets,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(offsets.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(offsets.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  
  ini <- log(offsets.x, 10)[indxStart]
  fin <- log(offsets.x, 10)[indxFinish]
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}

getRangePoly <- function(x, length.out=10){
  # offset>=0, scale>0, degree in positive integers
  
  offsets <- 10^seq(0,3,by=1)
  scales <-  10^seq(-3,3,by=1)
  degrees <- seq(1,3, by=1)
  
  res <- list(offset=offsets, scale=scales, degree=degrees)
  
  return(res)
}




myGpCreate <- function(learner, data){
  options <- gpOptions(approx=learner$optimizeParams$approx)
  options$kern$comp <- learner$hyperParams$non_data$kernel  
  options$numActive <- min(learner$optimizeParams$numActive, nrow(data$x))
  
  if(learner$optimizeParams$fixInducing){
    options$fixInducing <- TRUE
    wishInducing <- makeGrid(data$x, options$numActive)
    
    # obrain indices of nearest points to the inducing points we would like
    
    indxs <- knnx.index(data=data$x, query=wishInducing, k=1)

    options$fixIndices <- as.numeric(indxs)
  }
  
  model <- gpCreate(q=dim(data$x)[2], d=1, X=data$x, y=as.matrix(data$y), options)
  return(list(val=list(modelInit=model), seq=NULL))
}


######################################################################################################*
# optimization functions for learners' hyperparameters

optHP.CV <- function(learner, plot=FALSE, fac=1){
  #print("enters optHP.CV")
  
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    return(list())
  }
  
  numParams <- length(optimHyperParams)
  
  numEach <- lapply(optimHyperParams, function(el) length(el$seq))
  
  paramsList <- lapply(optimHyperParams, function(el) el$seq)
  names(paramsList) <- names(optimHyperParams)
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
  
  constructParams <- function (otherParams, paramsList) {
    pn = names(paramsList)#names(substitute(paramsList))[-1]
    ret = expand.grid(paramsList, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    params = lapply(1:nrow(ret), function(ind){ 
       aux <- as.list(ret[ind,])
       names(aux) <- colnames(ret)
       return(c(aux, otherParams))
      })
    paramNames = lapply(1:nrow(ret), function(ind) paste(pn, ret[ind, ], sep = "=", collapse = " "))
    names(params) = paramNames
    class(params) = "CVST.params"
    return(params)
  }
  
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  # not necessary coz CV is only done on data hyperparms??
  #params <- lapply(params, function(el) c(el, learner$hyperParams$non_data))
  #class(params) = "CVST.params"
  
  
  grid <- CV.parallel(learner, params, fac=fac)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
    
  # reshape back into one dimension per hyperparameter
  
  
  
  dimnms <- dimnames(grid)
  
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], numEach , dim(grid)[3])
  
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  
  dimnames(grid) <- dimnms
  
  # obtain best hyperparameters
  
  testGrid <- as.array(apply(grid, names(paramsList), function(el) el[learner$optimizeParams$testTrain,1])) #grid["test",,,1] #el["test",1]
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

CV.parallel.avg <- function(learner, params, fac=1, verbose=TRUE) {
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
      # f <- 1; p <- params[[1]]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- getSubset(trainData, -validationIndex)
      curTest <- getSubset(trainData, validationIndex)
      # either mean squared error or mean classification error
      
      
      learnerAux <- learner
      learnerAux$hyperParams$trainData <- curTrain
      
      nmsPars <- names(learnerAux$hyperParams$data$optimizable)
      for(pr in nmsPars){
        # pr <- nmsPars[1]
        learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
      }
      
      #learner <- learnerAux
      # learnerAux$learnParams$alpha
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


pred.CV <- function(learner, data){
  
  numFolds <- learner$optimizeParams$numFolds
  n <- getN(data)
  size <- floor(n / numFolds)
  
  res <- mcmapply(FUN=function(f){
    # f <- 1
    # print(f)
    last <- f*size
    if(f==numFolds) last <- n
    validationIndex <- seq((f-1)*size + 1, last)
    curTrain <- getSubset(data, -validationIndex)
    curTest <- getSubset(data, validationIndex)
    # either mean squared error or mean classification error
    
    
    learner1 <- learner
    learner1$hyperParams$trainData <- curTrain
    
    learner1 <- try(learner1$learn(learner=learner1, plot=FALSE))
    lambdaAux <- 1
    while(class(learner1)=="try-error"){
      print(paste("singular matrix for parameters = ", paste(paste(names(p), p, collapse="-") ,collapse=", "), "for fold = ", f))
      
      learner1 <- learner
      learner1$hyperParams$data$optimizable$lambda$val <- lambdaAux
      lambdaAux <- lambdaAux*10
      dummyTrain <- curTrain
      dummyTrain$x <- matrix(rnorm(nrow(curTrain$x)*ncol(curTrain$x), nrow(curTrain$x), ncol(curTrain$x)))
      dummyTrain$y <- rnorm(length(curTrain$y))
        
      learner1$hyperParams$trainData <- dummyTrain
      
      learner1 <- try(learner1$learn(learner1))
      
      if(class(learner1)!="try-error"){
        #print("exits!")
        predTrain <- learner1$predict(learner=learner1, data=curTrain)
        predTest <- learner1$predict(learner=learner1, data=curTest)
      
        nms <- names(predTrain)
        predTrain <- lapply(nms, function(el){
        # el <- nms[1]
        res <- predTrain[[el]]
        res[predTrain[[el]]>0 | predTrain[[el]]<=0] <- NA
        return(res)
      })
        names(predTrain) <- nms
      
        predTests <- lapply(nms, function(el){
        # el <- nms[1]
        res <- predTest[[el]]
        res[predTest[[el]]>0 | predTest[[el]]<=0] <- NA
        return(res)
      })
        names(predTest) <- nms
      
      
        res <- list(train=predTrain, test=predTest)
        return(res)
      }
    }
    
    predTrain <- learner1$predict(learner=learner1, data=curTrain)
    predTest <- learner1$predict(learner=learner1, data=curTest)
    
    
    res <- list(train=predTrain, test=predTest)
    
    
    return(res)
  }, f=1:numFolds, mc.cores=1, SIMPLIFY=FALSE)
  
  auxTrain <- lapply(res, function(el) el$train)
  auxTest <- lapply(res, function(el) el$test)
  nmsAux <- names(auxTrain[[1]])
  
  predTrain <- lapply(nmsAux, function(nm) do.call(rbind, lapply(auxTrain, function(el) el[[nm]])))
  predTest <- lapply(nmsAux, function(nm) do.call(rbind, lapply(auxTest, function(el) el[[nm]])))
  names(predTrain) <- names(predTest) <- nmsAux
  
  # plot(predTest$x, predTest$gy); ord <- order(predTest$x)
  # for(i in 1:ncol(predTest$gyh)) lines(predTest$x[ord], predTest$gyh[ord,i], col=i)
  
  # plot(predTrain$x, predTrain$gy); ord <- order(predTrain$x)
  # for(i in 1:ncol(predTrain$gyh)) lines(predTrain$x[ord], predTrain$gyh[ord,i], col=i)
  
  return(list(train=predTrain, test=predTest))
}


CV.parallel <- function(learner, params, fac=1, verbose=TRUE) {
  stopifnot(class(learner) == "emley.learner" && class(params) == "CVST.params")
  
  #print("enters CV parallel")
  
  numFolds <- learner$optimizeParams$numFolds
  trainData <- learner$hyperParams$trainData
  
  
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), names(params))
  
  
  
  pm <- proc.time()
  losses <- mcmapply(FUN=function(p){
    # p <- params[[1]]
    #print(p)
    learnerAux <- learner
    nmsPars <- names(learnerAux$hyperParams$data$optimizable)
    for(pr in nmsPars){
      # pr <- nmsPars[1]
      learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
    }
    
      learnerAux <- learnerAux$learn(learnerAux)
      preds <- pred.CV(learner=learnerAux, data=trainData)
    
    # x <- trainData$x; y <- trainData$y; ord <- order(x)
    # plot(x, y)
    # for(i in 1:ncol(preds$test$gyh)) lines(x[ord], preds$test$gyh[ord,i], col=i)
    
    lossesTrain <- sapply(names(learnerAux$optimizeParams$losses), function(func){
      do.call(func, list(learner=learnerAux, pred=preds$train))
    })
    
    lossesTest <- sapply(names(learnerAux$optimizeParams$losses), function(func){
      # func <- "cmem_L2_f"
      do.call(func, list(learner=learnerAux, pred=preds$test))
    })
    
    res <- cbind(lossesTrain, lossesTest)
    #print(res)
    
    return(res)
  }, params, mc.cores=1, SIMPLIFY="array")
  proc.time() - pm #
  
  losses <- aperm(losses, c(2,3,1))
  dimnames(losses) <- list(trainTest=c("train","test"), params=names(params), var=names(learner$optimizeParams$losses))
  
  #print("exits CV parallel")
  return(losses)
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


#############################################################################################################*
# Kernel Regression learners 
#############################################################################################################*

#############################################################################################################*
# learn and predict functions for different kernel regression learners


learn.vanilla <- function(learner, plot=FALSE) {
  return(learner)
}
predict.vanilla <- function(learner, data) {
  
  gy <- as.matrix(data$y)
  x <- data$x
  return(list(x=as.matrix(x), gy=gy, gyh=matrix(0, nrow(gy), ncol(gy))))
}


learn.krr <- function(learner, plot=FALSE) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")

  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  Kxs <- Kxs +matrix(1,N,N)
  N <- nrow(Kxs)
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  
  alpha <- solve(Matrix(Kxs + diag(lambda, N))) %*% y
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.krr <- function(learner, data) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  kxs <- kxs + matrix(1, nrow(kxs), ncol(kxs))
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}

learn.qhsic <- function (learner, plot=FALSE) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  
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
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, x, pars=parsXb)
  
  
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
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  trainData <- learner$hyperParams$trainData
  
  if(nrow(trainData$x) > learner$optimizeParams$maxPoints){
    trainData <- getSubset(trainData, learner$learnParams$indxInducing) 
  } else{
    indxs <- NA
  }
  
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}	


learn.kqr <- function(learner, plot=FALSE) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  N <- nrow(Kxs)
  lambda <- getHyperPar(learner, "lambda")*N
  taus <- getHyperPar(learner, "taus")  
  
  qrms <- lapply(c(3/N, taus, 1-3/N), try(function(t) kqr(x, y, tau = t , kernel = "rbfdot", kpar= parsXs, C=lambda, scaled=FALSE)))
  
  alphas <- sapply(qrms, function(qrm){
    if(class(qrm)=="try-error") alpha <- rep(NA, N) else alpha <- qrm@alpha
    return(alpha)
  }, simplify="matrix")
  colnames(alphas) <- c(3/N, taus, 1-3/N)
  
  learner$learnParams$alphas <- alphas
  
  bs <- sapply(qrms, function(qrm){
    if(class(qrm)=="try-error") b <- NA else b <- qrm@b
    return(b)
  })
  names(bs) <- c(3/N, taus, 1-3/N)
  
  learner$learnParams$bs <- bs
  
  return(learner)
}

predict.kqr <- function(learner, data) {
  
  
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- t(t(kxs %*% learner$learnParams$alpha)- learner$learnParams$b)
  
  
  gy <- data$y
  x <- data$x

  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=pred))
}

# helper function to estimate CDF by interpolation from, possibly non-monotonous quantiles
estimateCDF <- function(qs, ps){
  # remove NAs from these pairs
  indxNA <- which(!is.na(qs))
  qs <- qs[indxNA]
  ps <- ps[indxNA]
  
  # the cdf is non decreasing so we  delete such points - if we dont take appropriate measures
  # we could end up with only one point!! 
  # possible measures:
  # 1. use very smooth quantile corresponding to p=0 and p=1 such that they will almost always be q0<q1 and you wd at least have 2 pts
  # 2. correct q0 as q0 = min(qs) and q1 as q1 = max(qs)  before going on - this corresponds to a uniform in range of qs
  # 3. for training first and last quantiles use whole sample (we wd still need to do 1 or 2 but it might help)
  # in the spirit of letting the noise and causal direction to punish/reward the parameters i prefer option 2
  
  qs[1] <- min(qs)
  qs[length(qs)] <- max(qs)
  
  difqs <- diff(qs)
  negDif <- which(difqs<0)
  while(any(negDif)){
    qs[negDif+1] <- NA
    indxNA <- which(!is.na(qs))
    qs <- qs[indxNA]
    ps <- ps[indxNA]
    difqs <- diff(qs)
    negDif <- which(difqs<0)
  }
  #qs <- sort(qs)
  # itnerpolate to estimate cdf
  splf <- splinefun(qs, ps, method="monoH.FC") 
  # plot(seq(min(qs), max(qs), length.out=100), splf(seq(min(qs), max(qs), length.out=100)), type="l"); lines(qs, ps, type="p")
  return(list(qs=qs, CDF=splf))
}

# obtain predictions (liklihoods) and residuals (cdfs) of quantile regression model
resids.kqr <- function(learner, pred){
  ys <- pred$gy
  
  taus <- getHyperPar(learner, "taus")  
  preds <- pred$gyh
  # xs <- pred$x; plot(xs, ys); ord <- order(xs); for(i in 1:ncol(preds)) lines(xs[ord], preds[ord,i], col=i)
  # abline(v=xs[50], col="red")
  res <- sapply(1:length(ys), function(i){
    # i <- 50
    #print(paste("i: ", i))
    
    # CDF estimation
    # get tau quantiles and cdf pairs (q,p) from our quantile reg model
    qs_i <- preds[i, ]
    ps_i <- c(0,taus,1)
    CDFi <- estimateCDF(qs_i, ps_i)
    qs_i <- CDFi$qs
    CDFi <- CDFi$CDF
    
    
    
    # we are interested in evaluating the cdf and pdf py(.|x_i) evaluating at y_i so we get set pt <- y_i
    y_i <- ys[i]
    # we evaluate py(y_i|x_i) and Py(y_i|x_i) (pdf and cdf)
    lik <- CDFi(y_i, deriv=1)
    cdf <- CDFi(y_i, deriv=0)
    
    # for testing out of training sample y_i could fall out of support of py(y_i|x_i) in this case 
    # if there are some y_j for any j that fall within qi0-qi1 we simply take lik=min_j p(y_j|x_i) j such that qi0<=y_j<=qi1
    # else we look for nearby x_j (using kernel) and set lik= p(y_i|x_j) if qj0<=y_i<=qj1 and lik=min_k p(y_k|x_j) if
    # !(qj0<=y_i<=qj1) but there exists at least 1 yk such that qj0<=y_k<=qj1
    if(lik<=0 | findInterval(y_i, range(qs_i)) %in% c(0,2)){
      aux0 <- CDFi(ys, deriv=1)
      indxValid <- which(aux0>0 & findInterval(ys, range(qs_i))==1)
      if(length(indxValid)>0){
        indx <- which.min(aux0[indxValid])
        lik <- aux0[indxValid][indx]*0.5
      } else{
        print("entering no mans land!!!")
        lik <- 0
      }
    }
    if(cdf<0) cdf <- 0
    if(cdf>1) cdf <- 1
  
    return(c(lik, cdf))
  }, simplify="array")
  res <- t(res)
  colnames(res) <- c("lik","cdf")
  indxNA <- which(is.na(res[,"lik"]))
  res[indxNA,"lik"] <- 0#1e-6
  res[,"lik"] <- -log(res[,"lik"])
  colnames(res) <- c("pred","resid")
  return(res)
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
  Krg <- kernelMatrix(kernelRg$name, rs, pars=kernelRg$pars)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- Kxb%*%H
  Krgc <- Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc)) + lambda*t(alpha)%*%Kxs%*%alpha  
  return(res)
}
.hsicRegLossGrad <- function(alpha, y, Kxs, Kxb, kernelRg, lambda){
  rs <- y-Kxs%*%alpha
  Krg <- kernelMatrix(kernelRg$name, rs, pars=kernelRg$pars)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  ones <- matrix(1,N,1)
  aux <- Kxbc*Krg*(ones%*%t(rs)-rs%*%t(ones))
  gamma <- kernelRg$pars$sigma
  
  grad <- apply(diag(N),2, function(ek)t(ones)%*%(aux*(ones%*%t(ek)%*%Kxs-Kxs%*%ek%*%t(ones)))%*%ones)
  
  grad <- 2*gamma*grad  + 2*lambda*Kxs%*%alpha
  return(grad)
}


learn.hsic <- function (learner, plot=FALSE) {
  #print("enters hsic learn function")
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  parsRg <- getKernelPars(learner, kernelName="kernelRg")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, x, pars=parsXb)
  kernelRg <- list(name=learner$hyperParams$non_data$kernelXb$name, pars=parsRg)
  
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
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=learner$hyperParams$trainData$x, pars=parsXs)
  

  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}



learn.gptk <- function(learner, plot=FALSE) return(learner)
predict.gptk <- function(learner, data){
  aux <- gpPosteriorMeanVar(learner$hyperParams$data$optimizable$model$val, X=data$x)
  pred <- list(x=data$x, gy=data$y, gyh=aux)
  return(pred)
}

learn.logReg <- function(learner, plot=FALSE) {
  
  phix <- learner$makeFeature(learner, data=learner$hyperParams$trainData, "X")
  
  y <- learner$hyperParams$trainData$y
  
  dataMat <- as.data.frame(cbind(y, phix))
  colnames(dataMat) <- c("y", paste("phix", seq(ncol(phix)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  fit <- glm(formula=form, data=dataMat, family = binomial(link="logit"))
  
  
  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logReg <- function(learner, data) {
  
  phix <-learner$makeFeature(learner, data, "X") 
  
  y <- learner$hyperParams$trainData$y
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  
  pred <- predict(fit, newdata=dataMat, type="response")
  
  
  return(list(x=as.matrix(data$x), gy=as.matrix(data$y), gyh=as.matrix(pred)))
}


learn.logKReg <- function(learner, plot=FALSE) {
  
  phix <- learner$makeFeature(learner, data=learner$hyperParams$trainData, "X")
  
  y <- learner$hyperParams$trainData$y
  
  dataMat <- as.data.frame(cbind(y, phix))
  colnames(dataMat) <- c("y", paste("phix", seq(ncol(phix)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  fit <- glmnet(x=phix, y=y, family = "binomial", alpha = 0, lambda = lambda)
  #fit <-  glm(form, data=dataMat, family=binomial(link="logit"))  

  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logKReg <- function(learner, data) {
  
  phix <-learner$makeFeature(learner, data, "X") 

  y <- learner$hyperParams$trainData$y
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  #pred <- predict(fit, newdata=dataMat, type="response")
  pred <- predict(fit, newx=phix, s=lambda, type="response")
  
  return(list(x=as.matrix(data$x), gy=as.matrix(data$y), gyh=as.matrix(pred)))
}

learn.logRegInt <- function(learner, plot=FALSE) {
  
  phix <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  
  lambdas <- getHyperPar(learner, "lambda")*N
  
  numFolds <- learner$optimizeParams$numFolds
  
  cv.fit <- cv.glmnet(x=phix, y=y, family="binomial", alpha=0, lambda=lambdas, nfolds=numFolds)
  
  #plot(log(cv.fit$lambda,10), cv.fit$cvm)
  #plot(cv.fit)
  
  
  indx <- which.min(cv.fit$cvm)
  #fit <- cv.fit$glmnet.fit
  lambda <- cv.fit$lambda[indx]
  
  learner$hyperParams$data$non_optimizable$lambda$val <- lambda
  
  fit <- glmnet(x=phix, y=y, family = "binomial", alpha = 0, lambda = lambda)
  #fit <-  glm(form, data=dataMat, family=binomial(link="logit"))  
  
  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logRegInt <- function(learner, data) {
  
  phix <-data$x 
  
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  #pred <- predict(fit, newdata=dataMat, type="response")
  pred <- predict(fit, newx=phix, s=lambda, type="response")
  
  return(list(x=as.matrix(data$x), gy=as.matrix(data$y), gyh=as.matrix(pred)))
}


learn.logNNReg <- function(learner, plot=FALSE){
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  depth <- getHyperPar(learner, "depth")
  breadth <- getHyperPar(learner, "breadth") 
  
  dataMat <- as.data.frame(cbind(y, x))
  colnames(dataMat) <- c("y", paste("x", seq(ncol(x)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  
  maxTrys <- 10
  count <- 1
  minThresh <- 0.01
  weights <- NULL
  pm <- proc.time()
  while(count < maxTrys & is.null(weights)){
     print(count)
     nn <- neuralnet(formula=form, data=dataMat, hidden=rep(breadth, depth), 
                    linear.output=FALSE, threshold=minThresh, rep=1,
                    stepmax = 1e3, lifesign=c("none","minimal","full")[1], lifesign.step=20)
     
     
    weights <- nn$weights 
    minThresh <- minThresh*10
    count <- count + 1
  }
  proc.time() - pm #
  # plot(nn)
  
  learner$learnParams$model_nn <- nn
  
  phix <- learner$makeFeature(learner, learner$hyperParams$trainData, "X")
  
  dataMat <- as.data.frame(cbind(y, phix))
  colnames(dataMat) <- c("y", paste("phix", seq(ncol(phix)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  fit <- glm(formula=form, data=dataMat, family = binomial(link="logit"))
  
  learner$learnParams$model_log <- fit

  return(learner)  
} 
predict.logNNReg <- function(learner, data){
  phix <-learner$makeFeature(learner, data, "X") 
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model_log
  pred <- predict(fit, newdata=dataMat, type="response")
  
  
  return(list(x=as.matrix(data$x), gy=as.matrix(data$y), gyh=as.matrix(pred)))
}

# generic transformation and residual functions

NCE <- function(x, classifier, keepHyperParams=NULL, plot=FALSE){
  n <- nrow(x)
  p <- ncol(x)
  midPt <- round(n/2)
  indx.tr <- 1:midPt
  indx.te <- (midPt+1):n
  n_tr <- length(indx.tr)
  n_te <- length(indx.te)
  
  # make train classification data
  x_tr <- as.matrix(rbind(x[indx.tr,, drop=F], matrix(runif(n_tr*p), n_tr, p)))
  y_tr <- c(rep(1, n_tr), rep(0, n_tr)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr))
  trainClassData <- constructData(x_tr[smpl,,drop=F], y_tr[smpl])
  
  # make test classification data
  x_te <- as.matrix(rbind(x[indx.te,,drop=F], matrix(runif(n_te*p), n_te, p)))
  y_te <- c(rep(1, n_te), rep(0, n_te)) 
  set.seed(123)
  smpl <- sample(nrow(x_te))
  testClassData <- constructData(x_te[smpl, ,drop=F], y_te[smpl])
  
  # train classifier
  classifier <- setParams(learner=classifier, trainData=trainClassData)
  
  # getHyperPar(classifier, "sigma")
  keepHypers <- lapply(keepHyperParams, function(par) getHyperPar(classifier, par))
  
  # NOTE: need to check if giving each group of dists a different seed is better for 
  # KCDC classification as I observed for sinx example with discreete bins
  data <- constructData(x, c(y_tr, y_te))
  classifier <- classifier$learn(learner=classifier)
  
  phix <- classifier$makeFeature(learner=classifier, data, "X")
  
  meanPhix <- apply(phix, 2, mean)
  
  pred <- classifier$pred(classifier, testClassData)
  loss <- do.call(classifier$optimizeParams$losses[[1]], list(learner=classifier, pred=pred))
  
  res <- list(meanPhix=meanPhix, loss=loss, hyperPars=keepHypers)
  
  if(ncol(x)==1){
    xx <- seq(0, 1, length.out=100)
    pdfData <- constructData(as.matrix(xx), rep(0,100))
    predPdf <- classifier$pred(learner=classifier, data=pdfData)
    pdf <- predPdf$gyh/(1-predPdf$gyh)
    res <- c(res, list(pdf=pdf))
    numBins <- round(n/10)
    if(plot){
      hist(x, numBins, prob=T)
      lines(xx, pdf, col="red")
    }
  }
  
  return(res)
}


makePhi <- function(learner, data, var=c("X","Y")){
  var2 <- c("x","y")[match(var, c("X","Y"))]
  pars <- getKernelPars(learner, kernelName=paste("feature", var, sep=""))
  x <- data[[var2]]
  auxPars <- pars
  auxPars$x <- x
  phi_char <- learner$hyperParams$non_data[[paste("feature", var, sep="")]]$name
  phi_char <- strsplit(phi_char, split="T")[[1]][1]
  aux <- strsplit(phi_char, split="_")[[1]]
  if(aux[1] == "rff") phi_char <- "rff"
  # names(auxPars)
  # rff(x=auxPars$x, num_f=auxPars$num_f, seed=auxPars$seed, p_w=auxPars$p_w, map=auxPars$map, sigma=auxPars$sigma)
  phi  <- do.call(phi_char, auxPars)
  return(phi)
}

makeScoreFeat_bin <- function(learner, data){
  
}

makeNNFeats <- function(learner, data, var=NULL){
  nn <- learner$learnParams$model_nn
  pred_nn <- predict(nn, newdata=cbind(x=data$x) , all.units=T)
  phi_nn <- pred_nn[[length(pred_nn)-1]]
  phi_nn <- phi_nn[,1:(ncol(phi_nn)-1), drop=F]
  return(phi_nn)
}

resids.add <- function(learner, pred){
  preds <- pred$gyh
  resids <- as.matrix(pred$gyh - pred$gy)
  res <- cbind(preds, resids)
  colnames(res) <- c("pred","resid")
  return(res)
} 

#############################################################################################################*
# Loss functions for cross validation and evaluation of learners

sse <- function(learner, pred){
  #print("enters sse")
  resids <- learner$resids(learner, pred)[,"resid"]
  res <- sum((resids)^2)
  #print("exits sse")
  return(res)
}
mse <- function(learner, pred){
  #print("enters mse")
  resids <- learner$resids(learner, pred)[,"resid"]
  res <- mean((resids)^2)
  #print("exits mse")
  return(res)
}
rmse <- function(learner, pred){
  #print("enters rmse")
  resids <- learner$resids(learner, pred)[,"resid"]
  res <- sqrt(mean((resids)^2))
  #print("exits rmse")
  return(res)
}



negLogLik <- function(learner, pred){
  negLogLiks <- learner$resids(learner, pred)[,"pred"]
  return(sum(negLogLiks))
  
}
hsicLoss2 <- function(learner, pred){
  cdfs <- learner$resids(learner, pred)[,"resid"]
  res <- dhsic.test(pred$x, cdfs)$statistic
  return(res)
}
pinball <- function(learner, pred){
  taus <-  getHyperPar(learner, "taus")  
  preds <- pred$gyh
  pinballs <- sapply(1:(length(taus)+2), function(i) {
    # i <- 5
    loss(preds[,i], pred$gy, loss.type=list(name="pinball", tau=c(3/(nrow(preds)), taus, 1-3/(nrow(preds)))[i]))
  }, simplify="array")
  
  pinballs <- apply(pinballs, 2, function(col) col/sd(col, na.rm=T))
  pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
  return(pinball)
  
}


negCE <- function(learner, pred){
  res <- LogLoss(pred$gyh, pred$gy)
  return(res)
}

# only for kernel methods - evaluates regularizer t(alpha)%*%Kxs%*%alpha
regL <- function(learner, pred){
  
  alpha <- learner$learnParams$alpha
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, pred$x, pars=parsXs)
  
  return(as.numeric(t(alpha)%*%Kxs%*%alpha))
}

qhsicLoss <- function(learner, pred){
  
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, pred$x, pars=parsXb)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  
  res <- t(pred$gyh)%*%Kxbc%*%pred$gyh + t(pred$gy)%*%Kxbc%*%pred$gy - 2*t(pred$gyh)%*%Kxbc%*%pred$gyh
  return(res)
}
nqhsicLoss <- function(learner, pred){
  
  resids <- learner$resids(learner, pred)[,"resid"]
  #resids <- pred$gy-pred$gyh
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, pred$x, pars=parsXb)
  N <- nrow(Kxb)
  Kr <- kernelMatrix("kern_lin", pred$x, pars=list(offset=0))
  
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Krc <- H%*%Kr%*%H
  res <- sum(diag(Kxbc%*%Krc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krc%*%Krc))))
  return(res)
}
hsicLoss  <- function(learner, pred){
  
  #print("enters hsicLoss")
  
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, pred$x, pars=parsXb)
  
  
  resids <- learner$resids(learner, pred)[,"resid"] #pred$gy-pred$gyh
  
  parsRg <- getKernelPars(learner, kernelName="kernelRg")
  Krg <- kernelMatrix(learner$hyperParams$non_data$kernelRg$name, resids, pars=parsRg)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- Kxb%*%H
  Krgc <- Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc)) 
  #print("exits hsicLoss")
  return(res)
}
nhsicLoss  <- function(learner, pred){
  
  #print("enters hsicLoss")
  
  parsXb <- getKernelPars(learner, kernelName="kernelXb")
  Kxb <- kernelMatrix(learner$hyperParams$non_data$kernelXb$name, pred$x, pars=parsXb)
  
  resids <- learner$resids(learner, pred)[,"resid"]
  #resids <- pred$gy-pred$gyh
  parsRg <- getKernelPars(learner, kernelName="kernelRg")
  Krg <- kernelMatrix(learner$hyperParams$non_data$kernelRg$name, resids, pars=parsRg)
  N <- nrow(Kxb)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  Krgc <- H%*%Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc))/(sqrt(sum(diag(Kxbc%*%Kxbc)))*sqrt(sum(diag(Krgc%*%Krgc))))
  #print("exits hsicLoss")
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
  learnerAux$hyperParams$non_data$kernelXb <- list(name="kern_rbf", pars=c(sigma="beta"))
  learnerAux$hyperParams$non_data$kernelYhk <- list(name="kern_rbf", pars=c(sigma="kappa"))
  
  
  parsXb <- getKernelPars(learnerAux, kernelName="kernelXb")
  Kxb <- kernelMatrix(learnerAux$hyperParams$non_data$kernelXb$name, x[,indxSens], pars=parsXb)
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  
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
  learnerAux$hyperParams$non_data$kernelXb <- list(name="kern_rbf", pars=c(sigma="beta"))
  learnerAux$hyperParams$non_data$kernelYhk <- list(name="kern_rbf", pars=c(sigma="kappa"))
  
  parsXb <- getKernelPars(learnerAux, kernelName="kernelXb")
  Kxb <- kernelMatrix(learnerAux$hyperParams$non_data$kernelXb$name, x[,indxSens], pars=parsXb)
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  
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
  learnerAux$hyperParams$non_data$kernelYhk <- list(name="kern_rbf", pars=c(sigma="kappa"))
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  
  
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
  
  parsXb <- getKernelPars(learnerAux, kernelName="kernelXb")
  Kxb <- kernelMatrix(learnerAux$hyperParams$non_data$kernelXb$name, x[,indxSens], pars=parsXb)
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  
  
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
  
  parsXb <- getKernelPars(learnerAux, kernelName="kernelXb")
  Kxb <- kernelMatrix(learnerAux$hyperParams$non_data$kernelXb$name, x[,indxSens], pars=parsXb)
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  

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
  
  parsYhk <- getKernelPars(learnerAux, kernelName="kernelYhk")
  Kyhk <- kernelMatrix(learnerAux$hyperParams$non_data$kernelYhk$name, pred$gyh, pars=parsYhk)
  
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

#############################################################################################################*
# Kernel CMEM learners 
#############################################################################################################*

#############################################################################################################*
# learn and predict functions for different kernel cmem learners
learn.cmem_L2 <- function(learner, plot=FALSE) {
  
  parsX <- getKernelPars(learner, kernelName="kernelX")
  parsY <- getKernelPars(learner, kernelName="kernelY")
  
  x <- learner$hyperParams$trainData$x
  y <- as.matrix(learner$hyperParams$trainData$y)
  
  kernelX_char <- learner$hyperParams$non_data$kernelX$name
  kernelX_char <- strsplit(kernelX_char, split="T")[[1]][1]
  kernelY_char <- learner$hyperParams$non_data$kernelY$name
  kernelY_char <- strsplit(kernelY_char, split="T")[[1]][1]
  Lx  <- kernelMatrix(kernelX_char, x, pars=parsX) 
  Ky  <- kernelMatrix(kernelY_char, y, pars=parsY)
  Cks <- kernelGradNormMatrix(kernel=kernelX_char, x,  pars=parsX)
  
  
  n <- nrow(Lx)
  lambda <- getHyperPar(learner, "lambda") #*n

  n <- nrow(x)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  if(getHyperPar(learner, "centerLx")) Lx <- H%*%Lx%*%H
  if(getHyperPar(learner, "centerKy")) Ky <- H%*%Ky%*%Ky
  
  BlambdaAux <- try(solve(Lx+n*lambda*I))
  if(class(BlambdaAux)=="try-error") BlambdaAux <- I
  
  learner$learnParams$Lx <- Lx
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- BlambdaAux
  learner$learnParams$Cks <- Cks
  
  return(learner)
}
predict.cmem_L2 <- function(learner, data) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelX")
  parsYg <- getKernelPars(learner, kernelName="kernelY")
  
  trainData <- learner$hyperParams$trainData
  kernelX_char <- learner$hyperParams$non_data$kernelX$name
  kernelX_char <- strsplit(kernelX_char, split="T")[[1]][1]
  kernelY_char <- learner$hyperParams$non_data$kernelY$name
  kernelY_char <- strsplit(kernelY_char, split="T")[[1]][1]
  lx <- kernelMatrix(kernelX_char, x=data$x, y=trainData$x, pars=parsXs)
  ky <- kernelMatrix(kernelY_char, x=as.matrix(data$y), y=as.matrix(trainData$y), pars=parsYg)
  Ky <- kernelMatrix(kernelY_char, x=as.matrix(data$y), pars=parsYg)
  
  nTe <- nrow(as.matrix(data$y))
  nTr <- nrow(trainData$x)
  Ite <- diag(nTe)
  Itr <- diag(nTr)
  HTe <- Ite-matrix(1/nTe,nTe,nTe)
  HTr <- Itr-matrix(1/nTr,nTr,nTr)
  if(getHyperPar(learner, "centerLx")) lx <- HTe%*%lx%*%HTr
  if(getHyperPar(learner, "centerKy")) ky <- HTe%*%ky%*%HTr
  
  KTr  <- learner$learnParams$Ky
  
  Blambda <- learner$learnParams$Blambda
  Alambda <- Blambda%*%KTr%*%Blambda
  LAL <- lx%*%Alambda%*%t(lx)
  LBK <- lx%*%Blambda%*%t(ky) 
  
  return(list(dLAL=diag(LAL), dLBK=diag(LBK), dKy=diag(Ky)))
  #return(list(lx=lx, ky=ky, Ky=Ky))
}

learn.cmfm_L2 <- function(learner, plot=FALSE) {
  
  parsX <- getKernelPars(learner, kernelName="kernelX")
  
  
  x <- learner$hyperParams$trainData$x
  y <- as.matrix(learner$hyperParams$trainData$y)
  
  kernelX_char <- learner$hyperParams$non_data$kernelX$name
  kernelX_char <- strsplit(kernelX_char, split="T")[[1]][1]
  Lx  <- kernelMatrix(kernelX_char, x, pars=parsX) 
  Cks <- kernelGradNormMatrix(kernel=kernelX_char, x,  pars=parsX)
  
  n <- nrow(Lx)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  
  parsY <- getKernelPars(learner, kernelName="featureY")
  auxParsY <- parsY
  auxParsY$x <- y
  phiy_char <- learner$hyperParams$non_data$featureY$name
  phiy_char <- strsplit(phiy_char, split="T")[[1]][1]
  aux <- strsplit(phiy_char, split="_")[[1]]
  if(aux[1] == "rff") phiy_char <- "rff"
  phiy  <- do.call(phiy_char, auxParsY)
  
  if(getHyperPar(learner, "centerLx")) Lx <- H%*%Lx%*%H
  if(getHyperPar(learner, "centerKy")){
    #phiy <- H%*%phiy
    #  phiy <- apply(phiy,2,norml)
    phiy <- H%*%phiy
  }
  
  Ky <- phiy%*%t(phiy)
  lambda <- getHyperPar(learner, "lambda") #*n
  Blambda <- solve(Lx+n*lambda*I)
  
 
  alpha0 <- Blambda%*%phiy
  
  alphaF <- NULL
  if(FALSE){
  
  lossKCDC <- function(alphaVec, Lx, phiy, lambda, n){
    alpha <- matrix(alphaVec, n, ncol(phiy))
    aux <- Lx%*%alpha%*%t(alpha)%*%Lx
    res1 <- -2*sum(diag(Lx%*%alpha%*%t(phiy)))
    res2 <- sum(diag(aux)) 
    
    res3 <- n*lambda*sum(diag(aux*aux)) 
    #res3 <- lambda*n*diag(aux*aux)[3]
    #i<- 2; res3 <- n*diag(aux*aux)[i]
    
    res4 <- -lambda*sum(diag(aux))^2
    
    return(res1+res2+res3+res4)
  } 
  
  gradKCDC <- function(alphaVec, Lx, phiy, lambda, n){
    alpha <- matrix(alphaVec, n, ncol(phiy))
    aux <- Lx%*%alpha%*%t(alpha)%*%Lx
    
    #aux2 <- do.call(cbind, lapply(1:n, function(i) Lx[,i, drop=F] %*% Lx[i,, drop=F]))
    
    
    
    In <- diag(n)
    
    res1 <- -2*Lx%*%phiy # right!
    res2 <- 2*Lx%*%Lx%*%alpha #+ right! 
    
    #res3 <- 4*lambda*n*aux2%*%((alpha%*%t(alpha))%x%In)%*%t(aux2) # wrong! 
    #res3 <- 4*lambda*n*(Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha wrong!
    
    # i <- 2
    # Di <- rep(0, n)
    # Di[i] <- 1
    # Di <- diag(Di)
    # res3 <- 4*n*(Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha #right-ish!
    
    # res3 <- sapply(1:n, function(i){
    #   Di <- rep(0, n)
    #   Di[i] <- 1
    #   Di <- diag(Di)
    #   res <- (Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha
    # }, simplify="array")

    res3 <- sapply(1:n, function(i){
      res <- as.numeric(Lx[i,]%*%alpha%*%t(alpha)%*%Lx[,i])* (Lx[,i,drop=F]%*%Lx[i,,drop=F])%*%alpha
    }, simplify="array")
    
    res3 <- 4*n*lambda*apply(res3, c(1,2), sum)
    
    res4 <- -4*lambda*sum(diag(aux))*Lx%*%Lx%*%alpha #right !
    
    res <- res1+res2+res3+res4
    
    return(as.numeric(res))
  }
  
  #lossKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
  #length(gradKCDC(alphaVec=as.numeric(alpha0), Lx, phiy, lambda, n))
  
  #library(numDeriv)
  
  # alpha0 <- matrix(seq(12),3,4)
  # Lx <- matrix(rnorm(9),3,3)
  # Lx <- Lx+t(Lx)
  # n <- 3
  # phiy <- matrix(runif(12),3,4)
  
  # pm <- proc.time()
  # Jnum <- jacobian(lossKCDC, as.numeric(alpha0), Lx=Lx, phiy=phiy, lambda=lambda, n=n)
  # proc.time() - pm # 2 mins
  # Jact <- gradKCDC(alphaVec=as.numeric(alpha0), Lx=Lx, phiy=phiy, lambda=lambda, n=n)
  # length(Jact); dim(Jnum)
  # plot(Jact, Jnum); abline(a=0, b=1, col="red")
  
  # optimizar con esa inicialización pero mas iteraciones
   res <- lbfgs(lossKCDC, gradKCDC, vars=as.numeric(alpha0), Lx=Lx, phiy=phiy, lambda=lambda, n=n, invisible=1, max_iterations=100, epsilon=1e-5)
   alphaF <- matrix(res$par, n, ncol(phiy))
  # 
  # #lossKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
  # #lossKCDC(as.numeric(alphaF), Lx, phiy, lambda, n)
  # 
  # grad0 <- gradKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
  # gradF <- gradKCDC(as.numeric(alphaF), Lx, phiy, lambda, n)
  # plot(grad0, gradF)
  # t(grad0)%*%grad0
  # t(gradF)%*%gradF
  
  
  
  #plot(alpha0, alphaF); abline(a=0, b=1, col="red")
  #max(abs(alpha0-alphaF))
  
  #Blambda <- alphaF%*%t(phiy)%*%solve(Ky)  
  
  }
   
   if(plot){
     # plot and save
     # a) data x vs y
     # b) E[g1(y)|x_i] vs E[g2(y)|x_i]
     # c) E[g1(y)|x_i] and E[g2(y)|x] vs x (on top of data?)
     # d) alpha(x_i) for 5 i on grid
     #a)
     # plot(x,y)
     # b)
     
     ordx <- order(x)
     ordy <- order(y)
     xs <- x[ordx,1]
     ys <- y[ordx,1]
     phiys <- phiy[ordx,]
     
     
     #d)
     alphas <- t(Lx)%*%Blambda
     #alphas <-  try(solve(Ky, phiy%*%t(alphaF)%*%Lx))
     if(!class(alphas)=="try-error"){
      alphass <- alphas[ordx,]
      alphass <- alphass[,ordx]
     } else{
       alphass <- NA
     }
     
     
     
     # c)
     #EphiyGX <- Lx%*%alphaF
     EphiyGX <- alphas%*%phiy
     EphiyGXs <- EphiyGX[ordx,]
     
     repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_causaLearner/experiments/CMFM/results/toyPlots/"
     #repos <- "/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_causaLearner/experiments/CMFM/results/toyPlots/"
     #repos <- "/home/soulivanh/Documents/proyectos/CAUSALITY/toyPlots/"
     
     kernCombo <- paste(strsplit(learner$hyperParams$non_data$kernelX$name, "_")[[1]][2], strsplit(learner$hyperParams$non_data$featureY$name, "_")[[1]][2], sep="_")
     lambda <- log(learner$hyperParams$data$non_optimizable$lambda$val,10)
     
     dir(repos)
     dir.exists(repos)
     
     folder <- paste(repos, kernCombo, sep="")
     dir.exists(folder)
     dir.create(folder)
     
     files <- dir(folder)
     
     if(length(files) %in% c(0,1)){
       dataset <- 1
     } else{
       datasetsDone <- as.numeric(sapply(strsplit(files, "_"), function(el) el[length(el)-2]))
       #d1 <- duplicated(datasetsDone)
       #d2 <- rev(duplicated(rev(datasetsDone)))
       #d <- d1|d2
       
       #datasetsDone <- datasetsDone[which(d)]
       
       tabDataSetsDone <- table(datasetsDone)
       
       dataSetsStarted <- as.numeric(names(tabDataSetsDone))
       print(paste("dataset started:", dataSetsStarted))
       
       datasetsDone <- as.numeric(names(tabDataSetsDone))[which(tabDataSetsDone==18)]
       
       if(length(datasetsDone)>0){
         dataset <- max(datasetsDone)+1
       } else{
         dataset <- max(dataSetsStarted)
       }
       print(paste("dataset: ", dataset))
     }
     
     if(paste(paste(kernCombo, "l", lambda, "ds", dataset,"dir","yx", sep="_"), "png", sep=".") %in% files){
       dir <- "xy"
     } else{
       dir <- "yx"
     }
     
     
     Alambda <- Blambda%*%Ky%*%t(Blambda)
     #Alambda <- alphaF%*%t(alphaF)
     LAL <- Lx%*%Alambda%*%t(Lx)
     cmplxs <- diag(LAL)
     cmplxss <- cmplxs[ordx]
     
     # start to plot
     
     png(file=paste(folder,"/",paste(paste(kernCombo, "l", lambda, "ds", dataset,"dir",dir, sep="_"), "png", sep="."), sep=""))
     par(mfrow=c(3,3))
     plot(xs, ys)
     
     # rffs
     plot(y[ordy,1], phiy[ordy,1], type="l", col="red", main=paste(kernCombo, " lam: ", lambda))
     lines(y[ordy,1], phiy[ordy,2], col="blue")
     
     # b) 
     plot(phiy[,1], phiy[,2])
     
     # c)
     
     absRelErr1 <- sum(abs(phiys[,1]-EphiyGXs[,1])/abs(phiys[,1]))
     plot(xs, phiys[,1], main=paste("AbsRelErr: ", formatC(absRelErr1,digits=4)))
     lines(xs, EphiyGXs[,1], col="red")
     
     absRelErr2 <- sum(abs(phiys[,2]-EphiyGXs[,2])/abs(phiys[,2]))
     plot(xs, phiys[,2], main=paste("AbsRelErr: ", formatC(absRelErr2,digits=4)))
     lines(xs, EphiyGXs[,2], col="blue")
     
     dPhiy <- diag(phiys%*%t(phiys))
     dEphiy <- diag(EphiyGXs%*%t(EphiyGXs))
     # ||phi(y)|| vs ||phi(y)_hat||
     absRelErr3 <- sum(abs(dPhiy-dEphiy)/abs(dPhiy))
     plot(xs, dPhiy, main=paste("AbsRelErr: ", formatC(absRelErr3,digits=4)), ylim=range(dPhiy, dEphiy))
     lines(xs, dEphiy, col="red")
     
     # xs vs complexity
     plot(xs, cmplxss)
     
     # d)
     Lss <- Lx[ordx,]
     Lss <- Lss[,ordx]
     barplotAlpha <- function(xIndx, xs, alphass, Lss){
       # xIndx <- 50
       wdth <- min(diff(xs))/2
       cumSumSp <- xs/wdth+seq(0,length(xs)-1)*wdth+wdth*0.5
       sp <- cumSumSp[2:length(cumSumSp)]-cumSumSp[1:(length(cumSumSp)-1)]
       sp <- c(min(sp), sp)
       sp <- (sp-min(sp))/(max(sp)-min(sp))
       sp <- sp*100+2*wdth
       
       
       
       normAlphas <- alphass[xIndx,]*alphass[xIndx,]
       normAlphas <- normAlphas / sum(normAlphas)
       logNormAlphas <- log(normAlphas)
       entAlphass <- -sum(normAlphas*logNormAlphas)
       
       maxAlphaIndx <- which.max(normAlphas)
       dists <- Lss[xIndx, maxAlphaIndx]
       
       barplot(alphass[xIndx,], width=wdth, space=sp, col=c("red","blue"), 
               main=paste("Ent:", formatC(entAlphass, digits=4), "dist: ", formatC(dists, digits=4)))
       xsT <- cumsum(sp)*wdth+seq(0,length(xs)-1)*wdth+wdth*0.5
       abline(v=xsT[xIndx], col="red")
       abline(h=0, col="red")
       #abline(v=(x-min(xs))/(max(xs)-min(xs))*(max(xsT)-min(xsT))+min(xsT), col="green")
       
     }
     
     if(!is.na(alphass)){
      barplotAlpha(xIndx=10, xs, alphass, Lss)
      # barplotAlpha(50, xs, alphass, Lss)
      barplotAlpha(90, xs, alphass, Lss)
     }
     
     dev.off()
     par(mfrow=c(1,1))
   }
   
   
  learner$learnParams$Lx <- Lx
  learner$learnParams$phiy <- phiy
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- Blambda
  learner$learnParams$alpha <- alphaF
  learner$learnParams$Cks <- Cks
  
  return(learner)
}
predict.cmfm_L2 <- function(learner, data) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelX")
  parsYg <- getKernelPars(learner, kernelName="featureY")
  
  trainData <- learner$hyperParams$trainData
  kernelX_char <- learner$hyperParams$non_data$kernelX$name
  kernelX_char <- strsplit(kernelX_char, split="\\.")[[1]][1]
  lx <- kernelMatrix(kernelX_char, x=data$x, y=trainData$x, pars=parsXs)
  
  phiy_char <- learner$hyperParams$non_data$featureY$name
  phiy_char <- strsplit(phiy_char, split="T")[[1]][1]
  aux <- strsplit(phiy_char, split="_")[[1]]
  if(aux[1] == "rff") phiy_char <- "rff"
  
  y <- as.matrix(data$y)
  auxParsYg <- parsYg
  auxParsYg$x <- y
  
  phiy  <- do.call(phiy_char, auxParsYg)
  
  nTe <- nrow(y)
  nTr <- nrow(trainData$x)
  Ite <- diag(nTe)
  Itr <- diag(nTr)
  HTe <- Ite-matrix(1/nTe,nTe,nTe)
  HTr <- Itr-matrix(1/nTr,nTr,nTr)
  if(getHyperPar(learner, "centerLx")) lx <- HTe%*%lx%*%HTr
  if(getHyperPar(learner, "centerKy")){
    #phiy <- HTe%*%phiy
    phiy <- apply(phiy,2,norml)
    phiy <- HTe%*%phiy
  }
  
  # test vs test
  Ky <- phiy%*%t(phiy)
  # test vs train
  ky <- phiy%*%t(learner$learnParams$phiy)
  
  KTr  <- learner$learnParams$Ky
  
  Blambda <- learner$learnParams$Blambda
  alpha <- learner$learnParams$alpha
  
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
    LBK <- lx%*%alpha%*%t(phiy) 
    EphiyGX <- lx%*%alpha
  } else{
    Alambda <- Blambda%*%KTr%*%t(Blambda)
    LBK <- lx%*%Blambda%*%t(ky) 
    alphas <- lx%*%Blambda
    EphiyGX <- alphas%*%learner$learnParams$phiy
  }
  
  
  LAL <- lx%*%Alambda%*%t(lx)
  
  #LBK <- lx%*%Blambda%*%t(ky) 
  #LBK <- lx%*%Blambda%*%t(phiy%*%t(learner$learnParams$phiy)) 
  #LBK <- lx%*%Blambda%*%learner$learnParams$phiy%*%t(phiy) 
  
  
  #alphas <- lx%*%Blambda
  
  #EphiyGX <- alphas%*%learner$learnParams$phiy
  #EphiyGX <- lx%*%Blambda%*%learner$learnParams$phiy
  
  
  #plot(data$x[order(data$x)], phiy[order(data$x),1])
  #lines(data$x[order(data$x)], EphiyGX[order(data$x),1], col="red")
  
  
  dEphiy <- EphiyGX%*%t(EphiyGX)
  
  #plot(data$x[order(data$x)], diag(Ky)[order(data$x)])
  #lines(data$x[order(data$x)], diag(dEphiy)[order(data$x)], col="red")
  
  return(list(dLAL=diag(LAL), dLBK=diag(LBK), dKy=diag(Ky), dEphiy=diag(dEphiy)))
  #return(list(lx=lx, ky=ky, Ky=Ky))
}

learn.cmfm_bin <- function(learner, plot=FALSE, logReg=TRUE) {
  
  x <- learner$hyperParams$trainData$x
  y <- as.matrix(learner$hyperParams$trainData$y)
  
  numBins <- getHyperPar(learner, "numBins")
  xqs <- quantile(as.numeric(x), probs=seq(0,1, length.out=numBins+2))
  #xqs <- xqs[2:(length(xqs)-1)]

  
  n <- nrow(x)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  phiy <- learner$makeFeature(learner, data=learner$hyperParams$trainData, "Y")
  
  if(getHyperPar(learner, "centerKy")){
    phiyc <- H%*%phiy
    
  } else{
    phiyc <- phiy
  }
  
  Ky <- phiyc%*%t(phiyc)
  
  learner$learnParams$xqs <- xqs
  learner$learnParams$Ky <- Ky
  learner$learnParams$phiy <- phiyc

  # now train classifiers (one for each bin) 
  # to distinguish between p(y| x in bin) vs p_fake
  
  # wd be nice if I could sometimes call the learn function
  # without running the logistic regression which is to evaluate
  # loss but since I call it from within CV which is generic 
  # for any learner I would need to add a dummy parameter
  # to all learners, maybe later
  if(logReg){
    #prepare fake data for classifier training
    NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
    classifier <- NCE_learner$classifier
    fakeDist <-  NCE_learner$fakeDist
    fakeDistPars <- NCE_learner$fakeDistPars
    kappa <- NCE_learner$kappa
  
    fits <- lapply(1:(length(xqs)-1), function(i){
    
    indxq <- which(x >= xqs[i] & x<= xqs[i+1])
    xq <- x[indxq]  
    yq <- y[indxq] 
    nq <- length(indxq)
    num_fake <- nq*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[-indxq] for KDE non parametric dist
    if("x" %in% names(fakeDistPars)) fakeDistPars$x <- x[-indxq]
    fakeDistPars <- c(fakeDistPars, list(n=num_fake))
    y_fake <- as.matrix(do.call(fakeDist, fakeDistPars))
    
    trainDataAux <- learner$hyperParams$trainData
    trainDataAux$y <- y_fake
    logRegPhiy_fake <- learner$makeLogRegFeature(learner, data=trainDataAux, "Y")
    
    
    if(getHyperPar(learner, "centerKy")){
      Hfake <- I-matrix(apply(logRegPhiy, 2, mean), num_fake, ncol(logRegPhiy), byrow=T)
      logRegPhiy_fake <- Hfake%*%logRegPhiy
    }
    # make train classification data
    x_tr <- as.matrix(rbind(logRegPhiyc[indxq,,drop=F], logRegPhiy_fake))
    y_tr <- c(rep(1, nq), rep(0, num_fake)) 
    set.seed(12)
    smpl <- sample(nrow(x_tr)) #jumble them up
    x_tr <- x_tr[smpl,,drop=F]
    y_tr <- y_tr[smpl]
    trainClassData <- constructData(x_tr, y_tr)
  
    # train classifier
    classifier <- eval(parse(text=classifier))
    classifier <- setParams(learner=classifier, trainData=trainClassData)
    fit <- classifier$learn(learner=classifier)
    return(fit)
  })
    learner$learnParams$classifier <- fits
  }
    
  return(learner)
}
predict.cmfm_bin <- function(learner, data) {
  
  parsYg <- getKernelPars(learner, kernelName="featureY")
  
  
  phiy_char <- learner$hyperParams$non_data$featureY$name
  phiy_char <- strsplit(phiy_char, split="T")[[1]][1]
  aux <- strsplit(phiy_char, split="_")[[1]]
  if(aux[1] == "rff") phiy_char <- "rff"

  y <- as.matrix(data$y)
  x <- data$x
  
  auxParsYg <- parsYg
  auxParsYg$x <- y
  
  phiy  <- do.call(phiy_char, auxParsYg)
  
  nTe <- nrow(y)
  Ite <- diag(nTe)
  HTe <- Ite-matrix(1/nTe,nTe,nTe)
  
  if(getHyperPar(learner, "centerKy")){
    #phiy <- HTe%*%phiy
    #phiy <- apply(phiy,2,norml)
    phiy <- HTe%*%phiy
  }
  
  # test vs test
  Ky <- phiy%*%t(phiy)
  
  #prepare fake data for classifier training
  
  NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
  classifier <- learner$learnParams$classifier
  fakeDist <-  NCE_learner$fakeDist
  fakeDistPars <- NCE_learner$fakeDistPars
  kappa <- NCE_learner$kappa
  
  xqs <- learner$learnParams$xqs
  preds <- lapply(1:(length(xqs)-1), function(i){
    # i <- 1
    #print(i)
    indxq <- which(x >= xqs[i] & x<= xqs[i+1])
    xq <- x[indxq, ,drop=F]  
    yq <- y[indxq] 
    nq <- length(indxq)
    num_fake <- nq*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[-indxq] for KDE non parametric dist
    if("x" %in% names(fakeDistPars)) fakeDistPars$x <- x[-indxq]
    fakeDistPars <- c(fakeDistPars, list(n=num_fake))
    y_fake <- as.matrix(do.call(fakeDist, fakeDistPars))
    auxParsYg$x <- y_fake
    phiy_fake <- do.call(phiy_char, auxParsYg)
    if(getHyperPar(learner, "centerKy")){
      Hfake <- I-matrix(apply(phiy, 2, mean), num_fake, ncol(phiy), byrow=T)
      phiy_fake <- Hfake%*%phiy
    }
    # make train classification data
    x_te <- as.matrix(rbind(phiy[indxq,,drop=F], phiy_fake))
    y_te <- c(rep(1, nq), rep(0, num_fake)) 
    set.seed(12)
    testClassData <- constructData(x_te, y_te)
    
    # train classifier
    classifier1 <- classifier[[i]]
    pred <- classifier1$pred(learner=classifier1, testClassData)
    return(pred)
  })
  
  xs <- lapply(preds, function(el) el$x)
  xs <- do.call(rbind, xs)
  gy <- lapply(preds, function(el) el$gy)
  gy <- do.call(rbind, gy)
  gyh <- lapply(preds, function(el) el$gyh)
  gyh <- do.call(rbind, gyh)
  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(gyh)))
  
}

  
#############################################################################################################*
# Loss functions for cross validation and evaluation of learners

NCE_loss <- function(learner, pred){
  n <- nrow(x)
  p <- ncol(x)
  midPt <- round(n/2)
  indx.tr <- 1:midPt
  indx.te <- (midPt+1):n
  n_tr <- length(indx.tr)
  n_te <- length(indx.te)
  
  # make train classification data
  x_tr <- as.matrix(rbind(x[indx.tr,, drop=F], matrix(runif(n_tr*p), n_tr, p)))
  y_tr <- c(rep(1, n_tr), rep(0, n_tr)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr))
  trainClassData <- constructData(x_tr[smpl,,drop=F], y_tr[smpl])
  
  # make test classification data
  x_te <- as.matrix(rbind(x[indx.te,,drop=F], matrix(runif(n_te*p), n_te, p)))
  y_te <- c(rep(1, n_te), rep(0, n_te)) 
  set.seed(123)
  smpl <- sample(nrow(x_te))
  testClassData <- constructData(x_te[smpl, ,drop=F], y_te[smpl])
  
  # train classifier
  classifier <- setParams(learner=classifier, trainData=trainClassData)
  
  # getHyperPar(classifier, "sigma")
  keepHypers <- lapply(keepHyperParams, function(par) getHyperPar(classifier, par))
  
  # NOTE: need to check if giving each group of dists a different seed is better for 
  # KCDC classification as I observed for sinx example with discreete bins
  data <- constructData(x, c(y_tr, y_te))
  classifier <- classifier$learn(learner=classifier)
  
  phix <- classifier$makeFeature(learner=classifier, data, "X")
  
  meanPhix <- apply(phix, 2, mean)
  
  pred <- classifier$pred(classifier, testClassData)
  loss <- do.call(classifier$optimizeParams$losses[[1]], list(learner=classifier, pred=pred))
  
  res <- list(meanPhix=meanPhix, loss=loss, hyperPars=keepHypers)
  
  if(ncol(x)==1){
    xx <- seq(0, 1, length.out=100)
    pdfData <- constructData(as.matrix(xx), rep(0,100))
    predPdf <- classifier$pred(learner=classifier, data=pdfData)
    pdf <- predPdf$gyh/(1-predPdf$gyh)
    res <- c(res, list(pdf=pdf))
    numBins <- round(n/10)
    if(plot){
      hist(x, numBins, prob=T)
      lines(xx, pdf, col="red")
    }
  }
  
  return(res)
}


cmem_L2 <- function(learner, pred){
  
  #LTr  <- learner$learnParams$Lx
  #KTr  <- learner$learnParams$Ky
  
  #Blambda <- learner$learnParams$Blambda
  #Alambda <- Blambda%*%KTr%*%Blambda
  
  #LTeTr  <- pred$lx 
  #KTeTr  <- pred$ky
  #KTe <- pred$Ky 
  
  #LAL <- LTeTr%*%Alambda%*%t(LTeTr)
  #LBK <- LTeTr%*%Blambda%*%t(KTeTr)
  
  res <- mean(pred$dLAL + pred$dKy - 2*pred$dLBK)
  #res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}

cmem_L2_sd <- function(learner, pred){
  
  #LTr  <- learner$learnParams$Lx
  #KTr  <- learner$learnParams$Ky
  
  #Blambda <- learner$learnParams$Blambda
  #Alambda <- Blambda%*%KTr%*%Blambda
  
  #LTeTr  <- pred$lx 
  #KTeTr  <- pred$ky
  #KTe <- pred$Ky 
  
  #LAL <- LTeTr%*%Alambda%*%t(LTeTr)
  #LBK <- LTeTr%*%Blambda%*%t(KTeTr)
  
  res <- sd(pred$dLAL + pred$dKy - 2*pred$dLBK)
  #res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}


cmem_L2_rel <- function(learner, pred){
  
  #LTr  <- learner$learnParams$Lx
  #KTr  <- learner$learnParams$Ky
  
  #Blambda <- learner$learnParams$Blambda
  #Alambda <- Blambda%*%KTr%*%Blambda
  
  #LTeTr  <- pred$lx 
  #KTeTr  <- pred$ky
  #KTe <- pred$Ky 
  
  #LAL <- LTeTr%*%Alambda%*%t(LTeTr)
  #LBK <- LTeTr%*%Blambda%*%t(KTeTr)
  
  #dKTe <- diag(KTe)
  
  res <- pred$dLAL + pred$dKy - 2*pred$dLBK
  res <- sum(res/pred$dKy)
  #res <- diag(LAL) + dKTe - 2*diag(LBK)
  #res <- sum(res/dKTe)
  
  return(res)
}

# Training relative error wrt to norm of phi(y)

TNRE <-function(learner, pred){
  
  res <- sqrt(mean((pred$dKy-pred$dEphiy)^2)) # / pred$dKy^2
  return(res)
  
}

TNRSE <-function(learner, pred){
  
  res <- sd((pred$dKy-pred$dEphiy)^2) # / pred$dKy^2
  return(res)
  
}

numUnd <-function(learner, pred){
  
  res <- abs(0.5-sum(pred$dEphiy<pred$dKy)/length(pred$dKy))
  return(res)
  
}


#############################################################################################################*
# measure functions for different cmem learners

calcMsrs <- function(cmemLearner){
  msrs_char <- sapply(cmemLearner$msrs, function(el) el$func)
  msrs_val <- sapply(msrs_char, function(el){
    # el <- "KCDC"
     print(paste("measure: ", el))
    res <- do.call(el, list(learner=cmemLearner))
    })
  names(msrs_val) <- msrs_char
  return(msrs_val)
}

# Kernel Conditional Deviance for Causal inference (KCDC)

KCDC <- function(learner, pred=NULL){
  if(is.null(getHyperPar(learner, "numBins"))){
    res <- KCDC_cont(learner)
  } else{
    res <- KCDC_bins(learner)
  }
  return(res)
}
KCDC_cont <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
  } else{
    Blambda <- learner$learnParams$Blambda
    Alambda <- Blambda%*%K%*%t(Blambda)
  }
    
  LAL <- L%*%Alambda%*%t(L)
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  res <- (c/n) - (b/n)^2
  
  return(res)
}
KCDC_bins <- function(learner){
  xqs <- learner$learnParams$xqs
  Ky <- learner$learnParams$Ky
  meanKys <- sapply(1:(length(xqs)-1), function(i){
    # i <- 1
    #print(i)
    indxq <- which(x >= xqs[i] & x<= xqs[i+1])
    Ky_block <- Ky[indxq,]
    Ky_block <- Ky_block[,indxq]
    meanKy <- mean(Ky_block)
    return(meanKy)
  })
  
  # # check
  # phiy <- learner$learnParams$phiy
  # meanPhiy <- sapply(1:(length(xqs)-1), function(i){
  #   # i <- 1
  #   #print(i)
  #   indxq <- which(x >= xqs[i] & x<= xqs[i+1])
  #   phiy_block <- phiy[indxq,]
  #   meanPhiy <- apply(phiy_block, 2, mean)
  #   res <- sum(meanPhiy*meanPhiy)
  #   return(res)
  # })
  # print("check")
  # print(var(sqrt(meanPhiy)))
  # print(var(sqrt(meanKys)))
  
  res <- var(sqrt(meanKys))
  return(res)
}


KCDCrel <- function(learner){
  x <- learner$hyperParams$trainData$x
  n <- nrow(x)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCDCrel"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCDC(learner)
  
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCDC(learnerAux) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}

KCDCpval <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCDCpval"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCDC(learner)
  
  
  KCDC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCDC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCDC_dist_null>mesr)/learner$msrs[["KCDCpval"]]$pars$numPerms
  
  
  return(pval)
}


# Mean complexity of functions
MCX <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
  } else{
    Blambda <- learner$learnParams$Blambda
    Alambda <- Blambda%*%K%*%t(Blambda)
  }
  
  LAL <- L%*%Alambda%*%t(L)
  
  b <- sum(diag(LAL)^(0.5))
  res <- b/n 
  
  return(res)
}

# Scale sensitive  effective dimension of data

EDML <- function(learner){
  
  L  <- learner$learnParams$Lx
  #K  <- learner$learnParams$Ky
  n <- nrow(L)
  #Blambda <- learner$learnParams$Blambda
  #Alambda <- Blambda%*%K%*%t(Blambda)
  
  L.svd <- svd(L)
  #B.svd <- svd(Blambda)
  
  #B.eig <- B.svd$d/max(B.svd$d)
  L.eig <- L.svd$d/max(L.svd$d)
  
  lambda <- getHyperPar(learner, "lambda")
  
  res <- n-sum(L.eig/(L.eig+n*lambda))
  #res2 <- sum(B.eig)
  
  return(res)
}

EDMB <- function(learner){
  
  L  <- learner$learnParams$Lx
  #K  <- learner$learnParams$Ky
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    #TEMPORARY!!!!!
    Blambda <- alpha%*%t(alpha)
  } else{
    Blambda <- learner$learnParams$Blambda
  }
  Blambda <- learner$learnParams$Blambda
  #Alambda <- Blambda%*%K%*%t(Blambda)
  
  #L.svd <- svd(L)
  B.svd <- svd(Blambda)
  
  B.eig <- B.svd$d/max(B.svd$d)
  #L.eig <- L.svd$d/max(L.svd$d)
  
  #lambda <- getHyperPar(learner, "lambda")
  
  #res <- sum(L.eig/(L.eig+n*lambda))
  res <- sum(B.eig)
  
  return(res)
}


# Training relative error 

TRE <- function(learner){
  data <- learner$hyperParams$trainData
  pred <- predict.cmfm_L2(learner, data)
  res <- cmem_L2_rel(learner, pred)
  return(res)
}

# CV - training relative norm error

CVTeRNE <- function(learner){
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,])
  
}
  
CVTeRNSE <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "cmem_L2_sd" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optimHyperParams <- learnerAux$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learnerAux$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner=learnerAux, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,])
  
}

teErrNorm <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "TNRE" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optimHyperParams <- learnerAux$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learnerAux$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner=learnerAux, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,])
  
}

trErrNorm <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "TNRE" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optimHyperParams <- learnerAux$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learnerAux$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner=learnerAux, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["train",,])
  
}

teNumUnd <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "numUnd" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optimHyperParams <- learnerAux$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learnerAux$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner=learnerAux, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,])
  
}

trNumUnd <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "numUnd" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optimHyperParams <- learnerAux$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learnerAux$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner=learnerAux, params, fac=1)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["train",,])
  
}

CVTrRNE <- function(learner){
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list()
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  
  params <- c(paramsList, lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)) #, learner$hyperParams$non_data
  
  params <- do.call(constructParams, params)
  grid <- CV.parallel(learner, params, fac=fac)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], paramsList, dimnms["var"])
  names(dimnms) <- c("trainTest", names(paramsList), "var")
  dimnames(grid) <- dimnms
  
  return(grid["train",,])
  
}


# Entropy of alphas
EAL <- function(learner){

  L <- learner$learnParams$Lx
  K <- learner$learnParams$Ky
  
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    phiy <- learner$learnParams$phiy
    Alambda <- alpha%*%t(alpha)
    alphas <- try(solve(K, phiy%*%t(alpha)%*%L))
  } else{
    Blambda <- learner$learnParams$Blambda
    #Alambda <- Blambda%*%K%*%t(Blambda)
    alphas <- t(L)%*%Blambda
  }
  
  if(class(alphas)=="try-error") return(NA)
  
  normAlphas <- alphas*alphas
  normAlphas <- normAlphas / apply(normAlphas, 1, sum)
  logNormAlphas <- log(normAlphas)
  ents <- apply(-normAlphas*logNormAlphas, 1, sum)
  res <- mean(ents)
  return(res)
}

# dist to max alpha(x)
DAL <- function(learner){
  
  L <- learner$learnParams$Lx
  K <- learner$learnParams$Ky
  
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    phiy <- learner$learnParams$phiy
    Alambda <- alpha%*%t(alpha)
    alphas <- try(solve(K, phiy%*%t(alpha)%*%L))
    
  } else{
    Blambda <- learner$learnParams$Blambda
    #Alambda <- Blambda%*%K%*%t(Blambda)
    alphas <- t(L)%*%Blambda
  }
  
  if(class(alphas)=="try-error") return(NA)
  
  normAlphas <- alphas*alphas
  maxAlphaIndx <- apply(normAlphas, 1, which.max)
  indxMat <- matrix(c(1:n,maxAlphaIndx), n, 2)
  dists <- L[indxMat]
  
  res <- mean(dists)
  return(res)
}


# Kernel Conditional Relative Deviance for Causal inference (KCRDC)
KCRDC <- function(learner, pred=NULL){
  if(is.null(getHyperPar(learner, "numBins"))){
    res <- KCRDC_cont(learner)
  } else{
    res <- KCRDC_bins(learner)
  }
  return(res)
}
KCRDC_cont <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
  } else{
    Blambda <- learner$learnParams$Blambda
    Alambda <- Blambda%*%K%*%t(Blambda)
  }
  
  
  LAL <- L%*%Alambda%*%t(L)
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  res <- (c/b) - (b/n) 
  
  return(res)
}
KCRDC_bins <- function(learner){
  xqs <- learner$learnParams$xqs
  Ky <- learner$learnParams$Ky
  
  meanKys <- sapply(1:(length(xqs)-1), function(i){
    # i <- 1
    #print(i)
    indxq <- which(x >= xqs[i] & x<= xqs[i+1])
    Ky_block <- Ky[indxq,]
    Ky_block <- Ky_block[,indxq]
    meanKy <- mean(Ky_block)
    return(meanKy)
  })
  meanKys <- sqrt(meanKys)
  
  # phiy <- learner$learnParams$phiy
  # meanPhiy <- sapply(1:(length(xqs)-1), function(i){
  #   # i <- 1
  #   #print(i)
  #   indxq <- which(x >= xqs[i] & x<= xqs[i+1])
  #   phiy_block <- phiy[indxq,]
  #   meanPhiy <- apply(phiy_block, 2, mean)
  #   res <- sum(meanPhiy*meanPhiy)
  #   return(res)
  # })
  # print("check")
  # meanPhiy <- sqrt(meanPhiy)
  # print(var(meanPhiy)/mean(meanPhiy))
  # print(var(meanKys)/mean(meanKys))
  
  
  res <- var(meanKys)/mean(meanKys)
  return(res)
}


KCRDCpval <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCRDCpval"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCRDC(learner)
  
  
  KCRDC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCDC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCRDC_dist_null>mesr)/learner$msrs[["KCRDCpval"]]$pars$numPerms
  
  
  return(pval)
}


# Kernel Conditional Mean distance  for Causal inference (KCMC)
KCMC <- function(learner, pred=NULL){
  if(is.null(getHyperPar(learner, "numBins"))){
    res <- KCMC_cont(learner)
  } else{
    res <- KCMC_bins(learner)
  }
  return(res)
}
KCMC_cont <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
  } else{
    Blambda <- learner$learnParams$Blambda
    Alambda <- Blambda%*%K%*%t(Blambda)
  }
  
  LAL <- L%*%Alambda%*%t(L)
  
  f <- sum(LAL)
  c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  res <- (c/n) - (f/n^2)
  
  return(res)
}
KCMC_bins <- function(learner){
  xqs <- learner$learnParams$xqs
  Ky <- learner$learnParams$Ky
  
  meanKys <- sapply(1:(length(xqs)-1), function(i) sapply(1:(length(xqs)-1), 
                                                          function(j){
    # i <- 1; j <- 2
    #print(i)
    indxq_i <- which(x >= xqs[i] & x<= xqs[i+1])
    indxq_j <- which(x >= xqs[j] & x<= xqs[j+1])
    Ky_block <- Ky[indxq_i,]
    Ky_block <- Ky_block[,indxq_j]
    meanKy <- mean(Ky_block)
    return(meanKy)
  }, simplify="array"), simplify="array")
  
  # phiy <- learner$learnParams$phiy
  # distMeanPhiy <- sapply(1:(length(xqs)-1), function(i){
  #   # i <- 1
  #   #print(i)
  #   indxq <- which(x >= xqs[i] & x<= xqs[i+1])
  #   phiy_block <- phiy[indxq,]
  #   meanPhiy_block <- apply(phiy_block, 2, mean)
  #   meanPhiy <- apply(phiy, 2, mean)
  #   distToMean <- as.numeric(dist(rbind(meanPhiy, meanPhiy_block)))^2 
  #   
  #   return(distToMean)
  # })
  # print("check")
  # print(mean(distMeanPhiy))
  # print(mean(diag(meanKys)) - mean(meanKys))
  # 
  
  res <- mean(diag(meanKys)) - mean(meanKys)
  return(res)
}


KCMCpval <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCMCpval"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCMC(learner)
  
  
  KCMC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCMC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCMC_dist_null>mesr)/learner$msrs[["KCMCpval"]]$pars$numPerms
  
  
  return(pval)
}

# Kernel Conditional Sensitivity for Causal inference (KCSC)
KCSC <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky 
  
  n  <- nrow(L)
  Blambda <- learner$learnParams$Blambda
  Cks <- learner$learnParams$Cks
  
  
  # pm <- proc.time()
  # E2 <- kernelDataArray2(K, as.matrix(y))
  # proc.time() - pm # 0.319
  
  
  
  #plot(E, E2); abline(a=0, b=1, col="red")
  
  
  
  Chat <- apply(Cks, c(1,2), sum)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  reps <- nrow(Blambda)/n
  H <- lapply(1:reps, function(el) H)
  H <- do.call(cbind, H)
  H <- lapply(1:reps, function(el) H)
  H <- do.call(rbind, H)
  
  if(getHyperPar(learner, "centerLx")) Chat <- H%*%Chat%*%H
  
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    phiy <- learner$learnParams$phiy
    D <- t(alpha)%*%Chat%*%alpha
    E <- t(phiy)%*%phiy
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    pm <- proc.time()
    E <- K%*%K
    proc.time() - pm # 0.002
  }
  
  
  pm <- proc.time()
  res <- sum(D*E)
  proc.time() - pm #0.001
  
  #pm <- proc.time()
  #res2 <- sum(diag(K%*%D%*%K))
  #proc.time() - pm # 0.002
  
  return(res)
}

KCSCpval <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCSCpval"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCSC(learner)
  
  
  KCSC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCSC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCSC_dist_null>mesr)/learner$msrs[["KCSCpval"]]$pars$numPerms
  
  
  return(pval)
}

# Kernel Conditional Complexity for Causal inference (KCCC)

# complexity measure: entropy of absoulute value of components of sum of gradient accross all data points
KCCC_ent <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky 
  phiy <- learner$learnParams$phiy
  n  <- nrow(L)
  m <- ncol(phiy)
  Blambda <- learner$learnParams$Blambda
  Cks <- learner$learnParams$Cks
  
  #x <- learner$hyperParams$trainData$x
  #y <- learner$hyperParams$trainData$y
  #E <- K%*%K
  #E <- kernelDataArray2(K, as.matrix(y))
  
  Chat <- apply(Cks, c(1,2), sum)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  reps <- nrow(Blambda)/n
  H <- lapply(1:reps, function(el) H)
  H <- do.call(cbind, H)
  H <- lapply(1:reps, function(el) H)
  H <- do.call(rbind, H)
  if(getHyperPar(learner, "centerLx")) Chat <- H%*%Chat%*%H
  
  
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
    D <- t(alpha)%*%Chat%*%alpha
    S <- diag(D)
    T <- sum(Chat*Alambda)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    S <- diag(t(phiy)%*%D%*%phiy)
    T <- sum(D*K)
  }
  
  
  
  
  # n <- 100; m <- 10; phiy <- matrix(rnorm(n*m),n,m)
  
  #proc.time()
  #S <- sapply(1:n, function(i) sapply(1:n, function(j) D[i,j]*phiy[i,]*phiy[j,], simplify="array"), simplify="array")
  #S <- apply(S, 1, sum)
  #proc.time() -pm # 15 secs
  
  proc.time()
 
  proc.time() - pm #15.9
  
  #plot(S,S2); abline(a=0, b=1, col="red")
  
  Stilde <- S/T
  res <- sum(abs(Stilde)*log(abs(Stilde)))
  
  return(res)
}

KCCCpval_ent <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCCCpval_ent"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCCC_ent(learner)
  
  
  KCCC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCCC_ent(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCCC_dist_null>mesr)/learner$msrs[["KCCCpval_ent"]]$pars$numPerms
  
  
  return(pval)
}

# complexity measure: complexity as in (https://arxiv.org/pdf/1009.1498.pdf) of components of sum of gradient accross all data points
KCCC_comp <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky 
  phiy <- learner$learnParams$phiy
  n  <- nrow(L)
  m <- ncol(phiy)
  Blambda <- learner$learnParams$Blambda
  Cks <- learner$learnParams$Cks
  
  #x <- learner$hyperParams$trainData$x
  #y <- learner$hyperParams$trainData$y
  #E <- K%*%K
  #E <- kernelDataArray2(K, as.matrix(y))
  
  Chat <- apply(Cks, c(1,2), sum)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  reps <- nrow(Blambda)/n
  H <- lapply(1:reps, function(el) H)
  H <- do.call(cbind, H)
  H <- lapply(1:reps, function(el) H)
  H <- do.call(rbind, H)
  if(getHyperPar(learner, "centerLx")) Chat <- H%*%Chat%*%H
  
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
    D <- t(alpha)%*%Chat%*%alpha
    S <- diag(D)
    T <- sum(Chat*Alambda)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    S <- diag(t(phiy)%*%D%*%phiy)
    T <- sum(D*K)
  }
  
 
  
  # n <- 100; m <- 10; phiy <- matrix(rnorm(n*m),n,m)
  
  #S <- sapply(1:n, function(i) sapply(1:n, function(j) D[i,j]*phiy[i,]*phiy[j,], simplify="array"), simplify="array")
  #S <- apply(S, 1, sum)
  
  Stilde <- abs(S)
  Stilde <- Stilde/T
  ent <- sum(Stilde*log(Stilde))
  # disequilibrium
  
  diseq <- sum((Stilde-1/m)^2)
  res <- ent*diseq
  
  return(res)
}

KCCCpval_comp <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCCCpval_comp"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCCC_comp(learner)
  
  
  KCCC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCCC_comp(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCCC_dist_null>mesr)/learner$msrs[["KCCCpval_comp"]]$pars$numPerms
  
  
  return(pval)
}



HSIC_cmem <- function(learner){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  
  Ks <- vector("list", 2)
  Ks[[1]] <- K
  Ks[[2]] <- L
  
  res <- dhsic(K=Ks)$dHSIC
  
  return(res)
}


# Kernel Conditional Norm Sensitivity for Causal Inference(KCNSC)

KCNSC <- function(learner){
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky 
  n  <- nrow(L)
  Cks <- learner$learnParams$Cks
  Blambda <- learner$learnParams$Blambda
  
  alpha <- learner$learnParams$alpha
  if(!is.null(alpha)){
    Alambda <- alpha%*%t(alpha)
  } else{
    
    Alambda <- Blambda%*%K%*%t(Blambda)
  }
  
  #x <- learner$hyperParams$trainData$x
  
  
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  reps <- nrow(Blambda)/n
  H <- lapply(1:reps, function(el) H)
  H <- do.call(cbind, H)
  H <- lapply(1:reps, function(el) H)
  H <- do.call(rbind, H)
  if(getHyperPar(learner, "centerLx")) Cks <- sapply(1:dim(Cks)[3], function(k) H%*%Cks[,,k]%*%H, simplify="array")
  
  # pm <- proc.time()  
  # Cks2 <- apply(Cks, 3, function(Ck) H%*%Ck%*%H)
  # proc.time() - pm # 0.128 for n =100
  # 
  # pm <- proc.time()
  # Cks3 <- sapply(1:dim(Cks)[3], function(k) H%*%Cks[,,k]%*%H, simplify="array")
  # proc.time() - pm # 0.109 for n=100
  # 
 
  
  # smpl <- sample(prod(dim(Cks2)), size=100)
  # plot(as.numeric(Cks2)[smpl], as.numeric(Cks3)[smpl]); abline(a=0, b=1, col="red")
  
  Fis <- sapply(1:dim(Cks)[3], function(k) L%*%Alambda%*%Cks[,,k]%*%Alambda%*%t(L), simplify="array")
  Fi_ii <- sapply(1:n, function(i) Fis[i,i,i])
  res <- sum(Fi_ii)
  
  
  return(res)
  
}

KCNSCpval <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCNSCpval"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCNSC(learner)
  
  
  KCNSC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner
    Ky <- learner$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCNSC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCNSC_dist_null>mesr)/learner$msrs[["KCNSCpval"]]$pars$numPerms
  
  
  return(pval)
}



