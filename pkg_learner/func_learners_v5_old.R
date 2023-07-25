# learner package

print("func_learners_v5 pkg")

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
library(corpcor) # is.positive.definite

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
  # addParams <- list(makeKernel=makeKernel, makeFeature=makePhi, calcMsrs=calcMsrs, msrs=measureParams)
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
                      "hsicYhReg", "nhsicRegA", "hsicRegA", "hsicYhRegA", "corre", "cmem_L2_f","cmem_L2_k","cmem_L2_f_rel","cmem_L2_k_rel",
                      "negLogLik", "hsicLoss2", "pinball", "gauss_log_lik","negCE","PCE","MisCR","KCDC","KCRDC","KCMC","KCSC","KCNSC","KCCC_ent","KCCC_pca_ent")
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

# CV helper
constructParams <- function (otherParams, paramsLists) {
  
  params <- lapply(paramsLists, function(paramsList){
    pn = names(paramsList)#names(substitute(paramsList))[-1]
    ret = expand.grid(paramsList, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    params = lapply(1:nrow(ret), function(ind){ 
      aux <- as.list(ret[ind,])
      names(aux) <- colnames(ret)
      return(c(aux, otherParams))
    })
    paramNames = lapply(1:nrow(ret), function(ind) paste(pn, ret[ind, ], sep = "=", collapse = " "))
    names(params) = paramNames
    return(params)
  })
  params <- do.call("c", params)
  class(params) = "CVST.params"
  return(params)
}


# set heuristic based and optimization (CV, max-likelihood etc) parameters
setParams <- function(learner, trainData, plot=FALSE){
  #learner should be an object of class "emley.learner"
  stopifnot(class(learner)=="emley.learner")
  #stopifnot(class(trainData)=="CVST.data")
  
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
  learner$hyperParams$seed <- as.numeric(substr(abs(round(sum(trainData$x),4))*10^4,1,4))
  
  
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
  
  learner$hyperParams$data$grid <- dataOptimParams$grid
  learner$hyperParams$data$indxOptGrid <- dataOptimParams$indxOptGrid
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
      plot(predList[[1]]$x, predList[[1]]$gy, xlab="x",ylab="y", 
           ylim=range(predList[[1]]$gy,sapply(predList, function(el) range(el$gyh))))
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
  
  kernelChar <- hyperParams[[kernelName]]$val
  kernel_char <- sapply(strsplit(kernelChar, split="T"), function(el) el[1])
  kernel_char <- sapply(strsplit(kernel_char, "_"), function(el) el[2])
  xyChar <- substr(kernelName, nchar(kernelName), nchar(kernelName))
  aux <- strsplit(hyperParamsNames, "\\.")
  kernelParams <- lapply(kernel_char, function(krnl_char){
    # krnl_char <- kernel_char[1]
    indx <- which(sapply(aux, function(el) el[2]) == krnl_char & sapply(aux, function(el) el[3])==xyChar)
    kernelParamsFuncArgs <- sapply(aux[indx], function(el) el[1])
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
  })
  names(kernelParams) <- kernelChar
  return(list(kernelChar=kernelChar, kernelParams=kernelParams))
}  


getHyperPar <- function(learner, parName){
  hyperParams <- learner$getHyperParams(learner)
  hyperParamsNames <- names(hyperParams)
  indx <- which(hyperParamsNames==parName)
  if(length(indx)==0) return(NULL)
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

getNumBins <- function(learner, data){
  
  
  numX <- learner$hyperParams$data$optimizable$num.bin.X$length.out
  if(is.null(numX)) numX <- 5
  numY <- learner$hyperParams$data$optimizable$num.bin.Y$length.out
  if(is.null(numY)) numY <- 5
  
  num.binX.val <- 3
  num.binY.val <- 3
  num.binX.seq <- seq(3,length.out=numX)
  num.binY.seq <- seq(3,length.out=numY)
  
  return(list(val=list(num.bin.X=num.binX.val, num.bin.Y=num.binY.val), 
              seq=list(num.bin.X=num.binX.seq, num.bin.Y=num.binY.seq)))
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
  
  numX <- learner$hyperParams$data$optimizable$offset.quad.X$length.out
  if(is.null(numX)) numX <-10
  numY <- learner$hyperParams$data$optimizable$offset.quad.Y$length.out
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
  # x <- matrix(rnorm(100), 100, 1)
  sigma0 <- 1/quantile(as.numeric(dist(x)^2), 0.5)
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
  # estimate interpolation model for variance_of_kernel(sigma)
  #spl <- spline(log(sigmas[indx],10), vars[indx]) #just for plotting
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  # out of prior sequence of sigmas just remove NAs
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  # interpolate at those non-NA poing
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  # calculate index of sigma for which var is max
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  # our starting point is the point where the variance is closest to being
  # minPct (1%) of the max variance on the LEFT side of the max
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  # our finishing point is the point where the variance is the closest to being
  # minPct (1%) of the max variance on the RIGHT side of the max
  indxFinish <- min(100,indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct)))
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #abline(v=log(c(sigma0, sigma1),10), col="green")
  #sigmas.x[c(indxStart, indxMax.x, indxFinish)]
  #log(sigmas.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  #vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2 - we do again to get a closer estimation of the > 1% variance interval
  sigmas <- 10^seq(log(sigmas.x[indxStart],10), log(sigmas.x[indxFinish],10), length.out=10)
  minPctLeft <- 0.01
  minPctRight <- 0.25
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
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPctLeft))
  indxFinish <- min(100, indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPctRight)))
  
  sigmaMed <- 1/quantile(as.numeric(dist(x)^2), 0.5)
  
  
  ini <- log(sigmas.x, 10)[indxStart]
  fin <- log(sigmas.x, 10)[indxMax.x]
  #fin <- log(sigmas.x, 10)[indxFinish]
  #fin <- sigmaMed
  # based on script_analysis_lambdas_sigmas.R in pkg_learner/experiment/cmemLearners/
  # I have concluded that the max we need to go is sigma0=1/median(dist(xx)^2)
  seq <- 10^seq(ini, fin, length.out=length.out)
  
  #plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  #abline(v=log(seq,10), col="red")
  #abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  #abline(v=log(c(sigma0, sigma1),10), col="green")
  
  
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
  # x <- matrix(rnorm(100), 100, 1)
  
  offsets <- 10^seq(-16,16,by=1)
  minPct.left <- 0.01
  minPct.right <- 0.95
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
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct.left))
  indxFinish <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct.right))
  
  #plot(log(offsets,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(offsets.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(offsets.x,10)[c(indxStart, indxMax.x, indxFinish)], col=c("red","orange","red"))
  #offsets.x[c(indxStart, indxMax.x, indxFinish)]
  #log(offsets.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  #vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  offsets <- 10^seq(log(offsets.x[indxStart],10), log(offsets.x[indxFinish],10), length.out=10)
  minPct.left <- 0.01
  minPct.right <- 0.5
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
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct.left))
  indxFinish <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct.right))
  
  #plot(log(offsets,10), vars, ylab="var"); abline(h=0, col="red")
  #lines(spl, col="green")
  #lines(log(offsets.x, 10), vars.x, col="blue", type="l")
  #abline(h=maxVar*minPct, v=log(offsets.x,10)[c(indxStart, indxFinish)], col=c("red","orange"))
  
  ini <- log(offsets.x, 10)[indxStart]
  #fin <- log(offsets.x, 10)[indxFinish]
  fin <- log(offsets.x, 10)[indxMax.x]
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
  
  nms <- names(optimHyperParams)
  aux <- strsplit(nms, "\\.")
  indX <- which(nms=="kernelX")
  indY <- which(nms=="kernelY" | nms=="featureY")
  indXX <- which(sapply(aux, function(el) el[3])=="X")
  indYY <- which(sapply(aux, function(el) el[3])=="Y")
  indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
  
  if(length(indX)>0 & length(indY)>0){
    paramsListss <- lapply(optimHyperParams[[indX]]$seq, function(krnX){
      # krnX <- optimHyperParams[[indX]]$seq[1]
      indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
      
      paramsListss <- lapply(optimHyperParams[[indY]]$seq, function(krnY){
        # krnY <- optimHyperParams[[indY]]$seq[1]
        indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
        optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
    
        paramsListO <- lapply(optimHyperParamsO, function(el) el$seq)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListX <- lapply(optimHyperParamsX, function(el) el$seq)
        names(paramsListX) <- names(optimHyperParamsX)
        paramsListY <- lapply(optimHyperParamsY, function(el) el$seq)
        names(paramsListY) <- names(optimHyperParamsY)
        paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
        names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
        return(paramsList)
      })
      
      return(paramsListss)
      })
    paramsListss <- unlist(paramsListss, recursive=F)
  
  } else if(length(indX)>0){
    paramsListss <- lapply(optimHyperParams[[indX]]$seq, function(krnX){
      # krnX <- optimHyperParams[[indX]]$seq[1]
      indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
      
      optimHyperParamsO <- optimHyperParams[indO]
      optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
      paramsListO <- lapply(optimHyperParamsO, function(el) el$seq)
      names(paramsListO) <- names(optimHyperParamsO)
      paramsListX <- lapply(optimHyperParamsX, function(el) el$seq)
      names(paramsListX) <- names(optimHyperParamsX)
      paramsList <- c(krnX,  paramsListO, paramsListX)
      names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
      return(paramsList)
      })
      
  } else if(length(indY)>0){
    paramsListss <- lapply(optimHyperParams[[indY]]$seq, function(krnY){
        # krnY <- optimHyperParams[[indY]]$seq[1]
        indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
        
        paramsListO <- lapply(optimHyperParamsO, function(el) el$seq)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListY <- lapply(optimHyperParamsY, function(el) el$seq)
        names(paramsListY) <- names(optimHyperParamsY)
        paramsList <- c(krnY, paramsListO, paramsListY)
        names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
        return(paramsList)
      })
      
  } else{
      optimHyperParamsO <- optimHyperParams[indO]
      paramsListO <- lapply(optimHyperParamsO, function(el) el$seq)
      names(paramsListO) <- names(optimHyperParamsO)
      paramsList <- c(paramsListO)
  }
  
   otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
  
  
  params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
  
  # not necessary coz CV is only done on data hyperparms??
  #params <- lapply(params, function(el) c(el, learner$hyperParams$non_data))
  #class(params) = "CVST.params"
  
  
  grid <- CV.parallel(learner, params, fac=fac)
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  
  
  dimnms <- dimnames(grid)
  
  
  
  # obtain best hyperparameters
  
  mainLoss <- learner$optimizeParams$mainLoss
  testTrain <- learner$optimizeParams$testTrain
  
  testGrid <- grid[testTrain,,mainLoss]
  
  minTestGrid <- min(testGrid, na.rm=T)
  #we add a tiny bit of noise to get exactly one minimum
  testGrid <- testGrid+ rnorm(length(testGrid),mean=0, sd=max(minTestGrid,1e-10)*1e-10)
  minTestGrid <- min(testGrid, na.rm=T)
  optMinBool <- testGrid==minTestGrid
  
  # keep opt indices on grid to have easy access to chosen grid point later on
  
  indxOptGrid <- which.min(testGrid)
  
  
  
  
  opt <- params[[indxOptGrid]][intersect(nms, names(params[[indxOptGrid]]))]
  
  
  
  # check (only when lanbda and sigma are passed for CV)
  #grid["test",which.min(abs(as.numeric(dimnames(grid)$lambda)-opt$lambda)),which.min(abs(as.numeric(dimnames(grid)$sigma)-opt$sigma)),1]
  #min(grid["test",,,1])
  
  res <-  list(opt=opt, grid=grid, indxOptGrid=indxOptGrid)
  
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
      # learnerAux$learnParams$beta_f
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


pred.CV <- function(learner, data, numCores=1){
  
  numFolds <- learner$optimizeParams$numFolds
  n <- getN(data)
  size <- floor(n / numFolds)
  
  res <- mcmapply(FUN=function(f){
    # f <- 1
    #print(paste("fold", f))
    last <- f*size
    if(f==numFolds) last <- n
    validationIndex <- seq((f-1)*size + 1, last)
    curTrain <- getSubset(data, -validationIndex)
    curTest <- getSubset(data, validationIndex)
    # either mean squared error or mean classification error
    
    
    learner1 <- learner
    learner1$hyperParams$trainData <- curTrain
    
    #print("a")
    learner1 <- learner1$learn(learner=learner1, forLoss=T)
    #print("b")
    
    predTrain <- learner1$predict(learner=learner1, data=curTrain, forLoss=T)
    predTest <- learner1$predict(learner=learner1, data=curTest, forLoss=T)
    
    
    res <- list(train=predTrain, test=predTest)
    
    
    return(res)
  }, f=1:numFolds, mc.cores=numCores, SIMPLIFY=FALSE)
  
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


CV.parallel <- function(learner, params, fac=1, verbose=TRUE, numCoresFold=1) {
  stopifnot(class(learner) == "emley.learner" && class(params) == "CVST.params")
  
  #print("enters CV parallel")
  
  numFolds <- learner$optimizeParams$numFolds
  trainData <- learner$hyperParams$trainData
  
  
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), names(params))
  
  mc_cores <- min(40, detectCores()-4)
   # mc_cores <- 1
   #numCoresFold <- 1
  
  #count <<-0 
  pm <- proc.time()
  losses <- mcmapply(FUN=function(p){
    # p <- params[[16]]
    # print(p)
     #count <<- count + 1
     #print(count)
    learnerAux <- learner
    nmsPars <- names(learnerAux$hyperParams$data$optimizable)
    for(pr in nmsPars){
      # pr <- nmsPars[1]
      learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
    }
    
      # I dont actually know why im calling learn before goint into 
      # folds, unless I need something to calculate losses below?
      # gonna comment it out for now... no!! i know, its so i can
      # have KCDC, KCMC etc as loss functions since they use
      # learned stuff on all data
      learnerAux <- learnerAux$learn(learner=learnerAux, forLoss=F)
      preds <- pred.CV(learner=learnerAux, data=trainData, numCores = numCoresFold)
    
    # x <- trainData$x; y <- trainData$y; ord <- order(x)
    # plot(x, y)
    # for(i in 1:ncol(preds$test$gyh)) lines(x[ord], preds$test$gyh[ord,i], col=i)
    
    lossesTrain <- lapply(names(learnerAux$optimizeParams$losses), 
                          function(func){
                            # func <- "NCE"
                            # func <- 
                            #print(func)
                            do.call(func, list(learner=learnerAux, pred=preds$train))
                          })
    names(lossesTrain) <- names(learnerAux$optimizeParams$losses)
    lossesTrain <- unlist(lossesTrain)
    
    lossesTest <- lapply(names(learnerAux$optimizeParams$losses), function(func){
      # func <- "cmem_L2_f"
      do.call(func, list(learner=learnerAux, pred=preds$test))
    })
    names(lossesTest) <- names(learnerAux$optimizeParams$losses)
    lossesTest <- unlist(lossesTest)
    
    res <- cbind(lossesTrain, lossesTest)
    #print(res)
    
    return(res)
  }, params, mc.cores=mc_cores, SIMPLIFY="array")
  proc.time() - pm #
  
  losses <- aperm(losses, c(2,3,1))
  dimnames(losses) <- list(trainTest=c("train","test"), params=names(params), var=dimnames(losses)[[3]])
  
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


learn.vanilla <- function(learner, forLoss=F) {
  return(learner)
}
predict.vanilla <- function(learner, data, forLoss=F) {
  
  gy <- as.matrix(data$y)
  x <- data$x
  return(list(x=as.matrix(x), gy=gy, gyh=matrix(0, nrow(gy), ncol(gy))))
}


learn.krr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")

  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  N <- nrow(Kxs)
  
  lambda <- learner$hyperParams$data$optimizable$lambda$val*N
  
  alpha <- solve(Matrix(Kxs + diag(lambda, N))) %*% y
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.krr <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}

learn.qhsic <- function (learner, forLoss=F) {
  
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
predict.qhsic <- function (learner, data, forLoss=F){
  
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


learn.kqr <- function(learner, forLoss=F) {
  
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

predict.kqr <- function(learner, data, forLoss=F) {
  
  
  
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


learn.hsic <- function (learner, forLoss=F) {
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
predict.hsic <- function (learner, data, forLoss=F) {
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=learner$hyperParams$trainData$x, pars=parsXs)
  

  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}



learn.gptk <- function(learner) return(learner, forLoss=F)
predict.gptk <- function(learner, data, forLoss=F){
  aux <- gpPosteriorMeanVar(learner$hyperParams$data$optimizable$model$val, X=data$x)
  pred <- list(x=data$x, gy=data$y, gyh=aux)
  return(pred)
}

learn.logReg <- function(learner, forLoss=F) {
  
  phix <- learner$makeFeature(learner, data=learner$hyperParams$trainData, "X")
  
  y <- learner$hyperParams$trainData$y
  
  dataMat <- as.data.frame(cbind(y, phix))
  colnames(dataMat) <- c("y", paste("phix", seq(ncol(phix)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  fit <- glm(formula=form, data=dataMat, family = binomial(link="logit"))
  
  
  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logReg <- function(learner, data, forLoss=F) {
  
  phix <-learner$makeFeature(learner, data, "X") 
  
  y <- learner$hyperParams$trainData$y
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  
  pred <- predict(fit, newdata=dataMat, type="response")
  
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}


learn.logKReg <- function(learner, forLoss=F) {
  
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
predict.logKReg <- function(learner, data, forLoss=F) {
  
  phix <-learner$makeFeature(learner, data, "X") 

  y <- learner$hyperParams$trainData$y
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  #pred <- predict(fit, newdata=dataMat, type="response")
  pred <- predict(fit, newx=phix, s=lambda, type="response")
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}

learn.logRegInt <- function(learner, forLoss=F) {
  
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
predict.logRegInt <- function(learner, data, forLoss=F) {
  
  phix <-data$x 
  
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  #pred <- predict(fit, newdata=dataMat, type="response")
  pred <- predict(fit, newx=phix, s=lambda, type="response")
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}

learn.logRegInt0 <- function(learner, forLoss=F) {
  
  phix <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  dataMat <- as.data.frame(cbind(y, phix))
  colnames(dataMat) <- c("y", paste("phix", seq(ncol(phix)), sep="_"))
  form <- as.formula(paste("y~", paste(colnames(dataMat)[2:ncol(dataMat)], collapse="+")))
  
  #fit <- glmnet(x=phix, y=y, family = "binomial", alpha = 0, lambda = 0)
  fit <-  glm(form, data=dataMat, family=binomial(link="logit"))  
  
  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logRegInt0 <- function(learner, data, forLoss=F) {
  
  phix <-data$x 
  
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model
  
  pred <- predict(fit, newdata=dataMat, type="response")
  #pred <- predict(fit, newx=phix, s=0, type="response")
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}

learn.logRegInt3 <- function(learner, forLoss=F) {
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  lambdas <- getHyperPar(learner, "lambda")*N
  
  numFolds <- learner$optimizeParams$numFolds
  sigma0 <-  1/median(as.numeric(dist(x)^2))
  sigma0 <- min(sigma0, 1e10)
  sigma0 <- max(sigma0, 1e-10)
  
  learner$learnParams$sigma <- sigma0
  K <- kern_rbf(x=learner$hyperParams$trainData$x, sigma=sigma0)
  
  
  
  cv.fit <- cv.glmnet(x=K, y=y, family="binomial", alpha=0, lambda=lambdas, nfolds=numFolds)
  
  #plot(log(cv.fit$lambda,10), cv.fit$cvm)
  #plot(cv.fit)
  
  
  indx <- which.min(cv.fit$cvm)
  #fit <- cv.fit$glmnet.fit
  lambda <- cv.fit$lambda[indx]
  
  learner$hyperParams$data$non_optimizable$lambda$val <- lambda
  
  fit <- glmnet(x=K, y=y, family = "binomial", alpha = 0, lambda = lambda)
  #fit <-  glm(form, data=dataMat, family=binomial(link="logit"))  
  
  learner$learnParams$model <- fit
  
  return(learner)
}
predict.logRegInt3 <- function(learner, data, forLoss=F) {
  
  x <-data$x 
  
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  lambda <- getHyperPar(learner, "lambda")*N
  
  kxs <- kern_rbf(x=data$x, y=learner$hyperParams$trainData$x, sigma=learner$learnParams$sigma)
  
  
  fit <- learner$learnParams$model
  #pred <- predict(fit, newdata=dataMat, type="response")
  pred <- predict(fit, newx=kxs, s=lambda, type="response")
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}


learn.cmeClass <- function(learner, forLoss=F) return(learner)
predict.cmeClass <- function(learner, data, forLoss=F) return(data)

learn.logNNReg <- function(learner, forLoss=F){
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
predict.logNNReg <- function(learner, data, forLoss=F){
  phix <-learner$makeFeature(learner, data, "X") 
  
  dataMat <- as.data.frame(phix)
  colnames(dataMat) <- paste("phix", seq(ncol(phix)), sep="_")
  
  fit <- learner$learnParams$model_log
  pred <- predict(fit, newdata=dataMat, type="response")
  
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}

# generic transformation and residual functions

NCE <- function(x, classifier, keepHyperParams=NULL){
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
  featNm <- paste("feature", var, sep="")
  pars <- getKernelPars(learner, kernelName=featNm)
  phi_char <- pars$kernelChar
  phi_char <- sapply(strsplit(phi_char, split="T"), function(el) el[1])
  pars <- pars$kernelParams
  
  phi <- mcmapply(function(phi_char, pars){
    x <- data[[var2]]
    auxPars <- pars
    auxPars$x <- x
    aux <- strsplit(phi_char, split="_")[[1]]
    if(aux[1] == "rff") phi_char <- "rff"
    # names(auxPars)
    # rff(x=auxPars$x, num_f=auxPars$num_f, seed=auxPars$seed, p_w=auxPars$p_w, map=auxPars$map, sigma=auxPars$sigma)
    phi  <- do.call(phi_char, auxPars)
    return(phi)
  }, phi_char=phi_char, pars=pars, mc.cores=1, SIMPLIFY=F)
  
  phi <- do.call(cbind, phi)
  
  return(phi)
}


makeKernel <- function(learner, data1, data2=data1, var=c("X","Y"), grad=FALSE){
  var2 <- c("x","y")[match(var, c("X","Y"))]
  krnNm <- paste("kernel", var, sep="")
  pars <- getKernelPars(learner, kernelName=krnNm)
  kernel_char <- pars$kernelChar
  kernel_char <- sapply(strsplit(kernel_char, split="T"), function(el) el[1])
  pars <- pars$kernelParams
  
  Ks <- mcmapply(function(kernel_char, pars){
    x <- as.matrix(data1[[var2]])
    y <- as.matrix(data2[[var2]])
    res  <- list(K=kernelMatrix(kernel=kernel_char, x, y, pars=pars))
    if(grad){
      Cks <- kernelGradNormMatrix(kernel=kernel_char, x, y,  pars=pars)
      res <- c(res, list(Cks=Cks))
    }
    return(res)
  }, kernel_char=kernel_char, pars=pars, mc.cores=1, SIMPLIFY=F)
  
  if(length(Ks)>1){
    K <- sapply(Ks, function(el) el$K, simplify="array")
    K <- apply(K, c(1,2), sum)
    res <- list(K=K)
  } else{
    K <- Ks[[1]]$K
    res <- list(K=K)
  }
  if(grad){
    if(length(Ks)>1){
      Cks <- sapply(Ks, function(el) el$Cks, simplify="array")
      Cks <- apply(Cks, c(1,2,3), sum)
      res <- c(res, Cks=list(Cks=Cks))
    } else{
      Cks <- Ks[[1]]$Cks
      res <- c(res, Cks=list(Cks=Cks))
    }
  }
  return(res)
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
  nms <- names(pred)
  aux1 <- strsplit(nms, "\\.")
  nms1 <- sapply(aux1, function(el) paste(el[2], collapse="_"))
  indx_gyh <- which(nms1 == "gyh_class")
  indx_gy <- which(nms1 == "gy_class")
  indx <- sort(c(indx_gy, indx_gyh))
  aux1 <- aux1[indx]
  nms1 <- nms1[indx]
  
  nms2 <- sapply(aux1, function(el) el[1])
  
  res <- sapply(unique(nms2), function(featClass){
    # featClass <- unique(nms2)[1]
    indx <- which(nms2==featClass)
    nms <- paste(nms2[indx], c("gyh_class","gy_class"), sep=".")
    res <- LogLoss(pred[[nms[1]]], pred[[nms[2]]])  
  })
  
  
  return(res)
}

MisCR <- function(learner, pred, thrsh=0.5){
  nms <- names(pred)
  aux1 <- strsplit(nms, "\\.")
  nms1 <- sapply(aux1, function(el) paste(el[2], collapse="_"))
  indx_gyh <- which(nms1 == "gyh_class")
  indx_gy <- which(nms1 == "gy_class")
  indx <- sort(c(indx_gy, indx_gyh))
  aux1 <- aux1[indx]
  nms1 <- nms1[indx]
  
  nms2 <- sapply(aux1, function(el) el[1])
  
  res <- sapply(unique(nms2), function(featClass){
    # featClass <- unique(nms2)[1]
    indx <- which(nms2==featClass)
    nms <- paste(nms2[indx], c("gyh_class","gy_class"), sep=".")
    conTab <- contingencyTable(pred=(pred[[nms[1]]]>thrsh)*1, obs=pred[[nms[2]]])
    res <- misCR(conTab)
    names(res) <- NULL
    return(res)
  })
  
  
  return(res)
}

CCR <- function(learner, pred, thrsh=0.5){
  nms <- names(pred)
  aux1 <- strsplit(nms, "\\.")
  nms1 <- sapply(aux1, function(el) paste(el[2], collapse="_"))
  indx_gyh <- which(nms1 == "gyh_class")
  indx_gy <- which(nms1 == "gy_class")
  indx <- sort(c(indx_gy, indx_gyh))
  aux1 <- aux1[indx]
  nms1 <- nms1[indx]
  
  nms2 <- sapply(aux1, function(el) el[1])
  
  res <- sapply(unique(nms2), function(featClass){
    # featClass <- unique(nms2)[1]
    indx <- which(nms2==featClass)
    nms <- paste(nms2[indx], c("gyh_class","gy_class"), sep=".")
    conTab <- contingencyTable(pred=(pred[[nms[1]]]>thrsh)*1, obs=pred[[nms[2]]])
    res <- correctCR(conTab)
    names(res) <- NULL
    return(res)
  })
  
  
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
learnBlambda_L2 <- function(Lx, n, lambda){
  I <- diag(n)
  Blambda <- solve(Lx+n*lambda*I)
  return(Blambda)
}
learnBlambda_bin <- function(Lx, n, lambda){
  Blambda <- diag(1/diag(Lx%*%Lx))
  return(Blambda)
}
learnBlambda_KCMC <- function(Lx, n, lambda, gamma){
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  #Lx%*%(lambda*H)
  Blambda <- solve(Lx%*%(I+gamma*H)+n*lambda*I)
  return(Blambda)
}
learnBlambda_KCSC <- function(Lx, Cks, centerChat, n, lambda, gamma){
  I <- diag(n)
  Chat <- apply(Cks, c(1,2), sum)
  if(centerChat) Chat <- H%*%Chat%*%H
  Blambda <- solve(Lx%*%Lx+(gamma/n)*Chat+n*lambda*I)%*%Lx
  return(Blambda)
}


makeLogRegFeats_dep <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  
  # forLoss=F so that we don't try and assess a classifier we haven't built yet
  pred_real <- learner$predict(learner, data=data, forLoss=F)
  feats_class_real <- pred_real$gyh_k
  
  dataAux <- data
  
  dataAux$x <- fakes_x
  dataAux$y <- fakes_y
  pred_fake <- learner$predict(learner, data=dataAux, forLoss=F)
  feats_class_fake <- pred_fake$gy_k
  
  # make train classification data
  x_tr <- as.matrix(rbind(feats_class_real, feats_class_fake))
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  return(classData)  
}

makeLogRegFeats_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  # forLoss=F so that we don't try and assess a classifier we haven't built yet
  pred_real <- learner$predict(learner, data=data, forLoss=F)
  feats_class_real <- (pred_real$gy_k-pred_real$gyh_k)^2
  
  dataAux <- data
  
  dataAux$x <- fakes_x
  dataAux$y <- fakes_y
  pred_fake <- learner$predict(learner, data=dataAux, forLoss=F)
  pred_real_gy_k <- lapply(1:kappa, function(i) pred_real$gy_k)
  pred_real_gy_k <- do.call(rbind, pred_real_gy_k)
  feats_class_fake <- (pred_real_gy_k-pred_fake$gyh_k)^2
  
  # make train classification data
  x_tr <- as.matrix(rbind(feats_class_real, feats_class_fake))
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  return(classData)  
}
makeCME_cond_logRegfeats_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  
  probs <- diag(ky%*%Blambda%*%t(lx))
  #sum(probs <0 ); sum(probs>1)
  
  
  
  
  # make train classification data
  x_tr <- as.matrix(probs)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsNorm_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Alambda <- Blambda %*% Ky %*% Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  # krnNm <- "kernelX"
  # pars <- getKernelPars(learner, kernelName=krnNm)
  # kernel_char <- pars$kernelChar
  # kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  # pars <- pars$kernelParams
  # intlx <- mcmapply(function(ker, par){
  #   # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #   par$lim_max <- 0.5
  #   par$lim_min <- -0.5
  #   ker <- paste(ker, "nrml", sep="_")
  #   intlx  <- kernelMatrix(ker, trainData$x, dataWithFakes$x, pars=par)
  #   return(intlx)
  # }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  # if(dim(intlx)[3]==1) intlx <- intlx[,,1] else intlx <- apply(intlx, c(1,2), sum)
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
    

  probs <- diag(ky%*%Blambda%*%t(lx))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(ky%*%Blambda%*%intlx)
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  #sum(probsNorm < 0); sum(probs>1)
  
  
  if(FALSE){
    unNormProb <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
      # y <- 0;  x <- 0.1
      n <- length(x_tr)
      ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
      lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
      res <- ky%*%Blambda%*%lx
      return(as.numeric(res))
    }
    normProb2 <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
      # y <- 0;  x <- 0.1
      n <- length(x_tr)
      ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
      intlx <- kern_rbf_nrml(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax, lim_max=0.5, lim_min=-0.5)
      lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
      res <- (ky%*%Blambda%*%lx)/(ky%*%Blambda%*%intlx)
      return(as.numeric(res))
    }
    
    unNormProb <- function(xx, y, Blambda, y_tr, x_tr, offsety, offsetx){
      # y <- 0;  x <- 0.1
      n <- length(x_tr)
      ky <- kern_quad(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), offset=offsety)
      lx <- kern_quad(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), offset=offsetx)
      res <- ky%*%Blambda%*%lx
      return(as.numeric(res))
    }
    count <<- 0
    normProb2 <- function(xx, y, Blambda, y_tr, x_tr, offsety, offsetx){
      # y <- 0;  x <- 0.1
      count <<- count + 1
      print(count)
      n <- length(x_tr)
      ky <- kern_quad(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), offset=offsety)
      intlx <- kern_quad_nrml(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), offset=offsetx, lim_max=0.5, lim_min=-0.5)
      lx <- kern_quad(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), offset=offsetx)
      #print(dim(ky)); print(dim(lx)); print(dim(intlx))
      res <- (ky%*%Blambda%*%lx)/(ky%*%Blambda%*%intlx)
      return(as.numeric(res))
    }
    
    xx <- seq(-0.5,0.5,length.out=100)
    y0 <- dataWithFakes$y[1]
    sigmay <- learner$hyperParams$data$optimizable$sigma.rbf.Y$val
    sigmax <- learner$hyperParams$data$optimizable$sigma.rbf.X$val
    offsety <- learner$hyperParams$data$non_optimizable$offset.quad.Y$val
    offsetx <- learner$hyperParams$data$non_optimizable$offset.quad.X$val
    
    yy_unorm <- unNormProb(xx, y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, 
                           x_tr=learner$hyperParams$trainData$x, sigmay=sigmay, sigmax=sigmax)
    yy_norm <- normProb2(xx, y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, 
                         x_tr=learner$hyperParams$trainData$x, sigmay=sigmay, sigmax=sigmax)
    
    yy_unorm <- unNormProb(xx, y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, 
                           x_tr=learner$hyperParams$trainData$x, offsety=offsety, offsetx=offsetx)
    yy_norm <- normProb2(xx, y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, 
                         x_tr=learner$hyperParams$trainData$x, offsety=offsety, offsetx=offsetx)
    
    par(mfrow=c(1,2))
    plot(learner$hyperParams$trainData$x, learner$hyperParams$trainData$y); abline(h=y0, col="red") 
    plot(xx, yy_unorm, type="l", ylim=range(yy_unorm, yy_norm))
    lines(xx, yy_norm, type="l", col="red"); abline(h=0, col="blue")
    par(mfrow=c(1,1))
    
    # can I integrate normProb(x) to see if it integrates to one??
    
    sapply(c("Kronrod","Clenshaw","Simpson")[c(1,2,3)], function(meth) integral(normProb2, xmin=-0.5, xmax=0.5, method = meth, 
                                                                                 y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, x_tr=learner$hyperParams$trainData$x, sigmay=sigmay, sigmax=sigmax))
    
    sapply(c("Kronrod","Clenshaw","Simpson")[c(1,2,3)], function(meth) integral(normProb2, xmin=-0.5, xmax=0.5, method = meth, 
                                                                                y=y0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$y, x_tr=learner$hyperParams$trainData$x, offsety=offsety, offsetx=offsetx))
    
  }
  
  
  # make train classification data
  x_tr <- as.matrix(probsNorm)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsOdds_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Alambda <- Blambda %*% Ky %*% Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  # krnNm <- "kernelX"
  # pars <- getKernelPars(learner, kernelName=krnNm)
  # kernel_char <- pars$kernelChar
  # kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  # pars <- pars$kernelParams
  # intlx <- mcmapply(function(ker, par){
  #   # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #   par$lim_max <- 0.5
  #   par$lim_min <- -0.5
  #   ker <- paste(ker, "nrml", sep="_")
  #   intlx  <- kernelMatrix(ker, trainData$x, dataWithFakes$x, pars=par)
  #   return(intlx)
  # }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  # if(dim(intlx)[3]==1) intlx <- intlx[,,1] else intlx <- apply(intlx, c(1,2), sum)
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
  
  probs <- diag(ky%*%Blambda%*%t(lx))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(ky%*%Blambda%*%intlx)
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  probsNorm[which(probsNorm<=0)] <- 0.00001
  probsNorm[which(probsNorm>=1)] <- 0.99999
  #sum(probsNorm < 0); sum(probs>1)
  odds <- log(probsNorm/(1-probsNorm))
  
  
  # make train classification data
  x_tr <- as.matrix(odds)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsK_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  if(!is.positive.definite(Ky)){
    Ky <- make.positive.definite(Ky)
  }
  Alambda <- Blambda %*% Ky %*% Blambda
  if(!is.positive.definite(Alambda)){
    Alambda <- make.positive.definite(Alambda)
  }
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  # krnNm <- "kernelX"
  # pars <- getKernelPars(learner, kernelName=krnNm)
  # kernel_char <- pars$kernelChar
  # kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  # pars <- pars$kernelParams
  # intlx <- mcmapply(function(ker, par){
  #   # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #   par$lim_max <- 0.5
  #   par$lim_min <- -0.5
  #   ker <- paste(ker, "nrml", sep="_")
  #   intlx  <- kernelMatrix(ker, trainData$x, dataWithFakes$x, pars=par)
  #   return(intlx)
  # }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  # if(dim(intlx)[3]==1) intlx <- intlx[,,1] else intlx <- apply(intlx, c(1,2), sum)
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
  
  probs <- diag(ky%*%Blambda%*%t(lx))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(ky%*%Blambda%*%intlx)
  #print(summary(diag(lx%*%Alambda%*%t(lx))))
  
  #print(is.positive.definite(Blambda))
  #print(is.positive.definite(Ky))
  #print(is.positive.definite(Alambda))
  
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  
  
  # make train classification data
  #x_tr <- cbind(probs, const)
  x_tr <- as.matrix(probsNorm)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_feats_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  krnNm <- "kernelX"
  pars <- getKernelPars(learner, kernelName=krnNm)
  kernel_char <- pars$kernelChar
  kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  pars <- pars$kernelParams
  intlx <- mcmapply(function(ker, par){
    # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
    par$lim_max <- 0.5
    par$lim_min <- -0.5
    ker <- paste(ker, "nrml", sep="_")
    intlx  <- kernelMatrix(ker, trainData$x, dataWithFakes$x, pars=par)
    return(intlx)
  }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  if(dim(intlx)[3]==1) intlx <- intlx[,,1] else intlx <- apply(intlx, c(1,2), sum)
  
  
  probs <- diag(ky%*%Blambda%*%t(lx))
  #sum(probs <0 )
  const <- diag(ky%*%Blambda%*%intlx)
  #sum(const<0)
  probs <- probs/const
  #sum(probs < 0)  
  
  p_fake <- 1 # this assumes we have sampled fakes_x uniformly 
  # on [-0.5, 0.5] x ... x [-0.5, 0.5]
  class_probs <- probs/(probs + kappa*p_fake)
  #summary(class_probs)
  #sum(class_probs < 0)
  class_probs[which(class_probs < 0 | class_probs >1)] <- 0.5
  
  y_labs <- c(rep(1, n), rep(0, num_fake)) 
  
  # plot(class_probs, y_labs)
  
  res <- list(x_class=x, gy_class=y_labs, gyh_class=class_probs)
  return(res)
}


makeLogRegFeats_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  # forLoss=F so that we don't try and assess a classifier we haven't built yet
  pred_real <- learner$predict(learner, data=data, forLoss=F)
  feats_class_real <- (pred_real$gy_k-pred_real$gyh_k)^2
  
  dataAux <- data
  
  dataAux$x <- fakes_x
  dataAux$y <- fakes_y
  pred_fake <- learner$predict(learner, data=dataAux, forLoss=F)
  pred_real_gyh_k <- lapply(1:kappa, function(i) pred_real$gyh_k)
  pred_real_gyh_k <- do.call(rbind, pred_real_gyh_k)
  feats_class_fake <- (pred_fake$gy_k-pred_real_gyh_k)^2
  
  # make train classification data
  x_tr <- as.matrix(rbind(feats_class_real, feats_class_fake))
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  return(classData)  
}
makeCME_cond_logRegfeats_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  y_fake <- fakes_y
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_real_ext)
  y <- c(y_real, y_fake)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  
  
  probs <- diag(lx%*%Blambda%*%t(ky))
  #sum(probs <0 ); sum(probs>1)
  
  
  
  # make train classification data
  x_tr <- as.matrix(probs)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsNorm_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  y_fake <- fakes_y
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_real_ext)
  y <- c(y_real, y_fake)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Alambda <- Blambda %*% Ky %*% Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  #hyperParamsNames <- names(learner$getHyperParams(learner))
  # krnNm <- c("kernelY","featureY")[("featureY" %in% hyperParamsNames)*1+1]
  # if(krnNm == "kernelY"){
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   
  #   intky <- mcmapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     maxNumPars <- max(sapply(par, length))
  #     par$lim_max <- rep(0.5, maxNumPars)
  #     par$lim_min <- rep(-0.5, maxNumPars)
  #     ker <- paste(ker, "nrml", sep="_")
  #     intky  <- kernelMatrix(ker, as.matrix(trainData$y), dataWithFakes$y, pars=par)
  #     return(intky)
  #   }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  #   
  # } else if(krnNm == "featureY"){
  #   
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   phiy <- learner$learnParams$phiy #pred_withFakes$gy_f
  #   intky <- mapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     if(is.null(dim(y))) py <- 1 else py <- ncol(y)
  #     numf <- par$num_f
  #     parsPdf <- par[-which(names(par) %in% c("p_w","map","seed","num_f"))]
  #     parsPdf$n <- num_f*py
  #     set.seed(par$seed)
  #     w <- do.call(par$p_w, parsPdf)
  #     w <- matrix(w, py, numf)
  #     b <- matrix(runif(numf, 0, 2*pi), py, numf, byrow=T)
  #     intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w=w[,i], b=b[,i]))
  #     intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(phiy), ncol(phiy), byrow=T)
  #     intky <- phiy %*% t(intPhiy)
  #   }, ker=kernel_char, par=pars, SIMPLIFY="array")
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  # }  
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
  
  probs <- diag(lx%*%Blambda%*%t(ky))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(lx%*%Blambda%*%intky)
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  #sum(probsNorm < 0); sum(probs>1)
  
  # kernel normalization test
  if(FALSE){
    unNormProb <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
      # y <- 0;  x <- 0.1
      n <- length(x_tr)
      ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
      lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
      res <- ky%*%Blambda%*%lx
      return(as.numeric(res))
    }
    normProb2 <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
      # y <- 0;  x <- 0.1
      n <- length(x_tr)
      ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
      intlx <- kern_rbf_nrml(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax, lim_max=0.5, lim_min=-0.5)
      lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
      res <- (ky%*%Blambda%*%lx)/(ky%*%Blambda%*%intlx)
      return(as.numeric(res))
    }
    
    yy <- seq(-0.5,0.5,length.out=100)
    x0 <- dataWithFakes$x[1]
    sigmay <- learner$hyperParams$data$optimizable$sigma.rbf.Y$val
    sigmax <- learner$hyperParams$data$optimizable$sigma.rbf.X$val
    xx_unorm <- unNormProb(xx=yy, y=x0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$x, 
                           x_tr=learner$hyperParams$trainData$y, sigmay=sigmax, sigmax=sigmay)
    xx_norm <- normProb2(xx=yy, y=x0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$x, 
                         x_tr=learner$hyperParams$trainData$y, sigmay=sigmax, sigmax=sigmay)
    
    
    par(mfrow=c(1,2))
    plot(learner$hyperParams$trainData$x, learner$hyperParams$trainData$y); abline(v=x0, col="red") 
    plot(yy, xx_unorm, type="l", ylim=range(xx_unorm, xx_norm))
    lines(yy, xx_norm, type="l", col="red"); abline(h=0, col="blue")
    par(mfrow=c(1,1))
    
    # can I integrate normProb(x) to see if it integrates to one??
    
    sapply(c("Kronrod","Clenshaw","Simpson")[c(1,2,3)], function(meth) integral(normProb2, xmin=-0.5, xmax=0.5, method = meth, 
                                                                                y=x0, Blambda=Blambda, y_tr=learner$hyperParams$trainData$x, x_tr=learner$hyperParams$trainData$y, sigmay=sigmax, sigmax=sigmay))
    
  }
  
  # fourier feature normalization test
  if(FALSE){
    
    count <<- 0
    normProb <- function(yy, py, x, Blambda, Phiy, num_f, y_tr, x_tr, sigmay, sigmax){
      # yy <- 1
      count <<- count + 1
      print(count)
      yy <- matrix(yy,1,py)
      n <- nrow(x_tr)
      lx <- kern_rbf(x=x, y=x_tr, sigma=sigmax)
      phiy <- rff(x=yy, num_f, seed=1234, p_w="rnorm2", map="cos", sigma=sigmay)
      ky <- Phiy %*% t(phiy)
      
      
      parsPdf <- list(sigma=sigmay)
      numf <- num_f
      set.seed(1234)
      parsPdf$n <- num_f*py
      w <- do.call("rnorm2", parsPdf)
      w <- matrix(w, py, numf)
      #w[,1:4]
      b <- matrix(runif(num_f, 0, 2*pi), nrow(yy), numf, byrow=T)
      #print(dim(b))
      #print(b[,1:4])
      
      intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w[,i], b[,i]))
      intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(phiy), ncol(phiy), byrow=T)
      
      intky <- Phiy %*% t(intPhiy)
      
      res <- (lx%*%Blambda%*%ky)/(lx%*%Blambda%*%intky)
      return(as.numeric(res))
    }
    
    py <- 1
    num_f <- ncol(phiy)
    y_tr <- learner$hyperParams$trainData$y
    x_tr <- learner$hyperParams$trainData$x
    x0 <- dataWithFakes$x[5,,drop=F]
    sigmay <- learner$hyperParams$data$non_optimizable$sigma.rbf.Y$val
    sigmax <- learner$hyperParams$data$non_optimizable$sigma.rbf.X$val
    phiy <- learner$learnParams$phiy
    library(cubature)
    ints <- sapply(c("hcubature","pcubature", "cuhre", "divonne", "suave", "vegas"), 
                   function(meth) cubintegrate(normProb, lower=rep(-0.5,py), upper=rep(0.5,py),
                                               method = meth, maxEval=10^3,
                                               x=x0, py=py, Blambda=Blambda, Phiy=phiy, num_f=num_f, y_tr=y_tr, x_tr=x_tr, sigmax=sigmax, sigmay=sigmay))
    
    sapply(ints, function(el) el$integral)
    
  }
  
  # make train classification data
  x_tr <- as.matrix(probsNorm)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsOdds_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  y_fake <- fakes_y
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_real_ext)
  y <- c(y_real, y_fake)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Alambda <- Blambda %*% Ky %*% Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  # hyperParamsNames <- names(learner$getHyperParams(learner))
  # krnNm <- c("kernelY","featureY")[("featureY" %in% hyperParamsNames)*1+1]
  # if(krnNm == "kernelY"){
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   
  #   intky <- mcmapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     maxNumPars <- max(sapply(par, length))
  #     par$lim_max <- rep(0.5, maxNumPars)
  #     par$lim_min <- rep(-0.5, maxNumPars)
  #     ker <- paste(ker, "nrml", sep="_")
  #     intky  <- kernelMatrix(ker, as.matrix(trainData$y), dataWithFakes$y, pars=par)
  #     return(intky)
  #   }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  #   
  # } else if(krnNm == "featureY"){
  #   
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   phiy <- learner$learnParams$phiy #pred_withFakes$gy_f
  #   intky <- mapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     if(is.null(dim(y))) py <- 1 else py <- ncol(y)
  #     numf <- par$num_f
  #     parsPdf <- par[-which(names(par) %in% c("p_w","map","seed","num_f"))]
  #     parsPdf$n <- num_f*py
  #     set.seed(par$seed)
  #     w <- do.call(par$p_w, parsPdf)
  #     w <- matrix(w, py, numf)
  #     b <- matrix(runif(numf, 0, 2*pi), py, numf, byrow=T)
  #     intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w=w[,i], b=b[,i]))
  #     intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(phiy), ncol(phiy), byrow=T)
  #     intky <- phiy %*% t(intPhiy)
  #   }, ker=kernel_char, par=pars, SIMPLIFY="array")
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  # }
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
  
  probs <- diag(lx%*%Blambda%*%t(ky))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(lx%*%Blambda%*%intky)
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  probsNorm[which(probsNorm<=0)] <- 0.00001
  probsNorm[which(probsNorm>=1)] <- 0.99999
  #sum(probsNorm < 0); sum(probs>1)
  odds <- log(probsNorm/(1-probsNorm))
  
  
  # make train classification data
  x_tr <- as.matrix(odds)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_logRegfeatsK_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  y_fake <- fakes_y
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_real_ext)
  y <- c(y_real, y_fake)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Alambda <- Blambda %*% Ky %*% Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  # hyperParamsNames <- names(learner$getHyperParams(learner))
  # krnNm <- c("kernelY","featureY")[("featureY" %in% hyperParamsNames)*1+1]
  # if(krnNm == "kernelY"){
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   
  #   
  #   
  #   intky <- mcmapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     maxNumPars <- max(sapply(par, length))
  #     par$lim_max <- rep(0.5, maxNumPars)
  #     par$lim_min <- rep(-0.5, maxNumPars)
  #     ker <- paste(ker, "nrml", sep="_")
  #     intky  <- kernelMatrix(ker, as.matrix(trainData$y), dataWithFakes$y, pars=par)
  #     return(intky)
  #   }, ker=kernel_char, par=pars, SIMPLIFY = "array")
  #   
  #   
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  #   
  # } else if(krnNm == "featureY"){
  #   
  #   pars <- getKernelPars(learner, kernelName=krnNm)
  #   kernel_char <- pars$kernelChar
  #   kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
  #   pars <- pars$kernelParams
  #   phiy <- learner$learnParams$phiy #pred_withFakes$gy_f
  #   intky <- mapply(function(ker, par){
  #     # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
  #     if(is.null(dim(y))) py <- 1 else py <- ncol(y)
  #     numf <- par$num_f
  #     parsPdf <- par[-which(names(par) %in% c("p_w","map","seed","num_f"))]
  #     parsPdf$n <- num_f*py
  #     set.seed(par$seed)
  #     w <- do.call(par$p_w, parsPdf)
  #     w <- matrix(w, py, numf)
  #     b <- matrix(runif(numf, 0, 2*pi), py, numf, byrow=T)
  #     intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w=w[,i], b=b[,i]))
  #     intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(phiy), ncol(phiy), byrow=T)
  #     intky <- phiy %*% t(intPhiy)
  #   }, ker=kernel_char, par=pars, SIMPLIFY="array")
  #   if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  # }  
  
  dky_with_fakes <- pred_withFakes$gy_f_phiT
  
  probs <- diag(lx%*%Blambda%*%t(ky))
  #sum(probs <0 ); sum(probs>1)
  #const <- diag(lx%*%Blambda%*%intky)
  const <- sqrt(diag(lx%*%Alambda%*%t(lx)))*sqrt(as.numeric(dky_with_fakes))
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  #sum(const<0)
  
  

  
  # make train classification data
  # make train classification data
  #x_tr <- cbind(probs, const)
  x_tr <- as.matrix(probsNorm)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x_tr, y_tr)
  
  return(classData)  
  
}
makeCME_cond_feats_y <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  y_fake <- fakes_y
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_real_ext)
  y <- c(y_real, y_fake)
  dataWithFakes <- constructData(x, as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  hyperParamsNames <- names(learner$getHyperParams(learner))
  krnNm <- c("kernelY","featureY")[("featureY" %in% hyperParamsNames)*1+1]
  if(krnNm == "kernelY"){
    pars <- getKernelPars(learner, kernelName=krnNm)
    kernel_char <- pars$kernelChar
    kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
    pars <- pars$kernelParams
    
    intky <- mcmapply(function(ker, par){
      # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
      maxNumPars <- max(sapply(par, length))
      par$lim_max <- rep(0.5, maxNumPars)
      par$lim_min <- rep(-0.5, maxNumPars)
      ker <- paste(ker, "nrml", sep="_")
      intky  <- kernelMatrix(ker, as.matrix(trainData$y), dataWithFakes$y, pars=par)
      return(intky)
    }, ker=kernel_char, par=pars, SIMPLIFY = "array")
    if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
    
  } else if(krnNm == "featureY"){
    
    pars <- getKernelPars(learner, kernelName=krnNm)
    kernel_char <- pars$kernelChar
    kernel_char <- strsplit(kernel_char, split="T")[[1]][1]
    pars <- pars$kernelParams
    phiy <- learner$learnParams$phiy #pred_withFakes$gy_f
    intky <- mapply(function(ker, par){
      # i <- 1; ker <- kernel_char[i]; par <- pars[[i]]
      if(is.null(dim(y))) py <- 1 else py <- ncol(y)
      numf <- par$num_f
      parsPdf <- par[-which(names(par) %in% c("p_w","map","seed","num_f"))]
      parsPdf$n <- num_f*py
      set.seed(par$seed)
      w <- do.call(par$p_w, parsPdf)
      w <- matrix(w, py, numf)
      b <- matrix(runif(numf, 0, 2*pi), py, numf, byrow=T)
      intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w=w[,i], b=b[,i]))
      intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(phiy), ncol(phiy), byrow=T)
      intky <- phiy %*% t(intPhiy)
    }, ker=kernel_char, par=pars, SIMPLIFY="array")
    if(dim(intky)[3]==1) intky <- intky[,,1] else intky <- apply(intky, c(1,2), sum)
  }
  
  
  probs <- diag(lx%*%Blambda%*%t(ky))
  #sum(probs <0 )
  const <- diag(lx%*%Blambda%*%intky)
  #sum(const<0)
  probs <- probs/const
  #sum(probs < 0)  
  
  p_fake <- 1 # this assumes we have sampled fakes_x uniformly 
  # on [-0.5, 0.5] x ... x [-0.5, 0.5]
  class_probs <- probs/(probs + kappa*p_fake)
  #summary(class_probs)
  #sum(class_probs < 0)
  class_probs[which(class_probs < 0 | class_probs >1)] <- 0.5
  
  y_labs <- c(rep(1, n), rep(0, num_fake)) 
  
  # plot(class_probs, y_labs)
  
  res <- list(x_class=x, gy_class=y_labs, gyh_class=class_probs)
  return(res)
}


# learn and predict functions for different kernel cmem learners
learn.cmem <- function(learner, forLoss=F) {
  
  x <- learner$hyperParams$trainData$x
  y <- as.matrix(learner$hyperParams$trainData$y)
  n <- nrow(x)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  learnBlambda <- learner$learnParams$learnBlambda
  krnlX <- learner$makeKernel(learner, data1=learner$hyperParams$trainData, var="X", grad=(getHyperPar(learner, "gradX") | learnBlambda == "learnBlambda_KCSC"))
  Lx <- krnlX$K
  Cks <- krnlX$Cks
  krnlY <- learner$makeKernel(learner, data1=learner$hyperParams$trainData, var="Y", grad=F)
  Ky <- krnlY$K
  
  
  if(getHyperPar(learner, "centerLx")) Lx <- H%*%Lx%*%H
  if(getHyperPar(learner, "centerKy")) Ky <- H%*%Ky%*%H
  
  lambda <- getHyperPar(learner, "lambda") #*n
  gamma <- getHyperPar(learner, "gamma")
   
  # actual learning part
  if(learnBlambda == "learnBlambda_KCMC"){
    parsBlambda <- list(Lx=Lx,  n=n, lambda=lambda, gamma=gamma)
  } else if(learnBlambda == "learnBlambda_KCSC"){
    parsBlambda <- list(Lx=Lx, Cks=Cks,  centerChat=getHyperPar(learner, "centerLx"),  n=n, lambda=lambda, gamma=gamma)
  } else{
    parsBlambda <- list(Lx=Lx, n=n, lambda=lambda)
  }
  
  Blambda <- try(do.call(learnBlambda, parsBlambda))
  if(class(Blambda)=="try-error") Blambda <- matrix(0,n,n)
  
  beta_k <- Blambda %*% Ky
  alphas <- Blambda %*% Lx
  
  
  learner$learnParams$Lx <- Lx
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- Blambda
  learner$learnParams$Cks <- Cks
  learner$learnParams$beta_k <- beta_k
  learner$learnParams$alphas <- alphas
  
  if(forLoss & "gauss_log_lik" %in% names(learner$optimizeParams$losses)){
    resids <- Ky - Lx%*%beta_k
    cov <- (1/n)*t(resids)%*%resids
    covInv <- try(solve(cov + lambda*I))
    if(class(covInv)=="try-error") covInv <- matrix(0,n,n)
    cholCovInv <- try(chol(covInv, pivot=FALSE))
    if(class(cholCovInv)=="try-error") cholCovInv <- chol(covInv, pivot=TRUE)
    # chose above way coz its slightly faster
    # if(!isposdef(cov) & isposdef(cov, psd=T)){
    #   cholCovInv <- chol(covInv, pivot=TRUE)
    # } else{
    #   cholCovInv <- chol(covInv, pivot=FALSE)
    # }
    learner$learnParams$cholCovInv <- cholCovInv
    
  }
  
  if(forLoss & "negCE" %in% names(learner$optimizeParams$losses)){
    #prepare fake data for classifier training
    NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
    classifiers <- NCE_learner$classifier
    fakeDist_x <-  NCE_learner$fakeDist_x
    fakeDistPars_x <- NCE_learner$fakeDistPars_x
    fakeDist_y <-  NCE_learner$fakeDist_y
    fakeDistPars_y <- NCE_learner$fakeDistPars_y
    kappa <- NCE_learner$kappa
    featFuncs <- NCE_learner$featFunc
    
    
    num_fake <- n*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[smpl] for KDE non parametric dist
    if("y" %in% names(fakeDistPars_y)){
      # to sample we could use some sort of inverse of alpha(x) to give
      # weights to different samples
      smpl <- sample(y, n, replace=T)
      fakeDistPars_y$y <- y[smpl]
    }
    fakeDistPars_x <- c(fakeDistPars_x, list(n=num_fake))
    fakeDistPars_y <- c(fakeDistPars_y, list(n=num_fake))
    
    set.seed(learner$hyperParams$seed+1)
    fakes_x <- matrix(do.call(fakeDist_x, fakeDistPars_x), num_fake, ncol(x))
    fakes_y <- as.matrix(do.call(fakeDist_y, fakeDistPars_y))
    
    
    #trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
    fits <- mcmapply(function(featFunc, classifier){
      # i <- 1; featFunc <- featFuncs[i]; classifier <- classifiers[i]
      print("**************")
      print(paste(featFunc, classifier, sep="_"))
      trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
      
      # train classifier
      classifier <- eval(parse(text=classifier))
      classifier <- setParams(learner=classifier, trainData=trainClassData)
      fit <- classifier$learn(learner=classifier)
      
      return(fit)
    }, 
    featFunc=featFuncs, classifier=classifiers, mc.cores=1, SIMPLIFY=FALSE)
    
    names(fits) <- paste(featFuncs, classifiers, sep="-")
    learner$learnParams$classifier <- fits
  }
  
  return(learner)
}
predict.cmem <- function(learner, data, forLoss=F) {
  
 
  trainData <- learner$hyperParams$trainData
  
  krnlX <- learner$makeKernel(learner, data1=data, data2=trainData, var="X", grad=F)
  lx <- krnlX$K
  krnlY <- learner$makeKernel(learner, data1=data, data2=trainData, var="Y", grad=F)
  ky <- krnlY$K
  krnlY <- learner$makeKernel(learner, data1=data, data2=data, var="Y", grad=F)
  Ky <- krnlY$K
  
  nTe <- nrow(data$x)
  nTr <- nrow(trainData$x)
  Ite <- diag(nTe)
  Itr <- diag(nTr)
  HTe <- Ite-matrix(1/nTe,nTe,nTe)
  HTr <- Itr-matrix(1/nTr,nTr,nTr)
  if(getHyperPar(learner, "centerLx")) lx <- HTe%*%lx%*%HTr
  if(getHyperPar(learner, "centerKy")) ky <- HTe%*%ky%*%HTr
  if(getHyperPar(learner, "centerKy")) Ky <- HTe%*%Ky%*%HTe
  
  
  Blambda <- learner$learnParams$Blambda
  Ky_tr <- learner$learnParams$Ky
  beta_beta_ft <- Blambda %*% Ky_tr %*% t(Blambda)
  
  
  beta_k <- learner$learnParams$beta_k
  
  
  gy_f_phiT <- Ky #similarities between test and test
  gy_k <- ky      #similarities between test and train
  
  gyh_f_phiT <- lx%*%Blambda %*% t(ky) # gyh_f is nte x m, phiT is m x nte
  gyh_gyh_ft <- lx%*%beta_beta_ft%*%t(lx)
  gyh_k  <- lx%*%beta_k # nte x ntr
  
  residsPerObs <- diag(gy_f_phiT) + diag( lx%*%Blambda%*%(-2*t(ky) + Ky_tr%*%t(Blambda)%*%t(lx)))
  #sum(residsPerObs<0); sqrt(mean(residsPerObs))
  
  # All are equivalent but above version seems to be more accurate as it doesn't generate zeros
  #residsPerObs <- diag(gy_f_phiT) + diag( lx%*%(-2*Blambda %*% t(ky) + beta_beta_ft%*%t(lx)))
  #residsPerObs <- diag(gy_f_phiT) - 2*diag(gyh_f_phiT) + diag(gyh_gyh_ft)
  
  
  
  res <- list(x=data$x, y=data$y, lx=lx, gy_f_phiT=as.matrix(diag(gy_f_phiT)), 
              gyh_gyh_ft=as.matrix(diag(gyh_gyh_ft)), gy_k=gy_k, 
              gyh_f_phiT=as.matrix(diag(gyh_f_phiT)), gyh_k=gyh_k, residsPerObs=residsPerObs)
  
  if(forLoss & "gauss_log_lik" %in% names(learner$optimizeParams$losses)){
    resids <- ky - gyh_k
    lambda <- getHyperPar(learner, "lambda")
    logDetSigma_u <- nTr*log(nTr) + 2*nTr*log(lambda) + 
      2*as.numeric(determinant(Ky_tr, logarithm=TRUE)$modulus) +
      2*as.numeric(determinant(Blambda, logarithm=TRUE)$modulus)
    
    res1 <- (nTe^2 / 2)*log(2*pi) 
    res2 <- (nTe/2)*logDetSigma_u
    
    projResids <- resids %*% learner$learnParams$cholCovInv
    
    detConst <- res1 + res2
    resAux <- list(projResids=projResids, detConst=detConst)
    res <- c(res, resAux)
    
  }
  
  
  if(forLoss & "negCE" %in% names(learner$optimizeParams$losses)){
    NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
    classifiers <- learner$learnParams$classifier
    classifiers2 <- NCE_learner$classifier
    fakeDist_x <-  NCE_learner$fakeDist_x
    fakeDistPars_x <- NCE_learner$fakeDistPars_x
    fakeDist_y <-  NCE_learner$fakeDist_y
    fakeDistPars_y <- NCE_learner$fakeDistPars_y
    kappa <- NCE_learner$kappa
    featFuncs <- NCE_learner$featFunc
    
    y <- as.matrix(data$y)
    x <- data$x
    nTe <- nrow(x)
    num_fake <- nTe*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[-indxq] for KDE non parametric dist
    if("y" %in% names(fakeDistPars_y)){
      # to sample we could use some sort of inverse of alpha(x) to give
      # weights to different samples
      smpl <- sample(y, nTe, replace=T)
      fakeDistPars_y$y <- y[smpl]
    }
    fakeDistPars_x <- c(fakeDistPars_x, list(n=num_fake))
    fakeDistPars_y <- c(fakeDistPars_y, list(n=num_fake))
    set.seed(learner$hyperParams$seed)
    fakes_x <- matrix(do.call(fakeDist_x, fakeDistPars_x), num_fake, ncol(x))
    fakes_y <- as.matrix(do.call(fakeDist_y, fakeDistPars_y))
    
    preds_class <- mcmapply(function(featFunc, classifier){
      testClassData <- do.call(featFunc, list(learner=learner, data=data, fakes_x=fakes_x, fakes_y=fakes_y))
      
      
      # test classifier
      
      pred_class <- classifier$pred(learner=classifier, testClassData)
      
      resAux <- list(x_class=pred_class$x_class, gy_class=pred_class$gy_class, gyh_class=pred_class$gyh_class)
    }, featFunc=featFuncs, classifier=classifiers, mc.cores=1, SIMPLIFY=FALSE)
    names(preds_class) <- paste(featFuncs, classifiers2, sep="-")
    
    preds_class <- unlist(preds_class, recursive=F)
    
    res <- c(res, preds_class)
    
  }
  
  
  return(res)
  
}


learn.cmfm <- function(learner, forLoss=F) {
  
  x <- learner$hyperParams$trainData$x
  y <- as.matrix(learner$hyperParams$trainData$y)
  n <- nrow(y)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  learnBlambda <- learner$learnParams$learnBlambda
  krnlX <- learner$makeKernel(learner, data1=learner$hyperParams$trainData, var="X", grad=(getHyperPar(learner, "gradX")|learnBlambda == "learnBlambda_KCSC"))
  #Lx <- krnlX$K 
  Lx <- krnlX$K #+ matrix(1, n, n)
  Cks <- krnlX$Cks
  phiy <- learner$makeFeature(learner, data=learner$hyperParams$trainData, var="Y")
  
  
  if(getHyperPar(learner, "centerLx")) Lx <- H%*%Lx%*%H
  if(getHyperPar(learner, "centerKy")){
    #phiy <- H%*%phiy
    #  phiy <- apply(phiy,2,norml)
    phiy <- H%*%phiy
  }
  
  Ky <- phiy%*%t(phiy)
  lambda <- getHyperPar(learner, "lambda") #*n
  gamma <- getHyperPar(learner, "gamma")
  
  # actual learning part
  
  if(learnBlambda == "learnBlambda_KCMC"){
    parsBlambda <- list(Lx=Lx,  n=n, lambda=lambda, gamma=gamma)
  } else if(learnBlambda == "learnBlambda_KCSC"){
    parsBlambda <- list(Lx=Lx, Cks=Cks,  centerChat=getHyperPar(learner, "centerLx"),  n=n, lambda=lambda, gamma=gamma)
  } else{
    parsBlambda <- list(Lx=Lx, n=n, lambda=lambda)
  }
  Blambda <- try(do.call(learnBlambda, parsBlambda))
  if(class(Blambda)=="try-error"){
    print("errror no Blambda")
     Blambda <- matrix(0,n,n)
  }
  
  beta_f <- Blambda %*% phiy
  beta_k <- Blambda %*% Ky
  alphas <- Blambda %*% Lx
  
  learner$learnParams$Lx <- Lx
  learner$learnParams$phiy <- phiy
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- Blambda
  learner$learnParams$Cks <- Cks
  learner$learnParams$beta_f <- beta_f
  learner$learnParams$beta_k <- beta_k
  learner$learnParams$alphas <- alphas
  
  if(forLoss & "gauss_log_lik" %in% names(learner$optimizeParams$losses)){
    resids <- Ky - Lx%*%beta_k
    cov <- (1/n)*t(resids)%*%resids
    covInv <- try(solve(cov + lambda*I))
    if(class(covInv)=="try-error") covInv <- matrix(0,n,n)
    cholCovInv <- try(chol(covInv, pivot=FALSE))
    if(class(cholCovInv)=="try-error") cholCovInv <- chol(covInv, pivot=TRUE)
    # chose above way coz its slightly faster
    # if(!isposdef(cov) & isposdef(cov, psd=T)){
    #   cholCovInv <- chol(covInv, pivot=TRUE)
    # } else{
    #   cholCovInv <- chol(covInv, pivot=FALSE)
    # }
    learner$learnParams$cholCovInv <- cholCovInv
    
    
  }
  
  # if we need to calculate NCE loss then we train classifier. For example when:
  # 1) calling learn from within CV loop
  # 2) calling learn once a model has been trained
  # but not when...
  # calling learn just outside CV loop
  if(forLoss & "negCE" %in% names(learner$optimizeParams$losses)){
    #prepare fake data for classifier training
    NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
    classifiers <- NCE_learner$classifier
    fakeDist_x <-  NCE_learner$fakeDist_x
    fakeDistPars_x <- NCE_learner$fakeDistPars_x
    fakeDist_y <-  NCE_learner$fakeDist_y
    fakeDistPars_y <- NCE_learner$fakeDistPars_y
    kappa <- NCE_learner$kappa
    featFuncs <- NCE_learner$featFunc
    
  
    num_fake <- n*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[smpl] for KDE non parametric dist
    if("y" %in% names(fakeDistPars_y)){
      # to sample we could use some sort of inverse of alpha(x) to give
      # weights to different samples
      smpl <- sample(y, n, replace=T)
      fakeDistPars_y$y <- y[smpl]
    }
    fakeDistPars_x <- c(fakeDistPars_x, list(n=num_fake))
    fakeDistPars_y <- c(fakeDistPars_y, list(n=num_fake))
    
    set.seed(learner$hyperParams$seed+1)
    fakes_x <- matrix(do.call(fakeDist_x, fakeDistPars_x), num_fake, ncol(x))
    fakes_y <- as.matrix(do.call(fakeDist_y, fakeDistPars_y))
      
    
    #trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
    fits <- mcmapply(function(featFunc, classifier){
      # i <- 4; featFunc <- featFuncs[i]; classifier <- classifiers[i]
      #print("**************")
      #print(paste(featFunc, classifier, sep="_"))
      trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
    
      # train classifier
      classifier <- eval(parse(text=classifier))
      classifier <- setParams(learner=classifier, trainData=trainClassData)
      fit <- classifier$learn(learner=classifier)
      
      return(fit)
    }, 
    featFunc=featFuncs, classifier=classifiers, mc.cores=1, SIMPLIFY=FALSE)
    
    names(fits) <- paste(featFuncs, classifiers, sep="-")
    learner$learnParams$classifier <- fits
  }
  
  return(learner)
}
predict.cmfm <- function(learner, data, forLoss=F){
  
  trainData <- learner$hyperParams$trainData
  
  krnlX <- learner$makeKernel(learner, data1=data, data2=trainData, var="X", grad=F)
  lx <- krnlX$K
  lx <- lx #+ matrix(1, nrow(lx), ncol(lx))
  phiy <- learner$makeFeature(learner, data=data, var="Y")
  
  nTe <- nrow(data$x)
  nTr <- nrow(trainData$x)
  Ite <- diag(nTe)
  Itr <- diag(nTr)
  HTe <- Ite-matrix(1/nTe,nTe,nTe)
  HTr <- Itr-matrix(1/nTr,nTr,nTr)
  if(getHyperPar(learner, "centerLx")) lx <- HTe%*%lx%*%HTr
  if(getHyperPar(learner, "centerKy")){
    #phiy <- HTe%*%phiy
    #phiy <- apply(phiy,2,norml)
    phiy <- HTe%*%phiy
  }
  
  
  # test vs train
  ky <- phiy%*%t(learner$learnParams$phiy)
  # test vs test
  Ky <- phiy%*%t(phiy)
  
  Blambda <- learner$learnParams$Blambda
  Ky_tr <- learner$learnParams$Ky
  beta_beta_ft <- Blambda %*% Ky_tr %*% t(Blambda)
  
  
  beta_f <- learner$learnParams$beta_f
  beta_k <- learner$learnParams$beta_k
  
  
  gy_f <- phiy
  gy_f_phiT <- Ky #similarities between test and test
  gy_k <- ky      #similarities between test and train
  
  
  gyh_f      <- lx%*%beta_f 
  gyh_f_phiT <- lx%*%beta_f %*% t(phiy)
  gyh_gyh_ft <- lx%*%beta_beta_ft%*%t(lx)
  gyh_k      <- lx%*%beta_k
  
  #sum((gyh_f -gy_f)^2)
  #sum(diag(gy_f_phiT) - 2*diag(gyh_f_phiT) + diag(gyh_gyh_ft))
  # I chose to calculate resids here with below formula instead of in loss function
  # with above formula coz above formula, although equivalent was generarting negatives 
  # probably coz its computationaly less accurate. For cmem we have no choice but
  # we calculate ky directly so might not have same problem.
  
  resids <- (gyh_f-gy_f)^2
  residsPerObs <- apply(resids, 1, sum)
  
  res <- list(x=data$x, y=data$y, lx=lx, gy_f=gy_f, gy_f_phiT=as.matrix(diag(gy_f_phiT)), 
              gyh_gyh_ft=as.matrix(diag(gyh_gyh_ft)), gy_k=gy_k, gyh_f=gyh_f, 
              gyh_f_phiT=as.matrix(diag(gyh_f_phiT)), gyh_k=gyh_k, residsPerObs=residsPerObs)

  if(forLoss & "gauss_log_lik" %in% names(learner$optimizeParams$losses)){
    resids <- ky - gyh_k
    lambda <- getHyperPar(learner, "lambda")
    logDetSigma_u <- nTr*log(nTr) + 2*nTr*log(lambda) + 
      2*as.numeric(determinant(Ky_tr, logarithm=TRUE)$modulus) +
      2*as.numeric(determinant(Blambda, logarithm=TRUE)$modulus)
    
    res1 <- (nTe^2 / 2)*log(2*pi) 
    res2 <- (nTe/2)*logDetSigma_u
    
    projResids <- resids %*% learner$learnParams$cholCovInv
    
    detConst <- res1 + res2
    resAux <- list(projResids=projResids, detConst=detConst)
    res <- c(res, resAux)
  }
  
  
  if(forLoss & "negCE" %in% names(learner$optimizeParams$losses)){
    NCE_learner <- learner$hyperParams$data$non_optimizable$NCE_learner$val
    classifiers <- learner$learnParams$classifier
    classifiers2 <- NCE_learner$classifier
    fakeDist_x <-  NCE_learner$fakeDist_x
    fakeDistPars_x <- NCE_learner$fakeDistPars_x
    fakeDist_y <-  NCE_learner$fakeDist_y
    fakeDistPars_y <- NCE_learner$fakeDistPars_y
    kappa <- NCE_learner$kappa
    featFuncs <- NCE_learner$featFunc
    
    y <- as.matrix(data$y)
    x <- data$x
    nTe <- nrow(x)
    num_fake <- nTe*kappa
    # maybe we want a data dependent fakeDist? -> if x among parameters
    # add x = x[-indxq] for KDE non parametric dist
    if("y" %in% names(fakeDistPars_y)){
      # to sample we could use some sort of inverse of alpha(x) to give
      # weights to different samples
      smpl <- sample(y, nTe, replace=T)
      fakeDistPars_y$y <- y[smpl]
    }
    fakeDistPars_x <- c(fakeDistPars_x, list(n=num_fake))
    fakeDistPars_y <- c(fakeDistPars_y, list(n=num_fake))
    set.seed(learner$hyperParams$seed)
    fakes_x <- matrix(do.call(fakeDist_x, fakeDistPars_x), num_fake, ncol(x))
    fakes_y <- as.matrix(do.call(fakeDist_y, fakeDistPars_y))
    
    preds_class <- mcmapply(function(featFunc, classifier){
      testClassData <- do.call(featFunc, list(learner=learner, data=data, fakes_x=fakes_x, fakes_y=fakes_y))
    
    
      # test classifier
    
      pred_class <- classifier$pred(learner=classifier, testClassData)
     
      resAux <- list(x_class=pred_class$x_class, gy_class=pred_class$gy_class, gyh_class=pred_class$gyh_class)
    }, featFunc=featFuncs, classifier=classifiers, mc.cores=1, SIMPLIFY=FALSE)
    names(preds_class) <- paste(featFuncs, classifiers2, sep="-")
    
    preds_class <- unlist(preds_class, recursive=F)
    
    res <- c(res, preds_class)
    
  }
  
  return(res)
  
  
}


learn.cmfm_KCDC_prelim <- function(learner) {
  
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

    #lossKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
    #lossKCDC(as.numeric(alphaF), Lx, phiy, lambda, n)

    grad0 <- gradKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
    gradF <- gradKCDC(as.numeric(alphaF), Lx, phiy, lambda, n)
    plot(grad0, gradF)
    t(grad0)%*%grad0
    t(gradF)%*%gradF
    
    
    
    #plot(alpha0, alphaF); abline(a=0, b=1, col="red")
    #max(abs(alpha0-alphaF))
    
    #Blambda <- alphaF%*%t(phiy)%*%solve(Ky)  
    
  
  learner$learnParams$Lx <- Lx
  learner$learnParams$phiy <- phiy
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- Blambda
  learner$learnParams$alpha <- alphaF
  learner$learnParams$Cks <- Cks
  
  return(learner)
}

plotPredCMEM <- function(learner, pred, var=c("f","k"), indx=1:25){
  
  trainData <- learner$hyperParams$trainData
  pred_learn <- learner$predict(learner, data=trainData)
  
  gy_var <- paste("gy", var, sep="_")
  gyh_var <- paste("gyh", var, sep="_")
  
  df1 <- melt(pred_learn[[gy_var]])
  colnames(df1) <- c("obs","var","value")
  df1$type <- "obse"
  df1$mode <- "learn"
  df2 <- melt(pred_learn[[gyh_var]])
  colnames(df2) <- c("obs","var","value")
  df2$type <- "esti"
  df2$mode <- "learn"
  dfA <- rbind(df1, df2)
  dfA <- cast(dfA, obs+var+mode~type, value="value")
  df3 <- melt(pred_learn$x)
  colnames(df3) <- c("obs","var","value")
  df3$var2 <- "x"
  df4 <- melt(as.matrix(pred_learn$y))
  colnames(df4) <- c("obs","var","value")
  df4$var2 <- "y"
  dfA$x <- df3$value[match(dfA$obs, df3$obs)]
  dfA$y <- df4$value[match(dfA$obs, df4$obs)]
  
  gy_var <- paste("gy", var, sep="_")
  gyh_var <- paste("gyh", var, sep="_")
  df1 <- melt(pred[[gy_var]])
  colnames(df1) <- c("obs","var","value")
  df1$type <- "obse"
  df1$mode <- "pred"
  df2 <- melt(pred[[gyh_var]])
  colnames(df2) <- c("obs","var","value")
  df2$type <- "esti"
  df2$mode <- "pred"
  dfB <- rbind(df1, df2)
  dfB <- cast(dfB, obs+var+mode~type, value="value")
  df3 <- melt(pred$x)
  colnames(df3) <- c("obs","var","value")
  df3$var2 <- "x"
  df4 <- melt(as.matrix(pred$y))
  colnames(df4) <- c("obs","var","value")
  df4$var2 <- "y"
  dfB$x <- df3$value[match(dfB$obs, df3$obs)]
  dfB$y <- df4$value[match(dfB$obs, df4$obs)]
  
  df <- rbind(dfA, dfB)
  
  
  indxs <- which(df$var %in% indx)
  p <- ggplot(df[indxs,])
  p <- p + geom_point(aes(x=x, y=obse, col=mode), alpha=0.3, size=0.5)
  p <- p + geom_line(aes(x=x, y=esti, col=mode))
  p <- p + facet_wrap(var~., scales="free")
  p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
                 axis.text=element_blank(), axis.ticks=element_blank())
  p1 <- p + xlab("x") + ylab(gy_var)
  
  
  p <- ggplot(df[indxs,])
  p <- p + geom_point(aes(x=y, y=obse), col="steelblue", alpha=0.3, size=0.5)
  p <- p + geom_line(aes(x=y, y=obse), col="red")
  p <- p + facet_wrap(var~., scales="free")
  p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
                 axis.text=element_blank(), axis.ticks=element_blank())
  p2 <- p + xlab("y") + ylab(gy_var)
  
  print(p2)
  print(p1)
  
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


cmem_L2_f <- function(learner, pred){
  
  #  tr(phi %*% t(phi)) - 2 tr(phi_hat %*% phi) + tr(phi_hat %*% phi_hat)
  #res <- mean(pred$gy_f_phiT - 2*pred$gyh_f_phiT + pred$gyh_gyh_ft)
  #res <- sqrt(res)
  
  res <- sqrt(mean(pred$residsPerObs))
  
  return(res)
}

cmem_L2_f_rel <- function(learner, pred){
  
  res <- sqrt(mean(as.numeric(t(pred$residsPerObs))/pred$gy_f_phiT^2))
  
  
  return(res)
}

cmem_L2_k <- function(learner, pred){
  
  
  #  tr(K %*% t(K)) - 2 tr(K_hat %*% K) + tr(K_hat %*% K_hat)
  res <- mean(diag(pred$gy_k %*% t(pred$gy_k)) - 2*diag(pred$gyh_k %*% t(pred$gy_k)) + diag(pred$gyh_k  %*% t(pred$gyh_k)))
  #res <- sqrt(res)
  return(res)
}

cmem_L2_k_rel <- function(learner, pred){
  
  res <- diag(pred$gy_k %*% t(pred$gy_k)) - 2*diag(pred$gyh_k %*% t(pred$gy_k)) + diag(pred$gyh_k  %*% t(pred$gyh_k))
  #if(any(res<0)) print("error!")
  res <- mean(res/diag(pred$gy_k%*% t(pred$gy_k)))
  res <- sqrt(res)
  return(res)
}

# this version only for CV-avg
gauss_log_lik_avg <- function(learner, pred){
  
  n_tr <- nrow(learner$hyperParams$trainData$x)
  n_te <- nrow(pred$x)
  ky <- pred$gy_k # n_te x n_tr
  lx <- pred$lx
  Blambda <- learner$learnParams$Blambda
  Ky <- learner$learnParams$Ky
  Dlambda <- Ky %*% Blambda
  Itr <- diag(n_tr) 
  lambda <- getHyperPar(learner, "lambda") 
  Elambda <-  solve(Itr + n_tr*lambda*Dlambda %*% t(Dlambda))
  Elambda <- t(Dlambda) %*% Elambda %*% Dlambda
  
  # log |Sigma_u| = n_tr log n_tr + 2 n_tr log lambda + 2 log abs(|Ky|) -2 log abs(|Blambda|) 
  # log(|Ky %*% Ky|) = log( |Ky|^2) = 2*log(abs(|Ky|))
  # as.numeric(determinant(Ky%*%Ky, logarithm=TRUE)$modulus)
  # 2*as.numeric(determinant(Ky, logarithm=TRUE)$modulus)
  # as.numeric(determinant(Blambda%*%Blambda, logarithm=TRUE)$modulus)
  # 2*as.numeric(determinant(Blambda, logarithm=TRUE)$modulus)
  # https://math.stackexchange.com/questions/2244763/numerically-stable-determinant-of-matrix-product
  
  logDetSigma_u <- n_tr*log(n_tr) + 2*n_tr*log(lambda) + 
    2*as.numeric(determinant(Ky, logarithm=TRUE)$modulus) +
    2*as.numeric(determinant(Blambda, logarithm=TRUE)$modulus)
  
  
  res1 <- (n_te^2 / 2)*log(2*pi) 
  res2 <- (n_te/2)*logDetSigma_u 
  res3 <- (1/(2*lambda))*sum(diag(ky %*% t(ky) - 2*ky %*% Dlambda %*% t(lx) + lx%*% t(Dlambda) %*% Dlambda %*% lx))
  res4 <- (n_tr/2)*sum(diag(ky %*% Elambda %*% t(ky) - 2*ky %*% Elambda %*% Dlambda %*% t(lx) + lx%*% t(Dlambda) %*% Elambda %*% Dlambda %*% lx))
  
  res <- - res1 -res2 - res3 - res4
  
  return(res)
}

gauss_log_lik <- function(learner, pred){
  
  detConst <- pred$detConst
  projResids <- pred$projResids
  res <- sum(detConst) +sum(diag(projResids %*% t(projResids)))
  
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
  
  mcCores <- min(40, detectCores()-4)
  #mcCores <- 1
  msrs_val <- mcmapply(function(msr){
    # msr <- msrs_char[5]
     print(paste("measure: ", msr))
    res <- do.call(msr, list(learner=cmemLearner))
    }, msr=msrs_char, mc.cores=mcCores)
  
  
  names(msrs_val) <- msrs_char
  return(msrs_val)
}

# Kernel Conditional Deviance for Causal inference (KCDC)

KCDC <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
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
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
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
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    #TEMPORARY!!!!!
    Blambda <- beta_f%*%t(beta_f)
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
  pred <- predict.cmfm(learner, data)
  res <- cmem_L2_f_rel(learner, pred)
  return(res)
}

# CV - training relative norm error

PCEte <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  mainLossNm2 <- learner$optimizeParams$mainLoss
  mainLossNm <- strsplit(mainLossNm2, "\\.")[[1]][1] 
  featFunc <- learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
  classifier <- learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier
  if(length(featFunc)>1){
    featFunc <- "makeLogRegFeats_x"
    classifier <- "logRegInt1"
    learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  if(mainLossNm != "negCE") mainLossNm2 <- paste("negCE", paste(featFunc, classifier, sep="-") ,sep=".")
  
  if( mainLossNm2 %in% dimnames(grid)[[3]] & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var ==mainLossNm2)
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- -grid[indxMat]
  } else{
    
    lossFun <-  "negCE" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    
    res <- -grid["test",,mainLossNm2]
  }
    return(res)
}

PCEtr <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  mainLossNm2 <- learner$optimizeParams$mainLoss
  mainLossNm <- strsplit(mainLossNm2, "\\.")[[1]][1] 
  featFunc <- learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
  classifier <- learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier
  if(length(featFunc)>1){
    featFunc <- "makeLogRegFeats_x"
    classifier <- "logRegInt1"
    learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  if(mainLossNm != "negCE") mainLossNm2 <- paste("negCE", paste(featFunc, classifier, sep="-") ,sep=".")
  
  if(mainLossNm2 %in% dimnames(grid)[[3]] & !is.null(grid)){
    indxTrain <- which(dimnames(grid)$trainTest=="train")
    indxLoss <- which(dimnames(grid)$var == mainLossNm2)
    indxMat <- matrix(c(indxTrain, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- -grid[indxMat]
  } else{
    
    lossFun <-  "negCE" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    
    res <- -grid["train",,mainLossNm2]
  }
  return(res)
}

CCRte <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  mainLossNm2 <- learner$optimizeParams$mainLoss
  mainLossNm <- strsplit(mainLossNm2, "\\.")[[1]][1] 
  featFunc <- learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
  classifier <- learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier
  if(length(featFunc)>1){
    featFunc <- "makeLogRegFeats_x"
    classifier <- "logRegInt1"
    learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  if(mainLossNm != "MisCR") mainLossNm2 <- paste("MisCR", paste(featFunc, classifier, sep="-") ,sep=".")
  
  
  if(mainLossNm2 %in% dimnames(grid)[[3]] & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var == mainLossNm2)
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- 1-grid[indxMat]
  } else{
    mainLossNm2 <- paste("CCR", paste(featFunc, classifier, sep="-") ,sep=".")
    lossFun <-  c("CCR","negCE") 
    lossFunList <- lapply(lossFun, function(el) eval(parse(text=el)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    res <- grid["test",,mainLossNm2]
  }
  return(res)
  
}

CCRtr <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  mainLossNm2 <- learner$optimizeParams$mainLoss
  mainLossNm <- strsplit(mainLossNm2, "\\.")[[1]][1] 
  featFunc <- learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
  classifier <- learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier
  if(length(featFunc)>1){
    featFunc <- "makeLogRegFeats_x"
    classifier <- "logRegInt1"
    learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  if(mainLossNm != "MisCR") mainLossNm2 <- paste("MisCR", paste(featFunc, classifier, sep="-") ,sep=".")
  
  if(mainLossNm2 %in% dimnames(grid)[[3]] & !is.null(grid)){
    indxTrain <- which(dimnames(grid)$trainTest=="train")
    indxLoss <- which(dimnames(grid)$var ==mainLossNm2)
    indxMat <- matrix(c(indxTrain, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- 1-grid[indxMat]
  } else{
    
    mainLossNm2 <- paste("CCR", paste(featFunc, classifier, sep="-") ,sep=".")
    lossFun <-  c("CCR","negCE") 
    lossFunList <- lapply(lossFun, function(el) eval(parse(text=el)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    res <- grid["train",,mainLossNm2]
  }
  return(res)
  
}

L2_f_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("cmem_L2_f" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="cmem_L2_f")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "cmem_L2_f" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    
    res <- grid["test",,1]
  }
  return(res)
}

L2_f_tr <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("cmem_L2_f" %in% names(learner$optimizeParams$losses)& !is.null(grid)){
    indxTrain <- which(dimnames(grid)$trainTest=="train")
    indxLoss <- which(dimnames(grid)$var =="cmem_L2_f")
    indxMat <- matrix(c(indxTrain, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "cmem_L2_f" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    
    res <- grid["train",,1]
  }
  return(res)
}

gll_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("gauss_log_lik" %in% names(learner$optimizeParams$losses)& !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="gauss_log_lik")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "gauss_log_lik" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optimHyperParams <- learner$hyperParams$data$optimizable
    
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    res <- grid["test",,1]
  }
  return(res)
}

gll_tr <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("gauss_log_lik" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTrain <- which(dimnames(grid)$trainTest=="train")
    indxLoss <- which(dimnames(grid)$var =="gauss_log_lik")
    indxMat <- matrix(c(indxTrain, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "gauss_log_lik" 
    lossFunList <- list(eval(parse(text=lossFun)))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optimHyperParams <- learner$hyperParams$data$optimizable
    numParams <- length(optimHyperParams)
    if(numParams==0){
      paramsList <- list(dummy=1)
    } else{
      nms <- names(optimHyperParams)
      aux <- strsplit(nms, "\\.")
      indX <- which(nms=="kernelX")
      indY <- which(nms=="kernelY" | nms=="featureY")
      indXX <- which(sapply(aux, function(el) el[3])=="X")
      indYY <- which(sapply(aux, function(el) el[3])=="Y")
      indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
      
      if(length(indX)>0 & length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
            # krnY <- optimHyperParams[[indY]]$seq[1]
            indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
            optimHyperParamsO <- optimHyperParams[indO]
            optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
            optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
            
            paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
            names(paramsListO) <- names(optimHyperParamsO)
            paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
            names(paramsListX) <- names(optimHyperParamsX)
            paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
            names(paramsListY) <- names(optimHyperParamsY)
            paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
            names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
            return(paramsList)
          })
          
          return(paramsListss)
        })
        paramsListss <- unlist(paramsListss, recursive=F)
        
      } else if(length(indX)>0){
        paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
          # krnX <- optimHyperParams[[indX]]$val[1]
          indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
          
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsList <- c(krnX,  paramsListO, paramsListX)
          names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
          return(paramsList)
        })
        
      } else if(length(indY)>0){
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnY, paramsListO, paramsListY)
          names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
          return(paramsList)
        })
        
      } else{
        optimHyperParamsO <- optimHyperParams[indO]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsList <- c(paramsListO)
      }
      
      
    }
    
    otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
    #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
    
    
    params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
    
    
    numFolds <- learner$optimizeParams$numFolds
    grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
    
    if("fold" %in% names(dimnames(grid))){
      grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
    }
    # reshape back into one dimension per hyperparameter
    
    res <- grid["train",,1]
  }
  return(res)
}


CVTeRNE <- function(learner){
  
  lossFun <-  "cmem_L2_f_rel" 
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learner$optimizeParams$losses<- lossFunList 
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  numParams <- length(optimHyperParams)
  if(numParams==0){
    paramsList <- list(dummy=1)
  } else{
    nms <- names(optimHyperParams)
    aux <- strsplit(nms, "\\.")
    indX <- which(nms=="kernelX")
    indY <- which(nms=="kernelY" | nms=="featureY")
    indXX <- which(sapply(aux, function(el) el[3])=="X")
    indYY <- which(sapply(aux, function(el) el[3])=="Y")
    indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
    
    if(length(indX)>0 & length(indY)>0){
      paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
        # krnX <- optimHyperParams[[indX]]$val[1]
        indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
        
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
          names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
          return(paramsList)
        })
        
        return(paramsListss)
      })
      paramsListss <- unlist(paramsListss, recursive=F)
      
    } else if(length(indX)>0){
      paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
        # krnX <- optimHyperParams[[indX]]$val[1]
        indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
        
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
        names(paramsListX) <- names(optimHyperParamsX)
        paramsList <- c(krnX,  paramsListO, paramsListX)
        names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
        return(paramsList)
      })
      
    } else if(length(indY)>0){
      paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
        # krnY <- optimHyperParams[[indY]]$seq[1]
        indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
        
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
        names(paramsListY) <- names(optimHyperParamsY)
        paramsList <- c(krnY, paramsListO, paramsListY)
        names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
        return(paramsList)
      })
      
    } else{
      optimHyperParamsO <- optimHyperParams[indO]
      paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
      names(paramsListO) <- names(optimHyperParamsO)
      paramsList <- c(paramsListO)
    }
    
    
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
  
  
  params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
  
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  
  return(grid["test",,1])
  
}
  
CVTrRNE <- function(learner){
  
  lossFun <-  "cmem_L2_f_rel" 
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learner$optimizeParams$losses<- lossFunList 
  
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  numParams <- length(optimHyperParams)
  if(numParams==0){
    paramsList <- list(dummy=1)
  } else{
    nms <- names(optimHyperParams)
    aux <- strsplit(nms, "\\.")
    indX <- which(nms=="kernelX")
    indY <- which(nms=="kernelY" | nms=="featureY")
    indXX <- which(sapply(aux, function(el) el[3])=="X")
    indYY <- which(sapply(aux, function(el) el[3])=="Y")
    indO <- setdiff(1:numParams, c(indX, indY, indXX, indYY))
    
    if(length(indX)>0 & length(indY)>0){
      paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
        # krnX <- optimHyperParams[[indX]]$val[1]
        indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
        
        paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
          # krnY <- optimHyperParams[[indY]]$seq[1]
          indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
          optimHyperParamsO <- optimHyperParams[indO]
          optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
          optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
          
          paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
          names(paramsListO) <- names(optimHyperParamsO)
          paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
          names(paramsListX) <- names(optimHyperParamsX)
          paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
          names(paramsListY) <- names(optimHyperParamsY)
          paramsList <- c(krnX, krnY, paramsListO, paramsListX, paramsListY)
          names(paramsList)[1:2] <- names(optimHyperParams)[c(indX, indY)]
          return(paramsList)
        })
        
        return(paramsListss)
      })
      paramsListss <- unlist(paramsListss, recursive=F)
      
    } else if(length(indX)>0){
      paramsListss <- lapply(optimHyperParams[[indX]]$val, function(krnX){
        # krnX <- optimHyperParams[[indX]]$val[1]
        indKrnX <- which(sapply(strsplit(nms[indXX], "\\."), function(el) el[2])==strsplit(krnX, "_")[[1]][2])
        
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsX <- optimHyperParams[indXX[indKrnX]]
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListX <- lapply(optimHyperParamsX, function(el) el$val)
        names(paramsListX) <- names(optimHyperParamsX)
        paramsList <- c(krnX,  paramsListO, paramsListX)
        names(paramsList)[1] <- names(optimHyperParams)[c(indX)]
        return(paramsList)
      })
      
    } else if(length(indY)>0){
      paramsListss <- lapply(optimHyperParams[[indY]]$val, function(krnY){
        # krnY <- optimHyperParams[[indY]]$seq[1]
        indKrnY <- which(sapply(strsplit(nms[indYY], "\\."), function(el) el[2])==strsplit(krnY, "_")[[1]][2])
        optimHyperParamsO <- optimHyperParams[indO]
        optimHyperParamsY <- optimHyperParams[indYY[indKrnY]]
        
        paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
        names(paramsListO) <- names(optimHyperParamsO)
        paramsListY <- lapply(optimHyperParamsY, function(el) el$val)
        names(paramsListY) <- names(optimHyperParamsY)
        paramsList <- c(krnY, paramsListO, paramsListY)
        names(paramsList)[1] <- names(optimHyperParams)[c(indY)]
        return(paramsList)
      })
      
    } else{
      optimHyperParamsO <- optimHyperParams[indO]
      paramsListO <- lapply(optimHyperParamsO, function(el) el$val)
      names(paramsListO) <- names(optimHyperParamsO)
      paramsList <- c(paramsListO)
    }
    
    
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
  
  
  params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
  
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  
  return(grid["train",,1])
  
}


CVTeRNSE <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "cmem_L2_sd" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list(dummy=1)
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], list("opt"), dimnms["var"])
  names(dimnms) <- c("trainTest", "params", "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,1])
  
}

teErrNorm <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "TNRE" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses <- lossFunList 
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list(dummy=1)
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], list("opt"), dimnms["var"])
  names(dimnms) <- c("trainTest", "params", "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,1])
  
}

trErrNorm <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "TNRE" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list(dummy=1)
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
  if("fold" %in% names(dimnames(grid))){
    grid <- apply(grid, setdiff(names(dimnames(grid)), "fold"), mean) #, na.rm=T add if its ok for some folds to have computationally unfeasible matrices
  }
  
  # reshape back into one dimension per hyperparameter
  
  dimnms <- dimnames(grid)
  #log(as.numeric(sapply(strsplit(strsplit(names(which.min(grid["test",,])), " ")[[1]] ,"="), function(el) el[2])),10)
  #grid2 <- grid
  dim(grid) <- c(dim(grid)[1], 1 , dim(grid)[3])
  dimnms <- c(dimnms["trainTest"], list("opt"), dimnms["var"])
  names(dimnms) <- c("trainTest", "params", "var")
  dimnames(grid) <- dimnms
  
  return(grid["test",,1])
  
}

teNumUnd <- function(learner){
  learnerAux <- learner
  
  lossFun <-  "numUnd" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list(dummy=1)
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
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
  
  optimHyperParams <- learner$hyperParams$data$optimizable
  
  if(length(optimHyperParams)==0){
    paramsList <- list(dummy=1)
  } else{
    paramsList <- lapply(optimHyperParams, function(el) el$val)
    names(paramsList) <- names(optimHyperParams)
  }
  
  otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  params <- do.call(constructParams, list(otherParams=otherParams, paramsList=paramsList))
  
  
  numFolds <- learner$optimizeParams$numFolds
  grid <- CV.parallel(learner, params, fac=1, numCoresFold=min(numFolds, detectCores()-2))
  
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
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f) & is.null(learner$learnParams$alphas)){
    phiy <- learner$learnParams$phiy
    Alambda <- beta_f%*%t(beta_f)
    alphas <- try(solve(K, phiy%*%t(beta_f)%*%L))
  } else{
    #Blambda <- learner$learnParams$Blambda
    #Alambda <- Blambda%*%K%*%t(Blambda)
    alphas <- learner$learnParams$alphas
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
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f) & is.null(learner$learnParams$alpha)){
    phiy <- learner$learnParams$phiy
    Alambda <- beta_f%*%t(beta_f)
    alphas <- try(solve(K, phiy%*%t(beta_f)%*%L))
    
  } else{
    #Blambda <- learner$learnParams$Blambda
    #Alambda <- Blambda%*%K%*%t(Blambda)
    #alphas <- t(L)%*%Blambda
    alphas <- learner$learnParams$alphas
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
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
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
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
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
KCSC <- function(learner, pred=NULL){
  
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
  
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    phiy <- learner$learnParams$phiy
    D <- t(beta_f)%*%Chat%*%beta_f
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
KCCC_ent <- function(learner, pred=NULL){
  
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
  
  
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
    D <- t(beta_f)%*%Chat%*%beta_f
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
  
  
  #plot(S,S2); abline(a=0, b=1, col="red")
  
  Stilde <- S/T
  res <- sum(abs(Stilde)*log(abs(Stilde)))
  
  return(res)
}

KCCC_pca_ent <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  phiy <- learner$learnParams$phiy
  
  
  res.pca <- prcomp(phiy, scale = F)
  #library(factoextra)
  #fviz_eig(res.pca)
  #res.pca <- princomp(covmat=phiy%*%t(phiy), cor=F)
  
  phiy <- res.pca$x
  K <- phiy %*% t(phiy)
  
  
  
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
  
  
  
  D <- t(Blambda)%*%Chat%*%Blambda
  S <- diag(t(phiy)%*%D%*%phiy)
  T <- sum(D*K)
  
  
  
  
  
  # n <- 100; m <- 10; phiy <- matrix(rnorm(n*m),n,m)
  
  #proc.time()
  #S <- sapply(1:n, function(i) sapply(1:n, function(j) D[i,j]*phiy[i,]*phiy[j,], simplify="array"), simplify="array")
  #S <- apply(S, 1, sum)
  #proc.time() -pm # 15 secs
  
  
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
  
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
    D <- t(beta_f)%*%Chat%*%beta_f
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

KCNSC <- function(learner, pred=NULL){
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky 
  n  <- nrow(L)
  Cks <- learner$learnParams$Cks
  Blambda <- learner$learnParams$Blambda
  
  beta_f <- learner$learnParams$beta_f
  if(!is.null(beta_f)){
    Alambda <- beta_f%*%t(beta_f)
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



