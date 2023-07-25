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
#library(neuralnet) #learn.nnc
library(corpcor) # is.positive.definite
library(quadprog) # learn.lapSVM
library(vegan) # isomap in learn.lapSVM
library(pROC) # roc in myAUC
library(magrittr) # %>% operator
library(tensor) # learn.lnkrr
library(modi) #weighted.var
library(wCorr) # weightedCorr
library(pdfCluster) #kepdf
library(statmod) # gauss.quad.prob
library(quantregForest) # quantregForest in fqr
library(rvinecopulib) # bicop in cqr
library(qrnn) # mcqrnn.fit in nnqr


#############################################################################################################*
# learners - general learner functions 
#############################################################################################################*

constructData <- function (...){
  #stopifnot(is.list(x) || is.vector(x) || is.matrix(x))
  #stopifnot(is.list(y) || is.vector(y) || is.matrix(y))
  data = list( ...)
  class(data) = "CVST.data"
  return(data)
}

getSubset <- function (data, subset) {
  stopifnot(class(data) == "CVST.data")
  x = getX(data, subset)
  y = data$y[subset]
  if("s" %in% names(data)){
    s = data$s[subset,,drop=F]
    ret = constructData(x = x, y = y, s=s, indxU=subset)
  } else{
    ret = constructData(x = x, y = y, indxU=subset)  
  }
  
  return(ret)
}

getSubset <- function (data, subset) {
  stopifnot(class(data) == "CVST.data")
  
  dat <- list()
  if("bag" %in% names(data)){
    indx <- which(data$bag %in% abs(subset))
    signo <- sign(subset)[match(data$bag[indx], abs(subset))]
    subset2 <- signo*indx
    dat[["y"]] <- data[["y"]][subset]
    dat[["x"]] <- data[["x"]][subset2,,drop=FALSE]
    # re-index the bags
    dat[["bag"]] <- as.numeric(as.factor(data[["bag"]][subset2]))  
    if("M" %in% names(data)){
      M <- data$M
      M <- M[,subset]
      M <- M[subset,]
      dat[["M"]] <- M
    }
    
  } else{
    for(nm in names(data)){
      # nm <- "x"
      if(is.null(dim(data[[nm]]))){
        dat[[nm]] <- data[[nm]][subset]
      } else{
          dat[[nm]] <- data[[nm]][subset,,drop=FALSE]
      }
    }
  }
  
  dat$indxU <- subset
  
  class(dat) = "CVST.data"
  
  
  
  return(dat)
}


getN <- function (data){
  stopifnot(class(data) == "CVST.data")
  if("bag" %in% names(data)){
    N <- length(unique(data$bag))
  } else{
    if (is.list(data$x) || is.vector(data$x)) {
      N = length(data$x)
    }
    else {
      N = nrow(data$x)
    }
  }
  return(N)
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
                      "negLogLik", "hsicLoss2", "pinball", "hingeLoss","gauss_log_lik","negCE","myAUC","RMSE2","PCE","MisCR","KCDC","KCRDC","KCMC","KCSC","KCNSC","KCCC_ent","KCCC_pca_ent")
  #stopifnot( all(sapply(optimizeParams$losses, function(el) is.function(el))) && all(names(optimizeParams$losses) %in% validLossFuncs))
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
setParams <- function(learner, trainData, mc_cores=1, numCoresFold=1, plot=FALSE){
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
  # in case I want to pass parameters to optimizeSet func. For ex when optimizeSet = optHP.CV
  # and i want to do semi supervised learning I dont want to use CV.parallel the default param
  # instead i want to use CV.SSL.parallel
  parsOptimizeSet <- list(learner=learner, mc_cores=mc_cores, numCoresFold=numCoresFold, plot=plot)
  if("optimizeSetPars" %in% names(learner)){
    parsOptimizeSet <- c(parsOptimizeSet, learner$optimizeSetPars)
  }
  # predFunc <- parsOptimizeSet$predFunc
  # CVfunc <- parsOptimizeSet$CVfunc
  dataOptimParams <- do.call(learner$optimizeSet, parsOptimizeSet)
  
  #dataOptimParams$opt
  #dataOptimParams$grid["test",,,"negLogLik"]
  
  for(nm in paramsOpt){
    # nm <- paramsOpt[1]
    learner$hyperParams$data$optimizable[[nm]]$val <- dataOptimParams$opt[[nm]]
    
  }
  
  learner$hyperParams$data$grid <- dataOptimParams$grid
  learner$hyperParams$data$gridFold <- dataOptimParams$gridFold
  
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
           ylim=range(predList[[1]]$gy,sapply(predList, function(el) range(el$gyh,na.rm=T))))
      for(i in 1:length(predList)){
        # i <- 1
        xx <- predList[[i]]$x
        o <- order(xx)
        xx <- xx[o]
        yy <- predList[[i]]$gyh
        yy <- yy[o,,drop=FALSE]
        
        
        for(j in 1:ncol(yy)){
          # j <- 1
          lines(xx,yy[,j], col=j+1, lwd=2)
        }
      }
      #legend("bottomleft", legend=names(predList), lwd=2, col="red")
      
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

# this method should really be a method of kernel learners but since we are only using a 
# general pseudo-learner class we'll leave it outside for now
getKernelPars2 <- function(learner, kernelName){
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
  if(is.infinite(sigma0)) sigma0 <- 1
  
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

getFixedParams_rbf_dr <- function(learner, data){}

getFixedParams_tauQuad <- function(learner, data){
  m <- learner$hyperParams$data$non_optimizable$taus$length.out
  uw <- gauss.quad.prob(m)
  
  return(list(val=list(taus=uw$nodes)))
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
getFixedParams_rbf2b <- function(learner, data){
  
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
  
  
  numX <- learner$hyperParams$data$optimizable$sigma$length.out
  if(is.null(numX)) numX <-10
  
  rngX <- getRangeRbf2(x, length.out=numX)
  
  #rng <- range(rngX, rngY)
  #rngX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numX)
  #rngY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numY)
  
  return(list(val=list(sigma=sigmaX), seq=list(sigma=rngX)))
}
getFixedParams_rbf_indSig <- function(learner, data){
  
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
  
  
  sigmaX <- apply(x,2, function(col) 1/median(as.numeric(dist(unique(col))^2)))

  
  return(list(val=list(sigmas=sigmaX)))
}

wcorrRKHS_col <- function(sigmay, B, Lte_tr, ytr, Lte=NULL, yte=NULL, method="pearson"){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- Lte_tr%*%B
  
  if(is.null(yte)){
    Lte <- Lte_tr
    Ktr_te <- K 
  }else{
    Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
  }
  
  LBK <- LB%*%Ktr_te
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    #res <- cor(Lte_col,LBK_col, method="spearman")
    res <- weightedCorr(Lte_col,LBK_col, method=method, weights=Lte_col/sum(Lte_col))
    #print(res)
    return(res)
  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- match(sort(corrs,decreasing=T),corrs)[3]
  #plot(Lte[indx,], LBK[indx,], cex=Lte[indx,])
  #res <- mean(corrs)
  res <- corrs
  return(res)
}
getFixedParams_rbf3 <- function(learner, data){
  #print("enters getFixedParams_rbf3")
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
  
  sigmaxSeq <- learner$hyperParams$data$optimizable$sigma.rbf.X$seq
  sigmayL <- learner$hyperParams$data$optimizable$sigma.rbf.Y
  sigmaySeq <- sigmayL$seq
  lambdaSeq <- learner$hyperParams$data$optimizable$lambda$seq
  
  
  if(any(is.null(c(sigmaxSeq,lambdaSeq)))){
    
    pm <- proc.time()
    krr2_aux <- setParams(learner=eval(parse(text="krr2b")),trainData=data,mc_cores=round(detectCores()/2))
    proc.time() - pm 
    
    if(FALSE){
      krr2_aux <- krr2_aux$learn(krr2_aux)
      newData <- constructData(x=matrix(seq(-0.5,0.5, length.out=1000)), y=rnorm(1000))
      pred_aux <- krr2_aux$predict(krr2_aux, newData)
      plot(krr2_aux$hyperParams$trainData$x,krr2_aux$hyperParams$trainData$y)
      xx <- pred_aux$x
      yy <- pred_aux$gyh
      o <- order(xx)
      lines(xx[o],yy[o], col="red")
    }
    
    rmseTest <- krr2_aux$hyperParams$data$grid["test",,"rmse"]
    names(rmseTest) <- NULL
    df <- expand.grid(lambda=krr2_aux$hyperParams$data$optimizable$lambda$seq, sigma=krr2_aux$hyperParams$data$optimizable$sigma$seq)
    df$rmseTest <- rmseTest
    
    if(FALSE){
      v <- ggplot(df)
      v <- v + geom_raster(aes(x=log(lambda,10), y=log(sigma,10),fill = rmseTest), data=df) 
      v <- v + geom_contour(aes(x=log(lambda,10), y=log(sigma,10), z=rmseTest),colour = "white", bins = 10, data=df)
      print(v)
    }
    
    barLevel <- 1.1
    indxGood <- which(df$rmseTest/min(df$rmseTest) <barLevel)
    while(length(indxGood)==1){
      barLevel <- barLevel+0.1
      indxGood <- which(df$rmseTest/min(df$rmseTest) <barLevel)
    }
    lambdaRng <- df$lambda[indxGood]
    sigmaxRng <- df$sigma[indxGood]
    
    if(length(unique(sigmaxRng))==1){
      sigmaxRng <- unique(df$sigma[order(df$rmseTest, decreasing=FALSE)])[1:2]
    }
    
    lambdaRng <- range(lambdaRng)  
    sigmaxRng <- range(sigmaxRng)  
    numX <- learner$hyperParams$data$optimizable$sigma.rbf.X$length.out
    if(is.null(numX)) numX <-5
    numL <- learner$hyperParams$data$optimizable$lambda$length.out
    if(is.null(numL)) numL <-5
    lambdaSeq <- 10^seq(log(lambdaRng[1],10), log(lambdaRng[2],10), length.out=numL)
    sigmaxSeq <- 10^seq(log(sigmaxRng[1],10), log(sigmaxRng[2],10), length.out=numX)
    #lambdaSeq <- 10^seq(log(lambdaRng[1],10), log(lambdaRng[2],10), length.out=10)
    #sigmaxSeq <- 10^seq(log(sigmaxRng[1],10), log(sigmaxRng[2],10), length.out=10)
  }
  
  if(is.null(sigmaySeq) & !is.null(sigmayL)){
    
    
    
    xTe <- data$x[floor(n/2+1):n,,drop=F]
    xTr <- data$x[1:(floor(n/2)),,drop=F]
    yTe <- as.matrix(data$y[floor(n/2+1):n])
    yTr <- as.matrix(data$y[1:(floor(n/2))])
    nTr <- floor(n/2)
    I <- diag(nTr)  
    
    if(FALSE){
      corrsRKHS <- sapply(lambdaSeq, function(lam){ 
        sapply(sigmaxSeq, function(sigmax){
          # i <- 2; j <- 3; lam <- lambdaSeq[i]; sigmax <- sigmaxSeq[j]
          L <- kern_rbf(xTr, sigma=sigmax)
          Ltr_te <- kern_rbf(xTe,xTr, sigma=sigmax)
          Blambda <- solve(L+nTr*lam*I)
          corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda, Ltr_te, ytr=yTr, yte=yTe ))/(nTr^2)
          # plot(log(sigmasy,10), corrsRKHS_te)
          return(corrsRKHS_te)
          
        }, simplify="array")}, simplify="array")
      dim(corrsRKHS)
      dimnames(corrsRKHS) <- list(sigmaY=sigmasy, sigmaX=sigmaxSeq, lambda=lambdaSeq) 
      
      dfCorrRKHS <- melt(corrsRKHS)
      colnames(dfCorrRKHS)[4] <- "corrRKHS"
      
      dfCorrRKHS$logLambda <- round(log(dfCorrRKHS$lambda, 10),2)
      dfCorrRKHS$logSigmaX <- round(log(dfCorrRKHS$sigmaX,10),2)
      dfCorrRKHS <- dfCorrRKHS[which(!is.nan(dfCorrRKHS$corrRKHS)),]
      #dfCorrRKHS <- dfCorrRKHS[which(dfCorrRKHS$corrRKHS/max(dfCorrRKHS$corrRKHS) >0.25),]
      
      p <- ggplot(dfCorrRKHS)
      p <- p + geom_line(aes(x=log(sigmaY,10), y=corrRKHS))
      p <- p + geom_point(aes(x=log(sigmaY,10), y=corrRKHS))
      p <- p + facet_grid(logLambda~logSigmaX, scales="free")
      p
    }
    
    #i <- 0
    #j <- 0
    
    limsSigmaY <- sapply(lambdaSeq, function(lam){ 
      #i <<- i +1
      sapply(sigmaxSeq, function(sigmax){
        #j <<- j + 1
        #print(paste("i: ",i, " j: ",j))
        # i <- 1; j <- 1; lam <- lambdaSeq[i]; sigmax <- sigmaxSeq[j]
        
        L <- kern_rbf(xTr, sigma=sigmax)
        
        Ltr_te <- kern_rbf(xTe,xTr, sigma=sigmax)
        
        Blambda <- solve(L+nTr*lam*I)
        
        corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda, Ltr_te, ytr=yTr, yte=yTe ))/(nTr^2)
        
        # plot(log(sigmasy,10), corrsRKHS_te)
        if(!all(is.nan(corrsRKHS_te))){
          
          xx <- log(sigmasy,10)
          yy <- corrsRKHS_te
          # plot(xx, yy)
          # get the indices of the sequence of max consecutive neg difs
          aux <- diff(yy)<0
          
          maxConsecs <- sapply(1:length(aux), function(ini){
            indxConsec <- cumprod(aux[ini:length(aux)])
            res <- which(indxConsec==1)
            if(length(res)>0) {
              res <- res[length(res)]
            } else{
              res <- -1
            }
            return(res)
          })
          
          indxMax <- which.max(maxConsecs)
          indxConsec <- indxMax:(indxMax+maxConsecs[indxMax])
          xx <- xx[indxConsec]
          yy <- yy[indxConsec]
          # plot(xx, yy)
          corrsRKHS_te_sp <- splinefun(xx,yy, method="monoH.FC")
          
          xxx <- seq(min(xx),max(xx), length.out=1000)
          yyy <- corrsRKHS_te_sp(xxx)
          dyyy_dx <- corrsRKHS_te_sp(xxx,deriv=1)
          
          chngPhase <- which(dyyy_dx < -0.01)
          #print(paste("length(chngPhase): ", length(chngPhase)))
          if(length(chngPhase)<2){
            chngPhase <- which(dyyy_dx < quantile(dyyy_dx, probs=0.05))
          }
          #print("sort(dyyy_dx)")
          #print(sort(dyyy_dx))
          ini_chngPhase <- chngPhase[1]
          end_chngPhase <- chngPhase[length(chngPhase)]
          mark <- sum(corrsRKHS_te_sp(xxx[c(ini_chngPhase, end_chngPhase)]))/2
          indxFin <- which.min(abs(yyy-mark))
          
          #plot(xx, yy)
          #lines(xxx,yyy, col="red")
          #abline(v=xxx[c(ini_chngPhase,end_chngPhase)],col="blue")
          #abline(h=mark,col="green")
          #abline(v=xxx[indxFin],col="green")
          
          res <- 10^xxx[c(ini_chngPhase, indxFin)]
        } else{
          res <- c(NA,NA)
        }
        return(res)
        
      }, simplify="array")}, simplify="array")
    
    #dim(limsSigmaY)
    sigmayRng <- c(min(limsSigmaY[1,,], na.rm=T),max(limsSigmaY[2,,], na.rm=T))
    #print(paste("sigmayRng:", paste(sigmayRng, collapse="-")))
    numY <- learner$hyperParams$data$optimizable$sigma.rbf.Y$length.out
    if(is.null(numY)) numY <-5
    sigmaySeq <- 10^seq(log(sigmayRng[1],10), log(sigmayRng[2],10), length.out=numY)
  } 
  
  sigmaX <- 1/median(as.numeric(dist(unique(x))^2))
  sigmaY <- 1/median(as.numeric(dist(unique(y))^2))
  
  #print("exits getFixedParams_rbf3")
  return(list(val=list(sigma.rbf.X=sigmaX, sigma.rbf.Y=sigmaY, lambda=1), seq=list(sigma.rbf.X=sigmaxSeq, sigma.rbf.Y=sigmaySeq, lambda=lambdaSeq)))
}
getFixedParams_rbf4 <- function(learner, data){
  #print("enters getFixedParams_rbf4")
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
    data <- getSubset(data, 1:10000) 
    n <- 1000
  }
  
  x <- matrix(data$x, n, px)
  y <- matrix(data$y, n, py)
  
  
  densX <- kepdf(x, eval.points = x, kernel = "gaussian", bwtype = "adaptive")
  
  sigmaX <- 1/median(as.numeric(dist(unique(x))^2))
  sigmaY <- 1/median(as.numeric(dist(unique(y))^2))
  
  sigmaxSeq <- learner$hyperParams$data$optimizable$sigma.rbf.X$seq
  sigmayL <- learner$hyperParams$data$optimizable$sigma.rbf.Y
  sigmaySeq <- sigmayL$seq
  lambdaSeq <- learner$hyperParams$data$optimizable$lambda$seq
  
  #print(paste("sigmayL: ", sigmayL))
  #print(paste("names(sigmayL): ", names(sigmayL)))
  #if(!is.null(sigmaxSeq)) print(paste("log-sigmaxSeq: ", paste(round(log(sigmaxSeq,10),2),collapse=",")))
  #if(!is.null(sigmaySeq)) print(paste("log-sigmaySeq: ", paste(round(log(sigmaySeq,10),2),collapse=",")))
  #if(!is.null(lambdaSeq)) print(paste("log-lambdaSeq: ", paste(round(log(lambdaSeq,10),2),collapse=",")))
  
  if(any(is.null(c(sigmaxSeq,lambdaSeq)))){
    #print("lambda sigmax seq")  
    numX <- learner$hyperParams$data$optimizable$sigma.rbf.X$length.out
    if(is.null(numX)) numX <-5
    sigmaxSeq <- getRangeRbf2(x, length.out=numX)
    
    numL <- learner$hyperParams$data$optimizable$lambda$length.out
    if(is.null(numL)) numL <-5
    lambdaSeq <- 10^seq(-7, 1, length.out=numL)
    
  }
  
  if(is.null(sigmaySeq) & !is.null(sigmayL)){
    
    
    
    sigmasy <- 10^seq(log(sigmaY,10),16,length.out=30)
    nboots <- 100
    set.seed(123)
    smpl_tr <- sapply(1:nboots, function(i) sample(n, replace=T), simplify="matrix")
    smpl_te <- sapply(1:nboots, function(i) sample(setdiff(1:n, unique(smpl_tr[,i])), n+1, replace=T), simplify="matrix")
    
    xTe <- x #[51:100,,drop=F]
    xTr <- x #[1:50,,drop=F]
    yTe <- as.matrix(y)#[51:100,])
    yTr <- as.matrix(y)#[1:50,])
    nTr <- nrow(xTr)
    nTe <- nrow(xTe)
    n <- nrow(x) 
    Itr <- diag(nTr)  
    I <- diag(n)  
    
    parmTab <- expand.grid(lam=lambdaSeq, sigmax=sigmaxSeq)
    
    pm <- proc.time()
    corrsRKHS_boot <- mcmapply(FUN=function(j){ 
      # j <- 1
      lam <- parmTab$lam[j]
      sigmax <- parmTab$sigmax[j]
      print(paste("j: ", j, " lam:", lam, " sigmax", sigmax))
      Lte <- kern_rbf(xTe, sigma=sigmax)
      Ltr <- kern_rbf(xTr, sigma=sigmax)
      Lte_tr <- kern_rbf(xTe,xTr, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr*lam*Itr)    
      
      difYs <- sapply(sigmasy, function(sigmay){
        yAux <- sort(unique(yTr))
        yMin <- min(yAux[2:length(yAux)]-yAux[1:(length(yAux)-1)])*0.01
        yAux <- yTr + rnorm(length(yTr),mean=0, sd=yMin)
        
        K <- kern_rbf(yAux, sigma=sigmay)
        quantile(apply(K,2, function(col) length(unique(col))), 0.05)
      })/nTr
      
      difY1s <- sapply(sigmasy, function(sigmay){
        yAux <- sort(unique(yTr))
        yMin <- min(yAux[2:length(yAux)]-yAux[1:(length(yAux)-1)])*0.01
        yAux <- yTr + rnorm(length(yTr),mean=0, sd=yMin)
        K <- kern_rbf(yAux, sigma=sigmay)
        min(apply(K,2, function(col) sum(col!=1)))
      })/nTr
      
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:nboots, function(i){
        #i <- 1
        
        xTe_b <- x[smpl_te[,i],,drop=F]
        xTr_b <- x[smpl_tr[,i],,drop=F]
        yTe_b <- as.matrix(y[smpl_te[,i],])
        yTr_b <- as.matrix(y[smpl_tr[,i],])
        nTr_b <- nrow(xTr_b)
        nTe_b <- nrow(xTe_b)
        Itr_b <- diag(nTr_b)  
        
        Ltr <- kern_rbf(xTr_b, sigma=sigmax)
        Lte <- kern_rbf(xTe_b, sigma=sigmax)
        Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
        Blambda_tr <- solve(Ltr+nTr*lam*Itr)
        if(FALSE){
          LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
          plot(c(Lte), c(LBK)) 
          plot(c(Lte[21,]), c(LBK[21,])) 
          cor(c(Lte[21,]), c(LBK[21,]))
        }
        
        corrsRKHS_te <- sapply(sigmasy, function(sigmay) wcorrRKHS_col(sigmay, B=Blambda_tr, Lte_tr=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b ), simplify="array")
        
        return(corrsRKHS_te)
      }, simplify="array")
      proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
      dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nTe+1), sigmay=sigmasy, boot=1:nboots)
      
      # idea is to remove from consideration sigmay values where there
      # is vary little variation between distances 
      corrsRKHS_te_boot[,which(difYs<(quantile(difYs,0.1)+0.001)|difY1s<max(difY1s)),] <- NA
      
      return(corrsRKHS_te_boot)
      
    }, j=1:nrow(parmTab), SIMPLIFY="array", mc.cores=1)#round(detectCores()/2)
    proc.time() - pm #21 mins
    
    #dimnames(corrsRKHS_boot)[4:5] <- list(lambda=lambdaSeq, sigmax=sigmaxSeq)
    #names(dimnames(corrsRKHS_boot))[4:5] <- c("lambda","sigmax")
    dimnames(corrsRKHS_boot)[4] <- list(i=1:nrow(parmTab))
    names(dimnames(corrsRKHS_boot))[4] <- "i"
    
    df <- melt(corrsRKHS_boot)
    colnames(df)[5] <- "corr"
    df$lambda <- parmTab$lam[match(df$i, 1:nrow(parmTab))]
    df$sigmax <- parmTab$sigmax[match(df$i, 1:nrow(parmTab))]
    df <- df[which(!is.na(df$corr)),]
    
    tab <- cast(df, lambda+sigmax+sigmay~., value="corr", fun.aggregate="median")
    colnames(tab)[ncol(tab)] <- "corr"
    
    #print("dim(tab)")
    #print(dim(tab))
    #print("head(tab)")
    #print(head(tab))
    #print("tab corr")
    #print(tab)
    #print("summary(tab corr)")
    #print(summary(tab))
    
    indxOpt <- which(tab$corr > quantile(tab$corr, 0.9))
    
    tabOpt <- tab[indxOpt,]
    summary(tabOpt)
    
    
    #dim(limsSigmaY)
    sigmaxRng <- range(tabOpt$sigmax)
    sigmayRng <- range(tabOpt$sigmay)
    lambdaRng <- range(tabOpt$lambda)
    #print(paste("sigmayRng:", paste(sigmayRng, collapse="-")))
    numY <- learner$hyperParams$data$optimizable$sigma.rbf.Y$length.out
    if(is.null(numY)) numY <-5
    sigmaySeq <- 10^seq(log(sigmayRng[1],10), log(sigmayRng[2],10), length.out=numY)
    #sigmaxSeq <- 10^seq(log(sigmaxRng[1],10), log(sigmaxRng[2],10), length.out=numX)
    #lambdaSeq <- 10^seq(log(lambdaRng[1],10), log(lambdaRng[2],10), length.out=numL)
  } 
  
  #if(!is.null(sigmaxSeq)) print(paste("log-sigmaxSeq: ", paste(round(log(sigmaxSeq,10),2),collapse=",")))
  #if(!is.null(sigmaySeq)) print(paste("log-sigmaySeq: ", paste(round(log(sigmaySeq,10),2),collapse=",")))
  #if(!is.null(lambdaSeq)) print(paste("log-lambdaSeq: ", paste(round(log(lambdaSeq,10),2),collapse=",")))
  
  
  
  #print("exits getFixedParams_rbf4")
  return(list(val=list(sigma.rbf.X=sigmaX, sigma.rbf.Y=sigmaY, lambda=1, densX=densX), seq=list(sigma.rbf.X=sigmaxSeq, sigma.rbf.Y=sigmaySeq, lambda=lambdaSeq)))
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

getFixedParams_rbf_ln <- function(learner, data){
  lrnr <- learner$hyperParams$data$non_optimizable$krrRefModel$type
  lrnr <- eval(parse(text=lrnr))
  dataAux <- constructData(x=data$x, y=data$y)
  
  lrnr <- setParams(learner=lrnr, trainData=dataAux, mc_cores=floor(detectCores()/2))
  lrnr <- lrnr$learn(lrnr)
  
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
  
  
  
  sigmaX <- 1/quantile(as.numeric(dist(unique(x))^2), 0.5)
  
  
  numX <- learner$hyperParams$data$optimizable$sigmaX$length.out
  if(is.null(numX)) numX <-10
  
  #rngX <- getRangeRbf(cbind(x,y), length.out=numX)
  rngX <- getRangeRbf(x, length.out=numX)
  
  #rng <- range(rngX, rngY)
  #rngX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numX)
  #rngY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=numY)
  
  
  
  return(list(val=list(sigmaX=sigmaX, krrRefModel=lrnr), seq=list(sigmaX=rngX)))
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
  
  sigma1 <- 1/quantile(dists2, 0.99)  
  sigma2 <- min(1/quantile(dists2, 0.2),1e10)
  
  res <- 10^seq(log(sigma1,10), log(sigma2,10), length.out=length.out)
  #res <- c(sigma1, sigma2)
  
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


optHP.CV <- function(learner, plot=FALSE, fac=1, mc_cores=1, numCoresFold=1, CVfunc="CV.parallel", predFunc="pred.CV"){
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
  indXX <- which(sapply(aux, function(el) el[3])=="X" & "kernelX" %in% nms)
  indYY <- which(sapply(aux, function(el) el[3])=="Y" & sum(c("kernelY", "featureY") %in% nms)>0)
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
      paramsListss <- list(c(paramsListO))
  }
  
   otherParams <- lapply(learner$hyperParams$data$non_optimizable, function(el) el$val)
  #params <- c(paramsList, otherParams) #, learner$hyperParams$non_data
  
  
  params <- do.call(constructParams, list(otherParams=otherParams, paramsLists=paramsListss))
  
  # not necessary coz CV is only done on data hyperparms??
  #params <- lapply(params, function(el) c(el, learner$hyperParams$non_data))
  #class(params) = "CVST.params"
  
  if(is.null(learner$hyperParams$data$grid)){
    print("calculates grid")
    grid <- do.call(CVfunc, list(learner, params, fac=fac, mc_cores=mc_cores, numCoresFold=numCoresFold, predFunc=predFunc))
  } else{
    print("takes previously existing grid")
    grid <- learner$hyperParams$data$grid
  }
    
  fold <- FALSE
  if("fold" %in% names(dimnames(grid))){
    fold <- TRUE
    gridFold  <- grid
    indxFoldDim <- which("fold"==names(dimnames(grid)))
    newDim <- dim(grid)[-indxFoldDim]
    
    dmnms <- dimnames(grid)[-indxFoldDim]
    grid <- sapply(learner$optimizeParams$losses, function(loss){
      # loss <- learner$optimizeParams$losses[[1]]
      #print(paste("loss: ", loss$func))
      res <- apply(grid[,,,loss$func], setdiff(names(dimnames(grid[,,,loss$func])), "fold"), loss$aggFunc) 
    }, simplify="array")
    dim(grid) <- newDim
    dimnames(grid) <- dmnms
    
  }
  
  
  # reshape back into one dimension per hyperparameter
  
  
  
  dimnms <- dimnames(grid)
  
  
  
  # obtain best hyperparameters
  
  mainLoss <- learner$optimizeParams$mainLoss
  if(is.null(mainLoss)) mainLoss <- dimnames(grid)$var[1]
  testTrain <- learner$optimizeParams$testTrain
  
  
  
  #testGrid <- grid[testTrain,,mainLoss]
  testGrid <- adrop(grid[testTrain,,,drop=F],drop=1)
  
  #minTestGrid <- min(testGrid, na.rm=T)
  #we add a tiny bit of noise to get exactly one minimum
  #testGrid <- testGrid+ rnorm(length(testGrid),mean=0, sd=max(minTestGrid,1e-10)*1e-10)
  ##minTestGrid <- min(testGrid, na.rm=T)
  ##optMinBool <- testGrid==minTestGrid
  
  # keep opt indices on grid to have easy access to chosen grid point later on
  
  #indxOptGrid <- which.min(testGrid)
  optLossFunc <- learner$optimizeParams$optLossFunc
  if(is.null(optLossFunc)){
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLoss,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
  }  
  indxOptGrid <- do.call(optLossFunc, list(grid=testGrid))
  
  
  
  opt <- params[[indxOptGrid]][intersect(nms, names(params[[indxOptGrid]]))]
  
  
  
  # check (only when lanbda and sigma are passed for CV)
  #grid["test",which.min(abs(as.numeric(dimnames(grid)$lambda)-opt$lambda)),which.min(abs(as.numeric(dimnames(grid)$sigma)-opt$sigma)),1]
  #min(grid["test",,,1])
  
  res <-  list(opt=opt, grid=grid, indxOptGrid=indxOptGrid)
  if(fold) res <- c(res, list(gridFold=gridFold))
  
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


CV.parallel.avg <- function(learner, params, fac=1, verbose=TRUE, mc_cores=1, numCoresFold=1, predFunc="pred.CV") {
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
      
      
      
      lossesTrain <- sapply(learnerAux$optimizeParams$losses, function(func){
        do.call(func$func, list(learner=learnerAux, pred=predTrain))
      })
      
      lossesTest <- sapply(learnerAux$optimizeParams$losses, function(func){
        do.call(func$func, list(learner=learnerAux, pred=predTest))
      })
      
      
      
      res <- cbind(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=numCoresFold, SIMPLIFY="array")
    return(res)
  }, params, mc.cores=mc_cores, SIMPLIFY="array")
  
  
  
  losses <- aperm(losses, c(2,3,4,1))
  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=names(params), var=names(learner$optimizeParams$losses))
  
  #print("exits CV parallel")
  return(losses)
}


CV.boot <- function(learner, params, fac=1, verbose=TRUE, mc_cores=1, numCoresFold=1, predFunc="pred.CV") {
  stopifnot(class(learner) == "emley.learner" && class(params) == "CVST.params")
  
  #print("enters CV parallel")
  
  print("length(params)")
  print(length(params))
  print("params")
  print(names(params))
  print("num cores")
  print(mc_cores)
  print("num coresFold")
  print(numCoresFold)
  
  initOnly <- learner$optimizeParams$initOnly
  numFolds <- learner$optimizeParams$numFolds
  trainData <- learner$hyperParams$trainData
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), names(params))
  
  n <- length(trainData$y)
  size <- ceiling(n / numFolds)
  
  # set bootstrap samples here coz it might be easier for reproducibility 
  # reasons (setting seed outside parallel)
  
  if(numFolds==1 || initOnly){
    boot.smpls <- matrix(1:n, n, numFolds)
    
  } else{
    set.seed(1234)
    boot.smpls <- sapply(1:numFolds, function(f) sample(1:n, size=n, replace=T), simplify="array")
    
    
  }
  
  set.seed(learner$optimizeParams$seed)  
  seedParms <- sample(100000, nParams)
  
  pm <- proc.time()
  #count <- 0
  losses <- mcmapply(FUN=function(p, seedParm){
    # p <- params[[1]]
    #count <<- count + 1
    #print(paste("count: ", count))
    #print(names(p[1]))
    set.seed(seedParm)
    seedFolds <- sample(100000, numFolds)
    
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- params[[1]]
      # print(paste("fold: ", f))
      validationIndex <- boot.smpls[,f]
      
      # in this case we test on the same trainData as we have no choice:
      # we only have latent factor for 
      curTrain <- getSubset(trainData, validationIndex)
      #curTest <- getSubset(trainData, validationIndex)
      # either mean squared error or mean classification error
      
      
      learnerAux <- learner
      learnerAux$hyperParams$trainData <- curTrain
      
      nmsPars <- names(learnerAux$hyperParams$data$optimizable)
      for(pr in nmsPars){
        # pr <- nmsPars[1]
        learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
      }
      
      learnerAux <- learnerAux$learn(learner=learnerAux, seed=seedFolds[f], inits=FALSE)
      curTrain <- learnerAux$hyperParams$trainData
      
      predTrain <- learnerAux$predict(learner=learnerAux, data=curTrain)
      #predTest <- learnerAux$predict(learnerAux, data=curTest)
      
      
      
      lossesTrain <- sapply(learnerAux$optimizeParams$losses, function(func){
        #print(func)
        do.call(func$func, list(learner=learnerAux, pred=predTrain))
      })
      
      
      #lossesTest <- sapply(learnerAux$optimizeParams$losses, function(func){
      #  do.call(func$func, list(learner=learnerAux, pred=predTest))
      #})
      
      
      res <- cbind(lossesTrain, lossesTrain)
      
      return(res)
    }, f=1:numFolds, mc.cores=numCoresFold, SIMPLIFY="array")
    return(res)
  }, p=params, seedParm=seedParms, mc.cores=mc_cores, SIMPLIFY="array")
  pm <- proc.time()-pm
  print(paste("CV boot took ", round(pm["elapsed"]/60, 1), " mins"))
  
  
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
     #f <- 1
    #print(paste("fold", f))
    last <- f*size
    first <- (f-1)*size + 1
    if(f==numFolds){
      first <- n-size+1
      last <- n
    }
    validationIndex <- seq(first, last)
    curTrain <- getSubset(data, subset=-validationIndex)
    curTest <- getSubset(data, validationIndex)
    # either mean squared error or mean classification error
    
    
    learner1 <- learner
    learner1$hyperParams$trainData <- curTrain
    
    learner1 <- learner1$learn(learner=learner1, forLoss=T)
    
    predTrain <- learner1$predict(learner=learner1, data=curTrain, forLoss=T)
    predTest <- learner1$predict(learner=learner1, data=curTest, forLoss=T)
    
    
    res <- list(train=predTrain, test=predTest)
    
    
    return(res)
  }, f=1:numFolds, mc.cores=numCores, SIMPLIFY=FALSE)
  
  auxTrain <- lapply(res, function(el) el$train)
  auxTest <- lapply(res, function(el) el$test)
  nmsAux <- names(auxTrain[[1]])
  #c("lx","gy_k","gyh_k")
  predTrain <- lapply(nmsAux, function(nm){
    # nm <- nmsAux[3]
    #print(paste("nm: ", nm))
    res <- do.call(rbind, lapply(auxTrain, function(el) el[[nm]]))
    return(res)
    })
  predTest <- lapply(nmsAux, function(nm) do.call(rbind, lapply(auxTest, function(el) el[[nm]])))
  names(predTrain) <- names(predTest) <- nmsAux
  
  # plot(predTest$x, predTest$gy); ord <- order(predTest$x)
  # for(i in 1:ncol(predTest$gyh)) lines(predTest$x[ord], predTest$gyh[ord,i], col=i)
  
  # plot(predTrain$x, predTrain$gy); ord <- order(predTrain$x)
  # for(i in 1:ncol(predTrain$gyh)) lines(predTrain$x[ord], predTrain$gyh[ord,i], col=i)
  
  return(list(train=predTrain, test=predTest))
}

CV.parallel <- function(learner, params, fac=1, verbose=TRUE, mc_cores=1, numCoresFold=1, predFunc="pred.CV") {
  stopifnot(class(learner) == "emley.learner" && class(params) == "CVST.params")
  
  #print("enters CV parallel")
  
  numFolds <- learner$optimizeParams$numFolds
  trainData <- learner$hyperParams$trainData
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), names(params))
  
  #mc_cores <- min(40, detectCores()-4)
  #mc_cores <- 1
  #numCoresFold <- mc_cores
  
  #count <<-0 
  pm <- proc.time()
  losses <- mcmapply(FUN=function(p){
    # p <- params[[1]]
    # print(p)
    # count <<- count + 1
    # print(paste("count: ", count, " out of ", length(params)))
    learnerAux <- learner
    nmsPars <- names(learnerAux$hyperParams$data$optimizable)
    for(pr in nmsPars){
      # pr <- nmsPars[3]
      #print(paste(pr, ": ", p[[match(pr, names(p))]]))
      learnerAux$hyperParams$data$optimizable[[pr]]$val <- p[[match(pr, names(p))]]
    }
    
      # I dont actually know why im calling learn before goint into 
      # folds, unless I need something to calculate losses below?
      # gonna comment it out for now... no!! i know, its so i can
      # have KCDC, KCMC etc as loss functions since they use
      # learned stuff on all data
      #print("learn")
      pm <- proc.time()
      learnerAux <- learnerAux$learn(learner=learnerAux, forLoss=F)
      proc.time() - pm # 84 secs for 1000 pts, 3+2=5 taus, 
      parsPred <- list(learner=learnerAux, data=trainData, numCores = numCoresFold)
      #print("pred.CV")
      pm <- proc.time()
      preds <- do.call(predFunc, parsPred)
      proc.time() - pm # 2.26 mins for 4 folds 1000 pts
      
      #print("names(preds$test)")
    #print(names(preds$test))  
      
    # x <- trainData$x; y <- trainData$y; ord <- order(x)
    # plot(x, y)
    # for(i in 1:ncol(preds$test$gyh)) lines(x[ord], preds$test$gyh[ord,i], col=i)
    
    #print("loss")
    #print("test")
    lossesTest <- mcmapply(function(func,nm){
      # func <- learnerAux$optimizeParams$losses[[1]]
      # print(nm)
      do.call(func$func, list(learner=learnerAux, pred=preds$test))
    }, func=learnerAux$optimizeParams$losses, nm=names(learnerAux$optimizeParams$losses),mc.cores=numCoresFold, SIMPLIFY=FALSE)
    names(lossesTest) <- sapply(learnerAux$optimizeParams$losses, function(func) func$func)
    lossesTest <- unlist(lossesTest)
    
    #print("train")
    lossesTrain <- mcmapply(function(func,nm){
      # func <- learnerAux$optimizeParams$losses[[1]]
      #print(nm)
      do.call(func$func, list(learner=learnerAux, pred=preds$train))
    }, func=learnerAux$optimizeParams$losses, nm=names(learnerAux$optimizeParams$losses),mc.cores=numCoresFold, SIMPLIFY=FALSE)
    names(lossesTrain) <- sapply(learnerAux$optimizeParams$losses, function(func) func$func)
    lossesTrain <- unlist(lossesTrain)
    
    if(length(lossesTrain)==length(lossesTest)){
      res <- cbind(lossesTrain, lossesTest)
    } else{
      res <- cbind(lossesTest)
    }
    #print(res)
    #print(proc.time()-pm)
    
    return(res)
  }, params, mc.cores=mc_cores, SIMPLIFY="array")
  proc.time() - pm #
  
  losses <- aperm(losses, c(2,3,1))
  if(dim(losses)[1] == 1) trTest <- "test" else trTest <- c("train","test")
  dimnames(losses) <- list(trainTest=trTest, params=names(params), var=dimnames(losses)[[3]])
  
  #print("exits CV parallel")
  return(losses)
}

pred.SSL.CV <- function(learner, data, numCores=1){
  
  numFolds <- learner$optimizeParams$numFolds
  n <- getN(data)
  nUnlab <- length(data$indxU)
  nLab <- n - nUnlab
  rho <- nLab / n
  
  dataLab <- getSubset(data, -data$indxU)
  
  classif <- TRUE#FALSE
  if(length(unique(dataLab$y))< nLab*0.25) classif <- TRUE
  
  res <- mcmapply(FUN=function(f){
    #f <- 1
    #print(paste("fold", f))
    
    # if classification data do stratified sampling by class
    # in any case sample proportion rho  as unobserved
    if(classif){
      
      uniqueClasses <- unique(dataLab$y)
      
      indxCls <- lapply(uniqueClasses, function(cls){
        indxCls <- which(dataLab$y==cls)
        nCls <- length(indxCls)
        indx <- sample(nCls, round((1-rho)*nCls))
        return(indxCls[indx])
      })
      
      validationIndex <- unlist(indxCls)
      
    } else{
      
      validationIndex <- sample(nLab, round((1-rho)*nLab))
    }
    
    dataLab$indxU <- validationIndex
    curTrain <- getSubset(dataLab, -validationIndex)
    curTest <- getSubset(dataLab, validationIndex)
    # either mean squared error or mean classification error
    
    
    learner1 <- learner
    learner1$hyperParams$trainData <- dataLab
    
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


myGpOptimise <- function(learner, mc_cores=1, plot){
  
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
  return(list(x_class=as.matrix(x), gy_class=gy, gyh_class=as.matrix(x)))
}



learn.krr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")

  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  N <- nrow(Kxs)
  
  lambda <- getHyperPar(learner, "lambda")*N
  
  alpha <- try(solve(Matrix(Kxs + diag(lambda, N))) %*% y)
  if(class(alpha)=="try-error") alpha <- matrix(rnorm(N), N, 1)
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.krr <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}
predict.krr_class <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x_class=as.matrix(x), gy_class=as.matrix(gy), gyh_class=as.matrix(pred)))
}

learn.wkrr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  N <- nrow(Kxs)
  ws <- getHyperPar(learner, "ws")
  ws <- 1/ws
  ws <- ws/(sum(ws))*N
  
  lambda <- getHyperPar(learner, "lambda")*N
  
  alpha <- solve(Matrix(Kxs + diag(lambda*ws, N))) %*% y
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.wkrr <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}


# disribution regression 
# IMPORTANT: its assumed that the bag of the ys is 1,.., num_bags and that
# the bag variable corresponds to the bag label of the xs
learn.krr_dr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  bags <- learner$hyperParams$trainData$bag
  y <- learner$hyperParams$trainData$y
  
  
  dfbags <- combn(unique(bags),2)
  dfbags <- cbind(dfbags, matrix(rep(unique(bags), rep(2, length(unique(bags)))), 2, length(unique(bags))))
  N <- length(unique(bags))
  
  #count <- 0
  Kxs_vals <- apply(dfbags, 2, function(col){
    #print(paste("count: ", count, " out of ", ncol(dfbags)))
    #count <<- count + 1
    indx1 <- which(bags==col[1])
    indx2 <- which(bags==col[2])
    Kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x[indx1,,drop=F],x[indx2,,drop=F], pars=parsXs)
    res <- mean(Kx)
    return(res)
  })
  Kxs <- matrix(NA, N, N)
  Kxs[t(dfbags)] <- Kxs_vals
  Kxs[t(dfbags[c(2,1),])] <- Kxs_vals
  
  
  lambda <- getHyperPar(learner, "lambda")*N
  
  alpha <- solve(Matrix(Kxs + diag(lambda, N))) %*% y
  learner$learnParams$alpha <- alpha
  
  return(learner)
}
predict.krr_dr <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  bagsTr <- trainData$bag
  bagsPred <- data$bag
  
  dfbags <- expand.grid(bagPred=unique(bagsPred),bagTr=unique(bagsTr))
  Ntr <- length(unique(bagsTr))
  Npr <- length(unique(bagsPred))
  
  #count <- 0
  kxs_vals <- apply(as.matrix(dfbags), 1, function(col){
    # col <- as.matrix(dfbags)[1,]
    #print(paste("count: ", count))
    #count <<- count + 1
    indx1 <- which(bagsTr==col["bagTr"])
    indx2 <- which(bagsPred==col["bagPred"])
    kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, data$x[indx2,,drop=F],trainData$x[indx1,,drop=F], pars=parsXs)
    res <- mean(kx)
    return(res)
  })
  kxs <- matrix(NA, Npr, Ntr)
  kxs[as.matrix(dfbags)] <- kxs_vals
  
  
  pred <- kxs %*% learner$learnParams$alpha
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}

# kernel SVM classification

learn.ksvm <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  N <- length(y)
  
  Kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  Kx <- as.kernelMatrix(Kx)
  lambda <- getHyperPar(learner, "lambda")*N
  nu <- getHyperPar(learner, "nu")
  
  svmFit <- try(ksvm(x=Kx, y=factor(y), C=lambda, nu=nu, type=c("C-svc","nu-svc")[2], prob.model=T))
  if(class(svmFit)=="try-error") svmFit <- NA
  learner$learnParams$model <- svmFit
  
  return(learner)
}
predict.ksvm <- function(learner, data, forLoss=F) {
  
  
  fit <- learner$learnParams$model
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  
  if(is.na(fit)){
    pred <- runif(nrow(data$x))
  } else{
    kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
    kx <- as.kernelMatrix(kx)
    pred <- predict(fit, kx, type="probabilities")[,2]
  }
  
  return(list(x_class=as.matrix(data$x), gy_class=as.matrix(data$y), gyh_class=as.matrix(pred)))
}

# distribution kernel SVM classification

learn.svm_dr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  bags <- learner$hyperParams$trainData$bag
  y <- learner$hyperParams$trainData$y
  
  
  dfbags <- combn(unique(bags),2)
  dfbags <- cbind(dfbags, matrix(rep(unique(bags), rep(2, length(unique(bags)))), 2, length(unique(bags))))
  N <- length(unique(bags))
  
  
  #count <- 0
  Kxs_vals <- apply(dfbags, 2, function(col){
    #print(paste("count: ", count, " out of ", ncol(dfbags)))
    #count <<- count + 1
    indx1 <- which(bags==col[1])
    indx2 <- which(bags==col[2])
    Kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x[indx1,,drop=F],x[indx2,,drop=F], pars=parsXs)
    res <- mean(Kx)
    return(res)
  })
  Kxs <- matrix(NA, N, N)
  Kxs[t(dfbags)] <- Kxs_vals
  Kxs[t(dfbags[c(2,1),])] <- Kxs_vals
  I <- diag(N)
  H <- I-matrix(1/N,N,N)
  Kxs <- H %*% Kxs %*% H
  Kxs <- as.kernelMatrix(Kxs)
  
  
  nu <- getHyperPar(learner, "nu")
  set.seed(123)
  svmFit <- try(ksvm(x=Kxs, y=factor(y),  nu=nu, type=c("C-svc","nu-svc")[2], prob.model=T))
  if(class(svmFit)=="try-error") svmFit <- NA
  learner$learnParams$model <- svmFit
  learner$learnParams$Kx <- Kxs
  #table(predict(svmFit, Kxs), y)
  
  return(learner)
}
predict.svm_dr <- function(learner, data, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  bagsTr <- trainData$bag
  bagsPred <- data$bag
  
  dfbags <- expand.grid(bagPred=unique(bagsPred),bagTr=unique(bagsTr))
  Ntr <- length(unique(bagsTr))
  Npr <- length(unique(bagsPred))
  
  
  fit <- learner$learnParams$model
  Kx <- learner$learnParams$Kx
  
  if(is.na(fit)){
    pred <- runif(length(unique(data$bag)))
  } else{
    #count <- 0
    kxs_vals <- apply(as.matrix(dfbags), 1, function(col){
      # col <- as.matrix(dfbags)[1,]
      #print(paste("count: ", count))
      #count <<- count + 1
      indx1 <- which(bagsTr==col["bagTr"])
      indx2 <- which(bagsPred==col["bagPred"])
      kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, data$x[indx2,,drop=F],trainData$x[indx1,,drop=F], pars=parsXs)
      res <- mean(kx)
      return(res)
    })
    kxs <- matrix(NA, Npr, Ntr)
    kxs[as.matrix(dfbags)] <- kxs_vals
    Ipr <- diag(Npr)
    Hpr <- Ipr-matrix(1/Npr,Npr,Npr)
    Itr <- diag(Ntr)
    Htr <- Itr-matrix(1/Ntr,Ntr,Ntr)
    kxs <- Hpr %*% kxs %*% Htr
    kxs <- as.kernelMatrix(kxs)
    #plot(c(Kx), c(kxs))
    pred <- predict(fit, kxs, type="probabilities")[,2]
    pred2 <- as.numeric(as.character(predict(fit, kxs)))
  }
  
  gy <- data$y
  x <- data$x
  return(list(x_class=as.matrix(x), gy_class=as.matrix(gy), gyh_prob=as.matrix(pred), gyh_class=as.matrix(pred2)))
}


# distribution kernel SVM classification - rbf based on 
# pre calculculated Euclidean distances

calcBagKernTr <- function(trainData){
  x <- trainData$x
  bags <- trainData$bag
  dfbags <- combn(unique(bags),2)
  dfbags <- cbind(dfbags, matrix(rep(unique(bags), rep(2, length(unique(bags)))), 2, length(unique(bags))))
  N <- length(unique(bags))
  
  #count <- 0
  Mx_vals <- mcmapply(function(i){
    #print(paste("count: ", count, " out of ", ncol(dfbags)))
    #count <<- count + 1
    indx1 <- which(bags==dfbags[1,i])
    indx2 <- which(bags==dfbags[2,i])
    M <- matNorm2(x[indx1,,drop=F],x[indx2,,drop=F])
    return(M)
  }, i=1:ncol(dfbags), mc.cores=detectCores()/2, SIMPLIFY=F)
  
  #A <- list(a=matrix(seq(4),2,2), b=matrix(seq(5,8),2,2), 
  #          c=matrix(seq(9,12),2,2),d=matrix(seq(13,16),2,2))
  #B <- matrix(A, 2,2)
  
  numRows <- rep(c(table(bags)), N)
  numCols <- rep(c(table(bags)), rep(N, N))
  
  M <- lapply(1:(N^2), function(i) matrix(NA, numRows[i], numCols[i]))
  M <- matrix(M, N, N)
  M[t(dfbags)] <- Mx_vals[1:length(Mx_vals)]
  M[t(dfbags[c(2,1),])] <- Mx_vals
  
  
  # chk
  #i <- 40
  #j <- 7
  #plot(M[i,j][[1]], matNorm2(x[which(bags==i),,drop=F], x[which(bags==j),,drop=F]))
  
  
  return(M)
}
calcBagKernPr <- function(trainData, data){
  xTr <- trainData$x
  xPred <- data$x
  
  bagsTr <- trainData$bag
  bagsPred <- data$bag
  
  dfbags <- expand.grid(bagPred=unique(bagsPred),bagTr=unique(bagsTr))
  Ntr <- length(unique(bagsTr))
  Npr <- length(unique(bagsPred))
  
  mc_cores <- detectCores()/2
  # mc_cores <- 1
  Mx_vals <- mcmapply(function(i){
    # i <- 8
    #print(paste("i: ", i, " out of ", nrow(dfbags)))
    dfbags[i,]
    indx1 <- which(bagsPred==dfbags[i,1])
    indx2 <- which(bagsTr==dfbags[i,2])
    M <- matNorm2(xPred[indx1,,drop=F],xTr[indx2,,drop=F])
    return(M)
  }, i=1:nrow(dfbags), mc.cores=mc_cores, SIMPLIFY=F)
  
  numRows <- rep(c(table(bagsPred)), Ntr)
  numCols <- rep(c(table(bagsTr)), rep(Npr, Ntr))
  
  M <- lapply(1:(Ntr*Npr), function(i){ 
    # i <- 10001
    #print(i)
    matrix(NA, numRows[i], numCols[i])
  })
  M <- matrix(M, Npr, Ntr)
  M[as.matrix(dfbags)] <- Mx_vals
  
  # chk
  #i <- 40
  #j <- 7
  #plot(M[i,j][[1]],matNorm2(xPred[which(bagsPred==i),,drop=F], xTr[which(bagsTr==j),,drop=F]))
  
  return(M)
}


learn.svm_dr_M <- function(learner, forLoss=F) {
  
  sigma <- getHyperPar(learner, "sigma")
  
  M <- learner$hyperParams$trainData$M
  y <- learner$hyperParams$trainData$y
  
  Kxs <- apply(M, c(1,2), function(el) mean(exp(-sigma*el[[1]])))
  
  N <- nrow(M)
  I <- diag(N)
  H <- I-matrix(1/N,N,N)
  Kxs <- H %*% Kxs %*% H
  Kxs <- as.kernelMatrix(Kxs)
  
  
  nu <- getHyperPar(learner, "nu")
  set.seed(123)
  svmFit <- try(ksvm(x=Kxs, y=factor(y),  nu=nu, type=c("C-svc","nu-svc")[2], prob.model=T))
  if(class(svmFit)=="try-error") svmFit <- NA
  learner$learnParams$model <- svmFit
  learner$learnParams$Kx <- Kxs
  #table(predict(svmFit, Kxs), y)
  
  return(learner)
}
predict.svm_dr_M <- function(learner, data, forLoss=F) {
  
  sigma <- getHyperPar(learner, "sigma")
  
  
  M <- data$M
  Ntr <- ncol(M)
  Npr <- nrow(M)
  
  fit <- learner$learnParams$model
  Kx <- learner$learnParams$Kx
  
  if(is.na(fit)){
    pred <- runif(length(unique(data$bag)))
  } else{
    #count <- 0
    kxs <- apply(M, c(1,2), function(el) mean(exp(-sigma*el[[1]])))
    
    Ipr <- diag(Npr)
    Hpr <- Ipr-matrix(1/Npr,Npr,Npr)
    Itr <- diag(Ntr)
    Htr <- Itr-matrix(1/Ntr,Ntr,Ntr)
    kxs <- Hpr %*% kxs %*% Htr
    kxs <- as.kernelMatrix(kxs)
    #plot(c(Kx), c(kxs))
    pred <- predict(fit, kxs, type="probabilities")[,2]
    pred2 <- as.numeric(as.character(predict(fit, kxs)))
  }
  
  gy <- data$y
  x <- data$x
  return(list(x_class=as.matrix(x), gy_class=as.matrix(gy), gyh_prob=as.matrix(pred), gyh_class=as.matrix(pred2)))
}


# based on RSSL:::adjacency_knn
adjacency_knn <- function (D, k = 6){
  #Ds <- as.matrix(dist(X, method = distance))
  neighbours <- apply(D, 1, function(x) sort(x, index.return = TRUE)$ix[2:(k + 1)]) %>% as.integer
  adj <- as.matrix(Matrix::sparseMatrix(i = rep(1:nrow(D), each = k), j = neighbours, x = 1, dims = c(nrow(D), nrow(D))))
  adj <- (adj | t(adj)) * 1
  return(adj)
}

# learner es una lista en la que estan los parametros que son:
#   trainData:               una lista con x (en R^n), y (en {0,1,NA}^n) y indxU
#                            los indices "unlabelled".. cuando y_i = NA i debe estar en indxU
#                            pero puede ser que si i en indxU y_i = 0/1 porque no queremos 
#                            usar ese label por que estamos validando por ej.
#
#  parametros de kernel :    tipo de kernel, hiperparametro de kernel
#                             se llama a la funcion makeKernel
#
# Gtype                 :   como debemos construir el grafo  y la matriz
#                           asociada de similutedes W que queremos
#                           usar para propagar labels 
#                                     
# lambda1               :  parametro de regularización sobre laplaciano
# 
# lambda2               : parametro de regularizacion de suavidad
#
# normL                 : si se debe normalizar la matriz laplaciana L o no
learn.lapSSLkrr1 <- function(learner, forLoss=F) {
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  indxU <- learner$hyperParams$trainData$indxU
  if(any(indxU<0)) stop("only positive indices in this case")
  
  y[y==0] <- -1
  y[y==1] <- 1
  y[indxU] <- 0
  
  n <- nrow(x)
  
  
  krnlX <- learner$makeKernel(learner, data1=learner$hyperParams$trainData, var="X", grad=F)
  K <- krnlX$K
  
  Gtype <- learner$hyperParams$data$non_optimizable$Gtype
  
  # Gtype$val = "K" - usar W = K el mismo kernel 
  #             de las soluciones f= K %*% alpha
  #
  #           = "ad_euc" - usar use euclidean distance of k (Gtype$pars$adjacency_k) nearest neighbors
  #
  #           = "isomap" - usar algoritmo isomap para construir el grafo
  #            
  
  if(Gtype$val=="K"){
    W <- K
  } else if(Gtype$val=="adj"){
    # 1- K only valid for rbf
    W <- do.call("adjacency_knn", list(D=1-K, k=Gtype$pars$adjacency_k))
  } else if(Gtype$val=="adj_euc"){
    
    W1 <- RSSL:::adjacency_knn(x, distance="euclidean", k=Gtype$pars$adjacency_k)
    W2 <- as.matrix(dist(x, method="euclidean"))
    #W2 <- (1-matrix(norml(W2),n,n)+0.5)*(max(W2)-min(W2))+min(W2)
    W <- W1*W2
    
  } else if(Gtype$val=="isomap"){
    # x <- as.matrix(rnorm(100)); K <- kern_rbf(x=x, sigma=1);
    k <- Gtype$pars$adjacency_k
    valve <- TRUE
    while(valve){
      # 1- K only valid for rbf
      W <- try(as.matrix(isomapdist(dist=1-K,  k=k)))
      if(class(W)=="try-error") k <- k + 1 else valve <- FALSE
    }
    W <- max(W) - W + min(W)
    #plot(K,W)
  } else{
    stop(paste("Gtype: ",Gtype, " not supported"))
  }
  
  lambda1 <- getHyperPar(learner, "lambda1")
  lambda2 <- getHyperPar(learner, "lambda2")
  
  # compute normalized Laplacian matrix
  d <- rowSums(W)
  dSqrtInv <- 1 / sqrt(d)
  D <- diag(d) # degree matrix
  L <- D - W # Laplacian matrix
  
  normL <- getHyperPar(learner, "normL")
  if(normL) L  <- t(dSqrtInv*t(dSqrtInv*L)) # equivalent to, but faster than, 
  #                                      (diag(dSqrtInv) %*% L) %*% diag(dSqrtInv)
  
  #idxU <- is.na(y) # unobserved label index
  u <- length(indxU) # unobserved count
  l <- length(y) - u # observed count
  
  
  I <- diag(n) # identity matrix
  
  # precompute some matrix multiplications
  JK <- K
  JK[indxU,] <- 0 # equivalent to, but faster than, J%*%W where J <- diag(!idxU)
  
  if(all(lambda1==0)) {
    LK <- 0 # no Laplacian regularization
  } else {
    LK <- L%*%K
  }
  
  # Equation (8) in Belkin et al. (2006)
  toBeInv <- JK + lambda2 * l * I + ((lambda1 * l) / (u + l)^2) * LK
  alpha <- solve(toBeInv, y)
    
  learner$learnParams$alpha <- alpha
  learner$learnParams$K <- K
  learner$learnParams$W <- W
  learner$learnParams$L <- L
  

  return(learner)
}
predict.lapSSLkrr1 <- function(learner, data, forLoss=F) {
  
  y <- data$y
  y[y==0] <- -1
  y[y==1] <- 1
  
  trainData <- learner$hyperParams$trainData
  krnlX <- learner$makeKernel(learner, data1=data, data2=trainData, var="X", grad=F)
  k <- krnlX$K
  
  
  alpha <- learner$learnParams$alpha  
 #  bias <- learner$learnParams$bias 
  pred <- k %*% alpha  #+ bias
  
  return(list(x=as.matrix(data$x), gy=as.matrix(y), gyh=pred))
  
}



# based on LaplacianSVM from RSSL
learn.lapSVM  <- function(learner, forLoss=F){
  
  trainData <- learner$hyperParams$trainData
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  indxU <- learner$hyperParams$trainData$indxU
  if(any(indxU<0)) stop("only positive indices in this case")
  
  
  n <- nrow(x)
  u <- length(indxU)
  l <- n - u
  m <- ncol(x)
  
  y[y==0] <- -1
  y[y==1] <- 1
  y[indxU] <- 0
  indxMixToSep <- c(setdiff(1:n,  indxU), indxU)
  indxSeptoMix <- order(indxMixToSep)
  y <- y[indxMixToSep]
  x <- x[indxMixToSep,,drop=F]
  indxU <- which(y==0)
  trainData$x <- x
  trainData$y <- y
  trainData$indxU <- indxU
  
  #plot(x, x2); abline(a=0,b=1, col="red"); summary(as.numeric((x-x2)^2));all(x==x2)
  #all(y[1:l]==y2)
  
  krnlX <- learner$makeKernel(learner, data1=trainData, var="X", grad=F)
  K <- krnlX$K
  
  #plot(K, K2); abline(a=0,b=1, col="red"); summary(as.numeric((K-K2)^2));all(K==K2)
  
  
  lambda1 <- getHyperPar(learner, "lambda1")
  lambda2 <- getHyperPar(learner, "lambda2")
  eps <- getHyperPar(learner, "eps") #1e-09
  
  
  Gtype <- learner$hyperParams$data$non_optimizable$Gtype
  
  if(Gtype$val=="K"){
    W <- K
  } else if(Gtype$val=="adj"){
    # 1- K only valid for rbf
    W <- do.call("adjacency_knn", list(D=1-K, k=Gtype$pars$adjacency_k))
    #plot(W, W2); abline(a=0,b=1, col="red"); summary(as.numeric((W-W2)^2));all(W==W2)
  } else if(Gtype$val=="adj_euc"){
    
    #W <- RSSL:::adjacency_knn(x, distance="euclidean", k=Gtype$pars$adjacency_k)
    W <- dist(x, method="euclidean")
    W <- 1-matrix(norml(W),n,n)+0.5
    
  } else if(Gtype$val=="isomap"){
    # x <- as.matrix(rnorm(100)); K <- kern_rbf(x=x, sigma=1);
    k <- Gtype$pars$adjacency_k
    valve <- TRUE
    while(valve){
      # 1- K only valid for rbf
      W <- try(as.matrix(isomapdist(dist=1-K,  k=k)))
      if(class(W)=="try-error") k <- k + 1 else valve <- FALSE
    }
    W <- max(W) - W + min(W)
    #plot(K,W)
  } else{
    stop(paste("Gtype: ",Gtype, " not supported"))
  }
  
  d <- rowSums(W)
  dSqrtInv <- 1 / sqrt(d)
  D <- diag(d)
  L <- D - W
    
  normL <- getHyperPar(learner, "normL")
  if(normL) L  <- t(dSqrtInv*t(dSqrtInv*L)) # equivalent to, but faster than, 
  #                                      (diag(dSqrtInv) %*% L) %*% diag(dSqrtInv)
    
  #dim(L); dim(L2); plot(L, L2); abline(a=0,b=1, col="red"); summary(as.numeric((L-L2)^2));all(L==L2)
  
  I <- diag(n) # identity matrix
  
  # precompute some matrix multiplications
  JK <- K
  JK <- JK[-indxU,] #<- 0 # equivalent to, but faster than, J%*%W where J <- diag(!idxU)
  
  if(all(lambda1==0)) {
    LK <- 0 # no Laplacian regularization
  } else {
    LK <- L%*%K
  }
  
  Y <- diag(y) # n x n
  Y <- Y[,-indxU] # n x l
  #J <- cbind(diag(l), matrix(0, l, u)) # l x n
  toBeInv <- 2 * lambda2 * I + 2 * (lambda1/((l + u)^2)) * LK # n x n
  Qprime <- solve(toBeInv, Y) # n x l
  Q <- Y %*% JK %*% Qprime
  
  #dim(Q); dim(Q2); plot(Q[1:l,], Q2); abline(a=0, b=1, col="red"); all(Q[1:l,]==Q2)
  
  Amat <- diag(l)
  Amat <- t(rbind(y[-indxU], Amat, -Amat))
  bvec <- c(rep(0, l + 1), rep(-1/l, l))
  
  Qaux <- Q[1:l,] + diag(l) * eps
  
  
  if(is.positive.definite(Qaux)){
  
    beta <- try(solve.QP(Qaux , rep(1, l), Amat, bvec, meq = 1)$solution)
    if (any(is.nan(beta)) | class(beta)=="try-error"){
      #stop("Quadratic Programming problem returned: NaN. Try different hyperparameters?")
      beta <- matrix(0, l, 1)
    }
     
    
  } else{
    beta <- matrix(0, l, 1)
  }  
  alpha <- (Qprime %*% beta)
  
  #plot(alpha, alpha2); abline(a=0, b=1, col="red"); all(alpha==alpha2)
  
  if(all(alpha==0)){
    bias <- 0
  } else {
  
    SVs <- (beta > 1e-08/l) & (beta < 1/l - 1e-08/l)
    C <- 1/(2 * lambda2 * l)
    SVplus <- (alpha > 0.1) & (alpha < C - 0.1)
    SVmin <- (alpha < -0.1) & (alpha > -C + 0.1)
  
  
    bias <- median(as.numeric(K[(1:l)[SVs], ] %*% alpha) - y[1:l][SVs])
    bias <- 0
  
    if (sum(SVplus) > 0 && sum(SVmin) > 0) {
      bias <- -median(c(K[SVmin, ] %*% alpha + 1, K[SVplus, ] %*% alpha - 1))
    } else {
        if (sum(SVplus) > 0) {
          bias <- -median(K[SVplus, ] %*% alpha - 1)
        }
        if (sum(SVmin) > 0) {
          bias <- -median(K[SVmin, ] %*% alpha + 1)
        }
    }
    bias <- -median(K[(1:l)[SVs], ] %*% alpha - y[1:l][SVs])
  }
  if(is.na(bias)) bias <- 0
  learner$learnParams$alpha <- alpha[indxSeptoMix,,drop=F]
  learner$learnParams$bias <- bias
  K <- K[indxSeptoMix,]; K <- K[,indxSeptoMix]
  learner$learnParams$K <- K
  
  return(learner)
}
# based on selectMethod("decisionvalues","SVM") from RSSL
predict.lapSVM <- function (learner, data, forLoss=F){
  
  y <- data$y
  y[y==0] <- -1
  y[y==1] <- 1
  
  
  trainData <- learner$hyperParams$trainData  
  krnlX <- learner$makeKernel(learner, data1=data, data2=trainData, var="X", grad=F)
  k <- krnlX$K
  
  
  alpha <- learner$learnParams$alpha  
  bias <- learner$learnParams$bias 
  pred <- k %*% alpha  + bias
    
  return(list(x=as.matrix(data$x), gy=as.matrix(y), gyh=pred))
}

learn.qhsic <- function (learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  parsXb <- getKernelPars2(learner, kernelName="kernelXb")
  
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
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
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

# kernel quantile regression
learn.kqr <- function(learner, forLoss=F) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  N <- nrow(Kxs)
  lambda <- getHyperPar(learner, "lambda")*N
  taus <- getHyperPar(learner, "taus")  
  
  #tau <- 0.5
  #pm <- proc.time()
  #qrm <- kqr(x=x, y=matrix(y), tau = tau , kernel = "rbfdot", kpar= parsXs, C=lambda, scaled=FALSE)
  #proc.time() - pm 
  # 20 secs 1000 pts, tau <- 0.003, reduced=FALSE
  
  pm <- proc.time()
  qrms <- lapply(c(3/N, taus, 1-3/N), function(tau) try(kqr(x, y, tau = tau , kernel = "rbfdot", kpar= parsXs, C=lambda, scaled=FALSE)))
  proc.time() - pm #76 seconds
  
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
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  
  pred <- t(t(kxs %*% learner$learnParams$alpha)- learner$learnParams$b)
  
  
  gy <- data$y
  x <- data$x

  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=pred))
}


# copula (2d) quantile regression
learn.cqr <- function(learner, forLoss=F) {
  
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  n <- nrow(x)  
  
  u <- apply(cbind(x,y), 2, function(x) rank(x, ties.method = "random")/(n + 1))
  cop <- bicop(data = u,
               family_set = "tll",
               nonpar_method = "constant")
  
  
  learner$learnParams$cop <- cop
  

  return(learner)
}
predict.cqr <- function(learner, data, forLoss=F) {
  
  
  trainData <- learner$hyperParams$trainData
  n <- nrow(data$x)
  u <- rank(data$x, ties.method = "random")/(n + 1)
  
  N <- nrow(trainData$x)
  taus <- getHyperPar(learner, "taus")
  taus <- c(3/N, taus, 1-3/N)
  
  cop <- learner$learnParams$cop
  
  pred <- sapply(taus, function(tau) predict(object = cop, newdata = cbind(u, tau), what = "hinv1"))
  
  cond_q <- t(apply(pred, 1, function(row) quantile(trainData$y, row)))
  
  gy <- data$y
  x <- data$x
  
  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=cond_q))
}

# quantile forest regression
learn.fqr <- function(learner, forLoss=F){
  
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  
  nodesize <- getHyperPar(learner, "nodesize")
  sampsize <- getHyperPar(learner, "sampsize")
  
  model <-  quantregForest(x=as.matrix(x), y=as.matrix(y), nodesize=nodesize,sampsize=sampsize)
  
  
  learner$learnParams$model <- model
  
  
  return(learner)
}
predict.fqr <- function(learner, data, forLoss=F) {
  
  
  trainData <- learner$hyperParams$trainData
  
  N <- length(trainData$y)
  taus <- getHyperPar(learner, "taus")
  taus <- c(3/N, taus, 1-3/N)
  
  model <- learner$learnParams$model
  
  pred <- sapply(taus, function(tau) predict(model,  as.matrix(data$x), what=tau))
  
  
  
  gy <- data$y
  x <- data$x
  
  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=pred))
}

# quantile forest regression
learn.nnqr <- function(learner, forLoss=F) {
  
  
  x <- learner$hyperParams$trainData$x
  y <- learner$hyperParams$trainData$y
  
  n_hidden <- getHyperPar(learner, "n_hidden")
  n_hidden2 <- getHyperPar(learner, "n_hidden2")
  n_trials <- getHyperPar(learner, "n_trials")
  
  N <- nrow(x)
  taus <- getHyperPar(learner, "taus")
  taus <- c(3/N, taus, 1-3/N)
  
  iter_max <- learner$optimizeParams$iter_max
  print(paste("iter_max: ", iter_max))
  
  pm <- 0
  count <- 0
  stay <- TRUE
  while(stay){
    count <- count + 1
    print(paste("attempt no: ", count))
    pm <- proc.time()
    model <-  mcqrnn.fit(x=as.matrix(x), y=as.matrix(y), tau=taus,
                         n.hidden=n_hidden, n.hidden2=n_hidden2, n.trials=n_trials,
                         iter.max=iter_max)
    
    pred <- mcqrnn.predict(as.matrix(rnorm(10)), model)
    crit <- pred[,2:ncol(pred)]-pred[,1:(ncol(pred)-1)]
    print(crit[1:4,1:4])
    crit <- sum(crit)
    print(paste("crit: ",crit))
    print("pred")
    print(pred[1:4, 1:4])
    if(crit >1e-6) stay <- FALSE
    pm <- (proc.time()-pm)[3]
  }
  learner$learnParams$model <- model
  
  
  return(learner)
}
predict.nnqr <- function(learner, data, forLoss=F) {
  
  
  trainData <- learner$hyperParams$trainData
  
  
  model <- learner$learnParams$model
  
  
  pred <- mcqrnn.predict(as.matrix(data$x), model)
  
  
  
  gy <- data$y
  x <- data$x
  
  
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=pred))
}

# marginal quantile learner
learn.qr_marg <- function(learner, forLoss=F){
  y <- learner$hyperParams$trainData$y
  
  N <- length(y)
  taus <- getHyperPar(learner, "taus")
  taus <- c(3/N, taus, 1-3/N)
  qs <- quantile(y, probs=taus)
  learner$learnParams$qs <- qs
  return(learner)
}
predict.qr_marg <- function(learner, data, forLoss=F){
  gy <- data$y
  x <- data$x
  qs <- learner$learnParams$qs
  pred <- matrix(qs, length(gy), length(qs), byrow=T)
  
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
  # abline(v=xs[c(27,50)], col=c("red","blue"))
  
  
  res <- sapply(1:length(ys), function(i){
    # i <- 50
    #print(paste("i: ", i))
    
    # CDF estimation
    # get tau quantiles and cdf pairs (q,p) from our quantile reg model
    qs_i <- preds[i, ]
    ps_i <- c(0,taus,1)
    CDFi <- estimateCDF(qs=qs_i, ps=ps_i)
    qs_i <- CDFi$qs
    CDFi <- CDFi$CDF
    
    
    
    # we are interested in evaluating the cdf and pdf py(.|x_i) evaluating at y_i so we get set pt <- y_i
    y_i <- ys[i]
    grpi <- findInterval(y_i, qs_i)
    # we evaluate py(y_i|x_i) and Py(y_i|x_i) (pdf and cdf)
    lik <- CDFi(y_i, deriv=1)
    cdf <- CDFi(y_i, deriv=0)
    # abline(h=cdf, col="blue")
    # abline(v=y_i, col="blue")
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
  
    return(c(lik, cdf, grpi))
  }, simplify="array")
  res <- t(res)
  colnames(res) <- c("lik","cdf","grp")
  indxNA <- which(is.na(res[,"lik"]))
  res[indxNA,"lik"] <- 0#1e-6
  res[,"lik"] <- -log(res[,"lik"])
  colnames(res) <- c("pred","resid","grp")
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
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  parsXb <- getKernelPars2(learner, kernelName="kernelXb")
  parsRg <- getKernelPars2(learner, kernelName="kernelRg")
  
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
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=learner$hyperParams$trainData$x, pars=parsXs)
  

  pred <- kxs %*% learner$learnParams$alpha
  pred <- pred - mean(pred) + learner$learnParams$avgy
  gy <- data$y
  x <- data$x
  return(list(x=as.matrix(x), gy=as.matrix(gy), gyh=as.matrix(pred)))
}

# loss gand gradient functions for latent noise krr (lnkrr)

loss_latno_fhat <- function(fhat, y, alpha, lambda){
  res <- 0
  res <- res - 2*t(y)%*%fhat 
  res <- res + t(fhat) %*% fhat 
  res <- res + lambda*t(alpha)%*%fhat 
  
  return(res)
}
grad_latno_fhat <- function(fhat, y, alpha, lambda){
  res <- 0
  res <- res - 2*y 
  res <- res + 2* fhat 
  res <- res + lambda*alpha
  
  return(res)
}
jac_fhat_res_alphaOpt <- function(r, Kxa, y, lambda, sigma_r){
  n <- nrow(r)
  q <- ncol(r)
  I <- diag(n)
  B <- solve( Kxa + lambda*I)
  
  alpha <- B %*% y 
  
  
  df_dr1 <- sapply(1:n, function(j){
    # j <- 3
    aux <- t(Kxa[j,]*t(r[j,]- t(r))) 
    res <- -t(aux%*%alpha)
    df_drj <- t(sapply(1:n, function(i){
      # i <- 2; j <- 3
      if(i != j){
        res <- aux[,i]*alpha[i]
      } 
      return(res)
    }))
    #if(q==1) df_drj <- t(df_drj)
    return(df_drj)
  } ,simplify="array")
  df_dr1 <- df_dr1*sigma_r*2
  if(q==1) df_dr1 <- aperm(df_dr1, c(3,1,2))
  
  Kxa_prime <- sapply(1:n, function(j) 2*sigma_r*Kxa[,j]*t( t(r) - r[j,]), simplify="array")
  
  I <- array(0,c(n,q,n))
  
  B_Kxap <- tensor(B, Kxa_prime, 2, 1)
  Bprime <- sapply(1:n, function(j){
    #j <- 3
    res2 <- adrop(B_Kxap[,,j,drop=F],3) %o% B[,j]
    #res3 <-  B[,j] %o% t(adrop(B_Kxp[,,j,drop=F],3)) 
    res <- -(res2+aperm(res2, c(3,2,1)))
    
    return(res)
  }, simplify="array")
  
  
  dalpha_dr <- tensor(Bprime, y[,1], 3, 1)
  df_dr2 <- tensor(Kxa, dalpha_dr,2,1)
  
  df_dr <- 0
  df_dr <- df_dr + df_dr1
  df_dr <- df_dr + df_dr2
  
  dalpha_dr <- aperm(dalpha_dr, c(3,1,2))
  
  df_dr <- aperm(df_dr, c(3,1,2))
  return(list(df_dr=df_dr, dalpha_dr=dalpha_dr))
}
opt_latno_alpha <- function(Kxa, y, lambda){
  n <- nrow(Kxa)
  I <- diag(n)
  res <- solve( Kxa + lambda*I) %*% y
  return(res)
}
loss_latno_res_alphaOpt <- function(r, y, Kx, lambda, gamma, eta, L, ps, indxP){
  #print("enters loss_latno_res")
  n <- nrow(Kx)
  q <- length(r)/n
  r <- matrix(r, n, q)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  Kxc <- H%*%Kx%*%H
  
  sigma_rs  <- 1/quantile(as.numeric(dist(unique(r))^2), probs=ps)
  matNormsR <- matNorm2(r)
  Krs <- sapply(sigma_rs, function(sigma_r) exp(-sigma_r*matNormsR), simplify="array")
  dim(Krs) <- c(n, n, length(ps))
  sigma_r <- sigma_rs[indxP]
  Kr <- Krs[,,indxP]
  if(length(sigma_rs)==1){
    Kr2 <- Krs[,,1] 
  } else{
    Kr2 <- apply(Krs, c(1,2), sum)
  }
  #Kr <- kern_rbf(r, sigma=sigma_r)
  Kxa <- Kx*Kr
  alpha <- opt_latno_alpha(Kxa, y, lambda)
  fhat <- Kxa%*%alpha
  res <- 0
  res <- res + loss_latno_fhat(fhat, y, alpha, lambda) 
  res <- res + gamma*sum(diag(Kxc%*%Kr2))
  res <- res + eta*sum(diag(t(r)%*%r))
  #res <- res + L*( sum(diag(t(r)%*%r))/(n-1) - eta)^2 # equality restriction
  return(res)
}
grad_latno_res_alphaOpt <- function(r, y, Kx, lambda, gamma, eta, L, ps, indxP){
  n <- nrow(Kx)
  q <- length(r)/n
  r <- matrix(r, n, q)
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  Kxc <- H%*%Kx%*%H
  #sigma_r  <- 1/median(as.numeric(dist(unique(r))^2))
  #Kr <- kern_rbf(r, sigma=sigma_r)
  sigma_rs  <- 1/quantile(as.numeric(dist(unique(r))^2), probs=ps)
  matNormsR <- matNorm2(r)
  Krs <- sapply(sigma_rs, function(sigma_r) exp(-sigma_r*matNormsR), simplify="array")
  dim(Krs) <- c(n, n, length(ps))
  sigma_r <- sigma_rs[indxP]
  Kr <- Krs[,,indxP]
  if(length(sigma_rs)==1){
    Kr2 <- Krs[,,1]
  } else{
    Kr2 <- apply(Krs, c(1,2), sum)
  }
  Kxa <- Kx*Kr
  alpha <- opt_latno_alpha(Kxa, y, lambda)
  fhat <- Kxa%*%alpha
  
  grad_fhat <- grad_latno_fhat(fhat, y, alpha, lambda)
  
  J_fhat_res <- jac_fhat_res_alphaOpt(r, Kxa, y, lambda, sigma_r)
  
  dalpha_dr <- J_fhat_res$dalpha_dr
  
  J_fhat_res <- J_fhat_res$df_dr
  dL_dr <- tensor(grad_fhat[,1], J_fhat_res, 1, 2)
  
  dL_dr <- dL_dr + lambda*tensor(fhat[,1], dalpha_dr, 1, 2)
  
  # dL_dr2 <- sapply(1:n, function(j){ 
  #   # j <- 1
  #   df_drj <- t(sapply(1:n, function(i){
  #     # i <- 2; j <- 3
  #     if(i != j){
  #       res <- 2*sigma_r*Kx[i,j]*t( r[i,] - r[j,])*alpha[j]
  #     } else{
  #       res <- -2*sigma_r*t(t(Kx[i,]*t(r[j,]- t(r))) %*%alpha)
  #     }
  #     return(res)
  #   }))
  #   
  #   res <- t(grad_fhat) %*% df_drj
  #   return(res)
  # })
  # dL_dr2 <- t(dL_dr2)
  
  #plot(dL_dr, dL_dr2); abline(a=0, b=1, col="red")
  
  dHSIC_dr <- sapply(1:n, function(j){
    # j <- 1
    res <- apply(Kxc[,j]*Kr2[,j]*t(t(r)-r[j,]), 2, sum)
  })
  
  # don't know where the extra 2 factor comes from but
  # seems to match with numeric gradient:
  # it comes from careful analysis of double sum:
  # need three indices and careful when indices are same or different
  dHSIC_dr <- matrix(dHSIC_dr, q, n)
  dHSIC_dr <- 4*sigma_r*gamma*t(dHSIC_dr)
  
  res <- 0
  res <- res + dL_dr
  res <- res + dHSIC_dr
  res <- res + 2*eta*r
  #res <- res + 4*L*(sum(diag(t(r) %*% r))/(n-1) - eta)*r # equality restriction
  return(res)
}


learn.lnkrr <- function(learner, forLoss=F, mc_cores=1, seed=NULL, inits=TRUE) {
  
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  
  trainData <- learner$hyperParams$trainData
  y <- trainData$y
  x <- trainData$x
  
  Kx <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x, pars=parsXs)
  n <- nrow(Kx)
  q <- getHyperPar(learner, "q")
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  Kxc <- H%*%Kx%*%H
  
  lambda <- getHyperPar(learner, "lambda")*n
  eta <- getHyperPar(learner, "eta")
  gamma <- getHyperPar(learner, "gamma")
  sigLevel <- getHyperPar(learner, "sigLevel")
  mc_cores2 <- min(detectCores(), learner$optimizeParams$mc_cores)
  numCoresFold2 <- min(detectCores(), learner$optimizeParams$numCoresFold)
  
  if(is.null(gamma)){
    gammas1N <- getHyperPar(learner,"gammas1N")
    gammas2N <- getHyperPar(learner,"gammas2N")
    logGammaRange <- getHyperPar(learner,"logGammaRange")
    
    # copy learner, add gamma to optimizable and change
    # aggFunc to mean as there will only be one "bootstrap"
    gammas <- c(0,10^seq(logGammaRange[1],logGammaRange[2],length.out=gammas1N))
    learnerAux <- learner
    learnerAux$hyperParams$data$non_optimizable <- c(learnerAux$hyperParams$data$non_optimizable, learnerAux$hyperParams$data$optimizable) 
    learnerAux$hyperParams$data$optimizable <-list(gamma=list(val=NULL, seq=gammas))
    #learnerAux$optimizeParams$losses <- list(hsicTestRX=list(func="hsicTestRX", aggFunc="mean"))
    #learnerAux$optimizeParams$optLossFunc <- function(grid) which.max(grid[,"hsicTestRX"])
    learnerAux$optimizeParams$losses <- list(hsicTestRX=list(func="hsicTestRX", aggFunc="ksTest"),
                                             hsicTestResX=list(func="hsicTestResX", aggFunc="ksTest"),
                                             hsicTestResR=list(func="hsicTestResR", aggFunc="ksTest"),
                                             hsicTestResZ=list(func="hsicTestResZ", aggFunc="ksTest"),
                                             rmse=list(func="rmse", aggFunc="mean"))
    learnerAux$optimizeParams$optLossFunc <- function(grid) which.max(grid[,"hsicTestRX"]*grid[,"hsicTestResX"])
    learnerAux$optimizeParams$numFolds <- learnerAux$optimizeParams$numFoldsInt
    learnerAux$optimizeParams$seed <- learner$optimizeParams$seed+1
    
    
    learnerAux <- setParams(learner=learnerAux, trainData, mc_cores=mc_cores2, numCoresFold=numCoresFold2)
    #pvals <- learnerAux$hyperParams$data$grid["test",,"hsicTestRX"]
    pvalsRX <- learnerAux$hyperParams$data$grid["test",,"hsicTestRX"]
    pvalsResX <- learnerAux$hyperParams$data$grid["test",,"hsicTestResX"]
    pvalsResR <- learnerAux$hyperParams$data$grid["test",,"hsicTestResR"]
    pvalsResZ <- learnerAux$hyperParams$data$grid["test",,"hsicTestResZ"]
    pvals <- pvalsRX*pvalsResZ
    err <- learnerAux$hyperParams$data$grid["test",,"rmse"]
    
    gammasAux <- gammas
    gammasAux[which(gammasAux==0)] <- 10^(log(min(gammas[which(gammas!=0)]),10)-1)
    #plot(log(gammasAux,10), pvalsRX, type="b")
    #plot(log(gammasAux,10), pvalsResX, type="b")
    #plot(log(gammasAux,10), pvalsResR, type="b")
    #plot(log(gammasAux,10), pvalsResZ, type="b")
    #plot(log(gammasAux,10), pvals, type="b")
    #plot(log(gammasAux,10), err, type="b")
    
    indxMax <- which.max(pvals)
    if(indxMax==1){
      indxIni <- 1
      indxFin <- 2
    } else if(indxMax==length(gammas)){
      
      indxIni <- length(gammas)-1
      indxFin <- length(gammas)
      
    } else{
      indxIni <- indxMax-1
      indxFin <- indxMax+1
    }
    gammas <- 10^seq(log(gammasAux[indxIni],10),log(gammasAux[indxFin],10), length.out=gammas2N)
    #gammas <- c(0, gammas)
    
    
    learnerAux <- learner
    learnerAux$hyperParams$data$non_optimizable <- c(learnerAux$hyperParams$data$non_optimizable, learnerAux$hyperParams$data$optimizable) 
    learnerAux$hyperParams$data$optimizable <-list(gamma=list(val=NULL, seq=gammas))
    learnerAux$optimizeParams$losses <- list(hsicTestRX=list(func="hsicTestRX", aggFunc="ksTest"),
                                             hsicTestResX=list(func="hsicTestResX", aggFunc="ksTest"),
                                             hsicTestResR=list(func="hsicTestResR", aggFunc="ksTest"),
                                             hsicTestResZ=list(func="hsicTestResZ", aggFunc="ksTest"),
                                             rmse=list(func="rmse", aggFunc="mean"))
    learnerAux$optimizeParams$optLossFunc <- function(grid) which.max(grid[,"hsicTestResX"]*grid[,"hsicTestRX"])
    learnerAux$optimizeParams$numFolds <- learnerAux$optimizeParams$numFoldsInt
    learnerAux$optimizeParams$seed <- learner$optimizeParams$seed+2
    
    learnerAux <- setParams(learner=learnerAux, trainData, mc_cores=mc_cores2, numCoresFold=numCoresFold2)
    pvalsRX <- learnerAux$hyperParams$data$grid["test",,"hsicTestRX"]
    pvalsResX <- learnerAux$hyperParams$data$grid["test",,"hsicTestResX"]
    pvalsResR <- learnerAux$hyperParams$data$grid["test",,"hsicTestResR"]
    pvalsResZ <- learnerAux$hyperParams$data$grid["test",,"hsicTestResZ"]
    pvals <- pvalsRX*pvalsResZ
    err <- learnerAux$hyperParams$data$grid["test",,"rmse"]
    
    #plot(log(gammas,10), pvalsRX, type="b")
    #plot(log(gammas,10), pvalsResX, type="b")
    #plot(log(gammas,10), pvalsResR, type="b")
    #plot(log(gammas,10), pvalsResZ, type="b")
    #plot(log(gammas,10), pvals, type="b")
    #plot(log(gammas,10), err, type="b")
    indxOpt <- which.max(pvals)
    gamma <- gammas[indxOpt]
    ksPvalRX <- pvalsRX[indxOpt]
    ksPvalResX <- pvalsResX[indxOpt]
    ksPvalResR <- pvalsResR[indxOpt]
    ksPvalResZ <- pvalsResZ[indxOpt]
    rmseOpt <- err[indxOpt]
    learner$learnParams$gamma <- gamma
    learner$learnParams$ksPvalHsicRX <- ksPvalRX
    learner$learnParams$ksPvalHsicResX <- ksPvalResX
    learner$learnParams$ksPvalHsicResR <- ksPvalResR
    learner$learnParams$ksPvalHsicResZ <- ksPvalResZ
    learner$learnParams$rmseOpt <- rmseOpt
  } 
  
  #gamma <- gamma*(n^2)
  
  ps <- getHyperPar(learner, "ps")
  indxP <- getHyperPar(learner, "indxP")
  
  
  if(inits){
    learnerAux <- learner
    learnerAux$hyperParams$data$non_optimizable <- c(learnerAux$hyperParams$data$non_optimizable, learnerAux$hyperParams$data$optimizable) 
    learnerAux$hyperParams$data$optimizable <-list(gamma=list(val=NULL, seq=gamma))
    learnerAux$optimizeParams$losses <- list(hsicTestRX=list(func="hsicTestRX", aggFunc="mean"),
                                             hsicTestResX=list(func="hsicTestResX", aggFunc="mean"),
                                             hsicTestResR=list(func="hsicTestResR", aggFunc="mean"),
                                             hsicTestResZ=list(func="hsicTestResZ", aggFunc="mean"),
                                             rmse=list(func="rmse", aggFunc="mean"))
    learnerAux$optimizeParams$optLossFunc <- function(grid) which.max(grid[,"hsicTestRX"]*grid[,"hsicTestResZ"])
    learnerAux$optimizeParams$initOnly <- TRUE
    learnerAux$optimizeParams$numFolds <- learnerAux$optimizeParams$numFoldsInt
    learnerAux$optimizeParams$seed <- learner$optimizeParams$seed+3
    learnerAux <- setParams(learner=learnerAux, trainData, mc_cores=mc_cores2, numCoresFold=numCoresFold2)
    pvalsRX <- learnerAux$hyperParams$data$gridFold["test",,,"hsicTestRX"]
    pvalsResX <- learnerAux$hyperParams$data$gridFold["test",,,"hsicTestResX"]
    pvalsResR <- learnerAux$hyperParams$data$gridFold["test",,,"hsicTestResR"]
    pvalsResZ <- learnerAux$hyperParams$data$gridFold["test",,,"hsicTestResZ"]
    pvals <- pvalsRX*pvalsResZ
    err <- learnerAux$hyperParams$data$gridFold["test",,,"rmse"]
    
    numFolds <- learnerAux$optimizeParams$numFoldsInt
    #plot(1:numFolds, pvalsRX, type="p")
    #plot(1:numFolds, pvalsResX, type="p")
    #plot(1:numFolds, pvalsResR, type="p")
    #plot(1:numFolds, pvalsResZ, type="p")
    #plot(1:numFolds, pvals, type="p")
    #plot(1:numFolds, err, type="p")
    
    indxOpt <- which.max(pvals)
    # pvals[indxOpt]; err[indxOpt]
    set.seed(learnerAux$optimizeParams$seed)
    seedParms <- sample(100000,1)
    set.seed(seedParms)
    seedFold <- sample(100000, learnerAux$optimizeParams$numFoldsInt)
    seed <- seedFold[indxOpt]
  }
  
  set.seed(seed)
  r0 <- matrix(rnorm(n*q, sd=1), n, q)
  r0 <- apply(r0, 2, stdrize)
  # sd(r0); t(r0) %*% r0
  rhat <- lbfgs(loss_latno_res_alphaOpt, grad_latno_res_alphaOpt, vars=r0, 
                y=matrix(y,n,1),  Kx=Kx, lambda=lambda, gamma=gamma, eta=eta, L=L, ps=ps, indxP=indxP,
                max_iterations=100, invisible=1)
  rhat <- matrix(rhat$par, n, q)
  
  #print(paste("summary rhat: ", summary(rhat)))
  trainData <- constructData(x=x, y=y, r=rhat, rt=trainData$rt)
  learner$hyperParams$trainData <- trainData  
  
  
  
  sigma_r  <- 1/median(as.numeric(dist(unique(rhat))^2))
  Kr <- kern_rbf(rhat, sigma=sigma_r)  
  Kxa <- Kx*Kr
  alpha <- try(opt_latno_alpha(Kxa, y, lambda))
  
  if(class(alpha)=="try-error"){
    alpha <- matrix(0, n, 1)
  }
  
  learner$learnParams$alpha <- alpha
  learner$learnParams$sigma_r <- sigma_r
  
  # check we get the same loss as for optimal initialization
  #pred <- learner$predict(learner, data=learner$hyperParams$trainData)
  #pvals2 <- hsicTestResZ(learner, pred)*hsicTestRX(learner, pred)
  #err2 <- rmse(learner, pred)
  # pvals2; pvals[indxOpt]; err2; err[indxOpt]
  
  return(learner)
}
predict.lnkrr <- function(learner, data, forLoss=F) {
  parsXs <- getKernelPars2(learner, kernelName="kernelXs")
  sigma_r <- learner$learnParams$sigma_r
  trainData <- learner$hyperParams$trainData
  
  kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, x=data$x, y=trainData$x, pars=parsXs)
  krs <- kern_rbf(trainData$r, data$r, sigma=sigma_r)
  kxas <- kxs*krs
  pred <- kxas %*% learner$learnParams$alpha
  gy <- data$y
  res <- list(x=as.matrix(data$x), r=as.matrix(data$r), gy=as.matrix(gy), gyh=as.matrix(pred))
  return(res)
}


learn.gptk <- function(learner, forLoss=F) return(learner)
predict.gptk <- function(learner, data, forLoss=F){
  aux <- gpPosteriorMeanVar(learner$hyperParams$data$optimizable$model$val, X=data$x)
  pred <- list(x=data$x, gy=data$y, gyh=aux)
  return(pred)
}

learn.logReg <- function(learner, forLoss=F) {
  
  phix <- learner$makeFeature(learner, data=learner$hyperParams$trainData, var="X")
  
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
  
  #print("sigma")
  #print(sigma0)
  #print("K summary")
  #print(summary(as.numeric(K)))
  #print(K[1:5,1:5])
  
  cv.fit <- try(cv.glmnet(x=K, y=y, family="binomial", alpha=0, lambda=lambdas, nfolds=numFolds))
  if(class(cv.fit)=="try-error"){
    K <- K + matrix(rnorm(nrow(K)*ncol(K), mean=0, sd=1e-9), nrow(K), ncol(K))
    cv.fit <- try(cv.glmnet(x=K, y=y, family="binomial", alpha=0, lambda=lambdas, nfolds=numFolds))
  }
  
  #plot(log(cv.fit$lambda,10), cv.fit$cvm)
  #plot(cv.fit)
  
  
  indx <- which.min(cv.fit$cvm)
  #fit <- cv.fit$glmnet.fit
  lambda <- cv.fit$lambda[indx]
  
  #print("lambda")
  #print(lambda)
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
  trainClassData <- constructData(x=x_tr[smpl,,drop=F], y=y_tr[smpl])
  
  # make test classification data
  x_te <- as.matrix(rbind(x[indx.te,,drop=F], matrix(runif(n_te*p), n_te, p)))
  y_te <- c(rep(1, n_te), rep(0, n_te)) 
  set.seed(123)
  smpl <- sample(nrow(x_te))
  testClassData <- constructData(x=x_te[smpl, ,drop=F], y=y_te[smpl])
  
  # train classifier
  classifier <- setParams(learner=classifier, trainData=trainClassData)
  
  # getHyperPar(classifier, "sigma")
  keepHypers <- lapply(keepHyperParams, function(par) getHyperPar(classifier, par))
  
  # NOTE: need to check if giving each group of dists a different seed is better for 
  # KCDC classification as I observed for sinx example with discreete bins
  data <- constructData(x=x, y=c(y_tr, y_te))
  classifier <- classifier$learn(learner=classifier)
  
  phix <- classifier$makeFeature(learner=classifier, data, "X")
  
  meanPhix <- apply(phix, 2, mean)
  
  pred <- classifier$pred(classifier, testClassData)
  loss <- do.call(classifier$optimizeParams$losses[[1]]$func, list(learner=classifier, pred=pred))
  
  res <- list(meanPhix=meanPhix, loss=loss, hyperPars=keepHypers)
  
  if(ncol(x)==1){
    xx <- seq(0, 1, length.out=100)
    pdfData <- constructData(x=as.matrix(xx), y=rep(0,100))
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
  res <- cbind(preds, resids, 1)
  colnames(res) <- c("pred","resid","grp")
  return(res)
} 

resids.add_class <- function(learner, pred){
  preds <- pred$gyh
  resids <- as.matrix(pred$gyh_class - pred$gy_class)
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

# loss functions used in lnkrr
ksTest <- function(pvals){
  #print("enters hsicTest")
  
  pval <- ks.test(pvals, y="punif")$p.value
  #print("exits hsicTest")
  return(pval)
}
hsicTestRX <- function(learner, pred){
  #print("enters hsicTest")
  r <- pred$r
  x <- pred$x
  n <- nrow(x)
  ps <- getHyperPar(learner, "ps") #seq(0.1, 0.5, 0.1)
  sigma_rs  <- 1/quantile(as.numeric(dist(unique(r))^2), prob=ps)
  sigma_x <- 1/quantile(as.numeric(dist(unique(x))^2), prob=0.5)
  
  Kx <- kern_rbf(x, sigma=sigma_x)
  Krs <- sapply(sigma_rs, function(sigma_r){
    Kr <- kern_rbf(r, sigma=sigma_r)
    return(Kr)
  }, simplify="array")
  dim(Krs) <- c(n, n, length(ps))
  Kr <- apply(Krs, c(1,2), sum)/length(ps)
  
  Ks <- vector("list", 2)
  Ks[[1]] <- Kx
  Ks[[2]] <- Kr
  
  pval <- dhsic.test(K=Ks, method="gamma")$p.value  
  
  
  #print("exits hsicTest")
  return(pval)
}
hsicTestResX <- function(learner, pred){
  #print("enters hsicTest")
  
  pval <- dhsic.test(X=learner$resids(learner,pred)[,"resid"], Y=pred$x, method="gamma")$p.value  
  
  
  #print("exits hsicTest")
  return(pval)
}
hsicTestResR <- function(learner, pred){
  #print("enters hsicTest")
  
  pval <- dhsic.test(X=learner$resids(learner,pred)[,"resid"], Y=pred$r, method="gamma")$p.value  
  
  
  #print("exits hsicTest")
  return(pval)
}
hsicTestResZ <- function(learner, pred){
  #print("enters hsicTest")
  
  z <- cbind(pred$x, pred$r)
  pval <- dhsic.test(X=learner$resids(learner,pred)[,"resid"], Y=z, method="gamma")$p.value  
  
  
  #print("exits hsicTest")
  return(pval)
}

ksTestHsicRX <- function(learner, pred) learner$learnParams$ksPvalHsicRX
ksTestHsicResX <- function(learner, pred) learner$learnParams$ksPvalHsicResX
ksTestHsicResR <- function(learner, pred) learner$learnParams$ksPvalHsicResR
ksTestHsicResZ <- function(learner, pred) learner$learnParams$ksPvalHsicResZ
rmse2 <- function(learner, pred) learner$learnParams$rmseOpt


negLogLik <- function(learner, pred){
  negLogLiks <- learner$resids(learner, pred)[,"pred"]
  return(sum(negLogLiks))
  
}
hsicLoss2 <- function(learner, pred){
  cdfs <- learner$resids(learner, pred)[,"resid"]
  res <- dhsic.test(pred$x, cdfs, method="gamma")$statistic
  return(res)
}
pinball <- function(learner, pred){
  taus <-  getHyperPar(learner, "taus")  
  preds <- pred$gyh
  pinballs <- sapply(1:(length(taus)+2), function(i) {
    # i <- 1
    opera:::loss(preds[,i], pred$gy, loss.type=list(name="pinball", tau=c(3/(nrow(preds)), taus, 1-3/(nrow(preds)))[i]))
  }, simplify="array")
  
  pinballs <- apply(pinballs, 2, function(col) col/sd(col, na.rm=T))
  pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
  return(pinball)
  
}




RMSE2 <- function(learner, pred){
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
    res <- RMSE(pred[[nms[1]]], pred[[nms[2]]])  
  })
  
  
  return(res)
}

RMSE2 <- function(learner, pred){
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
    res <- -cor(pred[[nms[1]]], pred[[nms[2]]], method="spearman")  
  })
  
  
  return(res)
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

negCE3 <- function(learner, pred){
  
  res <- MLmetrics:::LogLoss(pred$gyh_prob, pred$gy_class)  
  
  return(res)
}

acc <- function(learner, pred){
  
  res <- 1-MLmetrics:::Accuracy(pred$gyh_class, pred$gy_class)  
  
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

myAUC <- function(learner, pred){
  ypred <- pred$gyh
  ytrue <- pred$gy
  rocRes <- roc(c(ytrue), c(ypred), direction = "<", algorithm = 1, quiet=TRUE)
  return(1-rocRes$auc)
}

negCE2 <- function(learner, pred){
  ypred <- pred$gyh
  ytrue <- pred$gy
  LogLoss(ypred, ytrue)
}

hingeLoss <- function(learner, pred){
  ypred <- pred$gyh
  ytrue <- pred$gy
  d <- 1 - ytrue * ypred
  d[d < 0] <- 0
  return(mean(d))
}

# only for kernel methods - evaluates regularizer t(alpha)%*%Kxs%*%alpha
regL <- function(learner, pred){
  
  alpha <- learner$learnParams$alpha
  
  parsXs <- getKernelPars(learner, kernelName="kernelXs")
  
  Kxs <- kernelMatrix(learner$hyperParams$non_data$kernelXs$name, pred$x, pars=parsXs)
  
  return(as.numeric(t(alpha)%*%Kxs%*%alpha))
}

qhsicLoss <- function(learner, pred){
  
  parsXb <- getKernelPars2(learner, kernelName="kernelXb")
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
  classData <- constructData(x=x_tr, y=y_tr)
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
  classData <- constructData(x=x_tr, y=y_tr)
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  #print("a.3")
  if(!is.positive.definite(Blambda)){
    Blambda <- make.positive.definite(Blambda)
  }
  Ky <- learner$learnParams$Ky
  #print("a.4")
  if(!is.positive.definite(Ky)){
    Ky <- make.positive.definite(Ky)
  }
  #print("a.5")
  Alambda <- Blambda %*% Ky %*% Blambda
  #print("a.6")
  if(!is.positive.definite(Alambda)){
    #print("a.7")
    #print(Alambda[1:5,1:5])
    #print(summary(as.numeric(Alambda)))
    Alambda <- try(make.positive.definite(Alambda))
    #print("a.8")
    if(class(Alambda)=="try-error") Alambda <- matrix(rnorm(nrow(Ky)*ncol(Ky)), nrow(Ky), ncol(Ky))
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
  #print(summary(as.numeric(dky_with_fakes)))
  
  #print(is.positive.definite(Blambda))
  #print(is.positive.definite(Ky))
  #print(is.positive.definite(Alambda))
  
  const <- diag(lx%*%Alambda%*%t(lx))*as.numeric(dky_with_fakes)
  const[which(const<0)] <- 0
  const <- sqrt(const)
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  if(all(probsNorm==0)) probsNorm <- runif(length(probsNorm), min=0, max=1e-9)
  #print(summary(probsNorm))
  
  # make train classification data
  #x_tr <- cbind(probs, const)
  x_tr <- as.matrix(probsNorm)
  y_tr <- c(rep(1, n), rep(0, num_fake)) 
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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

makeCME_cond_regfeatsK_x <- function(learner, data, fakes_x, fakes_y){
  
  n <- nrow(data$x)
  num_fake <- length(fakes_y)
  kappa <- num_fake/n
  
  
  x_real <- data$x
  y_real <- data$y
  x_fake <- fakes_x
  y_real_ext <- lapply(1:kappa, function(i) y_real)
  y_real_ext <- do.call(c, y_real_ext)
  x_real_ext <- lapply(1:kappa, function(i) x_real)
  x_real_ext <- do.call(rbind, x_real_ext)
  
  x <- rbind(x_real, x_fake)
  y <- c(y_real, y_real_ext)
  
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
  trainData <- learner$hyperParams$trainData
  Blambda <- learner$learnParams$Blambda
  #print("a.3")
  if(!is.positive.definite(Blambda)){
    Blambda <- make.positive.definite(Blambda)
  }
  Ky <- learner$learnParams$Ky
  #print("a.4")
  if(!is.positive.definite(Ky)){
    Ky <- make.positive.definite(Ky)
  }
  #print("a.5")
  Alambda <- Blambda %*% Ky %*% Blambda
  #print("a.6")
  if(!is.positive.definite(Alambda)){
    #print("a.7")
    #print(Alambda[1:5,1:5])
    #print(summary(as.numeric(Alambda)))
    Alambda <- try(make.positive.definite(Alambda))
    #print("a.8")
    if(class(Alambda)=="try-error") Alambda <- matrix(rnorm(nrow(Ky)*ncol(Ky)), nrow(Ky), ncol(Ky))
  }
  
  pred_withFakes <- learner$predict(learner, data=dataWithFakes, forLoss=F)
  ky <- pred_withFakes$gy_k
  lx <- pred_withFakes$lx
  
  data1_xy <- constructData(x=rbind(x_real, x_fake), y=c(y_real, y_real_ext))
  data2_xy <- constructData(x=rbind(x_real, x_real_ext), y=c(y_real, y_real_ext))
  
  deg_fakeness_xy <- makeKernel(learner, data1=data1_xy, data2=data2_xy, var="X", grad=FALSE)$K
  deg_fakeness_xx <- makeKernel(learner, data1=data2_xy, var="X", grad=FALSE)$K
  deg_fakeness_yy <- makeKernel(learner, data1=data1_xy, var="X", grad=FALSE)$K
  deg_fakeness_xy <- diag(deg_fakeness_xy)
  deg_fakeness_xx <- diag(deg_fakeness_xx)
  deg_fakeness_yy <- diag(deg_fakeness_yy)
  # Normalization not necessary (as long as rbf kernel used!) coz diagonal of rbf kernel always 1
  deg_fakeness <- deg_fakeness_xy/(sqrt(deg_fakeness_xx)*sqrt(deg_fakeness_yy))
  deg_fakeness <- deg_fakeness - mean(deg_fakeness)
  
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
  #print(summary(as.numeric(dky_with_fakes)))
  
  #print(is.positive.definite(Blambda))
  #print(is.positive.definite(Ky))
  #print(is.positive.definite(Alambda))
  
  const <- diag(lx%*%Alambda%*%t(lx))*as.numeric(dky_with_fakes)
  const[which(const<0)] <- 0
  const <- sqrt(const)
  #sum(const<0)
  probsNorm <- probs /const
  probsNorm[which(probs==0 | const==0)] <- 0
  if(all(probsNorm==0)) probsNorm <- runif(length(probsNorm), min=0, max=1e-9)
  #print(summary(probsNorm))
  
  # make train classification data
  #x_tr <- cbind(probs, const)
  x_tr <- as.matrix(probsNorm)
  #y_tr <- c(rep(1, n), rep(0, num_fake))
  y_tr <- deg_fakeness
  set.seed(12)
  smpl <- sample(nrow(x_tr)) #jumble them up
  x_tr <- x_tr[smpl,,drop=F]
  y_tr <- y_tr[smpl]
  classData <- constructData(x=x_tr, y=y_tr)
  # plot(x_tr, y_tr)
  
  return(classData)  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  classData <- constructData(x=x_tr, y=y_tr)
  
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
  dataWithFakes <- constructData(x=x, y=as.matrix(y))
  
  
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
  if(class(Blambda)=="try-error"){
    #print("try-error")
    Blambda <- matrix(0,n,n)
  }
  
  
  
  beta_k <- Blambda %*% Ky
  alphas <- Blambda %*% Lx
  
  
  learner$learnParams$Lx <- Lx
  learner$learnParams$Ky <- Ky
  learner$learnParams$Blambda <- Blambda
  learner$learnParams$Cks <- Cks
  learner$learnParams$beta_k <- beta_k
  learner$learnParams$alphas <- alphas
  
  #print("a.1")
  
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
  
  if(forLoss & any( c("RMSE2","negCE") %in% names(learner$optimizeParams$losses))){
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
    
    #print("a.2")
    
    #trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
    fits <- mcmapply(function(featFunc, classifier){
      # i <- 1; featFunc <- featFuncs[i]; classifier <- classifiers[i]
      #print("**************")
      #print(paste(featFunc, classifier, sep="_"))
      
      trainClassData <- do.call(featFunc, list(learner=learner, data=learner$hyperParams$trainData, fakes_x=fakes_x, fakes_y=fakes_y))
      
      
      
      # train classifier
      classifier <- eval(parse(text=classifier))
      classifier <- setParams(learner=classifier, trainData=trainClassData, mc_cores=1)
      
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
              gyh_f_phiT=as.matrix(diag(gyh_f_phiT)), gyh_k=gyh_k, residsPerObs=as.matrix(residsPerObs))
  
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
  
  
  if(forLoss & any( c("RMSE2","negCE") %in% names(learner$optimizeParams$losses))){
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
  Lx <- krnlX$K
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
  if(forLoss & any( c("RMSE2","negCE") %in% names(learner$optimizeParams$losses))){
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
  trainClassData <- constructData(x=x_tr[smpl,,drop=F], y=y_tr[smpl])
  
  # make test classification data
  x_te <- as.matrix(rbind(x[indx.te,,drop=F], matrix(runif(n_te*p), n_te, p)))
  y_te <- c(rep(1, n_te), rep(0, n_te)) 
  set.seed(123)
  smpl <- sample(nrow(x_te))
  testClassData <- constructData(x=x_te[smpl, ,drop=F], y=y_te[smpl])
  
  # train classifier
  classifier <- setParams(learner=classifier, trainData=trainClassData)
  
  # getHyperPar(classifier, "sigma")
  keepHypers <- lapply(keepHyperParams, function(par) getHyperPar(classifier, par))
  
  # NOTE: need to check if giving each group of dists a different seed is better for 
  # KCDC classification as I observed for sinx example with discreete bins
  data <- constructData(x=x, y=c(y_tr, y_te))
  classifier <- classifier$learn(learner=classifier)
  
  phix <- classifier$makeFeature(learner=classifier, data, "X")
  
  meanPhix <- apply(phix, 2, mean)
  
  pred <- classifier$pred(classifier, testClassData)
  loss <- do.call(classifier$optimizeParams$losses[[1]]$func, list(learner=classifier, pred=pred))
  
  res <- list(meanPhix=meanPhix, loss=loss, hyperPars=keepHypers)
  
  if(ncol(x)==1){
    xx <- seq(0, 1, length.out=100)
    pdfData <- constructData(x=as.matrix(xx), y=rep(0,100))
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

wcmem_L2_f <- function(learner, pred){
  
  #print("summary(pred$residsPerObs)")
  #print(summary(pred$residsPerObs))
  
  lx <- pred$lx
  aux <- apply(lx,2, function(col) sqrt(abs(weighted.mean(pred$residsPerObs, w=col))))
  #print("summary(aux)")
  #print(summary(aux))
  res <- quantile(aux, 0.9, na.rm=T)
  names(res) <- NULL
  
  return(res)
}

wcmem_L2_f_nonagg_loc <- function(learner, pred){
  
  #print("summary(pred$residsPerObs)")
  #print(summary(pred$residsPerObs))
  numFolds <- learner$optimizeParams$numFolds
  n <- getN(learner$hyperParams$trainData)
  size <- floor(n / numFolds)
  indxTrain <- lapply(1:numFolds, function(f){
    last <- f*size
    first <- (f-1)*size + 1
    if(f==numFolds){
      first <- n-size+1
      last <- n
    }
    res <- first:last
  })
  indxTrain <- unlist(indxTrain)
  numTrain <- length(indxTrain)
  indxTest <- lapply(1:numFolds, function(f){
    last <- f*size
    first <- (f-1)*size + 1
    if(f==numFolds){
      first <- n-size+1
      last <- n
    }
    res <- first:last
    res <- setdiff(1:n, res)
  })
  indxTest <- unlist(indxTest)
  numTest <- length(indxTest)  
  if(length(pred$residsPerObs)==length(indxTrain)){
    indx <- indxTrain
  } else{
    indx <- indxTest
  }
  
  #lx <- pred$lx
  lx <- learner$learnParams$Lx
  lx <- lx[indx,]
  res <- apply(lx,2, function(col) sqrt(abs(weighted.mean(pred$residsPerObs, w=col))))
  #print("summary(aux)")
  #print(summary(aux))
  #res <- quantile(aux, 0.9, na.rm=T)
  names(res) <- NULL
  names(res) <- paste("wcmem_L2_f_nonagg_loc", 1:length(res), sep="_")
  return(res)
}
wcmem_L2_f_nonagg_superloc <- function(learner, pred){
  
  #print("summary(pred$residsPerObs)")
  #print(summary(pred$residsPerObs))
  
  #lx <- pred$lx
  lx <- learner$learnParams$Lx
  res <- pred$residsPerObs
  #res <- apply(res, 2, mean)
  #print("summary(aux)")
  #print(summary(aux))
  #res <- quantile(aux, 0.9, na.rm=T)
  names(res) <- NULL
  names(res) <- paste("wcmem_L2_f_nonagg_superloc", 1:length(res), sep="_")
  return(res)
}

cmem_L2_f_rel <- function(learner, pred){
  
  res <- sqrt(mean(as.numeric(t(pred$residsPerObs))/pred$gy_f_phiT^2))
  
  
  return(res)
}

wcmem_L2_f_rel_nonagg_loc <- function(learner, pred){
  
  #print("summary(pred$residsPerObs)")
  #print(summary(pred$residsPerObs))
  
  #res <- sqrt(mean(as.numeric(t(pred$residsPerObs))/pred$gy_f_phiT^2))
  
  numFolds <- learner$optimizeParams$numFolds
  n <- getN(learner$hyperParams$trainData)
  size <- floor(n / numFolds)
  indxTrain <- lapply(1:numFolds, function(f){
    last <- f*size
    first <- (f-1)*size + 1
    if(f==numFolds){
      first <- n-size+1
      last <- n
    }
    res <- first:last
  })
  indxTrain <- unlist(indxTrain)
  numTrain <- length(indxTrain)
  indxTest <- lapply(1:numFolds, function(f){
    last <- f*size
    first <- (f-1)*size + 1
    if(f==numFolds){
      first <- n-size+1
      last <- n
    }
    res <- first:last
    res <- setdiff(1:n, res)
  })
  indxTest <- unlist(indxTest)
  numTest <- length(indxTest)  
  if(length(pred$residsPerObs)==length(indxTrain)){
    indx <- indxTrain
  } else{
    indx <- indxTest
  }
  
  #lx <- pred$lx
  lx <- learner$learnParams$Lx
  lx <- lx[indx,]
  
  res1 <- apply(lx,2, function(col) sqrt(abs(weighted.mean(pred$residsPerObs, w=col))))
  res2 <- apply(lx,2, function(col) sqrt(weighted.mean(pred$gy_f_phiT^2, w=col)))
  res <- res1/res2
  names(res) <- NULL
  #print("summary(aux)")
  #print(summary(aux))
  #res <- quantile(aux, 0.9, na.rm=T)
  #names(res) <- NULL
  names(res) <- paste("wcmem_L2_f_rel_nonagg_loc", 1:length(res), sep="_")
  return(res)
}
wcmem_L2_f_rel_nonagg_superloc <- function(learner, pred){
  
  #print("summary(pred$residsPerObs)")
  #print(summary(pred$residsPerObs))
  
  #lx <- pred$lx
  lx <- learner$learnParams$Lx
  res1 <- pred$residsPerObs
  #res1 <- apply(res1,2,mean)
  res2 <- pred$gy_f_phiT^2
  #res2 <- apply(res2,2,mean)
  res <- res1/res2
  #print("summary(aux)")
  #print(summary(aux))
  #res <- quantile(aux, 0.9, na.rm=T)
  names(res) <- NULL
  names(res) <- paste("wcmem_L2_f_rel_nonagg_superloc", 1:length(res), sep="_")
  return(res)
}

cmem_L2_k <- function(learner, pred){
  
  
  #  tr(K %*% t(K)) - 2 tr(K_hat %*% K) + tr(K_hat %*% K_hat)
  res <- mean(diag(pred$gy_k %*% t(pred$gy_k)) - 2*diag(pred$gyh_k %*% t(pred$gy_k)) + diag(pred$gyh_k  %*% t(pred$gyh_k)))
  res <- sqrt(res)
  return(res)
}

wcmem_L2_k <- function(learner, pred){
  
  lx <- pred$lx
  #  tr(K %*% t(K)) - 2 tr(K_hat %*% K) + tr(K_hat %*% K_hat)
  aux <- diag(pred$gy_k %*% t(pred$gy_k)) - 
    2*diag(pred$gyh_k %*% t(pred$gy_k)) + diag(pred$gyh_k  %*% t(pred$gyh_k))
  aux <- apply(lx, 2, function(col) weighted.mean(aux,w=col))
  aux <- abs(sqrt(aux))
  res <- quantile(aux, 0.9, na.rm=T)
  #res <- sqrt(res)
  names(res) <- NULL
  return(res)
}

corrRKHS <- function(learner, pred){
  
  
  LBK <- pred$gyh_k
  Lte <- pred$lx
  method <- "spearman"
  #print("summary(c(LBK))")
  #print(summary(c(LBK)))
  #print("summary(c(Lte))")
  #print(summary(c(Lte)))
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    #res <- cor(Lte_col,LBK_col, method="spearman")
    res <- cor(Lte_col,LBK_col, method=method)
    #print(res)
    return(-res)
  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- match(sort(corrs,decreasing=T),corrs)[3]
  #plot(Lte[indx,], LBK[indx,], cex=Lte[indx,])
  #res <- mean(corrs)
  #print("summary(corrs)")
  #print(summary(corrs))
  res <- quantile(corrs, 0.5, na.rm=T)
  names(res) <- NULL
  return(res)
}

wcorrRKHS <- function(learner, pred){
  
  
  LBK <- pred$gyh_k
  Lte <- pred$lx
  method <- "spearman"
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    #res <- cor(Lte_col,LBK_col, method="spearman")
    res <- weightedCorr(Lte_col,LBK_col, method=method, weights=Lte_col/sum(Lte_col))
    #print(res)
    return(-res)
  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- match(sort(corrs,decreasing=T),corrs)[3]
  #plot(Lte[indx,], LBK[indx,], cex=Lte[indx,])
  #res <- mean(corrs)
  res <- quantile(corrs, 0.5, na.rm=T)
  names(res) <- NULL
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
  
  #print(msrs_char)
  mcCores <- min(40, detectCores()-4)
  #mcCores <- 1
  msrs_val <- mcmapply(function(msr){
    # msr <- msrs_char[6]
     #print(paste("measure: ", msr))
    res <- do.call(msr, list(learner=cmemLearner))
    }, msr=msrs_char, mc.cores=mcCores, SIMPLIFY=FALSE)
  names(msrs_val) <- msrs_char
  nms <- lapply(msrs_char, function(el){
    if(length(msrs_val[[el]])>1){
      res <- paste(el, 1:length(msrs_val[[el]]) ,sep="_")
    } else{
      res <- el
    }
    
    })
  msrs_val <- unlist(msrs_val)
  nms <- unlist(nms)
  names(msrs_val) <- nms
  return(msrs_val)
}

# Kernel Conditional Deviance for Causal inference (KCDC)

# weighted according to density - trial
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  #densX <- getHyperPar(learner, "densX")
  densX <- learner$hyperParams$data$non_optimizable$densX$val
  weights <-  1/densX@estimate
  weights <- weights/sum(weights)
  res <- weighted.var(dLAL, w=weights)
  #names(res) <- NULL
  return(res)
}

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
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  
  b <- sum(dLAL^(0.5))
  c <- sum(dLAL)
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  return(res)
}

# weighted according to density (and locally) - trial
WKCDC <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  #densX <- getHyperPar(learner, "densX")
  densX <- learner$hyperParams$data$non_optimizable$densX$val
  weights <-  1/densX@estimate
  weights <- weights/sum(weights)
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res <- quantile(apply(L, 2, function(col) weighted.var(dLAL, w=col*weights)),0.9)
  names(res) <- NULL
  return(res)
}

# aggregated
WKCDC <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res <- quantile(apply(L, 2, function(col) weighted.var(dLAL, w=col)),0.9)
  names(res) <- NULL
  return(res)
}

WKCDC_nonagg_loc <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res <- apply(L, 2, function(col) weighted.var(dLAL, w=col))
  names(res) <- paste("WKCDC_nonagg_loc", 1:length(res), sep="_")
  return(res)
}

WKCDC_nonagg_superloc <- function(learner, pred=NULL){
  
  L  <- kern_rbf(learner$hyperParams$trainData$x, sigma=getHyperPar(learner,"sigma.rbf.X")/100)#learner$learnParams$Lx
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res <- apply(L, 2, function(col) weighted.var(dLAL, w=col))
  names(res) <- paste("WKCDC_nonagg_loc", 1:length(res), sep="_")
  return(res)
}


# mean instead of variance
WKCLC_nonagg_loc <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res <- apply(L, 2, function(col) weighted.mean(dLAL, w=col))
  names(res) <- paste("WKCLC_nonagg_loc", 1:length(res), sep="_")
  return(res)
}

WKCLC_nonagg_superloc <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  res <- dLAL^0.5
  names(res) <- paste("WKCLC_nonagg_superloc", 1:length(res), sep="_")
  return(res)
}


# coef var
WKCRDC_nonagg_loc <- function(learner, pred=NULL){
  
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res1 <- apply(L, 2, function(col) weighted.var(dLAL, w=col))
  res2 <- apply(L, 2, function(col) weighted.mean(dLAL, w=col))
  res <- res1/res2
  names(res) <- paste("WKCRDC_nonagg_loc", 1:length(res), sep="_")
  return(res)
}

WKCRDC_nonagg_superloc <- function(learner, pred=NULL){
  
  #L  <- learner$learnParams$Lx
  L  <- kern_rbf(learner$hyperParams$trainData$x, sigma=getHyperPar(learner,"sigma.rbf.X")/100)#learner$learnParams$Lx
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
  
  #b <- sum(diag(LAL)^(0.5))
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (b/n)^2
  # var(diag(LAL)^0.5)*(n-1)/n
  
  dLAL <- diag(LAL)
  dLAL[which(dLAL<0)] <- 0
  dLAL <- dLAL^0.5
  res1 <- apply(L, 2, function(col) weighted.var(dLAL, w=col))
  res2 <- apply(L, 2, function(col) weighted.mean(dLAL, w=col))
  res <- res1/res2
  names(res) <- paste("WKCRDC_nonagg_loc", 1:length(res), sep="_")
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

KCDCpval2 <- function(learner){
  #x <- learner$hyperParams$trainData$x
  n <- nrow(learner$learnParams$Ky)
  set.seed(12345)            
  rperms <- sapply(1:learner$msrs[["KCDCpval2"]]$pars$numPerms, function(i) sample(n))
  
  mesr <- KCDC(learner)
  
  trainDataU <- learner$hyperParams$trainData
  trainDataU$y <- trainDataU$y[rperms[,1]]
  mc_cores <- min(40, detectCores()-4)
  learner2 <- learner
  #learner2$hyperParams$data$optimizable$lambda$val <- NULL
  learner2$hyperParams$data$optimizable$lambda$val <- 1e8
    
  #learner2 <- setParams(learner=learner2, trainData=trainDataU, mc_cores=mc_cores)
  #learner2 <- learner2$learn(learner=learner2, forLoss=F)
  
  Ky <- learner2$learnParams$Ky[rperms[,1],]
  Ky <- Ky[,rperms[,1]]
  learner2$learnParams$Ky <- Ky
  
  KCDC_dist_null <- apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    learnerAux <- learner2
    Ky <- learner2$learnParams$Ky[rperm,]
    Ky <- Ky[,rperm]
    learnerAux$learnParams$Ky <- Ky
    res <- KCDC(learnerAux) 
    return(res)
  })
  
  pval <- 1-sum(KCDC_dist_null>mesr)/learner$msrs[["KCDCpval2"]]$pars$numPerms
  
  
  med.gam<-mean(KCDC_dist_null) ## sample mean 
  var.gam<-var(KCDC_dist_null) ## sample variance 
  l.est<-med.gam/var.gam ## lambda estimate (corresponds to rate) 
  a.est<-((med.gam)^2)/var.gam ## alfa estimate 
  xx <- seq(0,max(KCDC_dist_null)*1.5, length.out=100)
  yy <- dgamma(xx, shape=a.est, rate=l.est)
  hist(KCDC_dist_null, prob=T)
  lines(xx, yy, col="red")
  pval <- pgamma(mesr, shape=a.est, rate=l.est)
  
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

RMSEte <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  mainLossNm2 <- learner$optimizeParams$mainLoss
  mainLossNm <- strsplit(mainLossNm2, "\\.")[[1]][1] 
  featFunc <- learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
  classifier <- learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier
  if(length(featFunc)>1){
    featFunc <- "makeCME_cond_regfeatsK_x"
    classifier <- "krr1"
    learner$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    learner$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  if(mainLossNm != "RMSE2") mainLossNm2 <- paste("RMSE2", paste(featFunc, classifier, sep="-") ,sep=".")
  
  if( mainLossNm2 %in% dimnames(grid)[[3]] & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var ==mainLossNm2)
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- -grid[indxMat]
  } else{
    
    lossFun <-  "RMSE2" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- lapply(lossFun, function(el) list(func=el))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- lapply(lossFun, function(el) list(func=el))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    learner$optimizeParams$mainLoss <- mainLossNm2
    optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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

wL2_f_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_f" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_f")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_f" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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


wL2_f_te_nonagg_loc <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_f_nonagg_loc" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    print("from grid")
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_f_nonagg_loc")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_f_nonagg_loc" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- grid["test",,]
  }
  return(res)
}

wL2_f_te_nonagg_superloc <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_f_nonagg_superloc" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_f_nonagg_superloc")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_f_nonagg_superloc" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- grid["test",,]
  }
  return(res)
}

wL2_f_rel_nonagg_loc <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_f_rel_nonagg_loc" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_f_rel_nonagg_loc")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_f_rel_nonagg_loc" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- grid["test",,]
  }
  return(res)
}

wL2_f_rel_nonagg_superloc <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_f_rel_nonagg_superloc" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_f_rel_nonagg_superloc")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_f_rel_nonagg_superloc" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- grid["test",,]
  }
  return(res)
}

 



L2_k_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("cmem_L2_k" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="cmem_L2_k")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "cmem_L2_k" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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

wL2_k_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcmem_L2_k" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcmem_L2_k")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- grid[indxMat]
  } else{
    lossFun <-  "wcmem_L2_k" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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

corrRKHS_k_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("corrRKHS" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="corrRKHS")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- -grid[indxMat]
  } else{
    lossFun <-  "corrRKHS" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- -grid["test",,1]
  }
  return(res)
}

wcorrRKHS_k_te <- function(learner){
  grid <- learner$hyperParams$data$grid
  indxOptGrid <- learner$hyperParams$data$indxOptGrid
  if("wcorrRKHS" %in% names(learner$optimizeParams$losses) & !is.null(grid)){
    indxTest <- which(dimnames(grid)$trainTest=="test")
    indxLoss <- which(dimnames(grid)$var =="wcorrRKHS")
    indxMat <- matrix(c(indxTest, indxOptGrid, indxLoss),  1, length(dim(grid)))
    res <- -grid[indxMat]
  } else{
    lossFun <-  "wcorrRKHS" 
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    
    res <- -grid["test",,1]
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
    lossFunList <- list(list(func=lossFun))
    names(lossFunList) <- lossFun
    learner$optimizeParams$losses<- lossFunList 
    optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
    optLossFunc <- eval(parse(text=optLossFunc))
    learner$optimizeParams$optLossFunc <- optLossFunc
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learner$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learner$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
  
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses <- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
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
  lossFunList <- list(list(func=lossFun))
  names(lossFunList) <- lossFun
  
  learnerAux$optimizeParams$losses<- lossFunList 
  optLossFunc <- paste("function(grid) which.min(grid[,'",lossFun,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  learner$optimizeParams$optLossFunc <- optLossFunc
  
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


# weighted by density - trial
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
  
  #f <- sum(LAL)
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (f/n^2)
  
  #densX <- getHyperPar(learner, "densX")
  densX <- learner$hyperParams$data$non_optimizable$densX$val
  ws <-  1/densX@estimate
  ws <- ws/sum(ws)
  
  res <- weighted.mean(diag(LAL), w=ws) - t(ws)%*%LAL%*%ws
  
  # mean(diag(LAL))-mean(LAL)
  # mean(diag(LAL)-LAL)
  
  return(res)
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
  # mean(diag(LAL))-mean(LAL)
  # mean(diag(LAL)-LAL)
  
  return(res)
}


# weighted by density (and locally) - trial
WKCMC <- function(learner, pred=NULL){
  
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
  
  #f <- sum(LAL)
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (f/n^2)
  
  #densX <- getHyperPar(learner, "densX")
  densX <- learner$hyperParams$data$non_optimizable$densX$val
  ws2 <-  1/densX@estimate
  ws2 <- ws2/sum(ws2)
  
  res <- apply(L, 2, function(col){
    # col <- L[,1]
    ws <- col/sum(col)
    ws <- ws*ws2
    ws <- ws/sum(ws)
    res <- weighted.mean(diag(LAL), w=ws) - t(ws)%*%LAL%*%ws
    return(res)
  })
  res <- quantile(res,0.9)
  names(res) <- NULL
  # mean(diag(LAL))-mean(LAL)
  # mean(diag(LAL)-LAL)
  
  return(res)
}


# Weighted Kernel Conditional Mean distance  for Causal inference (KCMC)
WKCMC <- function(learner, pred=NULL){
  
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
  
  #f <- sum(LAL)
  #c <- sum(diag(LAL))
  
  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  #res <- (c/n) - (f/n^2)
  
  res <- apply(L, 2, function(col){
    # col <- L[,1]
    ws <- col/sum(col)
    res <- weighted.mean(diag(LAL), w=ws) - t(ws)%*%LAL%*%ws
    return(res)
    })
  res <- quantile(res,0.9)
  names(res) <- NULL
  # mean(diag(LAL))-mean(LAL)
  # mean(diag(LAL)-LAL)
  
  return(res)
}


# KIIM - kernel intrinsic invariance measure 

KIIM <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  K  <- learner$learnParams$Ky
  n <- nrow(L)
  Blambda <- learner$learnParams$Blambda
  I <- diag(n)
  H <- I-matrix(1/n,n,n)
  
  M <- K %*% Blambda %*% L %*% H %*% L %*% Blambda %*% K
  
  if(!is.positive.definite(M)){
    M <- try(make.positive.definite(M))
    #print("a.8")
    if(class(M)=="try-error") M <- matrix(rnorm(nrow(K)*ncol(K)), nrow(K), ncol(K))
  }
  #M <- make.positive.definite(M)
  
  eigM <- rev(eigen(M, only.values=T, symmetric=T)$values)
  
  if(sum(eigM)==0){
    cumeigM <- rep(0,length(eigM))
    Fin <- length(eigM)
  } else{
    cumeigM <- cumsum(eigM)/sum(eigM)*100
    Fin <- which(cumeigM>90)[1]
    if(Fin==length(cumeigM)) Fin <- length(cumeigM)-1
  }
  
  
  res <- sum(eigM[1:Fin])
  
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
    Alambda <- beta_f%*%t(beta_f)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    pm <- proc.time()
    E <- K%*%K
    Alambda <- Blambda%*%K%*%t(Blambda)
    proc.time() - pm # 0.002
  }
  
  
  pm <- proc.time()
  res <- sum(diag(Alambda*Chat))
  proc.time() - pm #0.001
  
  #pm <- proc.time()
  #res2 <- sum(diag(K%*%D%*%K))
  #proc.time() - pm # 0.002
  
  return(res)
}

KICSC <- function(learner, pred=NULL){
  
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
    Alambda <- beta_f%*%t(beta_f)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    pm <- proc.time()
    E <- K%*%K
    Alambda <- Blambda%*%K%*%t(Blambda)
    proc.time() - pm # 0.002
  }
  
  
  #pm <- proc.time()
  #res <- sum(diag(Alambda*Chat))
  #proc.time() - pm #0.001
  
  #M <- Alambda*Chat
  M <- K %*% Blambda %*% Chat %*% Blambda %*% K
  
  eigM <- rev(eigen(M, only.values=T, symmetric=T)$values)
  cumeigM <- cumsum(eigM)/sum(eigM)*100
  Fin <- which(cumeigM>90)[1]
  if(Fin==length(cumeigM)) Fin <- length(cumeigM)-1
  res <- sum(eigM[1:Fin])
  
  
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

# pca on phi - i did this for UAI by accident
KCCC_pca_ent <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  phiy <- learner$learnParams$phiy
  
  
  res.pca <- prcomp(phiy, scale = F)
  #library(factoextra)
  #fviz_eig(res.pca)
  res.pca <- princomp(covmat=phiy%*%t(phiy), cor=F)
  
  #phiy <- res.pca$x
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

# eigen values of R matrix (see UAI notation) : entropy of all eigen values
KCCC_pca_ent2 <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  phiy <- learner$learnParams$phiy
  
  
  #res.pca <- prcomp(phiy, scale = F)
  #library(factoextra)
  #fviz_eig(res.pca)
  #res.pca <- princomp(covmat=phiy%*%t(phiy), cor=F)
  
  #phiy <- res.pca$x
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
  
  M <- t(phiy)%*%D%*%phiy
  eigM <- eigen(M, only.values=T, symmetric=T)$values
  
  res <- sum(abs(eigM)*log(abs(eigM)))
  
  #res <- sum(eigM)
  
  
  #S <- diag(M)
  #T <- sum(D*K)
  
  
  
  
  
  # n <- 100; m <- 10; phiy <- matrix(rnorm(n*m),n,m)
  
  #proc.time()
  #S <- sapply(1:n, function(i) sapply(1:n, function(j) D[i,j]*phiy[i,]*phiy[j,], simplify="array"), simplify="array")
  #S <- apply(S, 1, sum)
  #proc.time() -pm # 15 secs
  
  
  #plot(S,S2); abline(a=0, b=1, col="red")
  
  #Stilde <- S/T
  
  
  
  #res <- sum(abs(Stilde)*log(abs(Stilde)))
  
  return(res)
}

# entropy on lowest 90% eigenvalues
KICCC_pca_ent <- function(learner, pred=NULL){
  
  L  <- learner$learnParams$Lx
  phiy <- learner$learnParams$phiy
  
  
  #res.pca <- prcomp(phiy, scale = F)
  #library(factoextra)
  #fviz_eig(res.pca)
  #res.pca <- princomp(covmat=phiy%*%t(phiy), cor=F)
  
  #phiy <- res.pca$x
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
  
  M <- t(phiy)%*%D%*%phiy
  eigM <- rev(eigen(M, only.values=T, symmetric=T)$values)
  cumeigM <- cumsum(eigM)/sum(eigM)*100
  Fin <- which(cumeigM>90)[1]
  if(Fin==length(cumeigM)) Fin <- length(cumeigM)-1
  
  res <- sum(abs(eigM[1:Fin])*log(abs(eigM[1:Fin])))
  
  #res <- sum(eigM[1:Fin])
  
  
  #S <- diag(M)
  #T <- sum(D*K)
  
  
  
  
  
  # n <- 100; m <- 10; phiy <- matrix(rnorm(n*m),n,m)
  
  #proc.time()
  #S <- sapply(1:n, function(i) sapply(1:n, function(j) D[i,j]*phiy[i,]*phiy[j,], simplify="array"), simplify="array")
  #S <- apply(S, 1, sum)
  #proc.time() -pm # 15 secs
  
  
  #plot(S,S2); abline(a=0, b=1, col="red")
  
  #Stilde <- S/T
  
  
  
  #res <- sum(abs(Stilde)*log(abs(Stilde)))
  
  return(res)
}

# entropy on phi R pht^T : entropy of all eigen values
KICCC_pca_ent2b <- function(learner, pred=NULL){
  
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
    Alambda <- beta_f%*%t(beta_f)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    pm <- proc.time()
    E <- K%*%K
    Alambda <- Blambda%*%K%*%t(Blambda)
    proc.time() - pm # 0.002
  }
  
  
  #pm <- proc.time()
  #res <- sum(diag(Alambda*Chat))
  #proc.time() - pm #0.001
  
  #M <- Alambda*Chat
  M <- K %*% Blambda %*% Chat %*% Blambda %*% K
  
  eigM <- rev(eigen(M, only.values=T, symmetric=T)$values)
  
  res <- sum(abs(eigM)*log(abs(eigM)))
  
  #cumeigM <- cumsum(eigM)/sum(eigM)*100
  #Fin <- which(cumeigM>90)[1]
  #if(Fin==length(cumeigM)) Fin <- length(cumeigM)-1
  #res <- sum(eigM[1:Fin])
  
  
  #pm <- proc.time()
  #res2 <- sum(diag(K%*%D%*%K))
  #proc.time() - pm # 0.002
  
  return(res)
}

# entropy on phi R pht^T : entropy on lowest 90% eigenvalues
KICCC_pca_entb <- function(learner, pred=NULL){
  
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
    Alambda <- beta_f%*%t(beta_f)
  } else{
    D <- t(Blambda)%*%Chat%*%Blambda
    pm <- proc.time()
    E <- K%*%K
    Alambda <- Blambda%*%K%*%t(Blambda)
    proc.time() - pm # 0.002
  }
  
  
  #pm <- proc.time()
  #res <- sum(diag(Alambda*Chat))
  #proc.time() - pm #0.001
  
  #M <- Alambda*Chat
  M <- K %*% Blambda %*% Chat %*% Blambda %*% K
  
  eigM <- rev(eigen(M, only.values=T, symmetric=T)$values)
  
  
  
  cumeigM <- cumsum(eigM)/sum(eigM)*100
  Fin <- which(cumeigM>90)[1]
  if(Fin==length(cumeigM)) Fin <- length(cumeigM)-1
  #res <- sum(eigM[1:Fin])
  
  res <- sum(abs(eigM[1:Fin])*log(abs(eigM[1:Fin])))
  
  #pm <- proc.time()
  #res2 <- sum(diag(K%*%D%*%K))
  #proc.time() - pm # 0.002
  
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
  
  res <- dhsic.test(K=Ks)
  #res <- res$dHSIC
  res <- res$p.value
  
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



