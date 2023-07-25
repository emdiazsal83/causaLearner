remove(list=ls())

library(mvtnorm)

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

n <- 100
w <- matrix(rnorm(n), n, 1)
sigma <- 0.7
alpha <- rnorm(n)
fofx <- function(x, w, alpha, sigma) kern_rbf_R(x, w, sigma=sigma)%*%alpha
onesN <- matrix(rep(1,n), n,1)
m <- 1000
onesM <- matrix(rep(1,m), m, 1)
x <- matrix(seq(-3,3, length.out=m),m,1)
X <- x %*% t(onesN)  
W <- onesM %*% t(w)
dim(X); dim(W); dim(kern_rbf_R(x, w, sigma=sigma))
X[1:3, 1:4]; W[1:3,1:4]
fprimeofx <- function(x, w, alpha, sigma){
  n <- nrow(w)
  m <- nrow(x)
  onesN <- matrix(rep(1,n), n, 1)
  onesM <- matrix(rep(1,m), m, 1)
  X <- x %*% t(onesN)  
  W <- onesM %*% t(w)
  #dim(X); dim(W); dim(kern_rbf_R(x, w, sigma=sigma))
  res <- 2*sigma*((kern_rbf_R(x, w, sigma=sigma)*(X-W))%*%alpha) 
  return(res)
}


n <- 100
w <- matrix(rnorm(n, mean=-4), n, 1)
sigma <- 0.7
alpha <- abs(rnorm(n))
x <- matrix(seq(-3,3, length.out=m),m,1)
y <- fofx(x, w, alpha, sigma)
yp <- fprimeofx(x, w, alpha, sigma)
plot(x, y, "l", ylim=range(y,yp))
lines(x, yp, col="red")

rKalpha_rbf <- function(x, tau, xInd=NULL, inducing=FALSE, sigma=NULL, minDR=NULL, maxDR=NULL, seed=NULL, mean_alpha=0){
  #if((ncol(x) !=length(sigma) & length(sigma)!=1) & (ncol(x) !=length(cutoff) & length(cutoff)!=1)) error("one sigma per column or one sigma for all columns!")
  #if(is.null(sigma) & is.null(cutoff)) error("one of sigma and error shd be provided")
  n <- nrow(x)
  if(inducing){
    x2 <- xInd
    removeDiag <- FALSE
  } else{
    x2 <- x
    removeDiag <- TRUE
  }
  
  if(any(is.na(sigma))){
    
    #print("minDR"); print(minDR)
    #print("maxDR"); print(maxDR)
    sigmaAux <- sapply(which(is.na(sigma)), function(i){
      # i <- 2
      # fit sigma
      # print(paste("i: ", i))
      ord <- 10
      sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        #print(paste("sd: ", sd1))
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F], sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
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
      
      sigmas1 <- (10^seq(DR.sig[1],DR.sig[2], length.out=20)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F],sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
        res <- var(distsX)
        return(res)
      })
      indxMaxVar <- which.max(varsHsics.sig)
      sigmaVar <- sigmas1[indxMaxVar]
      
      # obtain dynamic range of sigma
      spl <- spline(log(sigmas1,10), varsHsics.sig, n=100*length(sigmas1))
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
      
      #minDR <- 0.2
      #maxDR <- 0.4
      
      indx <- which(spl$x>=DR.sig[1] & spl$x<DR.sig[2])
      xx <- spl$x[indx]
      yy <- spl$y [indx]
      indxMax <- which.max(yy)
      xx2 <- rep(NA, length(xx))
      xx2[1:indxMax] <- yy[1:indxMax]/(2*yy[indxMax])
      xx2[(indxMax+1):length(xx)] <- 0.5+(0.5-yy[(indxMax+1):length(xx)]/(2*max(yy[(indxMax+1):length(xx)])))
      #plot(xx2, xx)
      
      splf <- splinefun(c(0,xx2,1), c(xx[1],xx,xx[length(xx)]), method="monoH.FC")
      print(c(minDR[i], maxDR[i]))
      rngSigs <- splf(c(minDR[i], maxDR[i]))
      if(abs(rngSigs[2]-rngSigs[1])<1e-3){
        rngSigs[1] <- mean(rngSigs)-1e-3
        rngSigs[2] <- mean(rngSigs)+1e-3
      }
      
      res <- runif(1, min=rngSigs[1], max=rngSigs[2])
      res <- 10^res
      
      if(FALSE){
        plot(log(sigmas1,10), varsHsics.sig, col="red", type="p")
        lines(spl)
        abline(v=DR.sig, col="purple")
        abline(v=log(sigmaVar,10), col="blue")
        abline(h=max(varsHsics.sig), col="green")
        abline(v=rngSigs, col="orange", type="p")
      }
      
      return(res)
    })
    sigma[which(is.na(sigma))] <- sigmaAux 
  }
  
  
  sigma2 <- matrix(0, ncol(x), ncol(x))
  diag(sigma2) <- 1/sigma
  x2 <- x %*% sigma2
  x2 <- x2/sqrt(2)
  
  
  if(inducing){
    
    xInd2 <- xInd %*% sigma2
    xInd2 <- xInd2/sqrt(2)
    K <- kern_rbf(x2, xInd2, sigma=1)
    #alternative K calculation
    if(FALSE){
      
      K_alt <- sapply(1:length(sigma), function(i) kern_rbf(x=x[,i,drop=F], y=xInd[,i,drop=F], sigma=1/(2*sigma[i]^2)), simplify="array")
      K_alt <- apply(K_alt, c(1,2), prod)
      smpl1 <- sample(length(K), 100)
      plot(c(K)[smpl1], c(K_alt)[smpl1]); abline(a=0, b=1, col="red")
    }
    
    
  } else{
    K <- kern_rbf(x2, sigma=1) 
    n <- nrow(K)
    I <- diag(n)
    #C <- K + tau^2*I
    
    
  }
  
  #alpha <- matrix(abs(rnorm(ncol(K))), ncol(K), 1)
  alpha <- matrix(rnorm(ncol(K), mean=mean_alpha), ncol(K), 1)
  y <- K %*% alpha 
  
  # we want cov(y, y') = K(y, y') 
  # say y=t(C) %*% x and y'= t(C') %*% x', x,x'~ N(0,I) 
  # then cov(y, y') = cov(t(C)%*%x, t(C')%*%x')
  #                 = 
  
  
  return(y)  
}



SIMdist_Kalpha_sigPar <- function(n,  tau, sig_RDx=NA, sig_RDy=NA, sig_Ex=NA, sig_Ey=NA, sig_x=NA,  sig_RDz=NA, sig_Ez=NA, 
                              min_RDx=NULL, min_RDy=NULL, min_Ex=NULL, min_Ey=NULL, min_x=NULL, min_RDz=NULL, min_Ez=NULL,
                              max_RDx=NULL, max_RDy=NULL, max_Ex=NULL, max_Ey=NULL, max_x=NULL, max_RDz=NULL, max_Ez=NULL,
                              inducing=FALSE, numInducing=n, inducing_mean_x=0,
                              sig_RDxInd=sig_RDx, sig_RDyInd=sig_RDy, sig_ExInd=sig_Ex, sig_EyInd=sig_Ey, sig_xInd=sig_x, sig_RDzInd=sig_RDz, sig_EzInd=sig_Ez, 
                              min_RDxInd=min_RDx, min_RDyInd=min_RDy, min_ExInd=min_Ex, min_EyInd=min_Ey, min_xInd=min_x, min_RDzInd=min_RDz, min_EzInd=min_Ez,
                              max_RDxInd=max_RDx, max_RDyInd=max_RDy, max_ExInd=max_Ex, max_EyInd=max_Ey, max_xInd=max_x, max_RDzInd=max_RDz, max_EzInd=max_Ez,
                              sig_nois_x=NULL, sig_nois_y=NULL,
                              addNoise=TRUE, 
                              seedDist=NULL, seedFun=NULL, seedNois=NULL){
  # n <- 1000
  # sig_RDz=NULL; sig_Ez=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  # min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL
  # max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL
  if(is.null(sig_RDz)!=is.null(sig_Ez) & is.null(min_RDz)!=is.null(min_Ez) ) stop("if there is a confounder both sig_RDz and sig_Ez shoud be given otherwise both should be NULL")
  set.seed(seedDist)
  seedDist2 <- sample(1:10000, 7)
  set.seed(NULL)
  set.seed(seedNois)
  seedNois2 <- sample(1:10000, 8)
  set.seed(NULL)
  set.seed(seedFun)
  seedFun2 <- sample(10000,2)
  set.seed(NULL)
  
  Ex <- RDdist_GP(n, tau, sigma=sig_RDx, minDR=min_RDx, maxDR=max_RDx, seedDist=seedDist2[1], seedFun=seedDist2[2])
  Ex <- apply(Ex, 2, stdrize)
  Ey <- RDdist_GP(n, tau, sigma=sig_RDy, minDR=min_RDy, maxDR=max_RDy, seedDist=seedNois2[1], seedFun=seedNois2[2])
  Ey <- apply(Ey, 2, stdrize)
  
  Ex2 <- Ex
  noisNms <- c("Ex")
  sig_Ex2 <- sig_Ex
  min_Ex2 <- min_Ex
  max_Ex2 <- max_Ex
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    #print("conf")
    Ez <- RDdist_GP(n, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[3], seedFun=seedNois2[4])
    Ez <- apply(Ez, 2, stdrize)
    Ex2 <- cbind(Ex2, Ez)
    noisNms <- c(noisNms,"Ez")
    sig_Ex2 <- c(sig_Ex2, sig_Ez)
    min_Ex2 <- c(min_Ex2, min_Ez)
    max_Ex2 <- c(max_Ex2, max_Ez)
  }
  nois <- cbind(Ex2, Ey)
  noisNms <- c(noisNms, "Ey")
  
  x <- rGP_rbf(Ex2, tau, sigma=sig_Ex2, minDR=min_Ex2, maxDR=max_Ex2, seed=seedDist2[3])
  x <- apply(x, 2, stdrize)
  x2 <- cbind(x, Ey)
  
  sig_x2 <- c(sig_x, sig_Ey)
  min_x2 <- c(min_x, min_Ey)
  max_x2 <- c(max_x, max_Ey)
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    x2 <- cbind(x2, Ez)
    sig_x2 <- c(sig_x2, sig_Ez)
    min_x2 <- c(min_x2, min_Ez)
    max_x2 <- c(max_x2, max_Ez)
  }
  
  if(inducing){
    #numInducing <- 10
    ExInd <- RDdist_GP(numInducing, tau, sigma=sig_RDxInd, minDR=min_RDxInd, maxDR=max_RDxInd, seedDist=seedDist2[4], seedFun=seedDist2[5])
    ExInd <- apply(ExInd, 2, stdrize)
    
    EyInd <- RDdist_GP(numInducing, tau, sigma=sig_RDyInd, minDR=min_RDyInd, maxDR=max_RDyInd, seedDist=seedNois2[5], seedFun=seedNois2[6])
    EyInd <- apply(EyInd, 2, stdrize)
    
    Ex2Ind <- ExInd
    noisNmsInd <- c("ExInd")
    sig_Ex2Ind <- sig_ExInd
    min_Ex2Ind <- min_ExInd
    max_Ex2Ind <- max_ExInd
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      #print("conf")
      EzInd <- RDdist_GP(numInducing, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[7], seedFun=seedNois2[8])
      EzInd <- apply(Ez, 2, stdrize)
      Ex2Ind <- cbind(Ex2Ind, EzInd)
      noisNmsInd <- c(noisNmsInd,"EzInd")
      sig_Ex2Ind <- c(sig_Ex2Ind, sig_EzInd)
      min_Ex2Ind <- c(min_Ex2Ind, min_EzInd)
      max_Ex2Ind <- c(max_Ex2Ind, max_EzInd)
    }
    noisInd <- cbind(Ex2Ind, EyInd)
    noisNmsInd <- c(noisNmsInd, "EyInd")
    
    xInd <- rGP_rbf(Ex2Ind, tau, sigma=sig_Ex2Ind, minDR=min_Ex2Ind, maxDR=max_Ex2Ind, seed=seedDist2[6])
    xInd <- apply(xInd, 2, stdrize)
    xInd <- xInd + inducing_mean_x
    x2Ind <- cbind(xInd, EyInd)
    
    sig_x2Ind <- c(sig_xInd, sig_EyInd)
    min_x2Ind <- c(min_xInd, min_EyInd)
    max_x2Ind <- c(max_xInd, max_EyInd)
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      x2Ind <- cbind(x2Ind, EzInd)
      sig_x2Ind <- c(sig_x2Ind, sig_EzInd)
      min_x2Ind <- c(min_x2Ind, min_EzInd)
      max_x2Ind <- c(max_x2Ind, max_EzInd)
    }
  } else{
    xInd2 <- NULL
  }
  
  
  #x=x2; tau; xInd=x2Ind; inducing=inducing; sigma=sig_x2; minDR=min_x2; maxDR=max_x2; seed=seedFun2[1]
  y <- rKalpha_rbf(x=x2, tau, mean_alpha=-inducing_mean_x,xInd=x2Ind, inducing=inducing, sigma=sig_x2, minDR=min_x2, maxDR=max_x2, seed=seedFun2[1])
  set.seed(NULL)
  
  y <- apply(y, 2, stdrize)
  
  # add noise
  if(addNoise){
    set.seed(seedDist2[7])
    addNoiseX <- rnorm(n, sd=sig_nois_x)
    x <- x + addNoiseX
    set.seed(NULL)
    addNoiseY <- rnorm(n, sd=sig_nois_y) 
    y <- y + addNoiseY
    nois <- cbind(nois, addNoiseX, addNoiseY)
    noisNms <- c(noisNms, "addNoiseX","addNoiseY")
  }
  #plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y, nois)
  colnames(res) <- c("x","y", noisNms)
  res <- list(res=res)
  if(inducing){
    resInd <- cbind(xInd, noisInd)
    colnames(resInd) <- c("xInd", noisNmsInd)
    res <- c(res, list(resInd=resInd))
  }
  
  return(res)
}


SIMdist_Kalpha_sigPar_wrapper <- function(q=100, n=1000, seed=NULL, seed_rep=1, tau=NULL, type=c("SIM", "SIMc","SIMG","SIMln"), 
                                      calcStats=FALSE, nms=NULL, ...){
  #q=100; n=1000; seed=NULL; tau=NULL; sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL; sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  type <- match.arg(type, choices=c("SIM", "SIMc","SIMG","SIMln"))
  
  # pars <- list()
  pars <- list(...)
  print("extra pars: ")
  print(names(pars))
  pars$n <- n
  #n <- 100
  #type <- "SIM"
  
  
  pars <- switch(type, 
                 SIM={
                   set.seed(seed) 
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                   
                   
      
                 },
                 
                 SIMc={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDz"]]) & is.null(pars[["min_RDz"]])) pars$sig_RDz <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_Ez"]]) & is.null(pars[["min_Ez"]])) pars$sig_Ez <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                 },
                 SIMG={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                 },
                 SIMln={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=1.5*200), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                   
                   
                 })
  
  print("names(pars)")
  print(names(pars))
  
  pars <- lapply(1:q, function(i) lapply(pars, function(el) if(length(el)==1) el else el[i]))
  
  
  # par <- pars[[1]]
  if(FALSE){
    sig_RDx=NA; sig_RDy=NA; sig_Ex=NA; sig_Ey=NA; sig_x=NA;  sig_RDz=NA; sig_Ez=NA; 
    min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL;
    max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL;
    inducing=FALSE; numInducing=n; inducing_mean_x=0
    sig_nois_x=NULL; sig_nois_y=NULL;
    addNoise=TRUE; 
    seedDist=NULL; seedFun=NULL; seedNois=NULL
    indx <- setdiff(1:length(par),grep("ind",names(pars)))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
    sig_RDxInd=sig_RDx; sig_RDyInd=sig_RDy; sig_ExInd=sig_Ex; sig_EyInd=sig_Ey; sig_xInd=sig_x; sig_RDzInd=sig_RDz; sig_EzInd=sig_Ez; 
    min_RDxInd=min_RDx; min_RDyInd=min_RDy; min_ExInd=min_Ex; min_EyInd=min_Ey; min_xInd=min_x; min_RDzInd=min_RDz; min_EzInd=min_Ez;
    max_RDxInd=max_RDx; max_RDyInd=max_RDy; max_ExInd=max_Ex; max_EyInd=max_Ey; max_xInd=max_x; max_RDzInd=max_RDz; max_EzInd=max_Ez;
    indx <- grep("ind",names(pars))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
  }
  # dat <- do.call("SIMdist_Kalpha_sigPar", par)
  
  
  # plot(dat$res[,"x"], dat$res[,"y"])
  
  if(!is.null(nms)) nms <- paste(type, 1:q, nms, sep=".") else nms <- paste(type, 1:q, sep=".")
  
  mc_cores <- detectCores()/2
  #mc_cores <- 1
  print("xs")
  pm <- proc.time()
  xs <- mcmapply(function(par, nm){
    # par <- pars[[64]]
    print(paste("nm: ", nm))
    res <- do.call("SIMdist_Kalpha_sigPar", par)
    return(res)
  }, par=pars, nm=nms, mc.cores=mc_cores, SIMPLIFY=FALSE)
  proc.time() - pm 
  # 1.241 secs for q=100, n=100 and 4 cores with confounder (3 gps)
  # 7.15 mins for q=100, n=1000 and 4 cores with confounder (3 gs)
  print("xs2")
  xs2 <- lapply(xs, function(el) el$res[,c("x","y")])
  print("ind_xs")
  ind_xs <- lapply(xs, function(el) el$resInd)
  print("dags,ns")
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs, function(el) el$res[,-which(colnames(el$res) %in% c("x","y"))])
  names(xs2) <- names(ns) <- names(dags) <- names(ind_xs) <- nms
  print("dagList")
  dataList <- list(xs=xs2, ns=ns, dags=dags, names=nms, ind_xs=ind_xs)
  
  if(calcStats){
    print("calc stats")
    pm <- proc.time()
    count <- 0
    descr_indices <- sapply(xs2, function(xs){
      count <<- count + 1
      # count <- 1;  xs <- xs2[[count]]
      print("*****************")
      print(count)
      pm <- proc.time()
      num_obs <- 100
      num_reps <- floor(nrow(xs))/num_obs
      pm <- proc.time()
      res <- mcmapply(function(j){
        #j <- 1
        #print(paste("j: ",j))
        
        # plot(xs[,"x"], xs[,"y"])
        # hist(xs[,"x"])
        # hist(xs[,"y"])
        smpl <- sample(1:nrow(xs), num_obs)
        df <- as.data.frame(xs[smpl,])
        mod <- lm(y~x,df)
        #plot(df$x, residuals(mod))
        #hist(residuals(mod))
        xs_2 <- apply(xs[smpl,], 2, norml)
        #plot(xs_2[,"x"], xs_2[,"y"])
        trDat <- constructData(x=as.matrix(xs_2[,"x"]), y=as.matrix(xs_2[,"y"]))
        krrAux <- setParams(learner=krr1, trainData=trDat)
        krrAux <- krrAux$learn(krrAux)
        predKrr <- pred.CV(krrAux, trDat)$test
        #o <- order(trDat$x); plot(trDat$x, trDat$y); lines(trDat$x[o], predKrr$gyh_class[o], col="red")
        residsLin <- residuals(mod)
        # plot(xs[,"x"], residsLin)
        residsNonLin <- krrAux$resids(krrAux, predKrr)[,"resid"]
        #plot(predKrr$x_class, residsNonLin)
        #plot(o, residsNonLin[o])
        res1 <- dhsic.test(X=df[,"x",drop=F],Y=matrix(residsLin,ncol=1))$p.value
        res2 <- dhsic.test(X=predKrr$x_class,Y=residsNonLin)$p.value
        #o <- order(predKrr$x_class)
        #res2b <- dhsic.test(X=o,Y=residsNonLin[o])$p.value
        res3 <- ks.test(x=residsLin, y="pnorm", mean=mean(residsLin), sd=sd(residsLin))$p.value
        res4 <- ks.test(x=residsNonLin, y="pnorm", mean=mean(residsNonLin), sd=sd(residsNonLin))$p.value
        res5 <- ks.test(x=xs_2[,"x"], y="pnorm", mean=mean(xs_2[,"x"]), sd=sd(xs_2[,"x"]))$p.value
        res6 <- ks.test(x=xs_2[,"x"], y="punif", min=min(xs_2[,"x"]), max=max(xs_2[,"x"]))$p.value
        res7 <- Shannon_KDP(xs_2[,"x"])
        
        modDens <- kepdf(xs_2[,"x"], eval.points = xs_2[,"x"], kernel = "gaussian", bwtype = "adaptive")
        #hist(xs_2[,"x"], prob=T); o <- order(modDens@eval.points); lines(modDens@eval.points[o], modDens@estimate[o], col="red")
        res8 <- max(modDens@estimate)-min(modDens@estimate)
        res <- c(res1, res2, res3, res4, res5, res6, res7, res8)
        names(res) <- c("lm_indep","add_indep","lm_gauss","add_gauss",
                        "cause_gauss","cause_unif", "cause_ent","cause_rngPdf")
        #print(res)
        proc.time() - pm 
        return(res)
      }, SIMPLIFY="array", j=1:num_reps, mc.cores=detectCores()/2)
      proc.time()- pm #37.5  with 1 core, 22.5 with 4 cores (100pts)
      res <- apply(res, 1, mean)
      return(res)
    })
    proc.time() - pm #34 mins
    descr_indices <- t(descr_indices)
    dataList <- c(dataList, list(descr_indices=descr_indices))
  }
  
  
  #plotPairsList(dataList)
  return(dataList)  
}




numAlpha <- 7
numPerAlpha <- 100
q <- numAlpha*numPerAlpha
alpha <- seq(numAlpha)

# sig_Ey
min_min_par <- 0.6; min_max_par <- 0.99
max_min_par <- 0.61; max_max_par <- 1
#min_min_par <- 0.8; min_max_par <- 0.85
#max_min_par <- 0.81; max_max_par <- 0.86
min_par <- runif(numPerAlpha, min_min_par, min_max_par)
min_par <- rep(min_par, numAlpha)
max_par <- min_par+0.01
min_par; max_par
#min_par <- rep(min_par, rep(numPerAlpha,  numAlpha))
#max_par <- rep(max_par, rep(numPerAlpha,  numAlpha))
# sig_x
min_min_par2 <- 0.35; min_max_par2 <- 0.41
max_min_par2 <- 0.36; max_max_par2 <- 0.42
max_par2 <- seq(max_max_par2, max_min_par2, length.out=numAlpha)
min_par2 <- seq(min_max_par2, min_min_par2, length.out=numAlpha)
min_par2; max_par2
min_par2 <- rep(min_par2, rep(numPerAlpha, numAlpha))
max_par2 <- rep(max_par2, rep(numPerAlpha, numAlpha))

# X~rgmamma(shape, scale); E[X]= shape*scale, V[X]=shape*scale^2
# sig_RDx, sig_Ex

Ex <- RDdist_GP(n=100, tau=1e-4, sigma=0.5); hist(Ex)
x <- rGP_rbf(Ex, tau=1e-4, sigma=30); hist(x)

Ex <- RDdist_GP(n=100, tau=1e-4, sigma=1e9); hist(Ex)
x <- rGP_rbf(Ex, tau=1e-4, sigma=1e9); hist(x); ks.test(x, y="")


set.seed(12)
min_shape <- 5; max_shape <- 1e6
seq_shape <- 10^seq(log(min_shape,10), log(max_shape,10), length.out=numAlpha)
min_scale <- 0.1; max_scale <- 1e-3
seq_scale <- 10^seq(log(min_scale,10), log(max_scale,10), length.out=numAlpha)
seq_shape*seq_scale
seq_shape*seq_scale^2
seq_shape <- rep(seq_shape, rep(numPerAlpha, numAlpha))
seq_scale <- rep(seq_scale, rep(numPerAlpha, numAlpha))
sigAlpha_RDx <-  rgamma(q, shape=seq_shape, scale=seq_scale)
#min_shape <- 2; max_shape <- 1e6 # real SIM
min_shape <- 10; max_shape <- 1e5
seq_shape <- 10^seq(log(min_shape,10), log(max_shape,10), length.out=numAlpha)
#min_scale <- 1.5; max_scale <- 1e-3 # real SIM
min_scale <- 0.15; max_scale <- 1e-3
seq_scale <- 10^seq(log(min_scale,10), log(max_scale,10), length.out=numAlpha)
seq_shape*seq_scale
seq_shape*seq_scale^2
seq_shape <- rep(seq_shape, rep(numPerAlpha, numAlpha))
seq_scale <- rep(seq_scale, rep(numPerAlpha, numAlpha))
sigAlpha_Ex <-  rgamma(q, shape=seq_shape, scale=seq_scale)





inducing_mean_x <- c(-2,-4,2,4)#-2#c(-10,-4,-2,-1.5,-1,-0.5,0)
#inducing_mean_x <- rep(inducing_mean_x, rep(numPerAlpha,  numAlpha))
inducing_mean_x <- sample(inducing_mean_x, q, replace=T)#rep(inducing_mean_x, q)
inducing_mean_x <- inducing_mean_x
signInd <- sample(c(-1,1), numPerAlpha, replace=T)
signInd <- rep(signInd, numAlpha)
inducing_mean_x <- inducing_mean_x*signInd

set.seed(1234) #set.seed(123)
seedDist <- sample(100000, numPerAlpha); seedDist <- rep(seedDist, numAlpha)
seedFun <- sample(100000, numPerAlpha) ; seedFun <- rep(seedFun, numAlpha)
seedNois <- sample(100000, numPerAlpha); seedNois <- rep(seedNois, numAlpha)
set.seed(NULL)

nms <- paste(rep(1:numPerAlpha, numAlpha), rep(alpha, rep(numPerAlpha, numAlpha)), sep="_")
seed <- 12345 #12
seed_rep <- numAlpha
addNoise <- TRUE
inducing <- TRUE
numInducing <- sample(c(100,100), numPerAlpha*numAlpha, replace=T)
calcStats <- FALSE
n <- 1000
#type <- "SIMG"
type <- "SIM"
#, sig_Ex=sig_Ex, sig_RDx=sig_RDx
# pars <- list(addNoise=addNoise, inducing=inducing, numInducing=numInducing, min_Ey=min_par, max_Ey=max_par,min_x=min_par2, max_x=max_par2, seedDist=seedDist, seedFun=seedFun, seedNois=seedNois, inducing_mean_x=inducing_mean_x)
#SIMdist_GP_sigPar_wrapper <- function(q=100, n=1000, seed=NULL, tau=NULL, type=c("SIM", "SIMc","SIMG","SIMln"), 
#                                      calcStats=FALSE, nms=NULL, ...)
# 
dataList <- SIMdist_Kalpha_sigPar_wrapper(q=q, n=n, seed=seed, seed_rep=seed_rep, type=type, calcStats=calcStats, nms=nms, addNoise=addNoise, inducing=inducing, numInducing=numInducing, min_x=min_par2, max_x=max_par2, seedDist=seedDist, seedFun=seedFun, seedNois=seedNois, inducing_mean_x=inducing_mean_x,min_Ey=min_par, max_Ey=max_par, sig_RDx = sigAlpha_RDx, sig_Ex=sigAlpha_Ex)
# version with no seed control so that every function is random
#dataList <- SIMdist_Kalpha_sigPar_wrapper(q=numPerAlpha*numAlpha, n=n, type=type, calcStats=calcStats, nms=nms, addNoise=addNoise, inducing=inducing, numInducing=numInducing, min_Ey=min_par, max_Ey=max_par,min_x=min_par2, max_x=max_par2, inducing_mean_x=inducing_mean_x)

plotPairsList_alphaProg(dataList, reps=c(1,11,21,31,41,51,61,71,81,91))
  

dataList$names
data <- melt(dataList$xs)
colnames(data) <- c("numData","xY","val","dataset")
data <- cast(data, dataset+numData~xY, value="val")
colnames(data) <- c("dataset","numData","x","y")
head(data)
aux <- strsplit(sapply(strsplit(data$dataset, "\\."), function(el) el[3]), "_")
data$alpha <- as.numeric(sapply(aux, function(el) el[2]))
data$rep <- as.numeric(sapply(aux, function(el) el[1]))
table(data$alpha)
table(data$rep)

myKS <- function(x){
  #print(mean(x))
  #print(sd(x))
  ks.test(x, y="pnorm", mean=0, sd=1)$p.value
} 
data2 <- as.data.frame(data)
data2$var <- data2$x
data2$dataset <- factor(data2$dataset)
# ks.test form normality
library(reshape)
table(data2$dataset)
ksTest <- reshape:::cast(data, dataset~., fun.aggregate="myKS", value="val")

ksTest <- aggregate(data2$x, by=list(data2$dataset), FUN="myKS")
aux <- strsplit(sapply(strsplit(as.character(ksTest$Group.1), "\\."), function(el) el[3]), "_")
ksTest$rep <- as.numeric(sapply(aux, function(el) el[1]))
ksTest$alpha <- as.numeric(sapply(aux, function(el) el[2]))
table(ksTest$rep); table(ksTest$alpha)
p <- ggplot(ksTest)
p <- p + geom_boxplot(aes(x=factor(alpha), y=x))
p

# histograms
p <- ggplot(data[which(data$rep<=10),])
p <- p + geom_histogram(aes(x=x,y=..density..), position="identity", bins=20, colour="light grey") 
p <- p + geom_density(aes(x=x,y=..density..))
p <- p + facet_grid(alpha~rep)
p

dataList$names
i <- 1
repAlphas <- paste(i,1:numAlpha, sep="_"); 
reps <- (1:numAlpha-1)*(numPerAlpha)+i; 
nms_i <- paste("SIM",repAlphas,reps,sep=".") 
par(mfrow=c(3,3)); for(j in 1:numAlpha) hist(dataList$xs[[nms_i[j]]][,1]); par(mfrow=c(1,1))
i <- i + 1

i <- 1
repAlphas <- paste(i,1:numAlpha, sep="_"); 
reps <- (1:numAlpha-1)*(numPerAlpha)+i; 
nms_i <- paste("SIM",repAlphas,reps,sep=".") 
par(mfrow=c(1,2))
hist(dataList$xs[[nms_i[1]]][,1]) 
hist(dataList$xs[[nms_i[numAlpha]]][,1]) 
par(mfrow=c(1,1))
i <- i + 1