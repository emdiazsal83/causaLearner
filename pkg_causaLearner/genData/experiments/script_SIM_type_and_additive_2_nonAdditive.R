remove(list=ls())

library(mvtnorm)

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

# simulate SIM, "methods and benchmarks" (Mooij et al 2016) style data

# first the sampling procedure "RD"

x <- c(1,2,3)
sigma <- 0.7
kern_rbf(x=matrix(x), sigma=1/(2*sigma^2))
exp(-matNorm2(matrix(x)*1/sigma)/2)

n <- 10; m <- 100000
sigma <- 3
tau <- 0.1
x <- matrix(rnorm(n, sd=10), n, 1)
K <- kern_rbf(x=x, sigma=sigma)
alpha <- matrix(rnorm(n*m, sd=tau), n, m)
x <- K %*% alpha
I <- diag(n)
x2 <- t(rmvnorm(m, mean=rep(0,n), sigma=K+tau^2*I)) 
dim(x2)
plot(c(cor(t(x))), c(cor(t(x2)))); abline(a=0, b=1, col="red")



rGP_rbf <- function(x, tau, sigma=NULL, minRD=NULL, maxRD=NULL,  seed=NULL, seedInd=NULL){
  
  if(!is.null(seedInd)){
    set.seed(seedInd)
    xInd <- matrix(rnorm( (n-1)*ncol(x)), n-1, ncol(x))
    set.seed(NULL)
    x2 <- xInd
    removeDiag <- FALSE
  } else{
    x2 <- x
    removeDiag <- TRUE
  }
  
  if(is.null(sigma)){
      
    sigma <- sapply(1:length(minRD), function(i){
      # i <- 1
      # fit sigma
      # print(paste("i: ", i))
      ord <- 10
      sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
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
      spl <- spline(log(sigmas1,10), varsHsics.sig, n=10*length(sigmas1))
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
    
      #cutoff <- 0.5
      maxVar <-max(varsHsics.sig)*maxRD[i]
      minVar <-max(varsHsics.sig)*minRD[i]
      sigs <- spl$x[which(spl$y>maxVar)]
      vars-sigs
      res <- runif(1, min=rngSigs[1], max=rngSigs[2])
      res <- 10^res
      
      if(FALSE){
        plot(log(sigmas1,10), varsHsics.sig, col="red", type="p")
        lines(spl)
        abline(v=DR.sig, col="purple")
        abline(v=log(sigmaVar,10), col="blue")
        abline(h=max(varsHsics.sig), col="green")
      }
      
      return(res)
    })
  }
  
  
  sigma2 <- matrix(0, ncol(x), ncol(x))
  diag(sigma2) <- 1/sigma
  x2 <- x %*% sigma2
  x2 <- x2/sqrt(2)
  n <- nrow(x)
  
  
  if(!is.null(seedInd)){
    
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
    svdK <- svd(K) 
    set.seed(seed)
    y <- rnorm(nrow(x)-1, sd=sqrt(svdK$d)+tau)
    set.seed(NULL)
    y < - matrix(y, length(y), 1)
    #print("dim(y)"); print(dim(y))
    #print("dim(svdK$u)"); print(dim(svdK$u))
    y <- svdK$u %*%  y
    
  } else{
    K <- kern_rbf(x2, sigma=1) 
    n <- nrow(K)
    I <- diag(n)
    C <- K + tau^2*I
    if(!is.positive.definite(C)) C <- make.positive.definite((C))
    C  <- chol(C)
    set.seed(seed)
    y <- rnorm(nrow(x))
    set.seed(NULL)
    y <- t(C) %*% matrix(y,n,1)
    
  }
  
  # we want cov(y, y') = K(y, y') 
  # say y=t(C) %*% x and y'= t(C') %*% x', x,x'~ N(0,I) 
  # then cov(y, y') = cov(t(C)%*%x, t(C')%*%x')
  #                 = 
  
  
  return(y)  
}


RDdist_GP <- function(n, tau, sigma=NULL, cutoff=NULL, seedDist=NULL, seedFun=NULL){
  # n <- 100
  # sigma  <- rgamma(1, shape=5, scale=0.1)
  # tau <- 1e-4
  set.seed(seedDist)
  x <- matrix(rnorm(n))
  smpl <-sample(nrow(x))
  set.seed(NULL)
  F <- rGP_rbf(x, tau, sigma, cutoff, seed=seedFun)
  eF <- exp(F)
  o <- order(x)
  G <- cumtrapz(x[o], eF[o])
  G <- matrix(G[smpl])
  return(G)
}


hist(RDdist_GP(sigma=1000, tau=1e-4, n=1000))



SIMdist_GP_sigPar <- function(n,  tau, sig_RDx=NULL, sig_RDy=NULL, sig_Ex=NULL, sig_Ey=NULL, sig_x=NULL, sig_nois_x=NULL, sig_nois_y=NULL, sig_RDz=NULL, sig_Ez=NULL, 
                                       cof_RDx=NULL, cof_RDy=NULL, cof_Ex=NULL, cof_Ey=NULL, cof_x=NULL,                                   cof_RDz=NULL, cof_Ez=NULL,
                                       addNoise=TRUE, seedDist=NULL, seedFun=NULL, seedNois=NULL){
  # n <- 1000
  # sig_RDz=NULL; sig_Ez=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  if(is.null(sig_RDz)!=is.null(sig_Ez) & is.null(cof_RDz)!=is.null(cof_Ez) ) stop("if there is a confounder both sig_RDz and sig_Ez shoud be given otherwise both should be NULL")
  set.seed(seedDist)
  seedDist2 <- sample(1:10000, 4)
  set.seed(NULL)
  Ex <- RDdist_GP(n, tau, sigma=sig_RDx, cutoff=cof_RDx, seedDist=seedDist2[1], seedFun=seedDist2[2])
  Ex <- apply(Ex, 2, stdrize)
  
  
  set.seed(seedNois)
  seedNois2 <- sample(1:10000, 4)
  set.seed(NULL)
  Ey <- RDdist_GP(n, tau, sigma=sig_RDy, cutoff=cof_RDy, seedDist=seedNois2[1], seedFun=seedNois2[2])
  Ey <- apply(Ey, 2, stdrize)
  
  Ex2 <- Ex
  noisNms <- c("Ex")
  sig_Ex2 <- sig_Ex
  cof_Ex2 <- cof_Ex
  
  if(!is.null(sig_RDz) | !is.null(cof_RDz)){
    #print("conf")
    
    Ez <- RDdist_GP(n, tau, sigma=sig_RDz, cutoff=cof_RDz, seedDist=seedNois2[3], seedFun=seedNois2[4])
    Ez <- apply(Ez, 2, stdrize)
    Ex2 <- cbind(Ex2, Ez)
    noisNms <- c(noisNms,"Ez")
    sig_Ex2 <- c(sig_Ex2, sig_Ez)
    cof_Ex2 <- c(cof_Ex2, cof_Ez)
  }
  nois <- cbind(Ex2, Ey)
  noisNms <- c(noisNms, "Ey")
  
  
  x <- rGP_rbf(Ex2, tau, sigma=sig_Ex2, cutoff=cof_Ex2, seed=seedDist2[3])
  x <- apply(x, 2, stdrize)
  x2 <- cbind(x, Ey)
  sig_x2 <- c(sig_x, sig_Ey)
  cof_x2 <- c(cof_x, cof_Ey)
  
  if(!is.null(sig_RDz) | !is.null(cof_RDz)){
    x2 <- cbind(x2, Ez)
    sig_x2 <- c(sig_x2, sig_Ez)
    cof_x2 <- c(cof_x2, cof_Ez)
  }
  
  set.seed(seedFun)
  seedFun2 <- sample(10000,2)
  #x=x2; tau; sigma=sig_x2; cutoff=cof_x2; seed=seedFun2[1]; seedInd=seedFun2[2]
  y <- rGP_rbf(x=x2, tau, sigma=sig_x2, cutoff=cof_x2, seed=seedFun2[1], seedInd=seedFun2[2])
  set.seed(NULL)
  
  y <- apply(y, 2, stdrize)
  
  # add noise
  if(addNoise){
    set.seed(seedDist2[4])
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
  return(res)
}

if(FALSE){
sig_RDx <- rgamma(1, shape=5, scale=0.1)
sig_RDy <- rgamma(1, shape=5, scale=0.1)
sig_RDz <- rgamma(1, shape=5, scale=0.1)
sig_Ex <- rgamma(1, shape=2, scale=1.5)
sig_x <- rgamma(1, shape=2, scale=3)
sig_Ey <- rgamma(1, shape=2, scale=15)
sig_Ez <- rgamma(1, shape=2, scale=15)
sig_nois_x <- rgamma(1, shape=2, scale=0.1)
sig_nois_y <- rgamma(1, shape=2, scale=0.1)
tau <- 1e-4
head(SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=NULL))
head(SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=sig_RDz, sig_Ez=sig_Ez))

# constant dist
hist(SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=NULL, seedDist=1234)[,"x"])
head(SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=NULL, seedDist=1234))

# constant dist & noise
head(SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=sig_RDz, sig_Ez=sig_Ez, seedDist=1234, seedNois=123))

# non constant dist & noise, but constant function
dat <- SIMdist_GP_sigPar(n,  tau, sig_RDx, sig_RDy, sig_Ex, sig_Ey, sig_x, sig_nois_x, sig_nois_y, sig_RDz=NULL, addNoise=TRUE,seedDist=NULL, seedFun=123, seedNois=2)
head(dat)
plot(dat[,"x"],dat[,"y"])
}

SIMdist_GP_sigPar_wrapper <- function(q=100, n=1000, seed=NULL, tau=NULL, type=c("SIM", "SIMc","SIMG","SIMln"),
                                      sig_RDx=NULL, sig_RDy=NULL, sig_RDz=NULL, sig_Ex=NULL, sig_Ey=NULL, sig_Ez=NULL, sig_x=NULL,
                                      cof_RDx=NULL, cof_RDy=NULL, cof_RDz=NULL, cof_Ex=NULL, cof_Ey=NULL, cof_Ez=NULL, cof_x=NULL,
                                      sig_nois_x=NULL, sig_nois_y=NULL, addNoise=TRUE, 
                                      seedDist=NULL, seedFun=NULL, seedNois=NULL, calcStats=FALSE, nms=NULL){
  #q=100; n=1000; seed=NULL; tau=NULL; sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL; sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  type <- match.arg(type, choices=c("SIM", "SIMc","SIMG","SIMln"))
  
  #n <- 100
  #type <- "SIM"
  pars <- switch(type, 
                 SIM={
                   set.seed(seed)
                   if(is.null(sig_RDx) & is.null(cof_RDx)) sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   if(is.null(sig_RDy) & is.null(cof_RDy)) sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   sig_RDz <- NULL
                   if(is.null(sig_Ex) & is.null(cof_Ex)) sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   if(is.null(sig_x) & is.null(cof_x)) sig_x <- rgamma(q, shape=2, scale=3)
                   if(is.null(sig_Ey) & is.null(cof_Ey)) sig_Ey <- rgamma(q, shape=2, scale=15)
                   sig_Ez <- NULL
                   if(is.null(sig_nois_x)) sig_nois_x <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(sig_nois_y)) sig_nois_y <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(tau) & is.null(tau)) tau <- 1e-4
                   if(is.null(addNoise)) addNoise=TRUE
                
                   pars <- list(n=n,  tau=tau, 
                                sig_RDx=sig_RDx, sig_RDy=sig_RDy, sig_Ex=sig_Ex, sig_Ey=sig_Ey, sig_x=sig_x, sig_RDz=sig_RDz, sig_Ez=sig_Ez,
                                cof_RDx=cof_RDx, cof_RDy=cof_RDy, cof_Ex=cof_Ex, cof_Ey=cof_Ey, cof_x=cof_x, cof_RDz=cof_RDz, cof_Ez=cof_Ez,
                                sig_nois_x=sig_nois_x, sig_nois_y=sig_nois_y, 
                                 addNoise=addNoise, 
                                seedDist=seedDist, seedFun=seedFun, seedNois=seedFun)
                   set.seed(NULL)
                   pars
                   },
                 SIMc={
                   set.seed(seed)
                   if(is.null(sig_RDx) & is.null(cof_RDx)) sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   if(is.null(sig_RDy) & is.null(cof_RDy)) sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   if(is.null(sig_RDz) & is.null(cof_RDz)) sig_RDz <- rgamma(q, shape=5, scale=0.1)
                   if(is.null(sig_Ex) & is.null(cof_Ex)) sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   if(is.null(sig_x) & is.null(cof_x)) sig_x <- rgamma(q, shape=2, scale=3)
                   if(is.null(sig_Ey) & is.null(cof_Ey)) sig_Ey <- rgamma(q, shape=2, scale=15)
                   if(is.null(sig_Ez) & is.null(cof_Ez)) sig_Ez <- rgamma(q, shape=2, scale=15)
                   if(is.null(sig_nois_x)) sig_nois_x <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(sig_nois_y)) sig_nois_y <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(tau)) tau <- 1e-4
                   if(is.null(addNoise)) addNoise=TRUE
                   pars <- list(n=n,  tau=tau, 
                                sig_RDx=sig_RDx, sig_RDy=sig_RDy, sig_Ex=sig_Ex, sig_Ey=sig_Ey, sig_x=sig_x, sig_RDz=sig_RDz, sig_Ez=sig_Ez,
                                cof_RDx=cof_RDx, cof_RDy=cof_RDy, cof_Ex=cof_Ex, cof_Ey=cof_Ey, cof_x=cof_x, cof_RDz=cof_RDz, cof_Ez=cof_Ez,
                                sig_nois_x=sig_nois_x, sig_nois_y=sig_nois_y, 
                                addNoise=addNoise, 
                                seedDist=seedDist, seedFun=seedFun, seedNois=seedFun)
                   set.seed(NULL)
                   pars
                   },
                 SIMG={
                   set.seed(seed)
                   if(is.null(sig_RDx) & is.null(cof_RDx)) sig_RDx <- rgamma(q, shape=1e6, scale=1e-3)
                   if(is.null(sig_RDy) & is.null(cof_RDy)) sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   sig_RDz <- NULL
                   if(is.null(sig_Ex) & is.null(cof_Ex)) sig_Ex <- rgamma(q, shape=1e6, scale=1e-3)
                   if(is.null(sig_x) & is.null(cof_x)) sig_x <- rgamma(q, shape=2, scale=3)
                   if(is.null(sig_Ey) & is.null(cof_Ey)) sig_Ey <- rgamma(q, shape=2, scale=15)
                   sig_Ez <- NULL
                   if(is.null(sig_nois_x)) sig_nois_x <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(sig_nois_y)) sig_nois_y <- rgamma(q, shape=2, scale=0.1)
                   if(is.null(tau)) tau <- 1e-4
                   if(is.null(addNoise)) addNoise=TRUE
                   pars <- list(n=n,  tau=tau, 
                                sig_RDx=sig_RDx, sig_RDy=sig_RDy, sig_Ex=sig_Ex, sig_Ey=sig_Ey, sig_x=sig_x, sig_RDz=sig_RDz, sig_Ez=sig_Ez,
                                cof_RDx=cof_RDx, cof_RDy=cof_RDy, cof_Ex=cof_Ex, cof_Ey=cof_Ey, cof_x=cof_x, cof_RDz=cof_RDz, cof_Ez=cof_Ez,
                                sig_nois_x=sig_nois_x, sig_nois_y=sig_nois_y, 
                                addNoise=addNoise, 
                                seedDist=seedDist, seedFun=seedFun, seedNois=seedFun)
                   set.seed(NULL)
                   pars
                 },
                 SIMln={
                   set.seed(seed)
                   if(is.null(sig_RDx) & is.null(cof_RDx)) sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   if(is.null(sig_RDy) & is.null(cof_RDy)) sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   sig_RDz <- NULL
                   if(is.null(sig_Ex) & is.null(cof_Ex)) sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   if(is.null(sig_x) & is.null(cof_x)) sig_x <- rgamma(q, shape=2, scale=3)
                   if(is.null(sig_Ey) & is.null(cof_Ey)) sig_Ey <- rgamma(q, shape=2, scale=1.5*200)
                   sig_Ez <- NULL
                   if(is.null(sig_nois_x)) sig_nois_x <- rgamma(q, shape=2, scale=0.01)
                   if(is.null(sig_nois_y)) sig_nois_y <- rgamma(q, shape=2, scale=0.01)
                   if(is.null(tau) ) tau <- 1e-4
                   if(is.null(addNoise)) addNoise=TRUE
                   pars <- list(n=n,  tau=tau, 
                                sig_RDx=sig_RDx, sig_RDy=sig_RDy, sig_Ex=sig_Ex, sig_Ey=sig_Ey, sig_x=sig_x, sig_RDz=sig_RDz, sig_Ez=sig_Ez,
                                cof_RDx=cof_RDx, cof_RDy=cof_RDy, cof_Ex=cof_Ex, cof_Ey=cof_Ey, cof_x=cof_x, cof_RDz=cof_RDz, cof_Ez=cof_Ez,
                                sig_nois_x=sig_nois_x, sig_nois_y=sig_nois_y, 
                                addNoise=addNoise, 
                                seedDist=seedDist, seedFun=seedFun, seedNois=seedFun)
                   set.seed(NULL)
                   pars
                 })
  
  pars <- lapply(1:q, function(i) lapply(pars, function(el) if(length(el)==1) el else el[i]))
  
  
  # par <- pars[[1]]
  # dat <- do.call("SIMdist_GP_sigPar", par)
  # n <- 100; tau <- 1e-4; sig_RDx <- 0.2303679; sig_RDy <- 0.3757372; sig_Ex <- 2.700093; sig_Ey <- NULL; sig_x <- NULL; sig_RDz <- NULL; sig_Ez <- NULL;  cof_RDx <- NULL; cof_RDy <- NULL; cof_Ex <- NULL; cof_Ey <- 0; cof_x <- 0; cof_RDz <- NULL; cof_Ez <- NULL; sig_nois_x <- 0.26479; sig_nois_y <- 0.1801949; addNoise <- FALSE; seedDist <- 51663; seedFun <- 32953; seedNois <- 32953   
  
  # plot(dat[,"x"], dat[,"y"])
  
  if(!is.null(nms)) nms <- paste(type, nms, 1:q, sep=".") else nms <- paste(type, 1:q, sep=".")
  
  mc_cores <- detectCores()/2
  # mc_cores <- 1
  pm <- proc.time()
  xs <- mcmapply(function(par, nm){
    print(paste("nm: ", nm))
    res <- do.call("SIMdist_GP_sigPar", par)
    return(res)
    }, par=pars, nm=nms, mc.cores=mc_cores, SIMPLIFY=FALSE)
  proc.time() - pm 
  # 1.241 secs for q=100, n=100 and 4 cores with confounder (3 gps)
  # 7.15 mins for q=100, n=1000 and 4 cores with confounder (3 gs)
  xs2 <- lapply(xs, function(el) el[,c("x","y")])
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs2, function(el) el[,-which(colnames(el) %in% c("x","y"))])
  names(xs2) <- names(ns) <- names(dags) <- nms
  dataList <- list(xs=xs2, ns=ns, dags=dags, names=nms)
  
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

#load real SIM data
load("/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/data/TCEPs/dag2-ME2-SIM-1000_sims.RData")
dataList_SIM <- list(xs=dataList$xs[1:100], noiss=dataList$noiss[1:100], dags=dataList$dags[1:100], names=names(dataList$xs[1:100]))
dataList_SIMc <- list(xs=dataList$xs[101:200], noiss=dataList$noiss[101:200], dags=dataList$dags[101:200], names=names(dataList$xs[101:200]))
dataList_SIMG <- list(xs=dataList$xs[201:300], noiss=dataList$noiss[201:300], dags=dataList$dags[201:300], names=names(dataList$xs[201:300]))
dataList_SIMln <- list(xs=dataList$xs[301:400], noiss=dataList$noiss[301:400], dags=dataList$dags[301:400], names=names(dataList$xs[301:400]))

plotPairsList(dataList_SIM)
plotPairsList(dataList_SIMc)
plotPairsList(dataList_SIMG)
plotPairsList(dataList_SIMln)


dataList_mySIM <- SIMdist_GP_sigPar_wrapper(q=100, n=100, type="SIM", calcStats=FALSE)
dataList_mySIMc <- SIMdist_GP_sigPar_wrapper(q=100, n=1000, type="SIMc", calcStats=FALSE)
dataList_mySIMG <- SIMdist_GP_sigPar_wrapper(q=100, n=1000, type="SIMG", calcStats=FALSE)
dataList_mySIMln <- SIMdist_GP_sigPar_wrapper(q=100, n=1000, type="SIMln", calcStats=FALSE)
plotPairsList(dataList_mySIM)
plotPairsList(dataList_mySIMc)
plotPairsList(dataList_mySIMG)
plotPairsList(dataList_mySIMln)


###############################################################3
# Additivity experiment

# sig_RD_y params
mean(rgamma(10000, shape=2, scale=15)); 2*15
sd(rgamma(10000, shape=2, scale=15)); sqrt(2*15^2)

# sig_x params
mean(rgamma(10000, shape=2, scale=3)); 2*3
sd(rgamma(10000, shape=2, scale=3)); sqrt(2*3^2)

q=100; n=100; seed=NULL; tau=NULL; type=c("SIM", "SIMc","SIMG","SIMln")
sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL
cof_RDx=NULL; cof_RDy=NULL; cof_RDz=NULL; cof_Ex=NULL; cof_Ey=NULL; cof_Ez=NULL; cof_x=NULL
sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE 
seedDist=NULL; seedFun=NULL; seedNois=NULL; calcStats=FALSE; nms=NULL

numAlpha <- 6
numPerAlpha <- 100
alpha <- seq(numAlpha)
min_cof_x <- 0
max_cof_x <- 0.2
#(sig_x <- 10^seq(min_sig_x, max_sig_x, length.out=numAlpha))
#(sig_Ey <- (10^seq(max_sig_x, min_sig_x, length.out=numAlpha))*10)
#sig_x <- rep(sig_x, rep(numPerAlpha, numAlpha))
#sig_Ey <- rep(sig_RDy, rep(numPerAlpha, numAlpha))
#cof_x <- seq(min_cof_x, max_cof_x, length.out=numAlpha)
cof_x <- rep(0.5, numAlpha)
cof_x <- rep(cof_x, rep(numPerAlpha,  numAlpha))
cof_Ey <- seq(min_cof_x, max_cof_x, length.out=numAlpha)
#cof_Ey <- rep(0, numAlpha)
cof_Ey <- rep(cof_Ey, rep(numPerAlpha,  numAlpha))
nms <- paste(rep(1:100, numAlpha), rep(alpha, rep(numPerAlpha, numAlpha)), sep="_")
seed <- 1234
set.seed(123)
seedDist <- sample(100000, numPerAlpha); seedDist <- rep(seedDist, numAlpha)
seedFun <- sample(100000, numPerAlpha) ; seedFun <- rep(seedFun, numAlpha)
seedNois <- sample(100000, numPerAlpha); seedNois <- rep(seedNois, numAlpha)
addNoise <- FALSE
dataList <- SIMdist_GP_sigPar_wrapper(q=numPerAlpha*numAlpha, n=100, type="SIM", calcStats=TRUE, nms=nms, cof_Ey =cof_Ey, cof_x=cof_x, seed=seed, seedDist=seedDist, seedFun=seedFun, seedNois=seedNois, addNoise=addNoise)

dataList$names
data <- melt(dataList$xs)
colnames(data) <- c("numData","xY","val","dataset")
data <- cast(data, dataset+numData~xY, value="val")
colnames(data) <- c("dataset","numData","x","y")
head(data)
aux <- strsplit( sapply(strsplit(data$dataset, "\\."), function(el) el[2]), "_")
data$alpha <- as.numeric(sapply(aux, function(el) el[2]))
data$rep <- as.numeric(sapply(aux, function(el) el[1]))
table(data$alpha)
table(data$rep)

# scatter
sizePts=0.5; shapePts=4; alph=0.2
p <- ggplot(data[which(data$rep<=10),])
p <- p + geom_point(aes(x, y), colour="steelblue4", alpha=alph,  size=sizePts, shape=shapePts)
p <- p + facet_grid(alpha~rep, scales="free")
#p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
#p <- p + theme(axis.text=element_blank(), axis.ticks=element_blank())
p <- p + theme(axis.title.x=element_blank(),  axis.title.y=element_blank())
p <- p + theme(legend.position= "none")
p <- p + theme(panel.background = element_blank())
#p <- p + scale_color_manual(values=clrs)
print(p)

# histograms
p <- ggplot(data[which(data$rep<=10),])
p <- p + geom_histogram(aes(x=x,y=..density..), position="identity", bins=20, colour="light grey") 
p <- p + geom_density(aes(x=x,y=..density..))
p <- p + facet_grid(alpha~rep)
p

plot(dataList$xs$SIM.1_0.1[,"x"], dataList$xs$SIM.1_100.501[,"x"])

# additivity
head(dataList$descr_indices)
dim(dataList$descr_indices)
df <- dataList$descr_indices
rownames(df) <- dataList$names
df <- melt(df)
colnames(df) <- c("dataset", "msr", "value")
aux <- strsplit( sapply(strsplit(as.character(df$dataset), "\\."), function(el) el[2]), "_")
df$alpha <- as.numeric(sapply(aux, function(el) el[2]))
df$rep <- as.numeric(sapply(aux, function(el) el[1]))
table(df$alpha)
table(df$rep)


p <- ggplot(df)
p <- p + geom_boxplot(aes(x=factor(alpha), y=value, colour=msr))
p <- p + facet_grid(msr~., scales="free")
p 

p <- ggplot(dfSegGen) 
p <- p + geom_density(aes(x=add_indep, fill=types),alpha = 0.1)
p <- p + facet_grid(types~., scales="free")
p

###################33
# kernel based simulation

RDdist <- function(W, tau, n){
  # n <- 100
  # W  <- rgamma(1, 5, 1/0.1)
  # tau <- 1e-4
  x <- rnorm(n)
  x <- sort(x)
  K <- kern_rbf(matrix(x), sigma=1/W)
  alpha <- rnorm(n, mean=0, sd=tau)
  F <- K %*% alpha
  #plot(sort(x), exp(F), type="l")
  G <- cumtrapz(sort(x), exp(F))
  G <- matrix(G[sample(nrow(G))])
  #hist(G)
  return(G)
}

RDdist_GP <- function(W, tau, n){
  # n <- 100
  # W  <- rgamma(1, shape=5, scale=0.1)
  # tau <- 1e-4
  x <- rnorm(n)
  I <- iden(n)
  K <- kern_rbf(matrix(x), sigma=1/(2*W^2))
  C <- K + tau^2*I
  if(!is.positive.definite(C)) C <- make.positive.definite((C))
  C  <- chol(C) 
  F <- t(C) %*% matrix(rnorm(n),n,1)
  #F <-  t(rmvnorm(1, mean=rep(0,n), sigma=solve(K+tau^2*I)))
  
  eF <- exp(F)
  o <- order(x)
  G <- cumtrapz(x[o], eF[o])
  G <- matrix(G[sample(nrow(G))])
  #hist(G)
  return(G)
}

SIMdist_GP <- function(n){
  # n <- 1000
  W_Ex <- rgamma(1, shape=5, scale=0.1); print(W_Ex)
  tau <- 1e-4
  Ex <- RDdist_GP(W_Ex,tau,n)
  hist(Ex)
  Ex <- apply(Ex, 2, stdrize)
  sigma <- rgamma(1, shape=2, scale=1.5)
  K <- kern_rbf(x=Ex, sigma=1/(2*sigma^2))
  I <- diag(n)
  C <- K + tau^2*I
  if(!is.positive.definite(C)) C <- make.positive.definite((C))
  C  <- chol(C) 
  x <- t(C) %*% matrix(rnorm(n),n,1)
  #x <-  t(rmvnorm(1, mean=rep(0,n), sigma=solve(K+tau^2*I))) 
  x <- apply(x, 2, stdrize)
  hist(x)
  W_Ey <- rgamma(1, shape=5, scale=0.1)
  Ey <- RDdist_GP(W_Ey, tau, n)
  # hist(Ey)
  Ey <- apply(Ey, 2, stdrize)
  sigma_x <- rgamma(1, shape=2, scale=3)
  sigma_Ey <- rgamma(1, shape=2, scale=15)
  Kx <- kern_rbf(x, sigma=1/(2*sigma_x^2))
  Ke <- kern_rbf(Ey, sigma=1/(2*sigma_Ey^2))
  K <- Kx*Ke
  
  if(FALSE){
    x2 <- cbind(x,Ey)
    sigma2 <- diag(c(1/sigma_x, 1/sigma_Ey))
    x2 <- x2 %*% sigma2
    x2 <- x2/sqrt(2)
    K2 <- kern_rbf(x2, sigma=1)
    smpl <- sample(length(K), 100)
    plot(c(K)[smpl], c(K2)[smpl]); abline(a=0, b=1, col="red")
  }
  C <- K + tau^2*I
  if(!is.positive.definite(C)) C <- make.positive.definite((C))
  C  <- chol(C) 
  y <- t(C) %*% matrix(rnorm(n),n,1)
  #y <-  t(rmvnorm(1, mean=rep(0,n), sigma=solve(K+tau^2*I))) 
  y <- apply(y, 2, stdrize)
  # add noise
  x <- x + rnorm(n, sd=rgamma(1, shape=2, scale=0.1))
  y <- y + rnorm(n, sd=rgamma(1, shape=2, scale=0.1))
  
  plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y)
  return(res)
}


SIMdist_kern <- function(n){
  # n <- 1000
  W_Ex <- rgamma(1, 5, rate=0.1); print(W_Ex)
  W_w <- rgamma(1, 5, rate=0.1)
  tau <- 1e-4
  Ex <- RDdist(W_Ex,tau,n)
  hist(Ex)
  Ex <- apply(Ex, 2, stdrize)
  Ex2 <- RDdist(W_Ex,tau,n)
  Ex2 <- apply(Ex2, 2, stdrize)
  w <- RDdist(W_w, tau,n)
  sigma <- rgamma(1, 2, rate=1.5)
  K <- kern_rbf(x=Ex, sigma=1/sigma)
  K2 <- kern_rbf(x=Ex2, y=w, sigma=1/sigma)
  alpha <- rnorm(n, sd=tau)
  x <- K %*% alpha
  
  x <- apply(x, 2, stdrize)
  x2 <- K2%*% alpha
  x2 <- apply(x2, 2, stdrize)
  hist(x)
  W_Ey <- rgamma(1, 5, 1/0.1)
  Ey <- RDdist(W_Ey, tau, n)
  # hist(Ey)
  Ey <- apply(Ey, 2, stdrize)
  Ey2 <- RDdist(W_Ey, tau, n)
  Ey2 <- apply(Ey, 2, stdrize)
  # hist(Ey)
  sigma_x <- rgamma(1, 2, 1/15)
  sigma_Ey <- rgamma(1, 2, 1/15)
  Kx <- kern_rbf(x, sigma=1/sigma_x)
  Ke <- kern_rbf(Ey, sigma=1/sigma_Ey)
  K <- Kx*Ke
  alpha <- rnorm(n, sd=tau)
  y <- K%*%alpha
  y <- apply(y, 2, stdrize)
  
  # add noise
  
  plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y)
  return(res)
}


