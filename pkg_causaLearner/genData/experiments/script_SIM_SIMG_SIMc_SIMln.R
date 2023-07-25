
library(pracma) # cumtrapz
library(parallel) # detectCores()
library(corpcor) #make.positive.definite
library(reshape)
library(ggplot2)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_causaLearner/utilities/func_dagStuff.R")
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")
source("./pkg_causaLearner/genData/func_getData_v2.R")


rGP_rbf <- function(x, tau, sigma){
  #if((ncol(x) !=length(sigma) & length(sigma)!=1) & (ncol(x) !=length(cutoff) & length(cutoff)!=1)) error("one sigma per column or one sigma for all columns!")
  #if(is.null(sigma) & is.null(cutoff)) error("one of sigma and error shd be provided")
  n <- nrow(x)
  
  x2 <- x
  removeDiag <- TRUE
  
  sigma2 <- matrix(0, ncol(x), ncol(x))
  diag(sigma2) <- 1/sigma
  x2 <- x %*% sigma2
  x2 <- x2/sqrt(2)
  
  K <- kern_rbf(x2, sigma=1) 
  n <- nrow(K)
  I <- diag(n)
  C <- K + tau^2*I
  if(!is.positive.definite(C)) C <- make.positive.definite((C))
  C  <- chol(C)
  y <- rnorm(nrow(x))
  y <- t(C) %*% matrix(y,n,1)
  
  # we want cov(y, y') = K(y, y') 
  # say y=t(C) %*% x and y'= t(C') %*% x', x,x'~ N(0,I) 
  # then cov(y, y') = cov(t(C)%*%x, t(C')%*%x')
  #                 = 
  
  
  return(y)  
}

RDdist_GP <- function(n, tau, sigma){
  # n <- 100
  # sigma  <- rgamma(1, shape=5, scale=0.1)
  # tau <- 1e-4
  x <- matrix(rnorm(n))
  smpl <-sample(nrow(x))
  F <- rGP_rbf(x, tau, sigma=sigma)
  eF <- exp(F)
  o <- order(x)
  G <- cumtrapz(x[o], eF[o])
  G <- matrix(G[smpl])
  return(G)
}

SIMdist_GP_sigPar <- function(n,  tau, sig_RDx=NA, sig_RDy=NA, 
                              sig_Ex=NA, sig_Ey=NA, sig_x=NA,  
                              sig_RDz=NA, sig_Ez=NA, 
                              sig_nois_x=NULL, sig_nois_y=NULL){
  # n <- 1000
  
  Ex <- RDdist_GP(n, tau, sigma=sig_RDx)
  Ex <- apply(Ex, 2, stdrize)
  Ey <- RDdist_GP(n, tau, sigma=sig_RDy)
  Ey <- apply(Ey, 2, stdrize)
  
  Ex2 <- Ex
  noisNms <- c("Ex")
  sig_Ex2 <- sig_Ex
  
  if(!is.na(sig_RDz)){
    #print("conf")
    Ez <- RDdist_GP(n, tau, sigma=sig_RDz)
    Ez <- apply(Ez, 2, stdrize)
    Ex2 <- cbind(Ex2, Ez)
    noisNms <- c(noisNms,"Ez")
    sig_Ex2 <- c(sig_Ex2, sig_Ez)
  }
  nois <- cbind(Ex2, Ey)
  noisNms <- c(noisNms, "Ey")
  
  x <- rGP_rbf(Ex2, tau, sigma=sig_Ex2)
  x <- apply(x, 2, stdrize)
  x2 <- cbind(x, Ey)
  
  sig_x2 <- c(sig_x, sig_Ey)
  
  if(!is.na(sig_RDz)){
    x2 <- cbind(x2, Ez)
    sig_x2 <- c(sig_x2, sig_Ez)
  }
  
  y <- rGP_rbf(x=x2, tau,  sigma=sig_x2)
  set.seed(NULL)
  
  y <- apply(y, 2, stdrize)
  
  # add noise
  
  addNoiseX <- rnorm(n, sd=sig_nois_x)
  x <- x + addNoiseX
  addNoiseY <- rnorm(n, sd=sig_nois_y) 
  y <- y + addNoiseY
  nois <- cbind(nois, addNoiseX, addNoiseY)
  noisNms <- c(noisNms, "addNoiseX","addNoiseY")
  
  #plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y, nois)
  colnames(res) <- c("x","y", noisNms)
  res <- list(res=res)
  
  return(res)
}

SIMdist_GP_sigPar_wrapper <- function(q=100, n=1000, 
                                      type=c("SIM", "SIMc","SIMG","SIMln"), 
                                      nms=NULL){
  #q=100; n=1000; seed=NULL; tau=NULL; sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL; sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  type <- match.arg(type, choices=c("SIM", "SIMc","SIMG","SIMln"))
  
  pars <- list()
  pars$n <- n
  #n <- 100
  #type <- "SIM"
  
  pars$tau <- 1e-4
  
  pars <- switch(type, 
                 SIM={
                   pars$sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDz <- NULL
                   pars$sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   pars$sig_x <- rgamma(q, shape=2, scale=3)
                   pars$sig_Ey <- rgamma(q, shape=2, scale=15)
                   pars$sig_Ez <- NULL
                   pars[["sig_nois_x"]] <- rgamma(q, shape=2, scale=0.1)
                   pars[["sig_nois_y"]] <- rgamma(q, shape=2, scale=0.1)
                   pars
                   
                 },
                 SIMc={
                   pars$sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDz <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   pars$sig_x <- rgamma(q, shape=2, scale=3)
                   pars$sig_Ey <- rgamma(q, shape=2, scale=15)
                   pars$sig_Ez <- rgamma(q, shape=2, scale=15)
                   pars[["sig_nois_x"]] <- rgamma(q, shape=2, scale=0.1)
                   pars[["sig_nois_y"]] <- rgamma(q, shape=2, scale=0.1)
                   pars
                 },
                 SIMG={
                   pars$sig_RDx <- rgamma(q, shape=1e6, scale=1e-3)
                   pars$sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDz <- NULL
                   pars$sig_Ex <- rgamma(q, shape=1e6, scale=1e-3)
                   pars$sig_x <- rgamma(q, shape=2, scale=3)
                   pars$sig_Ey <- rgamma(q, shape=2, scale=15)
                   pars$sig_Ez <- NULL
                   pars[["sig_nois_x"]] <- rgamma(q, shape=2, scale=0.1)
                   pars[["sig_nois_y"]] <- rgamma(q, shape=2, scale=0.1)
                   pars
                 },
                 SIMln={
                   pars$sig_RDx <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDy <- rgamma(q, shape=5, scale=0.1)
                   pars$sig_RDz <- NULL
                   pars$sig_Ex <- rgamma(q, shape=2, scale=1.5)
                   pars$sig_x <- rgamma(q, shape=2, scale=3)
                   pars$sig_Ey <- rgamma(q, shape=2, scale=1.5*200)
                   pars$sig_Ez <- NULL
                   pars[["sig_nois_x"]] <- rgamma(q, shape=2, scale=0.01)
                   pars[["sig_nois_y"]] <- rgamma(q, shape=2, scale=0.01)
                   pars
                 })
  
  print("names(pars)")
  print(names(pars))
  
  pars <- lapply(1:q, function(i) lapply(pars, function(el) if(length(el)==1) el else el[i]))
  
  if(!is.null(nms)) nms <- paste(type, 1:q, nms, sep=".") else nms <- paste(type, 1:q, sep=".")
  
  mc_cores <- detectCores()/2
  #mc_cores <- 1
  pm <- proc.time()
  xs <- mcmapply(function(par, nm){
    # i <- 64; par <- pars[[i]]; nm <- nms[i]
    print(paste("nm: ", nm))
    res <- do.call("SIMdist_GP_sigPar", par)
    return(res)
  }, par=pars, nm=nms, mc.cores=mc_cores, SIMPLIFY=FALSE)
  proc.time() - pm 
  # 1.241 secs for q=100, n=100 and 4 cores with confounder (3 gps)
  # 7.15 mins for q=100, n=1000 and 4 cores with confounder (3 gs)
  xs2 <- lapply(xs, function(el) el$res[,c("x","y")])
  
  xs3 <- lapply(xs, function(el){
    z1 <- el$res[,"Ey"]
    z2 <- el$res[,"addNoiseY"]
    res <- cbind(z1, z2)
    colnames(res) <- c("z1","z2")
    if("Ez" %in% colnames(el$res)) res <- cbind(res, z3=el$res[,"Ez"])
    res <- cbind(el$res[,c("x","y")], res)
    return(res)
  })
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs, function(el) el$res[,-which(colnames(el$res) %in% c("x","y"))])
  names(xs3) <- names(xs2) <- names(ns) <- names(dags) <- nms
  #dataList <- list(xs=xs2, ns=ns, zs=xs3, dags=dags, names=nms)
  dataList <- list(xs=xs3, noiss=ns, dags=dags, names=nms)
  
  #plotPairsList(dataList)
  return(dataList)  
}


# SIM

dataList_SIM <- SIMdist_GP_sigPar_wrapper(type="SIM")
#plotPairsList(dataList_SIM)
head(dataList_SIM$xs[[1]])

# SIMc

dataList_SIMc <- SIMdist_GP_sigPar_wrapper(type="SIMc")
#plotPairsList(dataList_SIMc)
head(dataList_SIMc$xs[[1]])

# SIMln

dataList_SIMln <- SIMdist_GP_sigPar_wrapper(type="SIMln")
#plotPairsList(dataList_SIMln)


# SIMG

dataList_SIMG <- SIMdist_GP_sigPar_wrapper(type="SIMG")
#plotPairsList(dataList_SIMG)
head(dataList_SIMG$xs[[1]])

# join

dataList <- dataJoin(dataListList=list(dataList_SIM, dataList_SIMc,  dataList_SIMG, dataList_SIMln))
nms <- substr(dataList$names,2, nchar(dataList$names))
dataList$names <- names(dataList$dags) <- names(dataList$xs) <- names(dataList$noiss) <- nms

# write 

library(jsonlite)
jsonData2 <- jsonlite:::toJSON(x=dataList)

folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/data/TCEPs/"
save(dataList, file=paste(folderSave, "dag2-ME2-SIM-1000_withZ_sims.RData",sep=""))
write(jsonData2, file=paste(folderSave, "dag2-ME2-SIM-1000_withZ_sims.json",sep=""))

i <- 52
plot(dataList$xs[[i]][,1], dataList$xs[[i]][,2])
dataList$names[i]
