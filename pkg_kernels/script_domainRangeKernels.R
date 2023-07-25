# explore the domain and range of different kernels 
library(GGally)
remove(list=ls())
server <- "optimus.uv.es"
user <- "emiliano"
repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_causaLearner/genData/func_getData_v2.R")

n <- 100
p <- 10
x <- matrix(runif(n*p), n, p)



getRangeRbf <- function(x, length.out=10){

  sigma0 <- 1/quantile(as.numeric(dist(x)^2), (length(x)-1)/length(x))
  sigma1 <- 1/quantile(as.numeric(dist(x)^2), 0.6)  
   
  sigmas <- 10^seq(-10,10,by=1)
  minPct <- 0.01
  Ks <- sapply(sigmas, function(sd) kern_rbf(x, sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))

  plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  lines(spl, col="green")
  lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  abline(v=log(c(sigma0, sigma1),10), col="green")
  sigmas.x[c(indxStart, indxMax.x, indxFinish)]
  log(sigmas.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  vars.x[c(indxStart, indxMax.x, indxFinish)]

  # round 2
  sigmas <- 10^seq(log(sigmas.x[indxStart],10), log(sigmas.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(sigmas, function(sd) kern_rbf(x, sigma=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  spl <- spline(log(sigmas[indx],10), vars[indx])
  splf <- splinefun(log(sigmas[indx],10), vars[indx])
  sigmas.x <- unique(sort(c(sigmas[indxMax],10^seq(min(log(sigmas[indx], 10)), max(log(sigmas[indx], 10)), length.out=100))))
  vars.x <- splf(log(sigmas.x, 10), deriv=0)
  indxMax.x <- which(sigmas.x == sigmas[indxMax]); sigmas.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))

  plot(log(sigmas,10), vars, ylab="var"); abline(h=0, col="red")
  lines(spl, col="green")
  lines(log(sigmas.x, 10), vars.x, col="blue", type="l")
  abline(h=maxVar*minPct, v=log(sigmas.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  abline(v=log(c(sigma0, sigma1),10), col="green")
  
  ini <- log(sigmas.x, 10)[indxStart]
  fin <- log(sigmas.x, 10)[indxFinish]
  seq <- 10^seq(ini, fin, length.out=length.out)
  return(seq)
}

getRangeLog <- function(x, length.out=10){
  
  
  degrees1 <- c(10^seq(-9,-1, length.out=round(length.out/2)-2), 0.5, 0.98)
  
  
  degrees2 <- round(seq(2,50, length.out=round(length.out/2)))
  
  degrees <- c(degrees1,degrees2)
  
  Ks <- sapply(degrees, function(dg) kern_log(x, degree=dg), simplify="array")
  num0s <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    sum(mat2==0, na.rm=T)/2
  })
  
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)
    return(res)
  })
  plot(degrees, log(vars, 10))
  plot(degrees, vars)
  
  deg1st0 <- degrees[which(num0s>0)[1]]
  degG0 <- degrees[which(vars>1e-10)[1]]
  plot(degrees, num0s, ylab="num0s") 
  abline(h=0, v=c(deg1st0, degG0), col="red")
  plot(log(degrees, 10), log(vars,10), ylab="num0s") 
  abline(h=0, v=log(c(deg1st0, degG0),10), col="red")
  
  
  degrees <- degrees[which(vars>1e-10 & num0s==0)]
  
  return(degrees)
}

getRangeQuad <- function(x, length.out=10){
 
  
  offsets <- seq(-10,10 ,length.out=21)
  
  Ks <- sapply(offsets, function(of) kern_quad(x, offset=of), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/abs(mean(as.numeric(mat2), na.rm=T))
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  
  plot(offsets, vars)
  
  ini <- which(vars>1)[1]
  fin <- which(vars>1)[sum(vars>1, na.rm=T)]
  
  offsets <- seq(offsets[ini], offsets[fin], length.out=length.out)
  return(offsets)
}

getRangeLaplace <- function(x, length.out=10){
  
  scale0 <- quantile(as.numeric(dist(x)^2), (length(x)-1)/length(x))
  scale1 <- quantile(as.numeric(dist(x)^2), 0.6)  
  
  scales <- 10^seq(-10,10,by=1)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    hist(mat2)
    diag(mat2) <- NA
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    return(res)
  })
  
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  lines(spl, col="green")
  lines(log(scales.x, 10), vars.x, col="blue", type="l")
  abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  abline(v=log(c(scale0, scale1),10), col="green")
  scales.x[c(indxStart, indxMax.x, indxFinish)]
  log(scales.x, 10)[c(indxStart, indxMax.x, indxFinish)]
  vars.x[c(indxStart, indxMax.x, indxFinish)]
  
  # round 2
  scales <- 10^seq(log(scales.x[indxStart],10), log(scales.x[indxFinish],10), length.out=10)
  minPct <- 0.01
  Ks <- sapply(scales, function(sd) kern_laplace(x, scale=sd), simplify="array")
  vars <- apply(Ks, 3, function(mat){
    # mat <- Ks[,,12]
    mat2 <- mat
    diag(mat2) <- NA
    hist(mat2)
    res <- var(as.numeric(mat2), na.rm=T)/mean(as.numeric(mat2), na.rm=T)
    #res <- var(as.numeric(mat2), na.rm=T)/(max(mat2, na.rm=T)-min(mat2, na.rm=T))
    return(res)
  })
  maxVar <- max(vars, na.rm=T)
  indxMax <- which.max(vars)
  indx <- which(!is.na(vars))
  spl <- spline(log(scales[indx],10), vars[indx])
  splf <- splinefun(log(scales[indx],10), vars[indx])
  scales.x <- unique(sort(c(scales[indxMax],10^seq(min(log(scales[indx], 10)), max(log(scales[indx], 10)), length.out=100))))
  vars.x <- splf(log(scales.x, 10), deriv=0)
  indxMax.x <- which(scales.x == scales[indxMax]); scales.x[indxMax.x]
  indxStart <- which.min(abs(vars.x[1:indxMax.x]-maxVar*minPct))
  indxFinish <- indxMax.x + which.min(abs(vars.x[(indxMax.x+1):length(vars.x)]-maxVar*minPct))
  
  plot(log(scales,10), vars, ylab="var"); abline(h=0, col="red")
  lines(spl, col="green")
  lines(log(scales.x, 10), vars.x, col="blue", type="l")
  abline(h=maxVar*minPct, v=log(scales.x,10)[c(indxStart, indxMax.x, indxFinish)], col="orange")
  abline(v=log(c(scale0, scale1),10), col="green")
  
  ini <- log(scales.x, 10)[indxStart]
  fin <- log(scales.x, 10)[indxFinish]
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

x <- matrix(runif(n*p), n, p)

p <- 4
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)
x <- x$x
data <- as.data.frame(x)
ggpairs(data, aes(alpha = 0.4))

#######################################################################################*
# RBF
#######################################################################################*


x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)$x
getRangeRbf(x)

#######################################################################################*
# LOG
#######################################################################################*

x <- matrix(runif(n*p), n, p)
x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)$x
getRangeLog(x)

#######################################################################################*
# quad
#######################################################################################*


x <- matrix(runif(n*p), n, p)
x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)$x
getRangeQuad(x)


#######################################################################################*
# Laplace
#######################################################################################*


x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)$x
getRangeLaplace(x)


#######################################################################################*
# polynomial
#######################################################################################*

x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)$x

getRangePoly(x)

