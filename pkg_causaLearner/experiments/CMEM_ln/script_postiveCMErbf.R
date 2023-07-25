# see if CME evaluated at real and fake data points is always positive for rbf kernel
remove(list=ls())
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
dir("./pkg_kernels/")

source("./pkg_kernels/func_kernel_pkg.R")

n <- 1000
x_true <- rnorm(n)
eps <- rnorm(n)
y_true <- x_true^2*eps
x_true <- (x_true-min(x_true))/(max(x_true)-min(x_true))
y_true <- (y_true-min(y_true))/(max(y_true)-min(y_true))
summary(x_true)
summary(y_true)
plot(x_true, y_true)
y_fake <- runif(n)

sigmax <- 50
sigmay <- 50

Ky_true <- kern_rbf(matrix(y_true, n, 1), sigma=sigmay)
Ky_fake <- kern_rbf(matrix(y_fake, n, 1), sigma=sigmay)
Lx <- kern_rbf(matrix(x_true, n, 1), sigma=sigmax)
lambda <- 1e-5
I <- diag(n)

Blambda <- Lx + n*lambda*I
Blambda <- solve(Blambda)

x_tr <- x_true
y_tr <- y_true
unNormProb <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
  # y <- 0;  x <- 0.1
  n <- length(x_tr)
  ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
  lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
  res <- ky%*%Blambda%*%lx
  return(as.numeric(res))
}
normProb <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
  # y <- 0;  x <- 0.1
  n <- length(x_tr)
  ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
  lx <- kern_rbf_nrml(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax, lim_max=1, lim_min=0)
  res <- ky%*%Blambda%*%lx
  return(as.numeric(res))
}

normProb2 <- function(xx, y, Blambda, y_tr, x_tr, sigmay, sigmax){
  # y <- 0;  x <- 0.1
  n <- length(x_tr)
  ky <- kern_rbf(x=matrix(y, 1, 1), y=matrix(y_tr, n, 1), sigma=sigmay)
  intlx <- kern_rbf_nrml(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax, lim_max=1, lim_min=0)
  lx <- kern_rbf(x=matrix(x_tr, n, 1), y=matrix(xx,length(xx),1), sigma=sigmax)
  res <- (ky%*%Blambda%*%lx)/(ky%*%Blambda%*%intlx)
  return(as.numeric(res))
}


xx <- seq(0,1,length.out=100)
y <- 0.5
yy_unorm <- unNormProb(xx, y=y, Blambda=Blambda, y_tr=y_true, x_tr=x_true, sigmay=sigmay, sigmax=sigmax)
yy_norm <- normProb2(xx, y=y, Blambda=Blambda, y_tr=y_true, x_tr=x_true, sigmay=sigmay, sigmax=sigmax)

summary(yy_unorm)
summary(yy_norm)

par(mfrow=c(1,2))
plot(x_true, y_true); abline(h=y, col="red") 
plot(xx, yy_unorm, type="l", ylim=range(yy_unorm, yy_norm))
lines(xx, yy_norm, type="l", col="red"); abline(h=0, col="blue")

# can I integrate normProb(x) to see if it integrates to one??
library(pracma)
sapply(c("Kronrod","Clenshaw","Simpson")[c(1,2,3)], function(meth) integral(normProb2, xmin=0, xmax=1, method = meth, 
         y=y, Blambda=Blambda, y_tr=y_true, x_tr=x_true, sigmay=sigmay, sigmax=sigmax))



summary(as.numeric(Lx %*% Blambda %*% Ky_true))
summary(as.numeric(Lx %*% Blambda %*% Ky_fake))
sum(as.numeric(Lx %*% Blambda %*% Ky_true)<0)/(n^2)
sum(as.numeric(Lx %*% Blambda %*% Ky_fake)<0)/(n^2)
indx <- which(as.numeric(Lx %*% Blambda %*% Ky_true)<0)
indxC <- which(!as.numeric(Lx %*% Blambda %*% Ky_true)<0)
sum(-as.numeric(Lx %*% Blambda %*% Ky_true)[indx]) / sum(as.numeric(Lx %*% Blambda %*% Ky_true)[indxC])
indx <- which(as.numeric(Lx %*% Blambda %*% Ky_fake)<0)
indxC <- which(!as.numeric(Lx %*% Blambda %*% Ky_fake)<0)
sum(-as.numeric(Lx %*% Blambda %*% Ky_fake)[indx]) / sum(as.numeric(Lx %*% Blambda %*% Ky_fake)[indxC])

summary(diag(Lx %*% Blambda %*% Ky_true))
summary(diag(Lx %*% Blambda %*% Ky_fake))
sum(diag(Lx %*% Blambda %*% Ky_true)<0)/(n^2)
sum(diag(Lx %*% Blambda %*% Ky_fake)<0)/(n^2)

