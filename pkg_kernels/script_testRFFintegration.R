# test integration of fourier features

# see if CME evaluated at real and fake data points is always positive for rbf kernel
remove(list=ls())
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
dir("./pkg_kernels/")

source("./pkg_kernels/func_kernel_pkg.R")

n <- 10
py <- 2
y_tr <- matrix(rnorm((n*py)),n,py)
px <- 3
x_tr <- matrix(rnorm((n*px)),n,px)
num_f <- 1000
sigmax <- 2.5
sigmay <- 0.7
Phiy <- rff(x=y_tr, num_f, seed=1234, p_w="rnorm2", map="cos", sigma=sigmay)

Lx <- kern_rbf(x_tr, sigma=sigmax)
lambda <- 1e-5
I <- diag(n)

Blambda <- Lx + n*lambda*I

n_xx <- 1
xx <- rnorm((n_xx*px))
y <- matrix(rnorm((1*py)),1,py)

normProb <- function(xx, px, y, Blambda, y_tr, x_tr, sigmay, sigmax){
  # y <- 0;  x <- 0.1
  #print("eval")
  xx <- matrix(xx,1,px)
  n <- nrow(x_tr)
  ky <- kern_rbf(x=y, y=y_tr, sigma=sigmay)
  intlx <- kern_rbf_nrml(x=x_tr, y=xx, sigma=sigmax, lim_max=0.5, lim_min=-0.5)
  lx <- kern_rbf(x=x_tr, y=xx, sigma=sigmax)
  res <- (ky%*%Blambda%*%lx)/(ky%*%Blambda%*%intlx)
  return(as.numeric(res))
}


library(cubature)
ints <- sapply(c("hcubature","pcubature", "cuhre", "divonne", "suave", "vegas"), 
       function(meth) cubintegrate(normProb, lower=rep(-0.5,px), upper=rep(0.5,px),
                                   method = meth, maxEval=10^3,
                                   y=y, px=px, Blambda=Blambda, y_tr=y_tr, x_tr=x_tr, sigmax=sigmax, sigmay=sigmay))

sapply(ints, function(el) el$integral)


RFF <- function(yy, py, w, b){
  yy <- matrix(yy, 1, py)
  res <- cos(yy %*% w + b)
  return(as.numeric(res))
}

for(py in 2:14){
#py <- 2
print("*****************")
print(py)


n_yy <- 1
yy <- rnorm((n_xx*py))
w <- rnorm(py)
b <- runif(1)


ints <- sapply(c("hcubature","pcubature", "cuhre", "divonne", "suave", "vegas"), 
               function(meth) cubintegrate(RFF, lower=rep(-0.5,py), upper=rep(0.5,py),
                                           method = meth, maxEval=10^3,
                                           py=py, w=w, b=b))



intToBitss <- function(x, maxInt) paste(rev(strsplit(substr(paste(as.integer(intToBits(x-1)), collapse=""),1, ceiling(log(maxInt,2))), split="")[[1]]), collapse="")

intAnalyticRFF <- function(py, w, b){
  ais <- sapply(1:(2^py), intToBitss, maxInt=2^py)
  ais <- sapply(strsplit(ais, split=""), function(el) c(1,-1)[as.numeric(el)+1])
  summandsOdd  <- sin(0.5*apply(ais*w, 2, sum)+b)
  summandsEven <- cos(0.5*apply(ais*w, 2, sum)+b)
  coefs <- (-1)^floor(py/2)*apply(ais,2,prod)
  summands <- switch(c("even","odd")[py%%2+1], even=summandsEven, odd=summandsOdd)
  res <- (1/prod(w))*sum(coefs*summands)
  return(res)
}

print(intAnalyticRFF(py, w, b))
print(sapply(ints, function(el) el$integral))
}


n <- 10
py <- 2
y_tr <- matrix(rnorm((n*py)),n,py)
px <- 3
x_tr <- matrix(rnorm((n*px)),n,px)
num_f <- 10
sigmax <- 2.5
sigmay <- 0.7
Phiy <- rff(x=y_tr, num_f, seed=1234, p_w="rnorm2", map="cos", sigma=sigmay)

Lx <- kern_rbf(x_tr, sigma=sigmax)
lambda <- 1e-5
I <- diag(n)

Blambda <- Lx + n*lambda*I

n_yy <- 1
yy <- rnorm((n_yy*py))
x <- matrix(rnorm((1*px)),1,px)

count <<- 0
normProb <- function(yy, py, x, Blambda, Phiy,  y_tr, x_tr, sigmay, sigmax){
  # y <- 0;  x <- 0.1
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


library(cubature)
ints <- sapply(c("hcubature","pcubature", "cuhre", "divonne", "suave", "vegas"), 
               function(meth) cubintegrate(normProb, lower=rep(-0.5,py), upper=rep(0.5,py),
                                           method = meth, maxEval=10^3,
                                           x=x, py=py, Blambda=Blambda, Phiy=Phiy, y_tr=y_tr, x_tr=x_tr, sigmax=sigmax, sigmay=sigmay))

sapply(ints, function(el) el$integral)


