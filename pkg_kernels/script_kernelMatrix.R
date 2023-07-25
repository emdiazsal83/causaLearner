# test kernelMatrix vs me calculating it in R or in C

remove(list=ls())

library(kernlab)
library(Matrix) # nearPD
library(matrixcalc) #is.positive.definite
library(microbenchmark)
library(reshape)
library(ggplot2)
library(Rcpp)

server <- "optimus.uv.es"
user <- "emiliano"
repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)

matNorm <- function(x, y=x){
  n <- nrow(x)
  m <- nrow(y)
  d <- ncol(x)
  onesN <- rep(1, n)
  onesM <- rep(1, m)
  xtx <- apply(x*x, 1, sum)
  yty <- apply(y*y, 1, sum)
  matNorms <- xtx%*%t(onesM) + onesN%*%t(yty) - 2*x%*%t(y)  
  return(matNorms)
}

# For different kernels go to:
# http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/

#################################################################################################*
# Vanilla / Linear 
#################################################################################################*

# depends on x %*% t(y)
kern_lin <- function(x, y=x, offset){
  K <- x%*%t(y)+offset
  return(K)
}

cppFunction('NumericMatrix kern_lin_Caux(NumericMatrix x, NumericMatrix y, double offset) {
            int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += x(i, k)*y(j,k);
            }
            out(i,j) = total+offset;
            }
            }
            return out;
            }')

kern_lin_C <- function(x, y=x, offset) kern_lin_Caux(x, y, offset)


#################################################################################################*
# RBF Kernel
#################################################################################################*

# depends on ||x - y||_2^2
kern_rbf <- function(x, y=x, sigma){
  matNorms <- matNorm(x, y)
  K <- exp(-sigma*matNorms)
  return(K)
}

cppFunction('NumericMatrix kern_rbf_Caux(NumericMatrix x, NumericMatrix y, double sigma) {
  int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += pow(x(i, k)-y(j,k),2.0);
            }
            out(i,j) = exp(-sigma*total);
            }
            }
            return out;
            }')

kern_rbf_C <- function(x, y=x, sigma) kern_rbf_Caux(x, y, sigma)

cppFunction('NumericMatrix kern_rbf_Caux2(NumericMatrix x, NumericMatrix y, double sigma) {
            int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total =  sum(pow(x(i,_)-y(j,_),2.0));
            
            out(i,j) = exp(-sigma*total);
            }
            }
            return out;
}')

kern_rbf_C2 <- function(x, y=x, sigma) kern_rbf_Caux2(x, y, sigma)

#################################################################################################*
# Laplace Kernel
#################################################################################################*

# depends on ||x - y||_2^2
kern_laplace <- function(x, y=x, sigma){
  matNorms <- matNorm(x, y)
  matNorms[matNorms<0] <- 0
  # isSymmetric(matNorms)
  # is.positive.definite(matNorms)
  # is.positive.definite(as.matrix(nearPD(matNorms, keepDiag=TRUE, maxit=10000)$mat))
  # plot(matNorms, as.matrix(nearPD(matNorms)$mat))
  K <- exp(-sigma*matNorms^0.5)
  return(K)
}

cppFunction('NumericMatrix kern_laplace_Caux(NumericMatrix x, NumericMatrix y, double sigma) {
  int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += pow(x(i, k)-y(j,k),2.0);
            }
            out(i,j) = exp(-sigma*pow(total,0.5));
            }
            }
            return out;
            }')

kern_laplace_C <- function(x, y=x, sigma) kern_laplace_Caux(x, y, sigma)


#################################################################################################*
# Polynomial Kernel
#################################################################################################*

# depends on x %*% t(y)
kern_poly <- function(x, y=x, degree, scale, offset){
  K <- (scale*(x%*%t(y))+offset)^degree
  return(K)
}

cppFunction('NumericMatrix kern_poly_Caux(NumericMatrix x, NumericMatrix y, double degree, double scale, double offset) {
            int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += x(i, k)*y(j,k);
            }
            out(i,j) = pow((scale*total)+offset, degree);
            }
            }
            return out;
            }')

kern_poly_C <- function(x, y=x, degree, scale, offset) kern_poly_Caux(x, y, degree, scale, offset)

#################################################################################################*
# Hyperbolic tan Kernel
#################################################################################################*

# depends on x %*% t(y)
kern_tanh <- function(x, y=x, scale, offset){
  K <- tanh(scale*(x%*%t(y))+offset)
  return(K)
}

cppFunction('NumericMatrix kern_tanh_Caux(NumericMatrix x, NumericMatrix y, double scale, double offset) {
  int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += x(i, k)*y(j,k);
            }
            out(i,j) = tanh(scale*total+offset);
            }
            }
            return out;
            }')

kern_tanh_C <- function(x, y=x, scale, offset) kern_tanh_Caux(x, y, scale, offset)

#################################################################################################*
# Log Kernel
#################################################################################################*

# depends on ||x - y||_2^2
kern_log <- function(x, y=x, degree){
  matNorms <- matNorm(x, y)
  K <- -log(matNorms^degree+1)
  return(K)
}

cppFunction('NumericMatrix kern_log_Caux(NumericMatrix x, NumericMatrix y, double degree) {
  int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += pow(x(i, k)-y(j,k),2.0);
            }
            out(i,j) = -log(pow(total, degree) + 1);
            }
            }
            return out;
            }')

kern_log_C <- function(x, y=x, degree) kern_log_Caux(x, y, degree)

#################################################################################################*
# Rational Quadratic Kernel
#################################################################################################*

# depends on ||x - y||_2^2
kern_quad <- function(x, y=x, offset){
  matNorms <- matNorm(x, y)
  K <- 1-matNorms/(matNorms+offset)
  return(K)
}

cppFunction('NumericMatrix kern_quad_Caux(NumericMatrix x, NumericMatrix y, double offset) {
  int n = x.nrow();
            int m = y.nrow();
            int d = x.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total = 0;
            for(int k = 0;  k < d; k++){
            total += pow(x(i, k)-y(j,k),2.0);
            }
            out(i,j) = 1-total/(total+offset);
            }
            }
            return out;
            }')

kern_quad_C <- function(x, y=x, offset) kern_quad_Caux(x, y, offset)

  



# tests

n <- 100
m <- 20
d <- 10
x <- matrix(rnorm(n*d), n, d)
y <- matrix(rnorm(m*d), m, d)

K1 <- kern_rbf(x, y, sigma=1)
K2 <- kern_rbf_C(x, y, sigma=1)
K3 <- kern_rbf_C2(x, y, sigma=1)

dim(K1); dim(K2); dim(K3)
plot(K1, K2); abline(a=0, b=1, col="red")
plot(K2, K3); abline(a=0, b=1, col="red")

# lin
offset0 <- 0
kernelX <- vanilladot()
Kx1 <- kernelMatrix(kernelX, x)
Kx2 <- kern_lin(x, offset=offset0)
Kx3 <- kern_lin_C(x, offset=offset0)
max((Kx1-Kx2)^2); max((Kx2-Kx3)^2) 
all(Kx1==Kx2); all(Kx2==Kx3) 

# rbf
sigma0 <- 1
kernelX <- rbfdot(sigma=sigma0)
Kx1 <- kernelMatrix(kernelX, x)
Kx2 <- kern_rbf(x, sigma=sigma0)
Kx3 <- kern_rbf_C(x, sigma=sigma0)
max((Kx1-Kx2)^2); max((Kx2-Kx3)^2) 
all(Kx1==Kx2); all(Kx2==Kx3) 

# Laplacian
sigma0 <- 1
kernelX <- laplacedot(sigma=sigma0)
Kx1 <- kernelMatrix(kernelX, x)
Kx2 <- kern_laplace(x, sigma=sigma0)
Kx3 <- kern_laplace_C(x, sigma=sigma0)
max((Kx1-Kx2)^2, na.rm=T); max((Kx2-Kx3)^2, na.rm=T) 
all(Kx1==Kx2); all(Kx3==Kx2) 
plot(Kx1, Kx2); abline(a=0, b=1, col="red")
plot(Kx3, Kx2); abline(a=0, b=1, col="red")

# polynomial
degree0 <- scale0 <- offset0 <- 1
kernelX <- polydot(degree=degree0, scale=scale0, offset=offset0)
Kx1 <- kernelMatrix(kernelX, x)
Kx2 <- kern_poly(x, degree=degree0, scale=scale0, offset=offset0)
Kx3 <- kern_poly_C(x, degree=degree0, scale=scale0, offset=offset0)
max((Kx1-Kx2)^2); max((Kx3-Kx2)^2) 
all(Kx1==Kx2); all(Kx3==Kx2)

# Hyperbolic tangent
scale0 <- offset0 <- 1
kernelX <- tanhdot(scale=scale0, offset=offset0)
Kx1 <- kernelMatrix(kernelX, x)
Kx2 <- kern_tanh(x, scale=scale0, offset=offset0)
Kx3 <- kern_tanh_C(x, scale=scale0, offset=offset0)
max((Kx1-Kx2)^2); max((Kx3-Kx2)^2) 
all(Kx1==Kx2); all(Kx3==Kx2)

# log

degree0  <- 2
Kx2 <- kern_log(x, degree=degree0)
Kx3 <- kern_log_C(x, degree=degree0)
max((Kx3-Kx2)^2) 
all(Kx3==Kx2)

# quadratic

offset0 <- 1
Kx2 <- kern_quad(x,  offset=offset0)
Kx3 <- kern_quad_C(x, offset=offset0)
max((Kx3-Kx2)^2) 
all(Kx3==Kx2)


n <- 2000
d <- 10
x <- matrix(rnorm(n*d), n, d)
# lin - R version best
microbenchmark(kern_lin(x, offset=offset0), kern_lin_C(x, offset=offset0), kernelMatrix(vanilladot(), x))
# rbf - C version best
microbenchmark(kern_rbf(x, sigma=sigma0), kern_rbf_C(x, sigma=sigma0), kern_rbf_C2(x, sigma=sigma0), kernelMatrix(rbfdot(sigma=sigma0), x))
# laplace - kernlab, then C version best
microbenchmark(kern_laplace(x, sigma=sigma0), kern_laplace_C(x, sigma=sigma0), kernelMatrix(laplacedot(sigma=sigma0), x))
# poly - R version best
microbenchmark(kern_poly(x, degree=degree0, scale=scale0, offset=offset0), kern_poly_C(x, degree=degree0, scale=scale0, offset=offset0), kernelMatrix(polydot(degree=degree0, scale=scale0, offset=offset0), x))
# tanh- R version best
microbenchmark(kern_tanh(x, scale=scale0, offset=offset0), kern_tanh_C(x, scale=scale0, offset=offset0), kernelMatrix(tanhdot(scale=scale0, offset=offset0), x))
# log - C version best
microbenchmark(kern_log(x, degree=degree0), kern_log_C(x, degree=degree0))
# quad - C version best
microbenchmark(kern_quad(x, offset=offset0), kern_quad_C(x, offset=offset0))

# It seems that those kernels that are dependent on ||x-y||_2^2 are faster in C++
# and those that depend on x%*%t(y) are faster in R.  

# for varying n 
ns <- round(seq(10, 2000, length.out=50))
d <- 10


res <- lapply(ns, function(n){
  # n <- 100
  print(paste("n: ", n))
  x <- matrix(rnorm(n*d), n, d)
  res_rbf <- cbind(kernel="rbf", melt(cast(as.data.frame(microbenchmark(kern_rbf(x, sigma=sigma0), kern_rbf_C(x, sigma=sigma0), kern_rbf_C2(x, sigma=sigma0), kernelMatrix(rbfdot(sigma=sigma0), x))), expr~., fun.aggregate="mean"))[,1:2])
  res_laplace <- cbind(kernel="laplace",melt(cast(as.data.frame(microbenchmark(kern_laplace(x, sigma=sigma0), kern_laplace_C(x, sigma=sigma0), kernelMatrix(laplacedot(sigma=sigma0), x))), expr~., fun.aggregate="mean"))[,1:2])
  res_poly <- cbind(kernel="poly",melt(cast(as.data.frame(microbenchmark(kern_poly(x, degree=degree0, scale=scale0, offset=offset0), kern_poly_C(x, degree=degree0, scale=scale0, offset=offset0), kernelMatrix(polydot(degree=degree0, scale=scale0, offset=offset0), x))), expr~., fun.aggregate="mean"))[,1:2])
  res_tanh <- cbind(kernel="tanh",melt(cast(as.data.frame(microbenchmark(kern_tanh(x, scale=scale0, offset=offset0), kern_tanh_C(x, scale=scale0, offset=offset0), kernelMatrix(tanhdot(scale=scale0, offset=offset0), x))), expr~., fun.aggregate="mean"))[,1:2])
  res_log <- cbind(kernel="log", melt(cast(as.data.frame(microbenchmark(kern_log(x, degree=degree0), kern_log_C(x, degree=degree0))), expr~., fun.aggregate="mean"))[,1:2])
  res_quad <- cbind(kernel="quad",melt(cast(as.data.frame(microbenchmark(kern_quad(x, offset=offset0), kern_quad_C(x, offset=offset0))), expr~., fun.aggregate="mean"))[,1:2])
  res <- rbind(res_rbf, res_laplace, res_poly, res_tanh, res_log, res_quad)
  res$n <- n
  return(res)
})
res <- do.call(rbind, res)

# based on ||x-y|| or t(x)%*% y ???

res$based <- c(rep("diff",4),rep("dot",2))[match(res$kernel, c("rbf","laplace", "quad", "log","tanh", "poly"))]

# kernelMatrix -> R (kernLab); kern_type -> R (my implementation);  
# kern_type_C -> C++; kern_type_C2 -> C++ vectorized

aux <- sapply(strsplit(as.character(res$expr), split="\\("), function(el) el[[1]])

aux2 <-strsplit(aux, split="_") 

res$implementation <- sapply(aux2, function(el){
  if(el[[1]]=="kernelMatrix"){
    res <- "kernlab"
  } else if(length(el)==2){
    res <- "myR"
  } else if(el[[3]]=="C"){
    res <- "myC"
  } else{
    res <- "myC2"
  }
})

#save(list=ls(), file="compareKernelImpls.RData")

p <- ggplot(res)
p <- p + geom_point(aes(x=n, y=value, colour=implementation))
p <- p + facet_grid(kernel~based, scales="free")
p



# IT SEEMS myC is best for kernels that depend on ||x-y||^2 and myR for kernels that
# depend on t(x) %*% y

