print("func_kernel pkg")
# kernel function package

library(Rcpp)
library(multicool) #multinom
library(rmutil) # rlaplace
library(pracma) #Norm
library(abind) #adrop


# Following kernels implemented so far:

# A. Based on ||x-y||^2_2

# 1. Rbf
# 2. Laplace
# 3. Rational quadratic
# 4. Log 

# B. Based on t(x) %*% y

# 1. Hyperbolic tan
# 2. Polynomial 

# C. other

# For different kernels go to:
# http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/


###################################*
# kernels based on ||x - y||^2_2
###################################*

######*
# R
######*
matNorm2 <- function(x, y=NULL){
  n <- nrow(x)
  d <- ncol(x)
  onesN <- rep(1, n)
  xtx <- apply(x*x, 1, sum)
  
  if(!is.null(y)){
    m <- nrow(y)
    onesM <- rep(1, m)
    yty <- apply(y*y, 1, sum)
  } else{
    onesM <- onesN
    yty <- xtx
    y <- x
  }
  
  matNorms <- xtx%*%t(onesM) + onesN%*%t(yty) - 2*x%*%t(y)  
  return(matNorms)
}

matNormP <- function(x, y=NULL, p){
  
  if(is.null(y)){
    y <- x
  }
  
  matNorms <- apply(x, 1, function(xCol) apply(y, 1, function(yCol) Norm(xCol-yCol, p=p)))
  
  return(t(matNorms))
}

# function for norms_ijk= (xk-xi)^t %*% (xk-xj)
# note if y=NULL an z=NULL then res[1,1,1]=0
matNorm2Data <- function(x, y=NULL, z=NULL){
  n <- nrow(x)
  d <- ncol(x)
  onesN <- rep(1, n)
  
  
  if(!is.null(y)){
    m <- nrow(y)
    onesM <- rep(1, m)
    xyt <- x%*%t(y)
    if(!is.null(z)){
      k <- nrow(z)
      xzt <- x%*%t(z)
      yzt <- y%*%t(z)
      zzt <- z%*%t(z)
    } else{
      k <- n
      xzt <- x%*%t(x)
      yzt <- y%*%t(x)
      zzt <- xzt
    }
  } else{
    onesM <- onesN
    xyt <- x%*%t(x)
    if(!is.null(z)){
      k <- nrow(z)
      xzt <- x%*%t(z)
      yzt <- xzt
      zzt <- z%*%t(z)
    } else{
      k <- n
      xzt <- x%*%t(x)
      yzt <- xzt
      zzt <- xzt
    }
  }
  
  Ik <- diag(k)
  
  
  matNorms <- sapply(1:k, function(i){
    ek <- Ik[,i]
    mat <- onesN%*%zzt[i,i,drop=FALSE]%*%t(onesM)- xzt%*%ek%*%t(onesM)- onesN%*%t(ek)%*%t(yzt)+ xyt
  }, simplify="array")
  
  return(matNorms)
}

# p <- 3; n <- 4; m <- 5; k <- 7; a <- matrix(seq(n*p), n, p, byrow=T); b <- matrix(seq(m*p), m, p, byrow=T); c <- matrix(seq(k*p), k, p, byrow=T) 
# d <- matNorm2Data(a, b, c)
# d[4,,]

# group 
# indicator <- c(1,0,1,1,0,1)
kern_bin_R <- function(x, y=x, num){
  # x <- rnorm(10); y=x; num <- 3
  featsX <- feat_bin(x, num)
  featsY <- feat_bin(y, num)
  K <- featsX %*% t(featsY)
  return(K)
}

# linear - uncentered
kern_van_R <- function(x, y=x){
  matNorms <- matNorm2(x, y)
  K <- -matNorms 
  return(K)
}

kern_van_gradNorm_R <- function(x, y=x, z=x){
  Cpart1 <- 1
  Cpart2 <- matNorm2Data(x, y, z)
  C <- 4*Cpart1*Cpart2
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}



# RBF
kern_rbf_R <- function(x, y=x, sigma){
  matNorms <- matNorm2(x, y)
  K <- exp(-sigma*matNorms)
  return(K)
}

kern_rbf_gradNorm_R <- function(x, y=x, z=x, sigma){
  
  Cpart1 <- kernelDataArray("kern_rbf", x=x, y=y, z=z, pars=list(sigma=sigma))
  Cpart2 <- matNorm2Data(x, y, z)
  C <- 4*sigma^2*Cpart1*Cpart2
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}

# normalized kernel matrix with respect to second (y) argument
kern_rbf_nrml_R <- function(x, y=x, sigma, lim_max, lim_min){
  # setwd("/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/")
  # source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")
  # n <- 100; p <- 10; x <- rnorm(n*p); x <- matrix(x, n, p); x <- apply(x, 2, norml)
  # n <- 110; y <- rnorm(n*p); y <- matrix(y, n, p); y <- apply(y, 2, norml)
  # sigma <- 10; lim_max <- 0.5; lim_min <- -0.5
  p <- ncol(x)
  mus <- x
  sigma2 <- 1/(2*sigma)
  nrmlCnst <- pnorm(lim_max, mean=mus, sd=sqrt(sigma2)) - pnorm(lim_min, mean=mus, sd=sqrt(sigma2)) 
  nrmlCnst <- apply(nrmlCnst, 1, prod)
  nrmlCnst <- nrmlCnst*((pi/(sigma))^(p/2))
  K <- kern_rbf_R(x, y, sigma)
  K <- K*nrmlCnst
  # K[1:5, 1:6]
  return(K)
}


# RBF - 2 sigmas
kern_rbf2_R <- function(x1, x2=x1, y1, y2=y1, sigma1, sigma2){
  
  K1 <- kern_rbf_R(x1, y1, sigma1)
  K2 <- kern_rbf_R(x2, y2, sigma2)
  K <- K1*K2
  return(K)
}

# RBF - mult sigmas
kern_rbf_indSig_R <- function(x, y=x, sigmas){
  # n <- 100; p <- 2; x <- y <- matrix(rnorm(n*p), n, p); sigmas <- c(1, 10)
  Ks <- mcmapply(FUN=function(i,sigma) kern_rbf_R(x[,i,drop=F], y[,i,drop=F], sigma), i=1:ncol(x), sigma=sigmas, mc.cores=1,SIMPLIFY="array")
  K <- apply(Ks, c(1,2), function(v) prod(v))
  return(K)
}


# Laplace, norm 1

kern_laplace1_R <- function(x, y=x, scale){
  matNorms <- matNormP(x, y, p=1)
  K <- exp(-(1/scale)*matNorms)
  return(K)
}

kern_laplace1_gradNorm_R <- function(x, y=x, z=x, scale){
  
  
  # ni <- 10
  # nj <- 20
  # nk <- 30
  # p <- 5
  # x <- matrix(rnorm(ni*p), ni, p)
  # y <- matrix(rnorm(nj*p), nj, p)
  # z <- matrix(rnorm(nk*p), nk, p)
  
  ni <- nrow(x)
  nj <- nrow(y)
  nk <- nrow(z)
  p <- ncol(x)
  
  
  # pm <- proc.time()
  Cpart1 <- kernelDataArray(kernel="kern_laplace1", x=x, y=y, z=z, pars=list(scale=scale))
  # proc.time() - pm  # 2 secs
  # 
  # pm <- proc.time()
  # Cpart2a <-  sapply(1:nrow(x), function(i) sapply(1:nrow(y), function(j) sapply(1:nrow(z), function(k) t(((z[k,]>x[i,])*1))%*%((z[k,]>y[j,])*1), simplify="array"), simplify="array"), simplify="array")
  # proc.time() - pm # 72 secs
  # 
  # # xi <- x[1,]; yi <- y[1,]; zi <- z[1,]
  # pm <- proc.time()
  # Cpart2b <-  apply(x, 1, function(xi) apply(y, 1, function(yi) apply(z, 1, function(zi) t(((zi>xi)*1))%*%((zi>yi)*1))))
  # dim(Cpart2b) <- c(nk, nj, ni) 
  # proc.time() - pm # 66 secs
  
  
  pm <- proc.time()
  signsX <- apply(x, 1, function(xi) apply(z, 1, function(zi) (zi>=xi)*1-(zi<xi)*1))
  signsY <- apply(y, 1, function(yi) apply(z, 1, function(zi) (zi>=yi)*1-(zi<yi)*1))
  dim(signsX) <- c(p, nk, ni)
  dim(signsY) <- c(p, nk, nj)
  
  #k <- 1
  Cpart2c <- sapply(1:nk, function(k) t(adrop(signsX[,k,,drop=F], drop=2))%*%adrop(signsY[,k,,drop=F], drop=2), simplify="array")
  #Cpart2c <- aperm(Cpart2c, c(3,2,1))
  proc.time() - pm #0.561
  
  
  # dim(Cpart2a)
  # dim(Cpart2b)
  # dim(Cpart2c)
  # summary(as.numeric(Cpart2a))
  # summary(as.numeric(Cpart2b))
  # summary(as.numeric(Cpart2c))
  # smpl <- sample(prod(dim(Cpart2)), size=100)  
  # plot(as.numeric(Cpart2a)[smpl], as.numeric(Cpart2b)[smpl]); abline(a=0, b=1, col="red")
  # plot(as.numeric(Cpart2a)[smpl], as.numeric(Cpart2c)[smpl]); abline(a=0, b=1, col="red")
  
  C <- 1/(scale^2)*Cpart1*Cpart2c
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}


# Laplace, norm 2
kern_laplace2_R <- function(x, y=x, scale){
  matNorms <- matNorm2(x, y)
  # argh! this is horrible
  matNorms[matNorms<0] <- 0
  K <- exp(-(1/scale)*matNorms^0.5)
  return(K)
}

kern_laplace2_gradNorm_R <- function(x, y=x, z=x, scale){
  
  m <- nrow(y)
  n <- nrow(x)
  onesM <- rep(1, m)
  onesN <- rep(1, n)
  
  Cpart1 <- kernelDataArray("kern_laplace2", x=x, y=y, z=z, pars=list(scale=scale))
  Cpart2 <- matNorm2Data(x, y, z)
  
  
  Cpart3 <- matNorm2Data(x=x, y=NULL, z=z)
  Cpart3 <- apply(Cpart3, 3, function(mat) onesM%*%t(diag(mat)))
  dim(Cpart3) <- dim(Cpart1)
  Cpart3 <- 1/Cpart3
  Cpart4 <- matNorm2Data(x=y, y=NULL, z=z)
  Cpart4 <- apply(Cpart4, 3, function(mat) diag(mat)%*%t(onesN))
  dim(Cpart4) <- dim(Cpart1)
  Cpart4 <- 1/Cpart4
  
  
  C <- 1/(scale^2)*Cpart1*Cpart2*(Cpart3*Cpart4)^0.5
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}

# Quadratic
kern_quad_R <- function(x, y=x, offset, degree=1){
  matNorms <- matNorm2(x, y)
  K <- 1-(matNorms^degree)/(matNorms^degree+offset)
  return(K)
}

kern_quad_gradNorm_R <- function(x, y=x, z=x, offset, degree=1){
  
  Cpart1 <- kernelDataArray("kern_quad", x=x, y=y, z=z, pars=list(offset=offset))
  Cpart1 <- Cpart1*Cpart1
  Cpart2 <- matNorm2Data(x, y, z)
  C <- (4/offset^2)*Cpart1*Cpart2
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}

# approx normalized kernel matrix with respect to second (y) argument
# using approx fourier features
kern_quad_nrml_R <- function(x, y=x, offset, lim_max, lim_min){
  # setwd("/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/")
  # source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")
  # n <- 100; p <- 10; x <- rnorm(n*p); x <- matrix(x, n, p); x <- apply(x, 2, norml)
  # n <- 110; y <- rnorm(n*p); y <- matrix(y, n, p); y <- apply(y, 2, norml)
  # offset <- 1.5; lim_max <- 0.5; lim_min <- -0.5

  
  py <- ncol(y)
  numf <- 1000
  parsPdf <- list(offset=offset)
  parsPdf$n <- num_f*py
  set.seed(1234)
  w <- do.call("rlaplace2", parsPdf)
  w <- matrix(w, py, numf)
  b <- runif(numf, 0, 2*pi)
  #b <- matrix(runif(numf, 0, 2*pi), py, numf, byrow=T)
  intPhiy <- sapply(1:ncol(w), function(i) intAnalyticRFF(py, w=w[,i], b=b[i]))
  intPhiy <- matrix(sqrt(2/num_f)*intPhiy, nrow(y), numf, byrow=T)
  phiy <- rff(x, numf, seed=1234, p_w="rlaplace2", map="cos", offset=offset)
  K <- phiy %*% t(intPhiy)
  # K[1:5,1:6]
  
  return(K)
}


# Log
kern_log_R <- function(x, y=x, degree){
  matNorms <- matNorm2(x, y)
  K <- -log(matNorms^degree+1)
  return(K)
}

kern_log_gradNorm_R <- function(x, y=x, z=x, degree){
  # don't change! the gradNorm log kernel is composed of quad kernel!
  # n <- 15; m <- 10; k <- 5; p <- 4
  # x <-matrix(rnorm(n*p), n, p); y <-matrix(rnorm(m*p), m, p); z <-matrix(rnorm(k*p), k, p)
  m <- nrow(y)
  n <- nrow(x)
  onesM <- rep(1, m)
  onesN <- rep(1, n)
  Cpart1 <- kernelDataArray("kern_quad", x=x, y=y, z=z, pars=list(offset=1, degree=degree))
  Cpart2 <- matNorm2Data(x, y, z)
  Cpart3 <- matNorm2Data(x=x, y=NULL, z=z)
  Cpart3 <- apply(Cpart3, 3, function(mat) onesM%*%t(diag(mat)))
  dim(Cpart3) <- dim(Cpart1)
  Cpart3 <- Cpart3^(degree-1)
  Cpart4 <- matNorm2Data(x=y, y=NULL, z=z)
  Cpart4 <- apply(Cpart4, 3, function(mat) diag(mat)%*%t(onesN))
  dim(Cpart4) <- dim(Cpart1)
  Cpart4 <- Cpart4^(degree-1)
  C <- 4*degree^2*Cpart1*Cpart2*Cpart3*Cpart4
  #Chat <- apply(C, c(1,2), sum)
  return(C)
}

######*
# C++
######*

# RBF
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

# RBF - 2 sigmas
kern_rbf2_CR <- function(x1, x2=x1, y1, y2=y1, sigma1, sigma2){
  K1 <- kern_rbf_C(x1, y1, sigma1)
  K2 <- kern_rbf_C(x2, y2, sigma2)
  K <- K1*K2
  return(K)
}

# RBF - 2 sigmas
cppFunction('NumericMatrix kern_rbf2_Caux(NumericMatrix x1, NumericMatrix x2, NumericMatrix y1, NumericMatrix y2, double sigma1, double sigma2) {
            int n = x1.nrow();
            int m = y1.nrow();
            int d1 = x1.ncol();
            int d2 = x2.ncol();
            NumericMatrix out(n, m);
            
            for (int i = 0; i < n; i++) {
            
            for (int j = 0; j < m; j++) {
            double total1 = 0;
            for(int k = 0;  k < d1; k++){
              total1 += pow(x1(i, k)-y1(j,k),2.0);
            }
            double total2 = 0;
            for(int k = 0;  k < d2; k++){
              total2 += pow(x2(i, k)-y2(j,k),2.0);
            }
            out(i,j) = exp(-sigma1*total1-sigma2*total2);
            }
            }
            return out;
            }')
kern_rbf2_C <- function(x1, x2=x1, y1, y2=y1, sigma1, sigma2) kern_rbf2_Caux(x1, x2, y1, y2, sigma1, sigma2)

# Laplace
cppFunction('NumericMatrix kern_laplace2_Caux(NumericMatrix x, NumericMatrix y, double scale) {
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
            out(i,j) = exp(-(1/scale)*pow(total,0.5));
            }
            }
            return out;
            }')
kern_laplace2_C <- function(x, y=x, scale) kern_laplace2_Caux(x, y, scale)

# Quadratic
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
kern_quad_C <- function(x, y=x, offset, degree=1) kern_quad_Caux(x, y, offset, degree)

# Log
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

  ###################################*
  # kernels based on t(x) %*% y
  ###################################*
  
  ######*
  # R
  ######*
  
  # linear
  kern_lin_R <- function(x, y=x, offset){
    K <- x%*%t(y)+offset
    return(K)
  }
  
  # tanh
  kern_tanh_R <- function(x, y=x, scale, offset){
    K <- tanh(scale*(x%*%t(y))+offset)
    return(K)
  }
  
  # polynomial
  kern_poly_R <- function(x, y=x, degree, scale, offset){
    K <- (scale*(x%*%t(y))+offset)^degree
    return(K)
  }
  
  kern_poly_gradNorm_R <- function(x, y=x, z=x, degree, scale, offset){
    # x <- matrix(rnorm(10,3),10,3); y=matrix(rnorm(11,3),11,3); z=matrix(rnorm(12,3),12,3); degree=2; scale=3; offset=4
    Cpart1 <- kernelDataArray("kern_poly", x=x, y=y, z=z, pars=list(degree=degree-1, scale=scale, offset=offset))
    Cpart2 <- sapply(1:nrow(z), function(k) x%*%t(y), simplify="array")
    C <- (scale*degree)^2*Cpart1*Cpart2
    #Chat <- apply(C, c(1,2), sum)
    return(C)
  }


######*
# C++
######*

# linear
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

# tanh
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

# polynomial
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

###################################*
# others
###################################*

###################################*
# Switch between R or C
###################################*

# RBF
kern_rbf <- function(x, y=x, sigma) kern_rbf_R(x, y, sigma)
kern_rbf_gradNorm <- function(x, y=x, z=x, sigma) kern_rbf_gradNorm_R(x, y, z, sigma)
kern_rbf2 <- function(x1, x2=x1, y1, y2=y1, sigma1, sigma2) kern_rbf2_R(x1, x2, y1, y2, sigma1, sigma2) # still have to test that this is the fastest vs kern_rbf2_R and kern_rbf2_RC
kern_rbf_nrml <- function(x, y=x, sigma, lim_max, lim_min) kern_rbf_nrml_R(x, y, sigma, lim_max, lim_min) 
kern_rbf_indSig <- function(x, y=x, sigmas) kern_rbf_indSig_R(x, y, sigmas)


# laplace 1
kern_laplace1 <- function(x, y=x, scale) kern_laplace1_R(x, y, scale)
kern_laplace1_gradNorm <- function(x, y=x, z=x, scale) kern_laplace1_gradNorm_R(x, y, z, scale)

# laplace 2
kern_laplace2 <- function(x, y=x, scale) kern_laplace2_C(x, y, scale)
kern_laplace2_gradNorm <- function(x, y=x, z=x, scale) kern_laplace2_gradNorm_R(x, y, z, scale)

# quadratic
kern_quad <- function(x, y=x, offset, degree=1) kern_quad_R(x, y, offset, degree)
kern_quad_gradNorm <- function(x, y=x, z=x, offset, degree=1) kern_quad_gradNorm_R(x, y, z, offset, degree)
kern_quad_nrml <- function(x, y=x, offset, lim_max, lim_min) kern_quad_nrml_R(x, y, offset, lim_max, lim_min) 


# log
kern_log <- function(x, y=x, degree) kern_log_R(x, y, degree)
kern_log_gradNorm <- function(x, y=x, z=x, degree) kern_log_gradNorm_R(x, y, z, degree)

# vanilla
kern_van <- function(x, y=x) kern_van_R(x, y)
kern_van_gradNorm <- function(x, y=x, z=x) kern_van_gradNorm_R(x, y, z) 

# group
kern_bin <- function(x, y=x, num) kern_bin_R(x, y, num)
kern_bin_gradNorm <- function(x, y=x, z=x, num) NA 


# linear
kern_lin <- function(x, y=x, offset) kern_lin_R(x, y, offset)


# tanh
kern_tanh <- function(x, y=x, scale, offset) kern_tanh_R(x, y, scale, offset)

# polynomial
kern_poly <- function(x, y=x, degree, scale, offset) kern_poly_R(x, y, degree, scale, offset)
kern_poly_gradNorm <- function(x, y=x, z=x, degree, scale, offset) kern_poly_gradNorm_R(x, y, z, degree, scale, offset)

kernelMatrix <- function(kernel, x, y=x, pars=list()){
  if(length(unique(sapply(pars, length)))>1) stop("all parameter vectors must be same length")
  if(length(pars)==0){
    pars$x <- x
    pars$y <- y
    K <- do.call(kernel, pars)
  } else{
    Ks <- sapply(1:length(pars[[1]]), function(i){
      # i <- 2
      parsAux <- lapply(pars, function(par) par[i])
      parsAux$x <- x
      parsAux$y <- y
      K <- do.call(kernel, parsAux)
      return(K)
    }, simplify="array")
    K <- apply(Ks, c(1,2), sum)
  }
  return(K)
}

kernelDataArray <- function(kernel, x, y=x, z=x, pars=list()){
  # kernel <- "kern_rbf"; x <- matrix(rnorm(400), 200, 2); pars=list(sigma=1)
  
  pars$x <- x
  pars$y <- z
  Knk <- do.call(kernel, pars)
  
  if(!is.null(y)){
    pars$x <- y
    pars$y <- z
    Kmk <- do.call(kernel, pars)
  } else{
    Kmk <- Knk
  }
  
  k <- nrow(z)
  n <- nrow(x)
  m <- nrow(y)
  Ik <- diag(k)
  onesN <- rep(1, n)
  onesM <- rep(1, m)
  
  Eks <- sapply(1:k, function(i){
   ek <- Ik[,i]
   mat <- (Knk%*%ek%*%t(onesM))*(onesN%*%t(ek)%*%t(Kmk))  
   return(mat)
  }, simplify="array")
  
  
  
  return(Eks)
}

kernelDataArray2 <- function(K, x){
  # kernel <- "kern_rbf"; x <- matrix(rnorm(400), 200, 2); pars=list(sigma=1)
  
  n <- nrow(x)
  Ik <- diag(n)
  onesN <- rep(1, n)
  
  Eks <- sapply(1:n, function(i){
    ek <- Ik[,i]
    mat <- (K%*%ek%*%t(onesN))*(onesN%*%t(ek)%*%K)  
    return(mat)
  }, simplify="array")
  
  
  E <- apply(Eks, c(1,2) , sum)
  
  return(E)
}

# this is the Cks matrix
kernelGradNormMatrix <- function(kernel, x, y=x, z=x, pars=list()){
  pars$x <- x
  pars$y <- y
  pars$z <- z
  K <- do.call(paste(kernel,"gradNorm", sep="_"), pars)
  return(K)
}

####################################################################################################*
# Random Fourier features alternative
####################################################################################################*

# integrate a random fourier feature with respect to x
# py-dimension of y data matrix (not phiy)
# for normalizing RFF-kernels to intepret as probability
# int (cos wy + b) = prod_{i=1}^p (1/w_i) * sum_{i=1}^{2^p} (-1)^i cos(0.5 w^T a_i + b) if p even
# int (cos wy + b) = prod_{i=1}^p (1/w_i) * sum_{i=1}^{2^p} (-1)^i sin(0.5 w^T a_i + b) if p odd
# where ai's are the 2^p ways of forming a vector of size p with 0's and 1's
intAnalyticRFF <- function(py, w, b){
  # py <- 2; w <- rnorm(2); b <- rnorm(1)
  # py <- 2; w <- c(0.47143516, -1.19097569); b <- c(1.43270697)
  ais <- sapply(1:(2^py), intToBitss, maxInt=2^py)
  ais <- sapply(strsplit(ais, split=""), function(el) c(1,-1)[as.numeric(el)+1], simplify="array")
  if(is.null(dim(ais))) ais <- matrix(ais, nrow=1)
  summandsOdd  <- sin(0.5*apply(ais*w, 2, sum)+b)
  summandsEven <- cos(0.5*apply(ais*w, 2, sum)+b)
  coefs <- (-1)^floor(py/2)*apply(ais,2,prod)
  summands <- switch(c("even","odd")[py%%2+1], even=summandsEven, odd=summandsOdd)
  res <- (1/prod(w))*sum(coefs*summands)
  return(res)
}

# ref: eric strobl's function from: https://rdrr.io/github/ericstrobl/RCIT/man/random_fourier_features.html
# and: https://pergamos.lib.uoa.gr/uoa/dl/frontend/file/lib/default/data/1324573/theFile
# and: https://www.cs.cmu.edu/~dsutherl/papers/rff_uai15.pdf
rff <- function(x, num_f, seed=NA, p_w, map=c("cos","cosSin","expComp","exp"), ...){
  # x <- rnorm(100); p_w <- "rnorm2"; num_f <- 4; seed <- 1234; sigma <- 1
  # pars <- list(sigma=sigma)
  if(is.na(seed)) seed <- NULL
  pars <- list(...)
  #print(pars)
  if(is.null(dim(x))) x <-  as.matrix(x)  
  
  
  n <- nrow(x)
  p <- ncol(x)
  
  parsPdf <- expand.grid(pars)
  numRKHS <- nrow(parsPdf)
  
  # total number of fourier features: one set of num_f per stacked phi
  numf <- num_f*numRKHS
  
  set.seed(seed)
  #parsPdf <- pars
  
  # for each RKHS we need to simulate num_f vectors w of size p
  parsPdf$n <- num_f*p
  
  
  
  # w <- matrix(seq(6), 3, 2)
  #w <- matrix( do.call(p_w, parsPdf), p, num_f)
  w <- apply(parsPdf, 1, function(row) do.call(p_w, as.list(row)))
  
  # we re-shape w so that rffs of diff RKHS are stacked (numf=num_f*RKHS)
  w <- matrix(w, p, numf)
  #print(dim(w))
  #print(w[,1:4])
  
  #b <- matrix(runif(num_f, 0, 2*pi), n, num_f, byrow=T)
  b <- matrix(runif(numf, 0, 2*pi), n, numf, byrow=T)
  #print(dim(b))
  #print(b[,1:4])
  
  im <- complex(real=0, imaginary=-1)
  res <- switch(map, 
                cos= sqrt(2/numf)*(cos(x%*%w+b)),
                cosSin = sqrt(1/numf)*cbind(sin(x%*%w), cos(x%*%w)),  
                exp = sqrt(1/numf)*(exp(im*x%*%w)))  
                
  return(res)
}

# for rbf kernel
rnorm2 <- function(n, sigma) sqrt(2*sigma)*rnorm(n=n, mean=0, sd=1)
# for quadratic kernel
# this is still not correct parametrization as phi %*% t(phi) not converging to Kquad as num_f -> Inf
rlaplace2 <- function(n, offset) rlaplace(n=n, m=0, s=1/(offset^0.5))
# for laplace kernel
rcauchy2 <- function(n, scale) rcauchy(n=n, location=0, scale=1/scale)




####################################################################################################*
# Feature maps
####################################################################################################*




# need tree implementation for this, a tree will be a list where first element is the node
# and the rest are trees (leafs are lists with only node element)

# return all paths from root to each leaf
getBranch <- function(tree, aux=as.character()){
  
  if(length(tree)==1){
    return(c(aux, tree$node))
  } else{
    
    res <- as.character()
    for(i in 2:(length(tree))){
      # i <- 2
      res <- c(res, getBranch(tree=tree[[i]], aux=c(aux,tree$node)))
    }
    return(res)
  }
  
}


# myTree <- list(node="R", list(node=0, list(node=0, list(node=2)), list(node=1, list(node=1)), list(node=2, list(node=0))), 
#                          list(node=1, list(node=0, list(node=1)), list(node=1, list(node=0))), 
#                          list(node=2, list(node=0, list(node=0))))
# 
# matrix(getBranch(myTree), ncol=4, byrow=T)

# Build a tree representing all discrete points on m-dimensional simplex equalling n. 
# i.e. The set S (k1,...,k2) in (N union 0)^m such that  k_1+...+k_m=n
buildDiscreteSimplexNtree <- function(m, n, node="R"){

  if(m==0){
    tree <- list()
    tree$node <- node
    return(tree)
  }
  
  if(m==1){
    tree <- list()
    tree$node <- node
    tree[[paste("br", 2, sep="_")]] <- buildDiscreteSimplexNtree(m=m-1, n=0, node=n)
    return(tree)
  } else{
    tree <- list()
    tree$node <- node
    for(i in 0:n){
      # i <- 0
      tree[[paste("br", i+1, sep="_")]] <- buildDiscreteSimplexNtree(m=m-1, n=n-i, node=i)
    }
    return(tree)
  }
  
}

#myTree2 <- buildDiscreteSimplexNtree(m=3, n=2)
#getBranch(myTree2)
#matrix(getBranch(myTree2), ncol=4, byrow=T)


# polynomial
feat_poly <- function(x, degree, scale, offset){
  # n <- 20; d <-3; x <- matrix(rnorm(n*d), n, d); degree <- 2; scale <-1; offset <- 2
  p <- ncol(x)
  # T = size of set S where
  # S = positive integers or zero k_1,...,k_{p+1} such that k_1+...+k_{p+1}=degree
  ks <- matrix(getBranch(buildDiscreteSimplexNtree(m=p+1, n=degree)), ncol=p+2, byrow=T) # a T by p+1 matrix
  ks <- matrix(as.numeric(ks[,2:ncol(ks)]), nrow=nrow(ks), ncol=p+1)
  # table(apply(ks, 1, sum))
  dataMat <- cbind(sqrt(scale)*x, sqrt(offset))
  feats <- apply(dataMat, 1, function(xi){
    # i <- 1; xi <- dataMat[i,]
    featRow <- apply(ks, 1, function(k){
      # j <- 2; k <- ks[j,]
      featj <- sqrt(multicool:::multinom(x=k, counts=TRUE))*prod((xi^k))
      return(featj)
     })
    return(featRow)
    })
  return(t(feats))
} 

#n <- 20; d <-1; x <- matrix(rnorm(n*d), n, d); degree <- 2; scale <-1; offset <- 0
#y <- feat_poly(x, degree, scale, offset)

feat_toy <- function(x){
  feats <- cbind(x[,1], x[,1]^2)
  return(feats)
}

feat_van <- function(x) x

feat_bin <- function(x, num){
  # x <- rnorm(10); num <- 3
  qs <- quantile(x, probs=seq(0,1,length.out=num+1))
    belong <- findInterval(x, qs, rightmost.closed=T)
  # convert to one-hot encoding
  feats <- t(sapply(belong, function(i){
    row <- rep(0,length(qs)-1)
    row[i] <- 1
    return(row)
  }))
  return(feats)
}
