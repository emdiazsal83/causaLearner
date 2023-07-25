remove(list=ls())
server <- "optimus.uv.es"
user <- "emiliano"
repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_causaLearner/genData/func_getData_v2.R")
library(numDeriv) # grad, hessian, jacobian

# Check that kernel gradient norm matrices are propertly calculated
# by comparing them to numeric gradients. 

# We will do this for following kernels

# vanilla
# rbf
# laplace1
# laplace2
# quadratic (cauchy)
# log
# polynomial 

n <- 200
p <- 100
x <- matrix(runif(n*p, min=2, max=6), n, p)

# first lets test our function for norms_ijk= (xk-xi)^t %*% (xk-xj)
# note if y=NULL an z=NULL then res[1,1,1]=0

pm <- proc.time()
norms <- matNorm2Data(x)
proc.time() - pm # 0.862 for (n,p)=(200,100)
dim(norms)
norms[1,1,1]

pm <- proc.time()
norms2 <- sapply(1:n, function(k) sapply(1:n, function(j) sapply(1:n, function(i) t(x[k,]-x[i,])%*%(x[k,]-x[j,])
                                                                 , simplify="array"), simplify="array"), simplify="array")
proc.time() - pm # 96 secs for (n,p)=(200,100)
norms2[1,1,1]

smpl <- sample(prod(dim(norms)), size=100)
plot(as.numeric(norms)[smpl], as.numeric(norms2)[smpl]); abline(a=0, b=1, col="red")


# vanilla - done

pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_van, x[k,,drop=F], y=x)
  Ck <- J%*%t(J)
  return(Ck)
}, simplify="array")
proc.time() - pm # 80 secs for (n,p)=(200,100)

pm <- proc.time()
Cks2 <- kern_van_gradNorm(x)
proc.time() - pm # 0.911 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")


# rbf - done

sigma <- 0.001
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_rbf, x[k,,drop=F], y=x, sigma=sigma)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm # 83 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_rbf_gradNorm(x, sigma=sigma)
proc.time() - pm # 1.922 secs for(n,p)=(200,100)

plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")


# laplace1 - 
# seems to work for some, but numeric gradient assigns alot of zeros with Richardson method
# with "simple" method it works for all - probably an artifact of numerical derivatives with Richardson method??

scale <- 13
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_laplace1, x[k,,drop=F], y=x, scale=scale, method=c("Richardson","simple","complex")[2])
  #J <- pracma:::jacobian(kern_laplace1, x0=x[k,,drop=F], y=x, scale=scale)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm #  45.5 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_laplace1_gradNorm(x, scale=scale)
proc.time() - pm #  2.9 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")

# laplace2 - done

scale <- 5
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_laplace2, x[k,,drop=F], y=x, scale=scale)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm #  8.7 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_laplace2_gradNorm(x, scale=scale)
proc.time() - pm #  4.6 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")


# quadratic/cauchy - done

offset <- 5
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_quad, x[k,,drop=F], y=x, offset=offset)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm # 82 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_quad_gradNorm_R(x, offset=offset)
proc.time() - pm #  1.8 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")

# log - done

degree <- 3
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_log, x[k,,drop=F], y=x, degree=degree)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm #  82 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_log_gradNorm(x, degree=degree)
proc.time() - pm #  4 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")

# polynomial  - done
offset <- 3
scale <- 5
degree <- 4
pm <- proc.time()
Cks <- sapply(1:n, function(k){
  J <- jacobian(kern_poly, x[k,,drop=F], y=x, offset=offset, scale=scale, degree=degree)
  Ck <- J%*%t(J)
  return(Ck)
})
proc.time() - pm #  13.3 secs for (n,p)=(200,100)
pm <- proc.time()
Cks2 <- kern_poly_gradNorm(x, offset=offset, scale=scale, degree=degree)
proc.time() - pm #  1.3 secs for (n,p)=(200,100)
plot(as.numeric(Cks)[smpl], as.numeric(Cks2)[smpl], col=as.numeric(sapply(1:n, function(k) diag(n), simplify="array")+1)[smpl]); abline(a=0, b=1, col="red")

