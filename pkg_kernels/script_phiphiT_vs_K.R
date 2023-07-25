# Check that features correspond to kernels
# for random  features phi, phi%*%t(phi)  should approximate  kernel K as we take more features

# w ~ norm     ->  K rbf
# w ~ laplace  ->  K quadratic/cauchy kernel
# w ~ cauchy   ->  K laplace1
# w ~ prod cauchy -> K laplace2

# for deterministic features phi, phi%*%t(phi) should equal kernel K

# polynomial feature -> polynomial kernel


# explore the domain and range of different kernels 
library(GGally)
remove(list=ls())
server <- "optimus.uv.es"
user <- "emiliano"
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_causaLearner/genData/func_getData_v2.R")

n <- 100
p <- 10

p <- 4
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
x <- simRandSEM(p, n, nodes, sigma=3, sigmaErr=1,  geU=function(y, nois, scale, constant) y)
x <- x$x
data <- as.data.frame(x)
ggpairs(data, aes(alpha = 0.4))

x <- matrix(runif(100*p), n, p)

#x <- matrix(seq(12), 4, 3, byrow=T)

# Exact: phi%*%t(phi) == K

# polynomial

offsets <- 10^seq(0,3,by=1)
scales <-  10^seq(-3,3,by=1)
degrees <- seq(1,3, by=1)



Ks_phi <- sapply(offsets, function(off) sapply(scales, function(scl) sapply(degrees, function(deg){
    # off <- offsets[1]; scl <- scales[1]; deg <- degrees[1]
    phi <- feat_poly(x, offset=off, scale=scl, degree=deg)
    K <- phi%*%t(phi)
    return(K)
  }, simplify="array"), simplify="array"), simplify="array")
Ks_kern <- sapply(offsets, function(off) sapply(scales, function(scl) sapply(degrees, function(deg) kern_poly(x, offset=off, scale=scl, degree=deg), 
                                                                                 simplify="array"), simplify="array"), simplify="array")


dim(Ks_phi); dim(Ks_kern)

plot(log(Ks_phi), log(Ks_kern))
abline(a=0, b=1, col="red")


# Approximate: phi%*%t(phi) =~ K

# w ~ norm     ->  K rbf


sigmas <- 10^seq(-10,10,by=1)  
num_f <- 1000  

Ks_phi <- sapply(sigmas, function(sd){
    # sd <- 1
    # pars <- list(sigma=sd)
    phi <- rff(x, num_f, seed=1234, p_w="rnorm2", map="cos",sigma=sd)
    K <- Real(phi%*%Conj(t(phi)))
    return(K)
  }, simplify="array")
Ks_kern <- sapply(sigmas, function(sd)  kern_rbf(x, sigma=sd), simplify="array")

dim(Ks_phi); dim(Ks_kern)

plot(Ks_phi, Ks_kern, col=rep(1:length(sigmas), rep(n^2, length(sigmas))))
legend("bottom", legend=log(sigmas,10), col=1:length(sigmas),  horiz=TRUE, cex=0.5, pch=1)
abline(a=0, b=1, col="red")

# w ~ laplace1  ->  K quadratic/cauchy kernel -> this one es not quite right 

offsets <- seq(1,10 ,length.out=10) 
offsets <- rep(2,10)
num_f <- 1000

rlaplace2 <- function(n, offset) rlaplace(n=n, m=0, s=1/(offset^0.5))
rlaplace2 <- function(n, offset) rlaplace(n=n, m=0, s=1/(offset^0.5))
kern_quad_R <- function(x, y=x, offset, degree=1){
  matNorms <- matNorm2(x, y)
  K <- 1-(matNorms^degree)/(matNorms^degree+offset)
  return(K)
}
kern_quad2 <- function(x, y=x, df){
  p <- ncol(x)
  matNorms <- matNorm2(x, y)
  K <- (1-(1/df)*matNorms)^((p+df)/2)
  return(K)
}


Ks_phi <- sapply(offsets, function(off){
  #off <- 1
  phi <- rff(x, num_f, seed=1234, p_w="rlaplace2", map="cos", offset=off)
  K <- Real(phi%*%Conj(t(phi)))
  return(K)
}, simplify="array")
Ks_kern <- sapply(offsets, function(off)  kern_quad(x, offset=off), simplify="array")

Ks_kern2 <- sapply(offsets, function(off)  kern_quad2(x, df=off), simplify="array")


dim(Ks_phi); dim(Ks_kern); dim(Ks_kern2)

plot(as.numeric(Ks_phi), as.numeric(Ks_kern), col=rep(1:length(offsets), rep(n^2, length(offsets))))
legend("bottom", legend=offsets, col=1:length(offsets),  horiz=TRUE, cex=0.5, pch=1)
abline(a=0, b=1, col="red")

plot(as.numeric(Ks_phi), as.numeric(Ks_kern2), col=rep(1:length(offsets), rep(n^2, length(offsets))))
legend("bottom", legend=offsets, col=1:length(offsets),  horiz=TRUE, cex=0.5, pch=1)
abline(a=0, b=1, col="blue", lwd=5)

# w ~ cauchy   ->  K laplace1

scales <- 10^seq(-10,10,by=1)
num_f <- 10000

Ks_phi <- sapply(scales, function(scl){
  #scl <- 1
  print(scl)
  phi <- rff(x, num_f, seed=1234, p_w="rcauchy2", map="exp", scale=scl)
  summary(as.numeric(phi))
  K <- Real(phi%*%Conj(t(phi)))
  print(sum(is.infinite(as.numeric(K))))
  return(K)
}, simplify="array")
pm <- proc.time()
Ks_kern <- sapply(scales, function(scl)  kern_laplace1(x, scale=scl), simplify="array")
proc.time() - pm # 3.13 mins for 21 1000 by 1000 matrices

dim(Ks_phi); dim(Ks_kern)

plot(Ks_phi, Ks_kern, col=rep(1:length(scales), rep(n^2, length(scales))))
legend("bottom", legend=log(scales,10), col=1:length(scales),  horiz=TRUE, cex=0.5, pch=1)
abline(a=0, b=1, col="red")
print("done")

# w ~ cauchy   ->  K laplace 2
