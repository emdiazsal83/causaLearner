# HSIC Regression

remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R")
source("./functions.R")



# simulate non-lin non gaussian data

set.seed(1)

n <- 100
x <- runif(n)
ny <- runif(n)
y <- x^2 + ny
nsF <- constructData(x,y)
nsB <- constructData(y,x)

nTest <- 1000
x <- runif(nTest)
ny <- runif(nTest)
y <- x^2 + ny
nsTestF <- constructData(x,y)
nsTestB <- constructData(y,x)

xx <- seq(min(x), max(x), length.out=100)
yy <- seq(min(y), max(y), length.out=100)
nsSeqF <- constructData(xx,yy)
nsSeqB <- constructData(yy,xx)

# kernel ridge regression
# forward model
# quasi HSIC Regression
# HSIC regression

############################################################################################################


optF.krr  <- optLambda(trainData=nsF, learner=krr,   numFolds=5, parallel=TRUE, plot=TRUE)
optF.krr$opt
optF.qhsic <- optLambda(trainData=nsF, learner=qhsic, numFolds=5, parallel=TRUE, plot=TRUE)
optF.qhsic$opt

pm <- proc.time()
optF.hsic <- optLambda(trainData=nsF, learner=hsic, numFolds=5, parallel=TRUE, plot=TRUE)
proc.time() - pm # 5.15 mins with 8 cores (2/5 parallelization)
optF.hsic$opt


pm <- proc.time()
mus.F <- findRelevantMus(nsF, qkric, numFolds=5, numMus=20)
opts.F <- optLambdas4Mus(trainData=nsF, learner=qkric, numFolds=5, mus.F, parallel=TRUE, plot=FALSE)
plot.optLambas4Mus(opts.F, nsF, qkric, numFolds=5)
proc.time() - pm # mins with 8 cores (2/5 parallelization)

# DONE 1. adapt to include diff kernel for residuals
# DONE 2a. check qhsic corresponds to hscic with lin kernel for residuals
# DONE (they don't, they have a diff estimatior) 2b. check hsic and qhsic correspond to those from dHSIC package
# DONE 2c. using dhsic see what "median" heuristic should correspond to 
# 2d. try null hypothesis version for  tuning sigmaX and sigmaR to get a uniform p-value dist
# 3. add p-value for qhsic and nqhsic(X,R) = qhsic(X,R)/(sqrt(qhsic(X,X))*sqrt(qhsic(R,R))) to optLambda graphs 
# 4. add p-value for hsic and nhsic(X,R) = hsic(X,R)/(sqrt(hsic(X,X))*sqrt(hsic(R,R))) to optLambda graphs 
# 5. for optLambdaMu
#					 a. find "interesting" zone for nqhsic/nhsic/pavalues ie where they are not just zero or one
#					 b. for each mu, find optimal reg parm lambda in terms of the test E^2 + mu*qhsic

# 6. explore derivatives of nqhsic and nhsic, does nqhsic have a closed form? build corresponding learners
# 7. for fair type loss add mu*() + (1-mu)*()  mu in (0,1)
# 8. check that krr and qhsic closed solutions are insensitive to centering of data -> it did affect CV sol for hsic centering the Kx matrix,.... what are we doing when we center Ks? is it just centering the data?



#optB.krr  <- optLambda(nsB, krr,   numFolds=5, parallel=TRUE, plot=TRUE)
#optB.krr$opt
#optB.qhsic <- optLambda(nsB, qhsic, numFolds=5, parallel=TRUE, plot=TRUE)
#optB.qhsic$opt

#pm <- proc.time()
#optB.hsic <- optLambda(nsB, hsic, numFolds=5, parallel=TRUE, plot=TRUE)
#proc.time()-pm # 2.6 hours with 8 cores (2/5 parallelization)
#optB.hsic$opt


pm <- proc.time()
mus.B <- findRelevantMus(nsB, qkric, numFolds=5, numMus=20)
opts.B <- optLambdas4Mus(trainData=nsB, learner=qkric, numFolds=5, mus.B, parallel=TRUE, plot=FALSE)
plot.optLambas4Mus(opts.B, nsB, qkric, numFolds=5)
proc.time() - pm # mins with 8 cores (2/5 parallelization)

 
 
 

mF.krr      <- krr$learn(nsF, optF.krr$opt[[1]])
mF.qhsic    <- qhsic$learn(nsF, optF.qhsic$opt[[1]])
mF.hsic    <- hsic$learn(nsF, optF.hsic$opt[[1]])

mB.krr      <- krr$learn(nsB, optB.krr$opt[[1]])
mB.qhsic    <- qhsic$learn(nsB, optB.qhsic$opt[[1]])
mB.hsic    <- hsic$learn(nsB, optB.hsic$opt[[1]])





#################################################################################
# check if qhsic ~ hsic(vanilla kernel for residuals)
param <- getSigmaR2(nsTestF)
qhsicLoss(mF.krr, nsTestF, krr, param)
hsicLoss(mF.krr, nsTestF, krr, param) 
#they're the same if we DONT center the Kernel matrices in hsicLoss!!! (dont do Kxc = Kx%*%H, Krc = Kr%*%H, instead Kxc =Kx, Krc=Kr)
# also if we DO center the kernel matrix in qhsicloss (do Kxc = H%*%K%*%H)
resF.krr <- nsTestF$y-krr$predict(mF.krr,nsTestF)

kernelX <- do.call(param$kernelX, param$psKernX)
kernelR <- do.call(param$kernelR, param$psKernR)


Kx <- kernelMatrix(kernelX, nsTestF$x)
Kr <- kernelMatrix(kernelR, resF.krr)
N <- nrow(Kx)
H <- diag(N)-matrix(1/N,N,N)
Kxc <- Kx%*%H
Krc <- Kr%*%H


Ks <- vector("list", 2)
Ks[[1]] <- Kx
Ks[[2]] <- Kr

(1/nrow(Kxc)^2)*hsicLoss(mB.krr, nsTestB, krr, param) 
dhsic(list(nsTestB$x, resB.krr), kernel=c("kernelX","kernelR")) #takes a while to evaluate gram matrices
dhsic(K=Ks)
# de hsic estimator from dhsic function not the same as the one we calculate!!! Got to read why theirs is better
dhsic.test(list(nsTestB$x, resB.krr), kernel=c("kernelX","kernelR")) #test statisitic is just dhsic*(n^2)

dhsic(K=list(Kx,Kr))$dHSIC/(sqrt(dhsic(K=list(Kx,Kx))$dHSIC)*sqrt(dhsic(K=list(Kr,Kr))$dHSIC))

x <- rnorm(1000)
y <- 5+sin(x)

kernelX <- do.call("rbfdot", list(sigma=1/median(as.numeric(dist(x)^2))))
kernelY <- do.call("rbfdot", list(sigma=1/median(as.numeric(dist(y)^2))))
Kx <- kernelMatrix(kernelX, x)
Ky <- kernelMatrix(kernelY, y)
dhsic(K=list(Kx,Ky))$dHSIC/(sqrt(dhsic(K=list(Kx,Kx))$dHSIC)*sqrt(dhsic(K=list(Ky,Ky))$dHSIC))


#

# C median function used
#double median_bandwidth_rcpp(NumericMatrix x, int n, int d){
#  int len = n;  # len is number of rows as long as its less than 1000
#  if(n > 1000){
#    len = 1000;
#  }
#  double xnorm = 0.0;
#  int lentot = len*(len+1)/2-len; # (n-1)*(n-1)/1 is number
#  int middle = lentot/2;
#  NumericVector bandvec(lentot);
#  int count = 0;
#  for(int i = 0; i < len; ++i){
#    int j = i+1;
#    while(j < len){
#      for(int l = 0; l < d; ++l){
#        xnorm += pow(x(i,l)-x(j,l), 2.0);
#      }
#      bandvec[count] = xnorm;
#      xnorm = 0.0;
#      ++j;
#      ++count;
#    }
#  }
#  NumericVector v = clone(bandvec);
#  std::nth_element(v.begin(), v.begin() + middle, v.end());
#  double bandwidth = v[middle];
#  bandwidth = sqrt(bandwidth*0.5);
#  return bandwidth;
#}

n <- 10
d <- 3
X <- matrix(rnorm(n*d, mean=3, sd=4),n,d)


lentot <- ((n-1)*n)/2 #always even!
length(as.numeric(dist(X)))


dists2 <- as.numeric(dist(X))^2
dists2 <- sort(dists2)
sqrt(0.5*median(dists2)) #its a bit diff, is it just numerical error?
sqrt(0.5*dists2[(lentot/2)])
sqrt(0.5*dists2[(lentot/2)+1])
dHSIC:::median_bandwidth_rcpp(X,n,d) #they use the mid+1 point as median
sigma <- sqrt(0.5*dists2[(lentot/2)+1])



K1 <- dHSIC:::gaussian_grammat_rcpp(X, sigma, n, d)
K2 <- kernelMatrix(rbfdot(1/(2*sigma^2)), X)
K1/K2
# so sigma = 1/(2*bw^2) = 1/(2*(sqrt(0.5*median(dists^2)))^2) 
# = 1/(median(dists^2))

#################################################################################




par(mfrow=c(2,2))
plot(nsTestF$x,nsTestF$y)
lines(nsSeqF$x, krr$predict(mF.krr, nsSeqF), col="red", lwd=2)
lines(nsSeqF$x, qhsic$predict(mF.qhsic, nsSeqF), col="blue", lwd=2)
lines(nsSeqF$x, hsic$predict(mF.hsic, nsSeqF), col="green", lwd=2)
plot(nsTestB$x,nsTestB$y)
lines(nsSeqB$x, krr$predict(mB.krr,nsSeqB), col="red", lwd=2)
lines(nsSeqB$x, qhsic$predict(mB.qhsic,nsSeqB), col="blue", lwd=2)
lines(nsSeqB$x, hsic$predict(mB.hsic,nsSeqB), col="green", lwd=2)

resF.krr <- nsTestF$y-krr$predict(mF.krr,nsTestF)
resF.qhsic <- nsTestF$y-qhsic$predict(mF.qhsic,nsTestF)
resF.hsic <- nsTestF$y-hsic$predict(mF.hsic,nsTestF)
plot(nsTestF$x, resF.krr, col="red", ylim=range(c(resF.krr, resF.qhsic, resF.hsic)))
lines(nsTestF$x, resF.qhsic, col="blue", type="p")
lines(nsTestF$x, resF.hsic, col="green", type="p")

resB.krr <- nsTestB$y-krr$predict(mB.krr,nsTestB)
resB.qhsic <- nsTestB$y-qhsic$predict(mB.qhsic,nsTestB)
resB.hsic <- nsTestB$y-hsic$predict(mB.hsic,nsTestB)
plot(nsTestB$x, resB.krr, col="red", ylim=range(c(resB.krr, resB.qhsic, resB.hsic)))
lines(nsTestB$x, resB.qhsic, col="blue",type="p")
lines(nsTestB$x, resB.hsic, col="green",type="p")

par(mfrow=c(2,2))
plot(mF.krr$alpha, mF.qhsic$alpha)
plot(mF.hsic$alpha, mF.qhsic$alpha)
plot(mB.krr$alpha, mB.qhsic$alpha)
plot(mB.hsic$alpha, mB.qhsic$alpha)



# norm of difference
# Forward
sqrt(t(mF.krr$alpha-mF.qhsic$alpha)%*%(mF.krr$alpha-mF.qhsic$alpha))
sqrt(t(mF.qhsic$alpha-mF.hsic$alpha)%*%(mF.qhsic$alpha-mF.hsic$alpha))
sqrt(t(mF.krr$alpha-mF.hsic$alpha)%*%(mF.krr$alpha-mF.hsic$alpha))
# Backward
sqrt(t(mB.krr$alpha-mB.qhsic$alpha)%*%(mB.krr$alpha-mB.qhsic$alpha))
sqrt(t(mB.qhsic$alpha-mB.hsic$alpha)%*%(mB.qhsic$alpha-mB.hsic$alpha))
sqrt(t(mB.krr$alpha-mB.hsic$alpha)%*%(mB.krr$alpha-mB.hsic$alpha))

# build a diagram showing, for i in seq(0,1,length.out=100), the following metrics for solutions
# alpha_i = i*   

#save(list=ls(), file="backForwExp1.RData")
#load( file="backForwExp1.RData")
