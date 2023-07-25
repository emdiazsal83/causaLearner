# HSIC Regression

remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R")
source("./functions.R")



# simulate non-lin non gaussian data

set.seed(1)

n <- 100
x <- runif(n)
ny <- runif(n)
y <- exp(x^2)*ny
nsF <- constructData(x,log(y))
nsB <- constructData(log(y),x)

nTest <- 1000
x <- runif(nTest)
ny <- runif(nTest)
y <- exp(x^2)*ny
nsTestF <- constructData(x,log(y))
nsTestB <- constructData(log(y),x)

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
proc.time() - pm # 9 hours with 8 cores (2/5 parallelization)
optF.hsic$opt




optB.krr  <- optLambda(nsB, krr,   numFolds=5, parallel=TRUE, plot=TRUE)
optB.krr$opt
optB.qhsic <- optLambda(nsB, qhsic, numFolds=5, parallel=TRUE, plot=TRUE)
optB.qhsic$opt

pm <- proc.time()
optB.hsic <- optLambda(nsB, hsic, numFolds=5, parallel=TRUE, plot=TRUE)
proc.time()-pm # 2.6 hours with 8 cores (2/5 parallelization)
optB.hsic$opt




mF.krr      <- krr$learn(nsF, optF.krr$opt[[1]])
mF.qhsic    <- qhsic$learn(nsF, optF.qhsic$opt[[1]])
#mF.hsic    <- hsic$learn(nsF, optF.hsic$opt[[1]])

mB.krr      <- krr$learn(nsB, optB.krr$opt[[1]])
mB.qhsic    <- qhsic$learn(nsB, optB.qhsic$opt[[1]])
#mB.hsic    <- hsic$learn(nsB, optB.hsic$opt[[1]])



#################################################################################


par(mfrow=c(2,2))
plot(nsTestF$x,nsTestF$y)
lines(nsSeqF$x, krr$predict(mF.krr, nsSeqF), col="red", lwd=2)
lines(nsSeqF$x, qhsic$predict(mF.qhsic, nsSeqF), col="blue", lwd=2)
#lines(nsSeqF$x, hsic$predict(mF.hsic, nsSeqF), col="green", lwd=2)
plot(nsTestB$x,nsTestB$y)
lines(nsSeqB$x, krr$predict(mB.krr,nsSeqB), col="red", lwd=2)
lines(nsSeqB$x, qhsic$predict(mB.qhsic,nsSeqB), col="blue", lwd=2)
#lines(nsSeqB$x, hsic$predict(mB.hsic,nsSeqB), col="green", lwd=2)

resF.krr <- nsTestF$y-krr$predict(mF.krr,nsTestF)
resF.qhsic <- nsTestF$y-qhsic$predict(mF.qhsic,nsTestF)
#resF.hsic <- nsTestF$y-hsic$predict(mF.hsic,nsTestF)
plot(nsTestF$x, resF.krr, col="red", ylim=range(c(resF.krr, resF.qhsic))) #, resF.hsic
lines(nsTestF$x, resF.qhsic, col="blue", type="p")
#lines(nsTestF$x, resF.hsic, col="green", type="p")

resB.krr <- nsTestB$y-krr$predict(mB.krr,nsTestB)
resB.qhsic <- nsTestB$y-qhsic$predict(mB.qhsic,nsTestB)
#resB.hsic <- nsTestB$y-hsic$predict(mB.hsic,nsTestB)
plot(nsTestB$x, resB.krr, col="red", ylim=range(c(resB.krr, resB.qhsic))) # , resB.hsic
lines(nsTestB$x, resB.qhsic, col="blue",type="p")
#lines(nsTestB$x, resB.hsic, col="green",type="p")

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
