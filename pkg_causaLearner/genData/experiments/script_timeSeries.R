# simulate non-linear time series

remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner")
print("loading causal learners functions")
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner/genData/func_getData_v2.R", echo=FALSE)

#library(kernlab)
#library(pcalg)
#library(unifDAG)
#library(abind)
#library(gRbase) # topoSort
#library(reshape2)
#library(ggplot2)


# x1 in R n x p1
# x2 in R n x p2
# y1 in R m x p1
# y2 in R m x p2

# Now I want to do for time series dag:
# p- number of processes
# C- a vector of cycles C = {c1,c2,..,cn}
# L- a vector of lags for for each of above cycles L={l(c1),...,l(cn)}
# m =  p*(t(C) %*%  L) number of nodes in DAG
# T = time steps to evolve entire system
# burnin- number of inital time steps to throw away

# tasks

# A - sim DAG as a (p x p ) x C x L array
# 1. sim instantaneous DAG
# 2. sim rest of DAG (time dependencies)
#     for c in C
#       for l in 1:l(c)
#           sim DG (with high prob along diagonal)

# B - simulate SEM
# 1. simulate noise N in R(T x p)
# 2. initialize one unit of dag
# 3. simulate functions - 1 per process
#   a) work out number of parameters with DAG 
#   b) simulate function f(sigma, z, alpha)
# 4. evolve system using dag structure
# 5. convert to time series object and plot



# desirables
# linearity->non-linearity
# multidimensional
# causal sufficiency
# stochastic -> deterministic
# measurement error




# vary non-linearity
x <- seq(-1.2,1.2,length.out=100)
sigma0 <- 1/median(as.numeric(dist(x)^2))
rng <- quantile(1/as.numeric(dist(x)^2), probs=c(0.1,0.9))
sigmas <- exp(seq(log(rng[1]), log(rng[2]), length.out=9))
dist <- "runif"
distPars <- list(min=-0.1, max=0.1)


par(mfrow=c(3,3))
seed <- sample(1000,1)
aux <- lapply(1:9, function(i){
  # i <- 1
  print(paste("i: ",i))
  
  set.seed(seed)
  mysimfx <- simfxH(numVars=1,sigma=sigmas[i], dist=dist, distPars=distPars)
  #names(mysimfx)
  
  fx <- writefxH(mysimfx)
  distPars2 <- c(distPars, n=length(x))
  no <- do.call(dist,distPars2)
  
  #y <- applySimfxH(x=matrix(x,length(x),1), simfx=mysimfx)
  y <- do.call("fx", list(x=matrix(x,length(x),1), n=no))
  
  plot(x,y, main=i) # , ylim=c(-1,1)
  abline(h=c(-1,1), v=c(-1,1), col="red")
})



# vary additivity
x <- seq(-1.2,1.2,length.out=100)

sigma0 <- 1/median(as.numeric(dist(x)^2))
rng <- quantile(1/as.numeric(dist(x)^2), probs=c(0.1,0.9))
sigmas <- exp(seq(log(rng[1]), log(rng[2]), length.out=9))
dist <- "runif"
distPars <- list(min=-0.5, max=0.5)


par(mfrow=c(3,3))
seed <- sample(1000,1)
aux <- lapply(1:9, function(i){
  # i <- 1
  print(paste("i: ",i))
  
  set.seed(seed)
  mysimfx <- simfxH(numVars=1,sigma=sigma0, sigmaErr=sigmas[i], dist=dist, distPars=distPars)
  #names(mysimfx)
  
  fx <- writefxH(mysimfx)
  distPars2 <- c(distPars, n=length(x))
  no <- do.call(dist,distPars2)
  #y <- applySimfxH(x=matrix(x,length(x),1), simfx=mysimfx)
  y <- do.call("fx", list(x=matrix(x,length(x),1), n=no))
  
  
  plot(x,y, main=i) # , ylim=c(-1,1)
  abline(h=c(-1,1), v=c(-1,1), col="red")
})

# vary additivity 2
x <- seq(-1.2,1.2,length.out=100)

sigma0 <- 1/median(as.numeric(dist(x)^2))
rng <- quantile(1/as.numeric(dist(x)^2), probs=c(0.1,0.9))
sigmas <- exp(seq(log(rng[1]), log(rng[2]), length.out=9))
dist <- "runif"
distPars <- list(min=-0.3, max=0.3)


par(mfrow=c(3,3))
seed <- sample(1000,1)
aux <- lapply(1:9, function(i){
  # i <- 1
  print(paste("i: ",i))
  
  set.seed(seed)
  mysimfx <- simfxH(numVars=1,sigma=sigma0, dist=dist, distPars=distPars, 
                    geU=eval(parse(text=paste("function(y,n) sign(y+n)*(abs(y+n))^",1+(i-1)/5, sep=""))))
  #names(mysimfx)
  fx <- writefxH(simfx=mysimfx)
  distPars2 <- c(distPars, n=length(x))
  no <- do.call(dist,distPars2)
  #y <- applySimfxH(x=matrix(x,length(x),1), simfx=mysimfx)
  y <- do.call("fx", list(x=matrix(x,length(x),1), n=no))
  
  plot(x,y, main=i) # , ylim=c(-1,1)
  abline(h=c(-1,1), v=c(-1,1), col="red")
})


p <- 3
C <- c(1,3,12)
L <- c(2,2,2)
pDiag <- 0.6
pOffDiag <- 0.1
timeDag <- simTimeDag(p, C, L, pDiag, pOffDiag, instantaneous=TRUE)

summaryDag <- apply(timeDag, c("from","to"), function(edge) any(edge==1)*1)

# one variable
simTimeDag(p=1, C=1, L=1, pDiag=1, pOffDiag=1, instantaneous=FALSE)
simTimeDag(p=1, C=2, L=2, pDiag=1, pOffDiag=1, instantaneous=FALSE)
simTimeDag(p=1, C=c(1,2), L=c(1,2), pDiag=1, pOffDiag=1, instantaneous=FALSE)

# more variables
simTimeDag(p=2, C=1, L=1, pDiag=0.2, pOffDiag=0.1, instantaneous=TRUE)
simTimeDag(p=2, C=c(1,3,6), L=c(2,1,2), pDiag=0.2, pOffDiag=0.1, instantaneous=FALSE)


createTimeDag(C=1, L=1, matList=list(matrix(0,2,2),matrix(c(1,0,0,0),2,2)))
createTimeDag(C=1, L=1, matList=list(matrix(c(0,1,0,0),2,2),matrix(c(1,0,0,1),2,2)))
createTimeDag(C=c(1,6), L=c(1,1), matList=list(matrix(c(0,1,0,0),2,2),matrix(c(1,0,0,1),2,2),matrix(c(1,0,0,1),2,2)))



p <- 3
C <- c(1,3,12)
L <- c(2,2,2)
pDiag <- 0.6
pOffDiag <- 0.1
timeDag <- simTimeDag(p, C, L, pDiag, pOffDiag, instantaneous=TRUE)
pctHidVars <- 1
sigma <- 1
sigmaErr<- 0
distNs<- "runif"
distNsPars <- list(min=-0.5, max=0.5)
funcs <- simFunsSem(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars)


funcsErr <- simFunsErr(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars)



# plot timeDAG, summaryDAG, timeSeries, timeSeries with error,  cross-hsic and cross-correlation plots according 
# to max(cycleLag), scatter plots according to cycleLag

# parameters
burnin <- 100
n <- 2000
p <- 3
C <- c(1,3,12)
L <- c(2,2,2)
pDiag <- 0.6
pOffDiag <- 0.1
pctHidVars <- 1
sigma <- 1
sigmaErr <- 0
distNs <- "runif"
distNsPars <- list(min=-10^-20, max=10^-20)
distErr <- "runif"
distErrPars <- list(min=-0.5, max=0.5)

# simulate
set.seed(7)  # interesting seeds:7, 
timeDag <- simTimeDag(p, C, L, pDiag, pOffDiag, instantaneous=TRUE)
funcsSem <- simFunsSem(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars, geU=function(y,n) y+n)
funcsErr <- simFunsErr(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars)
smTmSem <- simTimeSEM(n, burnin, timeDag, funcsSem, funcsErr, distNs, distNsPars, distErr, distErrPars)

names(smTmSem)


# timeDAG
timeDag
cycleLags <- dimnames(timeDag)$cycleLag
for(cycLg in rev(cycleLags)) plot(getGraph(timeDag[,,cycLg]), main=cycLg)

# summaryDAG
summaryDag <- apply(timeDag, c("from","to"), function(edge) any(edge==1)*1)
plot(getGraph(summaryDag))


# timeSeries, timeSeries with error
plot(as.ts(smTmSem$timeSeries), main="")
plot(as.ts(smTmSem$timeSeriesErr), main="")

# for each node plot node vs. dependencies

timeSeries <- smTmSem$timeSeries

head(timeSeries)

  

df <- getTimeDF(timeSeries=smTmSem$timeSeries, timeDag, funcsSem)
p <- ggplot(df)
p <- p + geom_point(aes(x=value, y=depVal, colour=nodeIndep), alpha=0.1)
p <- p + facet_grid(nodeDep~variable, scales="free")
p

pacf(as.ts(smTmSem$timeSeries))
acf(as.ts(smTmSem$timeSeries))

#####################################################################################################################
# a simple p=2, lag=1, additive noise example

# parameters
burnin <- 100
n <- 2000
p <- 2
C <- c(1)
L <- c(1)
sigma <- 1
sigmaErr <- 0
distNs <- "runif"
distNsPars <- list(min=-0.1, max=0.1)
distErr <- "runif"
distErrPars <- list(min=-0.1, max=0.1)

# simulate
(seed <- sample(1:1000, 1))
set.seed(seed) # interesting seeds: 
timeDag <- createTimeDag(C=1, L=1, matList=list(matrix(0,2,2),matrix(c(1,1,0,1),2,2)))
funcsSem <- simFunsSem(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars, geU=function(y,n) y+n)
funcsErr <- simFunsErr(timeDag, pctHidVars, sigma, sigmaErr, distNs, distNsPars)
smTmSem <- simTimeSEM(n, burnin, timeDag, funcsSem, funcsErr, distNs, distNsPars, distErr, distErrPars)
names(smTmSem)


# timeDAG
timeDag
cycleLags <- dimnames(timeDag)$cycleLag
for(cycLg in rev(cycleLags)) plot(getGraph(timeDag[,,cycLg]), main=cycLg)
# summaryDAG
summaryDag <- apply(timeDag, c("from","to"), function(edge) any(edge==1)*1)
plot(getGraph(summaryDag))


# timeSeries, timeSeries with error
summary(smTmSem$timeSeries)
plot(as.ts(smTmSem$timeSeries), main="")
plot(as.ts(smTmSem$timeSeriesErr), main="")

# for each node plot node vs. dependencies

timeSeries <- smTmSem$timeSeries
df <- getTimeDF(timeSeries=smTmSem$timeSeries, timeDag, funcsSem)


p <- ggplot(df)
p <- p + geom_point(aes(x=value, y=depVal, colour=nodeIndep), alpha=0.1)
p <- p + facet_grid(nodeDep~variable, scales="free")
p

pacf(as.ts(smTmSem$timeSeries))
acf(as.ts(smTmSem$timeSeries))
