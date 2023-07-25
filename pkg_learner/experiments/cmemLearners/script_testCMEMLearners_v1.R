# Kernel Deviance first approach

remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner")
print("loading causal learners functions")
source("./func_causaLearners_v1.R", echo=FALSE)

p <- 2
nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
nodes <- rep(list(nodesX),p)
dag <- matrix(c(0,0,1,0),2,2)
colnames(dag) <- rownames(dag) <- c("x","y")

set.seed(4)
simTest <- simRandAddSEM(p,1000, nodes, dag=dag)
plot(getGraph(simTest$dag))

X <- simTest$x
X <- apply(X, 2, stdrize)

apply(X, 2, mean)
apply(X, 2, sd)

x <- X[,"x"]
y <- X[,"y"]


plot(x,y)
plot(y,x)

KCDC <- function(x,y, lambda){
  
  sigma0 <- 1/median(as.numeric(dist(x)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))

 

  n <- length(x)
  L  <- kernelMatrix(kernelX, x) 
  K  <- kernelMatrix(kernelY, y)
  I <- diag(n)

  Blambda <- solve(L+n*lambda*I)
  Alambda <- Blambda%*%K%*%Blambda
  
  LAL <- L%*%Alambda%*%L

  #b <- sum(diag(LAL))
  #c <- sum(diag(LAL^2))
  
  b <- sum(diag(LAL)^(0.5))
  c <- sum(diag(LAL))
  
  res <- (c/n) - (b/n)^2
  
  return(res)
}

KCDCrel <- function(x, y, lambda, numPerms){
  
  sigma0 <- 1/median(as.numeric(dist(x)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(x)
  set.seed(12345)            
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  
  mesr <- KCDC(x, y, lambda)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y,y, lambda)))
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    
    res <- KCDC(x, y[rperm], lambda)#/(sqrt(KCDC(x,x, lambda))*sqrt(KCDC(y[rperm],y[rperm], lambda))) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}

KCSC <- function(x, y, lambda){
  
  thetaX <- 1/median(as.numeric(dist(x)^2))
  kernelX <- do.call("rbfdot", list(sigma=thetaX))
  thetaY <- 1/median(as.numeric(dist(y)^2))
  kernelY <- do.call("rbfdot", list(sigma=thetaY))
  n <- length(x)
  L  <- kernelMatrix(kernelX, x) 
  K  <- kernelMatrix(kernelY, y)
  ones <- rep(1,n)
  I <- diag(n)
  B <- solve(L+n*lambda*I)
  
  sigma2X <- 1/(4*thetaX)
  sigma2Y <- 1/(4*thetaY)
  muX <- (x%*%t(ones)+ones%*%t(x))/2
  muY <- (y%*%t(ones)+ones%*%t(y))/2
  
  aux1 <- (x%*%t(ones))^2+(ones%*%t(x))^2
  aux1 <- exp(-aux1*thetaX/2)
  aux1 <- thetaX*(sqrt(pi/(2*thetaX)))*aux1
  aux2 <- pnorm(max(x), mean=muX, sd=sqrt(sigma2X))-pnorm(min(x), mean=muX, sd=sqrt(sigma2X))
  Chat <- aux1*aux2
  
  aux1 <- (y%*%t(ones)-ones%*%t(y))^2
  aux1 <- exp(-aux1*thetaY/2)
  aux1 <- sqrt(pi/(2*thetaY))*aux1
  aux2 <- pnorm(max(y), mean=muY, sd=sqrt(sigma2Y))-pnorm(min(y), mean=muY, sd=sqrt(sigma2Y))
  E <- aux1*aux2
  
  res <- as.numeric(t(ones)%*%(B%*%Chat%*%B * E)%*%ones)
  return(res)
}

lambda <- 10^-5
numPerms <- 100

KCDC(x,y, lambda)
KCDC(y,x, lambda)

KCDCrel(x,y,lambda, numPerms)
KCDCrel(y,x,lambda, numPerms)


# set lambda parameter according to "Conditional mean embeddings as regressors" Grünewälder et al

vecRegLoss <- function(xTr,yTr, xTe=NULL, yTe=NULL, lambda){
  
  sigma0 <- 1/median(as.numeric(dist(xTr)^2))
  kernelX <- do.call("rbfdot", list(sigma=sigma0))
  sigma0 <- 1/median(as.numeric(dist(yTr)^2))
  kernelY <- do.call("rbfdot", list(sigma=sigma0))
  
  n <- length(xTr)
  
  LTr  <- kernelMatrix(kernelX, xTr) 
  KTr  <- kernelMatrix(kernelY, yTr)
  I <- diag(n)
  
  if(!is.null(xTe)){
    nTe <- length(xTe)
    LTrTe  <- kernelMatrix(kernelX, xTr, xTe) 
    KTrTe  <- kernelMatrix(kernelY, yTr, yTe)
    KTe <- kernelMatrix(kernelY, yTe)
  } else{
    nTe <- n
    LTrTe <- L
    KTrTe <- K
    KTe <- K
  }
  
  
  
  Blambda <- solve(LTr+n*lambda*I)
  Alambda <- Blambda%*%KTr%*%Blambda
  
  LAL <- t(LTrTe)%*%Alambda%*%LTrTe
  LBK <- t(LTrTe)%*%Blambda%*%KTrTe
  
  
  res <- sum(diag(LAL)) + sum(diag(KTe)) - 2*sum(diag(LBK))
  
  return(res)
}

vecRegLoss(x,y,lambda=lambda)


CV.parallel <- function(x, y, numFolds, lambdas, fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  trainData <- cbind(x=x,y=y)

  nParams <- length(lambdas)
  dimnames <- list(as.character(1:numFolds), lambdas)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    res <- mcmapply(FUN=function(f){
      # f <- 1; p <- lambdas[1]
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      lossesTrain <- vecRegLoss(xTr=curTrain[,"x"], yTr=curTrain[,"y"],  lambda=p)
      lossesTest <- vecRegLoss(xTr=curTrain[,"x"], yTr=curTrain[,"y"], xTe=curTest[,"x"], yTe=curTest[,"y"], lambda=p)
      res <- c(lossesTrain, lossesTest)
      
      return(res)
    }, f=1:numFolds, mc.cores=5, SIMPLIFY="array")
    return(res)
  }, p=lambdas, mc.cores=2, SIMPLIFY="array")
  

  dimnames(losses) <- list(trainTest=c("train","test"), fold=1:numFolds, params=lambdas)
  
  #print("exits CV parallel")
  return(losses)
}

lambdas <- 10^seq(-8,0, length.out=50)
numFolds <- 5

lossesXY <- CV.parallel(x, y, numFolds, lambdas)
lossesMeanXY <- apply(lossesXY, c("trainTest", "params"), mean)
lossesYX <- CV.parallel(y, x, numFolds, lambdas)
lossesMeanYX <- apply(lossesYX, c("trainTest", "params"), mean)

# X->Y
plot(log(lambdas,10), lossesMeanXY["test",], ylim=range(c(lossesMeanXY)), type="b", pch=1)
lines(log(lambdas,10), lossesMeanXY["train",], col="red", type="b", pch=1)
plot(rep(log(lambdas,10), rep(numFolds, length(lambdas))), lossesXY["test",,], ylim=c(range(lossesXY["test",,])), type="p", pch=1)
lines(log(lambdas,10), lossesMeanXY["test",], col="red", type="b", pch=1)
lossesMeanXY["test",]
which.min(lossesMeanXY["test",])

# Y->X
plot(log(lambdas,10), lossesMeanYX["test",], ylim=range(c(lossesMeanYX)), type="b", pch=1)
lines(log(lambdas,10), lossesMeanYX["train",], col="red", type="b", pch=1)
plot(rep(log(lambdas,10), rep(numFolds, length(lambdas))), lossesYX["test",,], ylim=c(range(lossesYX["test",,])), type="p", pch=1)
lines(log(lambdas,10), lossesMeanYX["test",], col="red", type="b", pch=1)
lossesMeanYX["test",]
which.min(lossesMeanYX["test",])


# KCDC dependence on lambda
lambdas <- 10^seq(-8,8,1)
resX <- sapply(lambdas, function(lambda) KCDC(x, y, lambda))
resY <- sapply(lambdas, function(lambda) KCDC(y, x, lambda))
plot(log(lambdas,10), resX, ylim=range(resX,resY))
lines(log(lambdas,10),resY, col="red",type="p")

#KCDCrel dependence on lambda
lambdas <- 10^seq(-8,8,1)
resX <- sapply(lambdas, function(lambda) KCDCrel(x, y, lambda, numPerms))
resY <- sapply(lambdas, function(lambda) KCDCrel(y, x, lambda, numPerms))
plot(log(lambdas,10), resX, ylim=range(resX,resY))
lines(log(lambdas,10),resY, col="red",type="p")

# Lets see for 100 datasets what KCDC and KCSC is 

q <- 100
set.seed(4)
(p <- rep(2, q))
(n <- rep(1000, q))
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- lapply(p, function(p) rep(list(nodes),p))


# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
set.seed(5)
dataTestList <- simRandAddSEMs(q, p, n, nodes, dag=dag)


unique(dataTestList$dags, MARGIN=3)

i <- 50
plot(dataTestList$xs[[i]])

lambda0 <- 10^-5
lambdas <- 10^seq(-8,0, length.out=50)
numFolds <- 5

pm <- proc.time()
msrs <- mcmapply(function(el, nm){
  # i <- 1; el <- dataTestList$xs[[i]]
  X <- apply(el, 2, stdrize)
  
  print(paste("name: ", nm))
  #print("head(X)")
  #print(head(X))
  #print(apply(X, 2, mean))
  #print(apply(X, 2, sd))
  
  x <- X[,1]
  y <- X[,2]
  
  lossesXY <- CV.parallel(x, y, numFolds, lambdas)
  lossesMeanXY <- apply(lossesXY, c("trainTest", "params"), mean)
  lossesYX <- CV.parallel(y, x, numFolds, lambdas)
  lossesMeanYX <- apply(lossesYX, c("trainTest", "params"), mean)
  
  lambdaXYOpt <- lambdas[which.min(lossesMeanXY["test",])]
  lambdaYXOpt <- lambdas[which.min(lossesMeanYX["test",])]
  
  
  
  #KCDC lambda fix
  KCDCxy <- KCDC(x,y, lambda0)
  KCDCyx <- KCDC(y,x, lambda0)
  KCDCq <- KCDCyx/KCDCxy
  
  #KCSC lambda fix
  KCSCxy <- KCSC(x,y, lambda0)
  KCSCyx <- KCSC(y,x, lambda0)
  KCSCq <- KCSCyx/KCSCxy
  
  #KCDC lambda opt
  KCDCxy <- KCDC(x,y, lambdaXYOpt)
  KCDCyx <- KCDC(y,x, lambdaYXOpt)
  KCDCq_opt <- KCDCyx/KCDCxy
  
  #KCSC lambda opt
  KCSCxy <- KCSC(x,y, lambdaXYOpt)
  KCSCyx <- KCSC(y,x, lambdaYXOpt)
  KCSCq_opt <- KCSCyx/KCSCxy
  
  #KCDC lambda optMin
  lambdaOpt <- min(lambdaXYOpt, lambdaYXOpt)
  KCDCxy <- KCDC(x,y, lambdaOpt)
  KCDCyx <- KCDC(y,x, lambdaOpt)
  KCDCq_optMin <- KCDCyx/KCDCxy
  
  #KCSC lambda optMin
  
  KCSCxy <- KCSC(x,y, lambdaOpt)
  KCSCyx <- KCSC(y,x, lambdaOpt)
  KCSCq_optMin <- KCSCyx/KCSCxy
  
  
  res <- c(KCDCq, KCSCq, KCDCq_opt, KCSCq_opt, KCDCq_optMin, KCSCq_optMin, lambdaXYOpt, lambdaYXOpt)
  names(res) <- c("KCDC","KCSC","KCDCopt","KCSCopt","KCDCoptMin","KCSCoptMin","lambdaXY","lambdaYX")
  return(res)
}, el=dataTestList$xs, nm=dataTestList$names, 
SIMPLIFY="array", mc.cores=6)
proc.time() - pm #

names(dimnames(msrs)) <- c("measure", "database")

msrsDB <- melt(msrs)

summary(msrsDB)

p <- ggplot(msrsDB)
p <- p + geom_boxplot(aes(x=measure, y=value))
p <- p + facet_wrap(.~measure, scales="free")
p

summary(t(msrs))
pctRight <- function(x){
  
  res <- sum(x>1)/length(x)*100
  return(res)
}
  
# KCDC
pctRight(x=msrs["KCDC",])
pctRight(msrs["KCDCopt",])
pctRight(msrs["KCDCoptMin",])

# KCSC
pctRight(x=msrs["KCSC",])
pctRight(msrs["KCSCopt",])
pctRight(msrs["KCSCoptMin",])

# lambda
plot(log(msrs["lambdaXY",], 10), log(msrs["lambdaYX",], 10))
abline(a=0, b=1, col="red")
sum(msrs["lambdaXY",]<msrs["lambdaYX",])

cast(msrsDB, measure~., value="value", fun.aggregate="pctRight")



