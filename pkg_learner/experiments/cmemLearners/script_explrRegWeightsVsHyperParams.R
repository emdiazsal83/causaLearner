# explore CONDITIONAL MEAN EMBEDDING - SPECIFICALLY HOW THE REGRESSION WEIGHTS DEPEND ON HYPERPARAMETERS

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"


repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)


source("./pkg_causaLearner/genData/func_getData_v2.R")
#library(dHSIC)
source("./pkg_dHSIC/dHSIC.R")
library(kpcalg)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_learner/func_learners_v3.R")
source("./pkg_learner/ob_cme_learner_v2.R", echo=FALSE)
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")
library(reshape)
library(ggplot2)

#########################################################################################################################################*
# Norm of mean embedding of a distribution
#########################################################################################################################################*

sigma <- 0.01
n <- 1000
x <- matrix(rnorm(n, mean=0, sd=1), n, 1)
K <- kern_rbf(x, sigma=sigma)
mu <- (1/n^2)*sum(K)


num <- 5
pis <- seq(0,1, length.out=num+1); pis <- pis[2:(length(pis)-1)] #sort(runif(num-1))
belong <- runif(n)
belong <- findInterval(belong, pis)+1
table(belong)
mus <- rnorm(num, mean=0, sd=2)
sds <- rexp(num, rate=2)
x <- rnorm(n, mean=mus[belong], sd=sds[belong])
dens <- density(x)
hist(x, prob=T, 100)
lines(dens, col="red")

x <- stdrize(x)
dens <- density(x)
hist(x, prob=T, 100)
lines(dens, col="red")
x <- matrix(x, n, 1)
K <- kern_rbf(x, sigma=sigma)
mu <- (1/n^2)*sum(K)


random_fourier_features <- function(x, w=NULL, b=NULL, num_f=NULL, sigma=NULL, seed=NULL){
  
  if (is.null(num_f)) num_f = 25
  
  if(is.null(dim(x))){
    x <-  as.matrix(x)  
  }
  
  
  r <- nrow(x)
  c <- ncol(x)
  
  if (is.null(sigma)){sigma=1}
  
  
  if (is.null(w)){
    
    
    set.seed(seed)
    w=(1/sigma)*matrix(rnorm(num_f*c), num_f, c)
    set.seed(seed)
    b=matrix(2*pi*runif(num_f),num_f,r)
  } #else if (nrow(w)<num_f){
  #w=rbind(w,(1/sigma)*matrix(rnorm((num_f-nrow(w))*c),num_f-nrow(w),c));
  #b=rbind(b,repmat(2*pi*runif(num_f-nrow(b)),1,r));
  #}
  
  feat = sqrt(2)*t(cos(w[1:num_f,1:c, drop=FALSE]%*%t(x)+b[1:num_f,]))
  #feat = sqrt(2)*t(log(1+exp(w[1:num_f,1:c]%*%t(x)+b[1:num_f,])));
  
  out=list(feat=feat,w=w,b=b);
  return(feat)
}

x <- rnorm(1000)
ft <- random_fourier_features(x, num_f=2)

pairs(ft)

# now lets do it for different number of normals in mixture

sigmaX <- 1/median(as.numeric(dist(unique(x))^2))


# mean embedding without kernel

library(moments)

nums <- round(seq(1,30, length.out=20))
sigma <- 10
numReps <- 1000
xReps <- sapply(nums, function(num){
  # num <- 1
  print(paste("num: ", num))
  xRep <- sapply(1:numReps, function(rep){
    # rep <- 5
    pis <- seq(0,1, length.out=num+1)
    pis <- sort(pis[2:(length(pis)-1)]) #sort(runif(num-1))
    
    if(num==1){
      belong <- rep(1, n)
    } else{
      belong <- runif(n)  
      belong <- findInterval(belong, pis)+1  
    }
    
    # table(belong)
    mus <- rnorm(num, mean=rnorm(num, sd=rexp(num, rate=2)), sd=10)
    sds <- rexp(num, rate=10)
    x <- rnorm(n, mean=mus[belong], sd=sds[belong])
    x <- norml(x) 
    #x <- stdrize(x)
    #print(1/median(as.numeric(dist(unique(x))^2)))
    #dens <- density(x)
    #hist(x, prob=T, 100)
    #lines(dens, col="red")
    #x <- matrix(x, n, 1)
    
    #mu <- apply(cbind(x^3, x^4), 2, mean)
    #mu <- apply(random_fourier_features(x, num_f=2), 2, mean)
    #mu <- apply(cbind(sin(x), cos(x)),2 , mean)
    #mu <- apply(cbind(kurtosis(x), skewness(x)),2 , mean)
    return(x)
  }, simplify="array")
  return(xRep)
}, simplify="array")

length(xReps)
summary(as.numeric(xReps))

dim(xReps)
dimnames(xReps) <- list(obs=1:n, rep=1:numReps, complexity=nums)

xDB <- melt(xReps)

indx <- which(xDB$rep %in% 996:1000 & xDB$complexity %in% c(1,3,10,30))
p <- ggplot(xDB[indx,]) 
p <- p + geom_histogram(aes(x=value, y=..density..), fill="white", alpha=0.1 ,colour="grey", fill="grey", bins=10)
p <- p + geom_density(aes(x=value, y=..density..), alpha=.2, fill="#FF6666")
p <- p + facet_grid(rep~complexity)
p <- p + theme(strip.text = element_text(face="bold", size=20,lineheight=5.0))
p

dim(xReps)
mus <- apply(xReps, c("rep","complexity"), function(vec) c(kurtosis(vec), skewness(vec)))
dim(mus)
dimnames(mus) <- list(feature=c("kurtosis","skewness"), rep=1:numReps, complexity=nums)

musDB <- melt(mus)
musDB <- cast(musDB, rep+complexity~feature, value="value")
colnames(musDB) <- c("rep", "cmplx", "f1", "f2")

indx <- which(musDB$cmplx %in% c(1,3,10,30))
p <- ggplot(musDB[indx,])
p <- p + geom_point(aes(x=f1, y=f2, colour=factor(cmplx)), size=0.8, shape=1)
p <- p + facet_grid(~cmplx)
p

musDB$norm <- with(musDB, sqrt(f1^2 + f2^2))


musDB2 <- cast(musDB, cmplx~., value="norm", fun.aggregate=mean)
colnames(musDB2) <- c("complexity", "norm")
p <- ggplot(musDB2)
p <- p + geom_point(aes(x=complexity, y=norm))
p

musDB2 <- cast(musDB, cmplx~., value="norm", fun.aggregate=sd)
colnames(musDB2) <- c("complexity", "norm")
p <- ggplot(musDB2)
p <- p + geom_point(aes(x=complexity, y=norm))
p

# now do it for x*n example

# in causal direction





# mean embedding complexity with kernel
nums <- round(seq(3,30, length.out=20))
sigma <- 10
numReps <- 1000
mus <- sapply(nums, function(num){
  # num <- 5
  print(paste("num: ", num))
  musRep <- sapply(1:numReps, function(rep){
    pis <- seq(0,1, length.out=num+1); pis <- pis[2:(length(pis)-1)] #sort(runif(num-1))
    belong <- runif(n)
    belong <- findInterval(belong, pis)+1
    mus <- rnorm(num, mean=0, sd=2)
    sds <- rexp(num, rate=2)
    x <- rnorm(n, mean=mus[belong], sd=sds[belong])
    x <- norml(x)
    #print(1/median(as.numeric(dist(unique(x))^2)))
    #dens <- density(x)
    #hist(x, prob=T, 100)
    #lines(dens, col="red")
    x <- matrix(x, n, 1)
    K <- kern_rbf(x, sigma=sigma)
    mu <- (1/n^2)*sum(K)
    return(mu)
  })
  mu <- mean(musRep)
  return(mu)
})

plot(nums, mus)

#########################################################################################################################################*
# Make graphs
#########################################################################################################################################*
# XtoY
edgL <- vector("list", length=3)
names(edgL) <- c("x","y","t")
edgL[["x"]] <- list(edges=c("y"))
XtoY <- graphNEL(nodes=c("x","y","t"), edgeL=edgL, edgemode="directed")
plot(XtoY)


#########################################################################################################################################*
# Define SEMS, simulate and plot(x,y)
#########################################################################################################################################*
N <- 250
#############################################################*
# XtoY
funcsXtoY <- list(fx=function(n) n, fy=function(x, n) sin(x)*n,ft=function(n) n)
semXtoY <- list(dag=XtoY, funcs=funcsXtoY, simPars=list(n=N, nodes=list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                                                                        y=list(dist="runif", pars=list(min=-1, max=1), a=1, b=1), 
                                                                        t=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1))))
set.seed(3)
simXtoY <- simSEM(semXtoY)

x <- simXtoY$x[,"x", drop=F]
y <- simXtoY$x[,"y"]

plot(x,y)

ord <- order(x[,1])
x <- x[ord,,drop=F]
y <- y[ord]

xy <- apply(cbind(x,y),2,norml)

trainDataXY <- constructData(xy[,1, drop=FALSE], xy[,2])

x <- trainDataXY$x[,1]
y <- trainDataXY$y
ord <- order(y)
x <- x[ord]
y <- y[ord]

trainDataYX <- constructData(as.matrix(y), x)

pm <- proc.time()
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_log_quad_L2_lambda_kernParsX, trainData=trainDataXY, plot=TRUE)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_log_quad_L2_lambda_kernParsX, trainData=trainDataYX, plot=TRUE)
proc.time() - pm #  41 secsfor 250 points

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$sigmaX$seq
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$sigmaX$seq

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$offsetX$seq
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$offsetX$seq

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$degreeX$seq
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$degreeX$seq

getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "lambda")

getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX")

getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "offsetX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "offsetX")

getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "degreeX")
getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_yx, "degreeX")



#cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$sigmaX$val <- getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX")*100

cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)

cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)

alphas <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learnParams$Blambda %*% cmem_rbf_rbf_L2_lambda_kernParsX_xy$learnParams$Lx

# write alpha(x) and mu(y|x)

alpha <- function(x, learner){
  data <- constructData(matrix(x,1,1), 1)
  
  parsXs <- getKernelPars(learner, kernelName="kernelX")
  trainData <- learner$hyperParams$trainData
  lx <- kernelMatrix(learner$hyperParams$non_data$kernelX$name, x=data$x, y=trainData$x, pars=parsXs)
  
  #pred <- learner$predict(learner, data)
  res <- learner$learnParams$Blambda %*% t(lx)
  return(as.numeric(res))
}

mu <- function(y, x, learner){
  data <- constructData(matrix(x,1,1), y)
  pred <- learner$predict(learner, data)
  alphas <- learner$learnParams$Blambda %*% t(pred$lx)
  res <- as.numeric(t(alphas) %*% t(pred$ky))
  return(res)
}

# plot alpha(x) and mu(y|x)

x <- trainDataXY$x
y <- trainDataXY$y



barplotAlpha <- function(x, xs, learner){
  
  wdth <- min(diff(xs))/2
  cumSumSp <- xs/wdth+seq(0,length(xs)-1)*wdth+wdth*0.5
  sp <- cumSumSp[2:length(cumSumSp)]-cumSumSp[1:(length(cumSumSp)-1)]
  sp <- c(min(sp), sp)
  sp <- (sp-min(sp))/(max(sp)-min(sp))
  sp <- sp*100+2*wdth
  barplot(alpha(x, learner), width=wdth, space=sp)
  xsT <- cumsum(sp)*wdth+seq(0,length(xs)-1)*wdth+wdth*0.5
  #abline(v=xsT, col="red")
  abline(h=0, col="red")
  abline(v=(x-min(xs))/(max(xs)-min(xs))*(max(xsT)-min(xsT))+min(xsT), col="green")
  
}

xs <- trainDataXY$x[,1]

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$sigmaX <- 1e4
cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$degreeX <- 0.98
cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$offsetX <- 1e-3

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$lambda$val <- 1e-8

cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)

cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)

barplotAlpha(0.5, xs, cmem_rbf_rbf_L2_lambda_kernParsX_xy)

xs <- trainDataXY$x[,1] #seq(min(x), max(x), length.out=10)
for(i in 1:length(xs)) barplotAlpha(xs[i], xs, cmem_rbf_rbf_L2_lambda_kernParsX_xy)

xs <- trainDataYX$x[,1] #seq(min(x), max(x), length.out=10)
for(i in 1:length(xs)) barplotAlpha(xs[i], xs, cmem_rbf_rbf_L2_lambda_kernParsX_yx)

# plot mu
ys <- sort(trainDataXY$y) #seq(min(y), max(y), length.out=100)
xs <- sort(trainDataXY$x[,1]) #seq(min(x), max(x), length.out=101)
grid <- sapply(xs, function(x) sapply(ys, function(y) mu(y,x, cmem_rbf_rbf_L2_lambda_kernParsX_xy), simplify="array") ,simplify="array")
dimnames(grid) <- list(y=ys, x=xs)

dimnames(grid)
plot(ys, grid[,1], type="l", xlim=range(ys), ylim=range(grid))
for(i in 1: ncol(grid)) lines(ys, grid[,i], type="l")

gridDB <- melt(grid)

v <- ggplot(gridDB, aes(x, y, z = value))
v <- v + geom_raster(aes(fill = value)) 
v <- v + geom_contour(colour = "white", bins = 10)
v

# use alpha(x) for kernel density estimation (of pdf) -> see notes on computational statistics
kernel <- c("gaussian", "epanechnikov", "rectangular",
            "triangular", "biweight",
            "cosine", "optcosine")[2]
# XY
ys <- trainDataXY$y 
xs <- trainDataXY$x[,1]
densXY <- sapply(1:length(xs), function(i){
  alfa <- alpha(xs[i], cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  #barplotAlpha(x, xs, cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  alfa <- abs(alfa) 
  alfa <- alfa/sum(alfa)
  #barplot(alfa)
  density(ys, weights=alfa, kernel=kernel)$y
}, simplify="array")


# XY
ys <- trainDataYX$y 
xs <- trainDataYX$x[,1]
densYX <- sapply(1:length(xs), function(i){
  alfa <- alpha(xs[i], cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #barplotAlpha(x, xs, cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  alfa <- abs(alfa) 
  alfa <- alfa/sum(alfa)
  #barplot(alfa)
  density(ys, weights=alfa, kernel=kernel)$y
}, simplify="array")



summary(apply(densXY*densXY, 2, sum))
summary(apply(densYX*densYX, 2, sum))

sd(apply(densXY*densXY, 2, sum))
sd(apply(densYX*densYX, 2, sum))

sum(apply(densXY,1, sd))
sum(apply(densYX,1, sd))


par(mfrow=c(2,2))
# XY
plot(trainDataXY$x[,1], trainDataXY$y)
dens <- density(trainDataXY$y) 
plot(dens$x, dens$y, main="", ylim=range(densXY), type="l")
numCurves <- 100
cols <- heat.colors(length(xs))
for(i in seq(1, length(xs), length.out=numCurves)) lines(dens$x, densXY[,i], col=cols[i])
lines(dens$x, dens$y)
# YX
plot(trainDataYX$x[,1], trainDataYX$y)
dens <- density(trainDataYX$y) 
plot(dens$x, dens$y, main="", ylim=range(densYX), type="l")
numCurves <- 100
cols <- heat.colors(length(xs))
for(i in seq(1, length(xs), length.out=numCurves)) lines(dens$x, densYX[,i], col=cols[i])
lines(dens$x, dens$y)

#############################################################################################################################*

# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
# random functions
q <- 100
n <- 200
set.seed(4)
p <- 2
(ps <- rep(p, q))
(ns <- rep(n, q))
nodes <- list(dist="runif", pars=list(min=-2, max=2), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- c("x","y")
nodess <- lapply(p, function(p) nodes)
dag <- matrix(c(0,0,1,0),2,2)
colnames(dag) <- rownames(dag) <- c("x","y")
set.seed(5)
dataTestList <- simRandSEMs(q, ps, ns, nodess, sigma=5, sigmaErr=0, dagMat=dag)
plotPairsList(dataTestList)

# chosen function
edgL <- vector("list", length=2)
names(edgL) <- c("x","y")
edgL[["x"]] <- list(edges=c("y"))
XtoY <- graphNEL(nodes=c("x","y"), edgeL=edgL, edgemode="directed")
plot(XtoY)
nodesUnif <-  list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                   y=list(dist="runif", pars=list(min=-1, max=1), a=1, b=1))
funcs <- list(fx=function(n) n, fy=function(x, n) sin(x)*n)
sem <- list(dag=XtoY, funcs=funcs, simPars=list(n=200, nodes=nodesUnif))
dataTestList <- simSEMs(q, sem)
plotPairsList(dataTestList)


pm <- proc.time()
msrs <- mcmapply(function(el, nm){
  # i <- 9; el <- dataTestList$xs[[i]]; nm <- dataTestList$names[i]
  X <- apply(el, 2, stdrize)
  
  print(paste("name: ", nm))
  #print("head(X)")
  #print(head(X))
  #print(apply(X, 2, mean))
  #print(apply(X, 2, sd))
  
  x <- X[,1]
  y <- X[,2]
  
  trainDataXY <- constructData(as.matrix(x), y)
  trainDataYX <- constructData(as.matrix(y), x)
  
  # train hyperparameters
  cmem_rbf_rbf_L2_none_xy <- setParams(learner=cmem_log_quad_L2_none, trainData=trainDataXY, plot=FALSE)
  cmem_rbf_rbf_L2_none_yx <- setParams(learner=cmem_log_quad_L2_none, trainDataYX, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$degreeX <- 0.98
  cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$optimizable$lambda$val <- 1e-8
  cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$degreeX <- 0.98
  cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$optimizable$lambda$val <- 1e-8
  
  #cmem_rbf_rbf_L2_lambda_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataYX, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataYX, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataXY, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataYX, plot=FALSE)
  
  # train learn parameters
  cmem_rbf_rbf_L2_none_xy <- cmem_rbf_rbf_L2_none$learn(cmem_rbf_rbf_L2_none_xy)
  cmem_rbf_rbf_L2_none_yx <- cmem_rbf_rbf_L2_none$learn(cmem_rbf_rbf_L2_none_yx)
  #cmem_rbf_rbf_L2_lambda_xy <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_xy)
  #cmem_rbf_rbf_L2_lambda_yx <- cmem_rbf_rbf_L2_lambda$learn(cmem_rbf_rbf_L2_lambda_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- cmem_rbf_rbf_L2_lambda_kernParsXY$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  # calculate measures
  msrsFixXY <- cmem_rbf_rbf_L2_none_xy$calcMsrs(cmem_rbf_rbf_L2_none_xy)
  msrsFixYX <- cmem_rbf_rbf_L2_none_yx$calcMsrs(cmem_rbf_rbf_L2_none_yx)
  #msrsOpt1XY <- cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
  #msrsOpt1YX <- cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
  #msrsOpt2XY <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  #msrsOpt2YX <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #msrsOpt3XY <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #msrsOpt3YX <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  
  #KCDC lambda fix
  KCDCq <- msrsFixYX["KCDC"]/msrsFixXY["KCDC"] 
  #KCDCrel lambda fix
  KCDCrelq <- msrsFixYX["KCDCrel"]/msrsFixXY["KCDCrel"]
  
  #KCDC lambda opt1
  #KCDCq_opt1 <- msrsOpt1YX["KCDC"]/msrsOpt1XY["KCDC"]
  #KCDCrel lambda opt1
  #KCDCrelq_opt1 <- msrsOpt1YX["KCDCrel"]/msrsOpt1XY["KCDCrel"]
  
  #KCDC lambda opt2
  #KCDCq_opt2 <- msrsOpt2YX["KCDC"]/msrsOpt2XY["KCDC"]
  #KCDCrel lambda opt2
  #KCDCrelq_opt2 <- msrsOpt2YX["KCDCrel"]/msrsOpt2XY["KCDCrel"]
  
  #KCDC lambda opt3
  #KCDCq_opt3 <- msrsOpt3YX["KCDC"]/msrsOpt3XY["KCDC"]
  #KCDCrel lambda opt3
  #KCDCrelq_opt3 <- msrsOpt3YX["KCDCrel"]/msrsOpt3XY["KCDCrel"]
  
  res <- c(KCDCq, KCDCrelq)
  names(res) <- c("KCDC","KCDCrel")
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2, KCDCq_opt3, KCDCrelq_opt3)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2","KCDCopt3","KCDCrelopt3")
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2")
  return(res)
}, el=dataTestList$xs, nm=dataTestList$names, 
SIMPLIFY="array", mc.cores=6)
proc.time() - pm 

names(dimnames(msrs)) <- c("measure", "database")
msrsDB <- melt(msrs)
pctRight <- function(x){
  res <- sum(x>1)/length(x)*100
  return(res)
}
cast(msrsDB, measure~., value="value", fun.aggregate="pctRight")

