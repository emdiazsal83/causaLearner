remove(list=ls())

setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner/")
source("./hypothesisScorers/pkg_ANM_hypothesisScorer/func_learners_v1.R")
source("./dataTreatments/func_dataTreatments.R")

set.seed(3)
n <- 100
x <- runif(n, 0, 2*pi)
ny <- rnorm(n, 0, 0.1)
y <- sin(x) + ny
y <- y - mean(y)

plot(x, y)

# import a data pair

folder <- "../../datos/pairs/"
fileName <- "pairmeta.txt"
fileFolder <- paste(folder, fileName, sep="")
meta <- read.csv(fileFolder, sep="", header=F)
# From README file
colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
head(meta)
# which one is cause and which one effect?
indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
length(indxUni)
pairs <- meta$pairNumber[indxUni]

i <- 4
pair <- pairs[i]
fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
fileFolder <- paste(folder, fileName, sep="")
x <- read.csv(fileFolder, sep="", header=F)
indx <- which(meta$pairNumber==pair)  
res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
# order cause first, effect second
x <- x[, res]
colnames(x) <- c("x","y")
x <- as.matrix(x)

x <- apply(x, 2,  stdrize)

apply(x, 2, mean)
apply(x, 2, sd)

y <- x[,"y"]
x <- x[,"x"]



plot(x, y)
trainData <- constructData(as.matrix(x), y)


# implement hsic regression in forward direction


# qhsic - cross validate lambda, rbf kernel, saturate heuristic for sigma, max var for beta

nonOptimizableParams <- list(sigma=list(val=NULL, type="sat"), beta=list(val=NULL, type="var"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=seq(-8,8, length.out=50), log=TRUE))
nonDataParams <- list(kernelXs="rbfSigma", kernelXb="rbfBeta", kernelYhk="rbfKappa", kernelRg="rbfGamma")
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, avgy=NULL, indxInducing=NULL)
optimizeParams <- list(losses=list(qhsicLoss=qhsicLoss, sse=sse, rmse=rmse, corre=corre), numFolds=5, maxPoints=500)
heuristicSet <- getFixedParams
optimizeSet <- optHP.CV
qhsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.qhsic, predict.qhsic)
qhsic <- setParams(learner=qhsic, trainData)
qhsic <- qhsic$learn(qhsic)
pred.qhsic <- qhsic$predict(qhsic, trainData)
plot.emeley.1D(predList=list(qhsic=pred.qhsic))


# KRR - cross validate lambda, rbf kernel, median heuristic for sigma 
nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=seq(-8,8, length.out=50), log=TRUE))
nonDataParams <- list(kernelXs="rbfSigma", kernelXb="rbfBeta", kernelYhk="rbfKappa", kernelRg="rbfGamma")
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), numFolds=5)
heuristicSet <- getFixedParams
optimizeSet <- optHP.CV
krr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr)
krr1 <- setParams(krr1, trainData)
krr1 <- krr1$learn(krr1)
pred.krr1 <- krr1$predict(krr1, trainData)
plot.emeley.1D(predList=list(krr1=pred.krr1))

# KRR - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list(beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=seq(-8,8, length.out=50), log=TRUE), sigma=list(val=NULL, seq=seq(-8,8, length.out=17), log=TRUE))
nonDataParams <- list(kernelXs="rbfSigma", kernelXb="rbfBeta", kernelYhk="rbfKappa", kernelRg="rbfGamma")
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), numFolds=5)
heuristicSet <- getFixedParams
optimizeSet <- optHP.CV
krr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr)
krr2 <- setParams(krr2, trainData)
krr2 <- krr2$learn(krr2)
pred.krr2 <- krr2$predict(krr2, trainData)
plot.emeley.1D(predList=list(krr2=pred.krr2))




# hsic - cross validate lambda, rbf kernel, max var heuristic for sigma, max var for beta and gamma, max_iterations = 20, initializations=30

nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=seq(-8,0, length.out=10), log=TRUE))
nonDataParams <- list(kernelXs="rbfSigma", kernelXb="rbfBeta", kernelYhk="rbfKappa", kernelRg="rbfGamma")
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, avgy=NULL)
optimizeParams <- list(losses=list(hsicLoss=hsicLoss, sse=sse, rmse=rmse, corre=corre), numFolds=5, max_iterations1=10, max_iterations2=10, num_init=2)
heuristicSet <- getFixedParams
optimizeSet <- optHP.CV
hsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.hsic, predict.hsic)
pm <- proc.time()
hsic <- setParams(hsic, trainData)
proc.time() - pm # 
hsic <- hsic$learn(hsic)
pred.hsic <- hsic$predict(hsic, trainData)
plot.emeley.1D(predList=list(hsic=pred.hsic))

# GP.gptk - ML lambda, sigma, n_max = 50 (i.e. use fitc with 50 inducing points)


nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), approx="fitc", numActive=100, fixInducing=FALSE, iters=20000)
heuristicSet <- myGpCreate
optimizeSet <- myGpOptimise
gptk1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.gptk, predict.gptk)
gptk1 <- setParams(gptk1, trainData)
gptk1 <- gptk1$learn(gptk1)
pred.gptk <- gptk1$predict(gptk1, data=trainData)
plot.emeley.1D(predList=list(gptk1=pred.gptk))

# GP.gptk - ML lambda, sigma, n_max = 50 (i.e. use fitc with 50 LINEARLY SPACED inducing points)


nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), approx="fitc",  numActive=100, fixInducing=TRUE, iters=2000)
heuristicSet <- myGpCreate
optimizeSet <- myGpOptimise
gptk2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.gptk, predict.gptk)
gptk2 <- setParams(gptk2, trainData)
gptk2 <- gptk2$learn(gptk2)
pred.gptk <- gptk2$predict(gptk2, data=trainData)
plot.emeley.1D(predList=list(gptk2=pred.gptk))


# GP.tgp - treed LLM gp

nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), approx="fitc", numActive=90, iters=2000)
heuristicSet <- myGpCreate
optimizeSet <- myGpOptimise
tgp <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.tgp, predict.tgp)
tgp <- setParams(tgp, trainData)
tgp <- gptk$learn(tgp)
pred.tgp <- tpg$predict(tgp, data=trainData)
plot.emeley.1D(predList=list(gptk=pred.tgp))

# GP.lagp - locally approximate gp, active-learning cohn algorithm

# compare all learners
plot.emeley.1D(predList=list(krr1=pred.krr1, krr2=pred.krr2, qhsic=pred.qhsic, hsic=pred.hsic, gptk=pred.gptk))
