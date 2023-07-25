remove(list=ls())


server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"



repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_causaLearner/", sep="")
repos <- paste("/home/emiliano/Drives/erc/causaLearner/", sep="")
repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"


dir(repos)
setwd(repos)


source("./pkg_learner/func_learners_v5.R")
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")
library(FNN) #knnx.index, in gptk2



# simple model y = sin(x) + n
set.seed(3)
n <- 100
x <- runif(n, 0, 2*pi)
ny <- rnorm(n, 0, 0.5)
y <- sin(x)+ny + 25
y <- y - mean(y)
plot(x, y)
x <- cbind(x,y)
colnames(x) <- c("x","y")


# import a data pair TCEP

folder <- "data/TCEPs/pairs/"
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
trainData <- constructData(x=as.matrix(x), y=y[sample(length(y))])
trainData <- constructData(x=as.matrix(x), y=y)
trainDataYX <- constructData(x=as.matrix(y), y=x)
plot(trainData$x, trainData$y)
plot(trainDataYX$x, trainDataYX$y)

# implement hsic regression in forward direction


# qhsic - cross validate lambda, rbf kernel, saturate heuristic for sigma, max var for beta

nonOptimizableParams <- list(sigma=list(val=NULL, type="sat"), beta=list(val=NULL, type="var"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, avgy=NULL, indxInducing=NULL)
optimizeParams <- list(losses=list(qhsicLoss=list(func="qhsicLoss"), 
                                   sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), numFolds=5, testTrain="test", maxPoints=500)
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
qhsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.qhsic, predict.qhsic, resids=resids.add)
qhsic <- setParams(learner=qhsic, trainData)
qhsic <- qhsic$learn(qhsic)
pred.qhsic <- qhsic$predict(qhsic, trainData)
plot.emeley.1D(predList=list(qhsic=pred.qhsic))
plot(pred.qhsic$x, qhsic$resids(qhsic, pred.qhsic)[,"resid"])
plot(pred.qhsic$x, qhsic$resids(qhsic, pred.CV(qhsic, trainData)$test)[,"resid"])


# KRR - cross validate lambda, rbf kernel, median heuristic for sigma 
nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
krr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr_class, resids=resids.add_class)
krr1 <- setParams(learner=krr1, trainData)
krr1$hyperParams$data$optimizable$lambda$seq
krr1$hyperParams$data$optimizable$lambda$val
krr1 <- krr1$learn(krr1)
pred.krr1 <- krr1$predict(krr1, trainDataYX)
plot(pred.krr1$x, pred.krr1$gy_class, ylim=range(pred.krr1$gy_class, pred.krr1$gyh_class))
lines(pred.krr1$x[order(pred.krr1$x)], pred.krr1$gyh_class[order(pred.krr1$x)], col="red", lwd=3)
plot.emeley.1D(predList=list(krr1=pred.krr1))
plot(pred.krr1$x, krr1$resids(krr1, pred.krr1)[,"resid"])
plot(pred.krr1$x, krr1$resids(krr1, pred.CV(krr1, trainData)$test)[,"resid"])

# weighted KRR - cross validate lambda, rbf kernel, median heuristic for sigma 
ws <- kern_rbf(x=trainData$x, sigma=1e-7)
ws <- ws[,8]
ws <- ws/sum(ws)
hist(ws)
nonOptimizableParams <- list(ws=list(val=ws),sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)))


nonOptimizableParams <- list(lambda=list(val=1),ws=list(val=ws),sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list()

nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
wkrr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.wkrr, predict.wkrr, resids=resids.add)
wkrr1 <- setParams(learner=wkrr1, trainData)
wkrr1$hyperParams$data$optimizable$lambda$seq
wkrr1$hyperParams$data$optimizable$lambda$val
wkrr1 <- wkrr1$learn(learner=wkrr1)
pred.wkrr1 <- wkrr1$predict(wkrr1, trainData)
plot(pred.wkrr1$x, pred.wkrr1$gy, ylim=range(pred.wkrr1$gy, pred.wkrr1$gyh))
lines(pred.wkrr1$x[order(pred.wkrr1$x)], pred.wkrr1$gyh[order(pred.wkrr1$x)], col="red", lwd=3)
plot.emeley.1D(predList=list(wkrr1=pred.wkrr1))
plot(pred.wkrr1$x, wkrr1$resids(wkrr1, pred.wkrr1)[,"resid"])
plot(pred.wkrr1$x, wkrr1$resids(wkrr1, pred.CV(wkrr1, trainData)$test)[,"resid"])


# KRR - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list(beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)), sigma=list(val=NULL, seq=seq(-8,8, length.out=17)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
krr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)
krr2 <- setParams(krr2, trainData)
krr2 <- krr2$learn(krr2)
pred.krr2 <- krr2$predict(krr2, trainData)
plot.emeley.1D(predList=list(krr2=pred.krr2))
plot(pred.krr2$x, krr2$resids(krr2, pred.krr2)[,"resid"])
plot(pred.krr2$x, krr2$resids(krr2, pred.CV(krr2, trainData)$test)[,"resid"])

# KRR - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list()
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)), sigma=list(val=NULL, length.out=17))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf2b"
optimizeSet <- "optHP.CV"
krr2b <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)
krr2b <- setParams(krr2b, trainData)
krr2b <- krr2b$learn(krr2b)
pred.krr2b <- krr2b$predict(krr2b, trainData)
plot.emeley.1D(predList=list(krr2b=pred.krr2b))
plot(pred.krr2b$x, krr2$resids(krr2b, pred.krr2b)[,"resid"])
plot(pred.krr2b$x, krr2b$resids(krr2b, pred.CV(krr2b, trainData)$test)[,"resid"])


# KRR - rbf kernel, cross validate lambda and sigma - find range of sigmaX 
nonOptimizableParams <- list()
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,0, length.out=50)), sigma.rbf.X=list(val=NULL, seq=NULL, length.out=5))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
krr3 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)
krr3 <- setParams(learner=krr3, trainData)
krr3 <- krr3$learn(krr3)
pred.krr3 <- krr3$predict(krr3, trainData)
plot.emeley.1D(predList=list(krr3=pred.krr3))
plot(pred.krr3$x, qhsic$resids(krr3, pred.krr3)[,"resid"])
plot(pred.krr3$x, qhsic$resids(krr3, pred.CV(krr3, trainData)$test)[,"resid"])

# hsic - cross validate lambda, rbf kernel, max var heuristic for sigma, max var for beta and gamma, max_iterations = 20, initializations=30

nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,0, length.out=10)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, avgy=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test", max_iterations1=10, max_iterations2=100, num_init=2)
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
hsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.hsic, predict.hsic, resids=resids.add)
pm <- proc.time()
hsic <- setParams(hsic, trainData, mc_cores=4)
proc.time() - pm 
# 7 mins for 10-100 iterations, 2 initializations, 4 cores on my lap
hsic <- hsic$learn(hsic)
pred.hsic <- hsic$predict(hsic, trainData)
plot.emeley.1D(predList=list(hsic=pred.hsic))
plot(pred.hsic$x, qhsic$resids(hsic, pred.hsic)[,"resid"])
plot(pred.hsic$x, qhsic$resids(hsic, pred.CV(hsic, trainData)$test)[,"resid"])

# quantile kernel regression learner - cross validate nothing
nonOptimizableParams <- list(taus=list(val=list(c(0.05, 0.25, 0.5, 0.75, 0.95))), lambda=list(val=1), sigmaX=list(val=12791.43))
optimizableParams <- list()
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigmaX")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)
kqr1 <- setParams(learner=kqr1, trainData)
kqr1 <- kqr1$learn(learner=kqr1)
pred.kqr1 <- kqr1$predict(learner=kqr1, data=trainData)
plot.emeley.1D(predList=list(kqr1=pred.kqr1))
plot(pred.kqr1$x, kqr1$resids(kqr1, pred.kqr1)[,"resid"])
plot(pred.kqr1$x, kqr1$resids(kqr1, pred.CV(kqr1, trainData)$test)[,"resid"])
liks <- kqr1$resids(kqr1, pred.CV(kqr1, trainData)$test)[,"pred"]
negLogLik(learner=kqr1, pred.CV(kqr1, trainData)$test)
sum(liks)

# quantile kernel regression learner - cross validate lambda, fixed sigma val
nonOptimizableParams <- list(taus=list(val=list(c(0.25, 0.5, 0.75))), sigmaX=list(val=1))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.1, 1,10, 100, 1000)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigmaX")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(negLogLik=list(func=negLogLik), 
                                   hsicLoss2=list(func=hsicLoss2), 
                                   pinball=list(func=pinball)), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)
pm <- proc.time()
kqr2 <- setParams(learner=kqr2, trainData, mc_cores=4)
proc.time() - pm 
# 10 mins for 1000 pts, 3 taus, 5 lambdas on one core
# 7 mins for 1000 pts, 3 taus, 5 lambdas on one core
getHyperPar(kqr2, "lambda")
kqr2 <- kqr2$learn(learner=kqr2)
pred.kqr2 <- kqr2$predict(learner=kqr2, data=trainData)
plot.emeley.1D(predList=list(kqr2=pred.kqr2))
plot(pred.kqr2$x, kqr2$resids(kqr2, pred.kqr2)[,"resid"])
plot(pred.kqr2$x, kqr2$resids(kqr2, pred.CV(kqr2, trainData)$test)[,"resid"])

# quantile regression learner - cross validate lambda, data-driven sigma val
nonOptimizableParams <- list(taus=list(val=list(c(0.05, 0.25, 0.5, 0.75, 0.95))), sigma.rbf.X=list(val=NULL))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.1, 1,10, 100, 1000)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(negLogLik=list(func=negLogLik), 
                                   hsicLoss2=list(func=hsicLoss2), 
                                   pinball=list(func=pinball)), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr3 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)
kqr3 <- setParams(learner=kqr3, trainData)
getHyperPar(kqr3, "lambda")
getHyperPar(kqr3, "sigma.rbf.X")
kqr3 <- kqr3$learn(learner=kqr3)
pred.kqr3 <- kqr3$predict(learner=kqr3, data=trainData)
plot.emeley.1D(predList=list(kqr3=pred.kqr3))
plot(pred.kqr3$x, kqr2$resids(kqr3, pred.kqr3)[,"resid"])
plot(pred.kqr3$x, kqr2$resids(kqr3, pred.CV(learner=kqr3, data=trainData)$test)[,"resid"])

# quantile regression learner - cross validate lambda, sigma fixed seq
nonOptimizableParams <- list(taus=list(val=list(c(0.05, 0.25, 0.5, 0.75, 0.95))))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.1, 1,10, 100, 1000)), sigma.rbf.X=list(val=NULL, seq=c(0.1, 1, 10)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(negLogLik=list(func=negLogLik), 
                                   hsicLoss2=list(func=hsicLoss2), 
                                   pinball=list(func=pinball)), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr4 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)
kqr4 <- setParams(learner=kqr4, trainData)
getHyperPar(kqr4, "lambda")
getHyperPar(kqr4, "sigma.rbf.X")
kqr4 <- kqr4$learn(learner=kqr4)
pred.kqr4 <- kqr4$predict(learner=kqr4, data=trainData)
plot.emeley.1D(predList=list(kqr4=pred.kqr4))
plot(pred.kqr4$x, kqr4$resids(kqr4, pred.kqr4)[,"resid"])
plot(pred.kqr4$x, kqr4$resids(kqr4, pred.CV(kqr4, trainData)$test)[,"resid"])

# quantile regression learner - cross validate lambda, sigma data-driven seq
nonOptimizableParams <- list(taus=list(val=list(c(0.25, 0.5, 0.75))))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.1, 1,10, 100, 1000)), sigma.rbf.X=list(val=NULL, seq=NULL, length.out=5))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(negLogLik=list(func=negLogLik), 
                                   hsicLoss2=list(func=hsicLoss2), 
                                   pinball=list(func=pinball)), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr5 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)
kqr5 <- setParams(learner=kqr5, trainData)
getHyperPar(kqr5, "lambda")
getHyperPar(kqr5, "sigma.rbf.X")
kqr5 <- kqr5$learn(learner=kqr5)
pred.kqr5 <- kqr5$predict(learner=kqr5, data=trainData)
plot.emeley.1D(predList=list(kqr5=pred.kqr5))
plot(pred.kqr5$x, kqr5$resids(kqr5, pred.kqr5)[,"resid"])
plot(pred.kqr5$x, kqr5$resids(kqr5, pred.CV(kqr5, trainData)$test)[,"resid"])

# marginal quantile learner to learn likelihood of marginal
nonOptimizableParams <- list(taus=list(val=NULL, length.out=5))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(qs=NULL)
optimizeParams <- list(losses=list(), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_tauQuad"
optimizeSet <- "optHP.CV"
qr_marg <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.qr_marg, predict.qr_marg, resids=resids.kqr)
qr_marg <- setParams(learner=qr_marg, trainData)
qr_marg <- qr_marg$learn(learner=qr_marg)
pred.qr_marg <- qr_marg$predict(learner=qr_marg, data=trainData)
head(qr_marg$resids(qr_marg, pred.qr_marg))

# quantile copula regression learner - cross validate nothing
nonOptimizableParams <- list(taus=list(val=NULL, length.out=5))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(cop=NULL)
optimizeParams <- list(losses=list(), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_tauQuad"
optimizeSet <- "optHP.CV"
cqr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cqr, predict.cqr, resids=resids.kqr)
cqr1 <- setParams(learner=cqr1, trainData)
cqr1 <- cqr1$learn(learner=cqr1)
pred.cqr1 <- cqr1$predict(learner=cqr1, data=trainData)
plot.emeley.1D(predList=list(cqr1=pred.cqr1))
plot(pred.cqr1$x, cqr1$resids(cqr1, pred.cqr1)[,"resid"])
dhsic.test(pred.cqr1$x, cqr1$resids(cqr1, pred.cqr1)[,"resid"], method="gamma")$p.value
plot(pred.cqr1$x, cqr1$resids(cqr1, pred.CV(cqr1, trainData)$test)[,"resid"])
liks <- cqr1$resids(cqr1, pred.CV(cqr1, trainData)$test)[,"pred"]
negLogLik(learner=cqr1, pred.CV(cqr1, trainData)$test)
sum(liks)

# quantile forest regression learner - cross validate nothing
nonOptimizableParams <- list(taus=list(val=NULL, length.out=5), 
                             nodesize=list(val=10), 
                             sampsize=list(val=50))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_tauQuad"
optimizeSet <- "optHP.CV"
fqr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.fqr, predict.fqr, resids=resids.kqr)
fqr1 <- setParams(learner=fqr1, trainData)
fqr1 <- fqr1$learn(learner=fqr1)
pred.fqr1 <- fqr1$predict(learner=fqr1, data=trainData)
plot.emeley.1D(predList=list(fqr1=pred.fqr1))
plot(pred.fqr1$x, fqr1$resids(fqr1, pred.fqr1)[,"resid"])
dhsic.test(pred.fqr1$x, cqr1$resids(fqr1, pred.fqr1)[,"resid"], method="gamma")$p.value
plot(pred.fqr1$x, cqr1$resids(fqr1, pred.CV(fqr1, trainData)$test)[,"resid"])
liks <- fqr1$resids(fqr1, pred.CV(fqr1, trainData)$test)[,"pred"]
negLogLik(learner=fqr1, pred.CV(fqr1, trainData)$test)
sum(liks)


# quantile NN regression learner - cross validate nothing
nonOptimizableParams <- list(taus=list(val=NULL, length.out=5), 
                             n_hidden=list(val=10), 
                             n_hidden2=list(val=50),
                             n_trials=list(val=1))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(), numFolds=4, testTrain="test", iter_max=1000)
heuristicSet <- "getFixedParams_tauQuad"
optimizeSet <- "optHP.CV"
nnqr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.nnqr, predict.nnqr, resids=resids.kqr)
nnqr1 <- setParams(learner=nnqr1, trainData)
pm <- proc.time()
nnqr1 <- nnqr1$learn(learner=nnqr1)
proc.time() - pm
pred.nnqr1 <- nnqr1$predict(learner=nnqr1, data=trainData)
plot.emeley.1D(predList=list(nnqr1=pred.nnqr1))
plot(pred.nnqr1$x, pred.nnqr1$gy)
o <- order(pred.nnqr1$x)
lines(pred.nnqr1$x[o,], pred.nnqr1$gyh[o,5], col="red", lwd=3)
res <- pred.nnqr1$gy[o,]-pred.nnqr1$gyh[o,5]

library(ald)
tau <- nnqr1$hyperParams$data$non_optimizable$taus$val[4]
xx <- seq(min(res)-mean(res), max(res)-mean(res), length.out=100)
trainDataAux <- constructData(x=matrix(xx,length(xx),1), y=rnorm(100))
predq <- nnqr1$predict(learner=nnqr1, data=trainDataAux)
pars <- mleALD(res)
yy <- dALD(xx, mu=pars$par[1], sigma=pars$par[2], p=pars$par[3])
#yy <- tau*(1-tau)*exp(-qs)

hist(res, prob=T, ylim=c(0,max(yy)))
lines(xx, yy, col="red")



plot(pred.nnqr1$x, nnqr1$resids(nnqr1, pred.nnqr1)[,"resid"])
dhsic.test(pred.nnqr1$x, nnqr1$resids(learner=nnqr1, pred=pred.nnqr1)[,"resid"], method="gamma")$p.value
plot(pred.nnqr1$x, nnqr1$resids(nnqr1, pred.CV(nnqr1, trainData)$test)[,"resid"])
liks <- fqr1$resids(nnqr1, pred.CV(nnqr1, trainData)$test)[,"pred"]
negLogLik(learner=nnqr1, pred.CV(nnqr1, trainData)$test)
sum(liks)



# GP.gptk - ML lambda, sigma, n_max = 50 (i.e. use fitc with 50 inducing points)

nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), approx="fitc", numActive=100, fixInducing=FALSE, iters=20000)
heuristicSet <- "myGpCreate"
optimizeSet <- "myGpOptimise"
gptk1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.gptk, predict.gptk, resids=resids.add)
gptk1 <- setParams(learner=gptk1, trainData)
gptk1 <- gptk1$learn(gptk1)
pred.gptk <- gptk1$predict(gptk1, data=trainData)
plot.emeley.1D(predList=list(gptk1=pred.gptk))
plot(pred.qhsic$x, qhsic$resids(qhsic, pred.qhsic)[,"resid"])
plot(pred.qhsic$x, qhsic$resids(qhsic, pred.CV(qhsic, trainData)$test)[,"resid"])

# GP.gptk - ML lambda, sigma, n_max = 50 (i.e. use fitc with 50 LINEARLY SPACED inducing points)

nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), approx="fitc",  numActive=100, fixInducing=TRUE, iters=2000)
heuristicSet <- "myGpCreate"
optimizeSet <- "myGpOptimise"
gptk2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.gptk, predict.gptk, resids=resids.add)
gptk2 <- setParams(gptk2, trainData)
gptk2 <- gptk2$learn(gptk2)
pred.gptk <- gptk2$predict(gptk2, data=trainData)
plot.emeley.1D(predList=list(gptk2=pred.gptk))

wishInducing <- makeGrid(trainData$x, 100)
indxs <- FNN:::knnx.index(data=trainData$x, query=wishInducing, k=1, algorithm="brute")


# GP.tgp - treed LLM gp

nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=sse, rmse=rmse, corre=corre), approx="fitc", numActive=90, iters=2000)
heuristicSet <- "myGpCreate"
optimizeSet <- "myGpOptimise"
tgp <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.tgp, predict.tgp)
tgp <- setParams(tgp, trainData)
tgp <- gptk$learn(tgp)
pred.tgp <- tpg$predict(tgp, data=trainData)
plot.emeley.1D(predList=list(gptk=pred.tgp))

# GP.lagp - locally approximate gp, active-learning cohn algorithm

# compare all learners
plot.emeley.1D(predList=list(krr1=pred.krr1, krr2=pred.krr2, qhsic=pred.qhsic, hsic=pred.hsic, gptk=pred.gptk))


# CLASSIFICATION

set.seed(3)
n <- 100
x <- runif(n, 0, 2*pi)
x2 <- runif(n, 0, 2*pi)
ny <- rnorm(n, 0, 0.1)
ny2 <- rnorm(n, 0, 0.1)
y <- sin(x) + ny
y2 <- sin(x2) + ny2
y <- (y<mean(y))*1
y2 <- (y2<mean(y2))*1
plot(x, y)
trainData <- constructData(x=as.matrix(x)-mean(x), y=y)
testData <- constructData(x=as.matrix(x2)-mean(x2), y=y2)

sig <- 0.01
K <- kern_rbf(trainData$x, sigma=sig)
K <- as.kernelMatrix(K)
svmFit <- ksvm(x=K, y=factor(y), C=0.01, nu=-0.01, type=c("C-svc","nu-svc")[2], prob.model=T)
Ktest <- kern_rbf(testData$x,trainData$x, sigma=sig)
Ktest <- as.kernelMatrix(Ktest)
preds <- predict(svmFit, Ktest)
predsProb <- predict(svmFit, Ktest, type="probabilities")
head(predsProb)
table(preds, y2)
plot(x2, y2, col=c("red","green")[(preds==y2)*1+1])
MLmetrics:::LogLoss(predsProb[,2], y2)
MLmetrics:::LogLoss(as.numeric(as.character(preds)), y2)


## Demo of the plot function
x <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
y <- matrix(c(rep(1,60),rep(-1,60)))

svp <- ksvm(x,y,type="C-svc")
plot(svp,data=x)


### Use kernelMatrix
K <- as.kernelMatrix(crossprod(t(x)))

svp2 <- ksvm(K, y, type="C-svc")

svp2



indxU <- sample(length(trainData$y), round(0.3*length(trainData$y)))
trainDataSSL <- constructData(x=as.matrix(x)-mean(x), y=y, indxU=indxU)


# calculate phi features externally
phix <- rff(x=as.matrix(x), num_f=1000, seed=1234, p_w="rnorm2", map="cos", sigma=10)
trainData2 <- constructData(x=phix, y=y)

#source("./pkg_learner/func_learners_v3.R")

# logistic regression with vanilla feature - cross validate nout  
nonOptimizableParams <- list(featureX=list(val="feat_van"))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=10, testTrain="test")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
logReg1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logReg, predict.logReg, resids=resids.add, makeFeature=makePhi)
logReg1 <- setParams(learner=logReg1, trainData)
getHyperPar(logReg1, "num_f")
getHyperPar(logReg1, "sigma")
getHyperPar(logReg1, "lambda")
pm <- proc.time()
logReg1 <- logReg1$learn(learner=logReg1)
proc.time() - pm #
pred.logReg1 <- logReg1$predict(learner=logReg1, data=trainData)
plot.emeley.1D(predList=list(logReg1=pred.logReg1))
plot(pred.logReg1$gy_class, pred.logReg1$gyh_class)
negCE(learner=logReg1, pred=pred.logReg1)

# logistic regression with rbf RFF - cross validate lambda internally - phi calculated externally 
nonOptimizableParams <- list(lambda=list(val=10^seq(-4,1, length.out=6)))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func=negCE)), numFolds=10, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
logRegInt1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logRegInt, predict.logRegInt, resids=resids.add)
logRegInt1 <- setParams(learner=logRegInt1, trainData2)
getHyperPar(logRegInt1, "lambda")
pm <- proc.time()
logRegInt1 <- logRegInt1$learn(learner=logRegInt1)
proc.time() - pm #
getHyperPar(logRegInt1, "lambda")
pred.logRegInt1 <- logRegInt1$predict(learner=logRegInt1, data=trainData2)
plot.emeley.1D(predList=list(logRegInt1=pred.logRegInt1))
plot(pred.logRegInt1$gy_class, pred.logRegInt1$gyh_class)
negCE(logRegInt1, pred.logRegInt1)


# logistic regression with rbf RFF - cross validate lambda  
nonOptimizableParams <- list(sigma=list(val=10), num_f=list(val=1000), p_w=list(val="rnorm2"), seed=list(val=1234), map=list(val="cos"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-4,1, length.out=6)))
nonDataParams <- list(featureX=list(name="rff_rbf", pars=c(sigma="sigma", p_w="p_w", seed="seed", map="map", num_f="num_f")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func=negCE)), numFolds=10, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
logKReg1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logKReg, predict.logKReg, resids=resids.add, makeFeature=makePhi)
logKReg1 <- setParams(learner=logKReg1, trainData)
getHyperPar(logKReg1, "lambda")
getHyperPar(logKReg1, "num_f")
getHyperPar(logKReg1, "sigma")
pm <- proc.time()
logKReg1 <- logKReg1$learn(learner=logKReg1)
proc.time() - pm #
pred.logKReg1 <- logKReg1$predict(learner=logKReg1, data=trainData)
plot.emeley.1D(predList=list(logKReg1=pred.logKReg1))
plot(pred.logKReg1$gy, pred.logKReg1$gyh)
negCE(logKReg1, pred.logKReg1)

# logistic regression with rbf RFF - cross validate lambda and sigma  features 
nonOptimizableParams <- list(p_w=list(val="rnorm2"), seed=list(val=1234), map=list(val="cos"), num_f=list(val=1000))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-4,1, length.out=6)),sigma=list(val=NULL, seq=10^seq(-4, 5, length.out=10)))
nonDataParams <- list(featureX=list(name="rff_rbf", pars=c(sigma="sigma", p_w="p_w", seed="seed", map="map", num_f="num_f")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=10, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
logKReg2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logKReg, predict.logKReg, resids=resids.add, makeFeature=makePhi)
logKReg2 <- setParams(learner=logKReg2, trainData)
getHyperPar(logKReg2, "num_f")
getHyperPar(logKReg2, "sigma")
getHyperPar(logKReg2, "lambda")
pm <- proc.time()
logKReg2 <- logKReg2$learn(learner=logKReg2)
proc.time() - pm #
pred.logKReg2 <- logKReg2$predict(learner=logKReg2, data=trainData)
plot.emeley.1D(predList=list(logKReg2=pred.logKReg2))
plot(pred.logKReg2$gy, pred.logKReg2$gyh)
negCE(logKReg2, pred.logKReg2)


# logistic regression with rbf RFF - cross validate lambda, rff rbf features  = direct sum of rbf featurs 
nonOptimizableParams <- list(sigma=list(val=c(1,10,100,1000,10000)), p_w=list(val="rnorm2"), seed=list(val=1234), map=list(val="cos"), num_f=list(val=200))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-4,1, length.out=6)))
nonDataParams <- list(featureX=list(name="rff_rbf", pars=c(sigma="sigma", p_w="p_w", seed="seed", map="map", num_f="num_f")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=10, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
logKReg3 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logKReg, predict.logKReg, resids=resids.add, makeFeature=makePhi)
logKReg3 <- setParams(learner=logKReg3, trainData)
getHyperPar(logKReg3, "num_f")
getHyperPar(logKReg3, "sigma")
getHyperPar(logKReg3, "lambda")
pm <- proc.time()
logKReg3 <- logKReg3$learn(learner=logKReg3)
proc.time() - pm #
pred.logKReg3 <- logKReg3$predict(learner=logKReg3, data=trainData)
plot.emeley.1D(predList=list(logKReg3=pred.logKReg3))
plot(pred.logKReg3$gy, pred.logKReg3$gyh)
negCE(logKReg3, pred.logKReg3)


# logistic regression with NN features - cross validate depth and breadth 
nonOptimizableParams <- list()
optimizableParams <- list(depth=list(val=NULL, seq=c(5,10,50)), breadth=list(val=NULL, seq=c(5,10)))
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model_nn=NULL, model_log=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
logNNReg1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logNNReg, predict.logNNReg, resids=resids.add, makeFeature=makeNNFeats)
logNNReg1 <- setParams(learner=logNNReg1, trainData)
getHyperPar(logNNReg1, "depth")
getHyperPar(logNNReg1, "breadth")
pm <- proc.time()
logNNReg1 <- logNNReg1$learn(learner=logNNReg1)
proc.time() - pm #
pred.logNNReg1 <- logNNReg1$predict(learner=logNNReg1, data=trainData)
plot.emeley.1D(predList=list(logNNReg1=pred.logNNReg1))
plot(pred.logNNReg1$gy, pred.logNNReg1$gyh)
negCE(logNNReg1, pred.logNNReg1)


# kernel SVM with rbf  - cross validate lambda, sigma and nu  features 
nonOptimizableParams <- list()
optimizableParams <- list(nu=list(val=NULL, seq=seq(0,1, length.out=6)),lambda=list(val=NULL, seq=c(0.001,0.01,0.1)), sigma.rbf.X=list(val=NULL, seq=c(0.001,0.01, 0.1)))
#nonOptimizableParams <- list(nu=list(val=0.2), lambda=list(val=0.0001), sigma.rbf.X=list(val=0.01))
#optimizableParams <- list()
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE3=list(func="negCE3")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
ksvm1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.ksvm, predict.ksvm, resids=resids.add)
pm <- proc.time()
ksvm1 <- setParams(learner=ksvm1, trainData)
proc.time() - pm #
getHyperPar(ksvm1, "lambda")
getHyperPar(ksvm1, "sigma.rbf.X")
getHyperPar(ksvm1, "nu")
pm <- proc.time()
ksvm1 <- ksvm1$learn(learner=ksvm1)
proc.time() - pm # 8 secs for 1000 pts
ksvm1$learnParams$model
pred.ksvm1 <- ksvm1$predict(learner=ksvm1, data=testData)
plot.emeley.1D(predList=list(ksvm1=pred.ksvm1))
plot(pred.ksvm1$x_class, pred.ksvm1$gy_class)
plot(pred.ksvm1$x_class, pred.ksvm1$gyh_class)
plot(pred.ksvm1$gyh_class,pred.ksvm1$gy_class)
negCE3(learner=ksvm1, pred=pred.ksvm1)
negCE3(learner=ksvm1, pred=ksvm1$predict(learner=ksvm1, data=trainData))
LogLoss(pred.ksvm1$gyh_class, pred.ksvm1$gy_class)




predData <- getSubset(trainDataSSL, trainDataSSL$indxU)
trainData <- getSubset(trainDataSSL, -trainDataSSL$indxU)

# laplacian KRR for semi super vised learning
#nonOptimizableParams <- list(Gtype=list(val="isomap", pars=list(adjacency_k=6)), normL=list(val=TRUE), kernelX=list(val="kern_rbf"), lambda1=list(val=0.01), lambda2=list(val=0.02), sigma.rbf.X=list(val=1))
#optimizableParams <- list()
nonOptimizableParams <- list(Gtype=list(val="isomap", pars=list(adjacency_k=6)), normL=list(val=TRUE)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
optimizableParams <- list(kernelX=list(val=NULL, seq=c("kern_rbf")), 
                          lambda1=list(val=NULL, seq=c(0,10^seq(-4,1, length.out=6))),#6
                          lambda2=list(val=NULL, seq=c(10^seq(-4,1, length.out=6))),#6
                          sigma.rbf.X=list(val=NULL, seq=c(10^seq(-4, 5, length.out=10))))#10

nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, K=NULL)
optimizeParams <- list(losses=list(myAUC=list(func="myAUC"), 
                                   rmse=list(func="rmse")), numFolds=5, testTrain="test", mainLoss="rmse")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
lapSSLkrr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSSLkrr1, predict.lapSSLkrr1, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)
lapSSLkrr1 <- setParams(learner=lapSSLkrr1, trainData=trainDataSSL)

getHyperPar(lapSSLkrr1, "kernelX")
getHyperPar(lapSSLkrr1, "sigma.rbf.X")
getHyperPar(lapSSLkrr1, "lambda1")
getHyperPar(lapSSLkrr1, "lambda2")
getHyperPar(lapSSLkrr1, "Gtype")
pm <- proc.time()
lapSSLkrr1 <- lapSSLkrr1$learn(learner=lapSSLkrr1)
proc.time() - pm #
pred.lapSSLkrr1 <- lapSSLkrr1$predict(learner=lapSSLkrr1, data=predData)
train.lapSSLkrr1 <- lapSSLkrr1$predict(learner=lapSSLkrr1, data=trainData)
dim(pred.lapSSLkrr1$gyh)
dim(train.lapSSLkrr1$gyh)
head(pred.lapSSLkrr1$gyh)
head(train.lapSSLkrr1$gyh)
head(pred.lapSSLkrr1$gy)
head(train.lapSSLkrr1$gy)

plot.emeley.1D(predList=list(lapSSLkrr1=pred.lapSSLkrr1))
plot.emeley.1D(predList=list(lapSSLkrr1=train.lapSSLkrr1))
plot(pred.lapSSLkrr1$gy, pred.lapSSLkrr1$gyh)
plot(train.lapSSLkrr1$gy, train.lapSSLkrr1$gyh)
1-myAUC(learner=lapSSLkrr1, pred=pred.lapSSLkrr1)
1-myAUC(learner=lapSSLkrr1, pred=train.lapSSLkrr1)
rmse(lapSSLkrr1, pred.lapSSLkrr1)
rmse(lapSSLkrr1, train.lapSSLkrr1)





predData <- getSubset(trainDataSSL, trainDataSSL$indxU)
trainData <- getSubset(trainDataSSL, -trainDataSSL$indxU)

# laplacian SVM for semi super vised learning
#nonOptimizableParams <- list(Gtype=list(val="adj", pars=list(adjacency_k=6)), normL=list(val=TRUE), kernelX=list(val="kern_rbf"), lambda1=list(val=10), lambda2=list(val=1e-4), sigma.rbf.X=list(val=2), eps=list(val=1e-09))
#optimizableParams <- list()
nonOptimizableParams <- list(Gtype=list(val="isomap", pars=list(adjacency_k=6)), normL=list(val=TRUE), eps=list(val=1e-08)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
optimizableParams <- list(kernelX=list(val=NULL, seq=c("kern_rbf")), 
                          lambda1=list(val=NULL, seq=c(0,10^seq(-4,1, length.out=6))),#6
                          lambda2=list(val=NULL, seq=c(10^seq(-4,1, length.out=6))),#6
                          sigma.rbf.X=list(val=NULL, seq=c(10^seq(-4, 5, length.out=10))))#10

nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, K=NULL, bias=NULL)
optimizeParams <- list(losses=list(myAUC=list(func="myAUC"), 
                                   hingeLoss=list(func="hingeLoss")), numFolds=5, testTrain="test", mainLoss="hingeLoss")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
lapSVM1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSVM, predict.lapSVM, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)
lapSVM1 <- setParams(learner=lapSVM1, trainData=trainDataSSL)

getHyperPar(lapSVM1, "kernelX")
getHyperPar(lapSVM1, "sigma.rbf.X")
getHyperPar(lapSVM1, "lambda1")
getHyperPar(lapSVM1, "lambda2")
getHyperPar(lapSVM1, "Gtype")
pm <- proc.time()
lapSVM1 <- lapSVM1$learn(learner=lapSVM1)
proc.time() - pm #
pred.lapSVM1 <- lapSVM1$predict(learner=lapSVM1, data=predData)
train.lapSVM1 <- lapSVM1$predict(learner=lapSVM1, data=trainData)
dim(pred.lapSVM1$gyh)
dim(train.lapSVM1$gyh)
head(pred.lapSVM1$gyh)
head(train.lapSVM1$gyh)
head(pred.lapSVM1$gy)
head(train.lapSVM1$gy)

plot.emeley.1D(predList=list(lapSVM1=pred.lapSVM1))
plot.emeley.1D(predList=list(lapSVM1=train.lapSVM1))
plot(pred.lapSVM1$gy, pred.lapSVM1$gyh)
plot(train.lapSVM1$gy, train.lapSVM1$gyh)
1-myAUC(learner=lapSVM1, pred=pred.lapSVM1)
1-myAUC(learner=lapSVM1, pred=train.lapSVM1)
hingeLoss(lapSVM1, pred.lapSVM1)
hingeLoss(lapSVM1, train.lapSVM1)

# check against LaplacianSVM for when sigma=1, lambda1=lambda2=0.01

myLaplacianSVM <- function (X, y, X_u = NULL, lambda = 1, gamma = 1, scale = TRUE, 
                            sigma = 1, adjacency_distance = "euclidean", 
                            adjacency_k = 6, normalized_laplacian = FALSE, eps = 1e-09) 
{
  ModelVariables <- RSSL:::PreProcessing(X = X, y = y, X_u = X_u, 
                                  scale = scale, intercept = FALSE, x_center = TRUE)
  #X <- ModelVariables$X
  #X_u <- ModelVariables$X_u
  Y <- ModelVariables$Y[, 1, drop = FALSE]
  scaling <- ModelVariables$scaling
  classnames <- ModelVariables$classnames
  modelform <- ModelVariables$modelform
  if (length(classnames) != 2) 
    stop("Dataset does not contain 2 classes")
  y <- as.numeric((Y * 2) - 1)
  y2 <<- y
  l <- nrow(X)
  m <- ncol(X)
  u <- nrow(X_u)
  n <- l + u
  if(TRUE){#(inherits(kernel, "kernel")) {
    Xtrain <- rbind(X, X_u)
    x2 <<- Xtrain
    K <- kern_rbf(x=Xtrain, sigma=sigma)#kernlab:::kernelMatrix(kernel, Xtrain, Xtrain)
    K2 <<- K
    W <- adjacency_knn(D=1-K, k = adjacency_k)
    W2 <<- W
    d <- rowSums(W)
    L <- diag(d) - W
    if (normalized_laplacian) {
      L <- diag(1/sqrt(d)) %*% L %*% diag(1/sqrt(d))
    }
    L2 <<-L
    Y <- diag(y)
    Y2 <<- Y
    J <- cbind(diag(l), matrix(0, l, u))
    J2 <<- J
    Qprime <- solve(2 * lambda * diag(l + u) + 2 * (gamma/((l + u)^2)) * L %*% K, t(J) %*% Y)
    Qprime2 <<- Qprime
    Q <- Y %*% J %*% K %*% Qprime
    Q2 <<- Q
    Amat <- diag(nrow(X))
    Amat <- t(rbind(y, Amat, -Amat))
    bvec <- c(rep(0, nrow(X) + 1), rep(-1/l, nrow(X)))
    beta <- solve.QP(Q + diag(l) * eps, rep(1, l), Amat, bvec, meq = 1)$solution
    beta2 <<- beta
    if (any(is.nan(beta))) 
      stop("Quadratic Programming problem returned: NaN. Try different hyperparameters?")
    alpha <- (Qprime %*% beta)
    alpha2 <<- alpha
    SVs <- (beta > 1e-08/l) & (beta < 1/l - 1e-08/l)
    C <- 1/(2 * lambda * l)
    SVplus <- (alpha > 0.1) & (alpha < C - 0.1)
    SVmin <- (alpha < -0.1) & (alpha > -C + 0.1)
    bias <- median(K[(1:l)[SVs], ] %*% alpha - y[SVs])
    bias <- 0
    if (sum(SVplus) > 0 && sum(SVmin) > 0) {
      print("enters a")
      bias <- -median(c(K[SVmin, ] %*% alpha + 1, K[SVplus, ] %*% alpha - 1))
      print(bias)
    }
    else {
      if (sum(SVplus) > 0) {
        print("enters b")
        bias <- -median(K[SVplus, ] %*% alpha - 1)
        print(bias)
      }
      if (sum(SVmin) > 0) {
        print("enters c")
        bias <- -median(K[SVmin, ] %*% alpha + 1)
        print(bias)
      }
    }
    print("d")
    bias <- -median(K[(1:l)[SVs], ] %*% alpha - y[SVs])
    print(bias)
  }
  else {
    stop("No appropriate kernel function from kernlab supplied. See, for instance, the help of vanilladot()")
  }
  return(new("SVM", alpha = as.numeric(alpha), bias = bias, 
             Xtrain = Xtrain, kernel = kernel, scaling = scaling, 
             modelform = modelform, classnames = classnames, intercept = FALSE, 
             name = "Laplacian SVM"))
}

x <- trainDataSSL$x
y <- trainDataSSL$y
y <- 1*(!y) 
indxU <- trainDataSSL$indxU
class_lap <- myLaplacianSVM(X=as.matrix(x[-indxU]), y=as.factor(y[-indxU]), X_u=as.matrix(x[indxU]),
                            sigma=2,
                            lambda=0.02,gamma=0.01,
                            normalized_laplacian = TRUE,
                            scale=FALSE)
indxMixToSep <- c(setdiff(1:length(y), indxU), indxU)
plot(class_lap@alpha, lapSVM1$learnParams$alpha[indxMixToSep,drop=F])
abline(a=0, b=1, col="red")
all(class_lap@alpha==lapSVM1$learnParams$alpha[indxMixToSep,drop=F])



# distribution regression data y = f(p(x)) + n

nb <- 1000
nb2 <- 300
n <- 20000
x <- rnorm(n*nb)
x2 <- rnorm(n*nb2)
bag <- sample(1:nb, nb*n, replace=T)
bag2 <- sample(1:nb2, nb2*n, replace=T)
mean_x_bag <- aggregate(x, list(bag), mean)[,2]
sd_x_bag <- aggregate(x, list(bag), sd)[,2]
mean_x_bag2 <- aggregate(x2, list(bag2), mean)[,2]
sd_x_bag2 <- aggregate(x2, list(bag2), sd)[,2]
y <- 2*mean_x_bag - 5*sd_x_bag + rnorm(nb,sd=0.01)
y2 <- 2*mean_x_bag2 - 5*sd_x_bag2 + rnorm(nb2,sd=0.01)
mean(y); mean(y2)
y <- y-mean(y)
y2 <- y2-mean(y2)
plot(mean_x_bag, y)
plot(mean_x_bag2, y2)
plot(sd_x_bag, y)
plot(sd_x_bag2, y2)
trainData <- constructData(x=as.matrix(x), bag=bag, y=y, mean_x_bag=mean_x_bag, sd_x_bag=sd_x_bag)
testData <- constructData(x=as.matrix(x2), bag=bag2, y=y2, mean_x_bag=mean_x_bag2, sd_x_bag=sd_x_bag2)
trainData2 <- trainData



# KRR distribution regression - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list()
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,1, length.out=5)), sigma=list(val=NULL, seq=seq(-5,4, length.out=5)))
#nonOptimizableParams <- list(lambda=list(val=0.001), sigma=list(val=1))
#optimizableParams <- list()
#nonOptimizableParams <- list(sigma=list(val=1))
#optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-2,0, 1)))

nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func=sse), 
                                   rmse=list(func=rmse), 
                                   corre=list(func=corre)), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf_dr"
optimizeSet <- "optHP.CV"
krr_dr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr_dr, predict.krr_dr, resids=resids.add)
pm <- proc.time()
krr_dr2 <- setParams(learner=krr_dr2, trainData, mc_cores=1)
proc.time() - pm
# grid of 25 params, 20000 pts per bag, 1000 bags, 4 outer param cores parallel
getHyperPar(krr_dr2, "lambda")
krr_dr2 <- krr_dr2$learn(learner=krr_dr2)
pred.krr_dr2 <- krr_dr2$predict(learner=krr_dr2, data=trainData)
pred.krr_dr2_test <- krr_dr2$predict(learner=krr_dr2, data=testData)
plot(pred.krr_dr2$gy, pred.krr_dr2$gyh); abline(a=0, b=1, col="red")
plot(pred.krr_dr2_test$gy, pred.krr_dr2_test$gyh); abline(a=0, b=1, col="red")

plot(trainData$mean_x_bag, trainData$y)
o <- order(trainData$mean_x_bag)
lines(trainData$mean_x_bag[o], pred.krr_dr2$gyh[o], col="red")
plot(trainData$sd_x_bag, trainData$y)
o <- order(trainData$sd_x_bag)
lines(trainData$sd_x_bag[o], pred.krr_dr2$gyh[o], col="red")


# distribution CLASSIFICATION data y = f(p(x)) + n

nb <- 100
nb2 <- 101
n <- 10
x <- rnorm(n*nb)
x2 <- rnorm(n*nb2)
bag <- sample(1:nb, nb*n, replace=T)
bag2 <- sample(1:nb2, nb2*n, replace=T)
mean_x_bag <- aggregate(x, list(bag), mean)[,2]
sd_x_bag <- aggregate(x, list(bag), sd)[,2]
mean_x_bag2 <- aggregate(x2, list(bag2), mean)[,2]
sd_x_bag2 <- aggregate(x2, list(bag2), sd)[,2]
y <- 2*mean_x_bag - 5*sd_x_bag + rnorm(nb,sd=0.01)
y2 <- 2*mean_x_bag2 - 5*sd_x_bag2 + rnorm(nb2,sd=0.01)
mean(y); mean(y2)
y <- y-mean(y)
y2 <- y2-mean(y2)
y <- (y<mean(y))*1
y2 <- (y2<mean(y2))*1 
table(y)
plot(mean_x_bag, y)
plot(mean_x_bag2, y2)
plot(sd_x_bag, y)
plot(sd_x_bag2, y2)
trainData <- constructData(x=as.matrix(x), bag=bag, y=y, mean_x_bag=mean_x_bag, sd_x_bag=sd_x_bag)
testData <- constructData(x=as.matrix(x2), bag=bag2, y=y2, mean_x_bag=mean_x_bag2, sd_x_bag=sd_x_bag2)
trainData2 <- trainData

fit <- ksvm(x=cbind(mean_x_bag, sd_x_bag), y=factor(y), C=1, nu=0.2, type=c("C-svc","nu-svc")[2], prob.model=T)
(tab <- table(y2, predict(fit, cbind(mean_x_bag2, sd_x_bag2))))
sum(tab[c(2,3)])/sum(tab)
error(fit)
error.ksvm
showMethods(classes="ksvm")
getMethod(error,"ksvm")


dfbags <- combn(unique(bag),2)
dfbags <- cbind(dfbags, matrix(rep(unique(bag), rep(2, length(unique(bag)))), 2, length(unique(bag))))
N <- length(unique(bag))
#count <- 0
Kxs_vals <- apply(dfbags, 2, function(col){
  #print(paste("count: ", count, " out of ", ncol(dfbags)))
  #count <<- count + 1
  indx1 <- which(bag==col[1])
  indx2 <- which(bag==col[2])
  Kx <- kern_rbf(matrix(x[indx1]),matrix(x[indx2]), sigma=1)
  res <- mean(Kx)
  return(res)
})
Kxs <- matrix(NA, N, N)
Kxs[t(dfbags)] <- Kxs_vals
Kxs[t(dfbags[c(2,1),])] <- Kxs_vals
I <- diag(N)
H <- I-matrix(1/N,N,N)
Kxs <- H %*% Kxs %*% H
Kxs <- as.kernelMatrix(Kxs)
fit2 <- ksvm(x=Kxs, y=factor(y), C=1, nu=0.2, type=c("C-svc","nu-svc")[2], prob.model=T)

dfbags <- expand.grid(bagPred=unique(bag2),bagTr=unique(bag))
Npr <- length(unique(bag2))
kxs_vals <- apply(as.matrix(dfbags), 1, function(col){
    # col <- as.matrix(dfbags)[1,]
    #print(paste("count: ", count))
    #count <<- count + 1
    indx1 <- which(bag==col["bagTr"])
    indx2 <- which(bag2==col["bagPred"])
    kx <- kern_rbf(matrix(x2[indx2]), matrix(x[indx1]), sigma=1)
    res <- mean(kx)
    return(res)
  })
kxs <- matrix(NA, Npr, N)
kxs[as.matrix(dfbags)] <- kxs_vals
Ipr <- diag(Npr)
Hpr <- Ipr-matrix(1/Npr,Npr,Npr)
kxs <- Hpr %*% kxs %*% H
kxs <- as.kernelMatrix(kxs)
table(y2, predict(fit2, kxs))



# SVM distribution CLASSIFICATION - rbf kernel, cross validate lambda, nu and sigma
nonOptimizableParams <- list(nu=list(val=0.2))
optimizableParams <- list( sigma=list(val=NULL, seq=10^seq(-4,0,1)))#nu=list(val=NULL, seq=c(0.2, 0.3, 0.4)),

#nonOptimizableParams <- list(nu=list(val=0.2), sigma=list(val=1))
#optimizableParams <- list()
#nonOptimizableParams <- list(sigma=list(val=1))
#optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-2,0, 1)))

nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(acc=list(func="acc"),negCE3=list(func="negCE3")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf_dr"
optimizeSet <- "optHP.CV"
svm_dr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.svm_dr, predict.svm_dr, resids=resids.add)
pm <- proc.time()
svm_dr2 <- setParams(learner=svm_dr2, trainData, mc_cores=1)
proc.time() - pm
svm_dr2$hyperParams$data$grid[,,"acc"]
# 26 mins grid of 3*5=15 params, 20 pts per bag, 100 bags, 1 outer param cores parallel
# 14 mins grid of 3*5=15 params, 10 pts per bag, 100 bags, 1 outer param cores parallel
#  mins grid of 3*5=15 params, 10 pts per bag, 100 bags, 4 outer param cores parallel
getHyperPar(svm_dr2, "sigma")
getHyperPar(svm_dr2, "nu")
svm_dr2 <- svm_dr2$learn(learner=svm_dr2)
svm_dr2$learnParams$model
pred.svm_dr2 <- svm_dr2$predict(learner=svm_dr2, data=trainData)
pred.svm_dr2_test <- svm_dr2$predict(learner=svm_dr2, data=testData)

negCE3(learner=svm_dr2, pred=pred.svm_dr2)
negCE3(learner=svm_dr2, pred=pred.svm_dr2_test)
(tabTr <- table(pred.svm_dr2$gy_class, pred.svm_dr2$gyh_class))
(tabPr <- table(pred.svm_dr2_test$gy_class, pred.svm_dr2_test$gyh_class))

svm_dr2$learnParams$model
sum(tabTr[c(1,4)])/sum(tabTr)
Accuracy(as.numeric(pred.svm_dr2$gyh_class), pred.svm_dr2$gy_class)
1-acc(svm_dr2, pred.svm_dr2)
Accuracy(as.numeric(pred.svm_dr2_test$gyh_class), pred.svm_dr2_test$gy_class)
1-acc(svm_dr2, pred.svm_dr2_test)
svm_dr2$hyperParams$data$grid[,,"acc"]

# SVM distribution CLASSIFICATION - rbf kernel, cross validate lambda, nu and sigma
# eucledean distances calculated between all points first



trainData$M <- calcBagKernTr(trainData)
testData$M <- calcBagKernPr(trainData, data=testData)

nonOptimizableParams <- list(nu=list(val=0.2))
optimizableParams <- list( sigma=list(val=NULL, seq=10^seq(-4,0,1)))#nu=list(val=NULL, seq=c(0.2, 0.3, 0.4)),

#nonOptimizableParams <- list(nu=list(val=0.2), sigma=list(val=1))
#optimizableParams <- list()
#nonOptimizableParams <- list(sigma=list(val=1))
#optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-2,0, 1)))

nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(acc=list(func="acc"),negCE3=list(func="negCE3")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf_dr"
optimizeSet <- "optHP.CV"
svm_dr_M <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.svm_dr_M, predict.svm_dr_M, resids=resids.add)
pm <- proc.time()
svm_dr_M <- setParams(learner=svm_dr_M, trainData, mc_cores=1)
proc.time() - pm
# grid of 3^3=27 params, 100 pts per bag, 20 bags, 1 outer param cores parallel
getHyperPar(svm_dr_M, "sigma")
getHyperPar(svm_dr_M, "nu")
svm_dr_M <- svm_dr_M$learn(learner=svm_dr_M)
svm_dr_M$learnParams$model
pred.svm_dr_M <- svm_dr_M$predict(learner=svm_dr_M, data=trainData)
pred.svm_dr_M_test <- svm_dr_M$predict(learner=svm_dr_M, data=testData)

negCE3(learner=svm_dr_M, pred=pred.svm_dr_M)
negCE3(learner=svm_dr_M, pred=pred.svm_dr_M_test)
(tabTr <- table(pred.svm_dr_M$gy_class, pred.svm_dr_M$gyh_class))
(tabPr <- table(pred.svm_dr_M_test$gy_class, pred.svm_dr_M_test$gyh_class))

svm_dr_M$learnParams$model
sum(tabTr[c(1,4)])/sum(tabTr)
Accuracy(as.numeric(pred.svm_dr_M$gyh_class), pred.svm_dr_M$gy_class)
1-acc(svm_dr_M, pred.svm_dr_M)
Accuracy(as.numeric(pred.svm_dr_M_test$gyh_class), pred.svm_dr_M_test$gy_class)
1-acc(svm_dr_M, pred.svm_dr_M_test)
svm_dr_M$hyperParams$data$grid[,,"acc"]


