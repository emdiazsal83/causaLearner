# construct learners

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

# Vanilla - cross validate lambda, rbf kernel, median heuristic for sigma 
nonOptimizableParams <- list()
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(), testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
vanilla <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.vanilla, predict.vanilla, resids=resids.add)

# KRR - arbitrary lambda, rbf kernel, median heuristic for sigma 
nonOptimizableParams <- list(lambda=list(val=0.1),sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list()
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(rmse=list(func="rmse")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
krr0 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr_class, resids=resids.add_class)


# KRR - cross validate lambda, rbf kernel, median heuristic for sigma 
nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-6,-1, length.out=6)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(rmse=list(func="rmse")), numFolds=5, testTrain="test", optLossFunc=function(grid) which.min(grid[,"rmse"]))
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
krr1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr_class, resids=resids.add_class)

# weighted KRR - cross validate lambda, rbf kernel, median heuristic for sigma 
ws <- rep(1, 1000)
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

# KRR - cross validate lambda, rbf kernel, one sigma per var, median heuristic for sigma 
nonOptimizableParams <- list(sigmas=list(val=NULL))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-6,-1, length.out=6)))
nonDataParams <- list(kernelXs=list(name="kern_rbf_indSig", pars=c(sigmas="sigmas")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(rmse=list(func="rmse")), numFolds=5, testTrain="test", optLossFunc=function(grid) which.min(grid[,"rmse"]))
heuristicSet <- "getFixedParams_rbf_indSig"
optimizeSet <- "optHP.CV"
krrIndSig <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)


# KRR - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list(beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)), sigma=list(val=NULL, seq=10^seq(-8,8, length.out=17)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), numFolds=5, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
krr2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)

# KRR - rbf kernel, cross validate lambda and sigma
nonOptimizableParams <- list()
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)), sigma=list(val=NULL, length.out=17))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL)
optimizeParams <- list(losses=list(sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), numFolds=5, testTrain="test", mainLoss="rmse")
heuristicSet <- "getFixedParams_rbf2b"
optimizeSet <- "optHP.CV"
krr2b <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.krr, predict.krr, resids=resids.add)


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
                                   corre=list(func="corre")), numFolds=5, maxPoints=100, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
qhsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.qhsic, predict.qhsic, resids=resids.add)

# hsic - cross validate lambda, rbf kernel, max var heuristic for sigma, max var for beta and gamma, max_iterations = 20, initializations=30

nonOptimizableParams <- list(sigma=list(val=NULL, type="med"), beta=list(val=NULL, type="med"), kappa=list(val=NULL, type="med"), gamma=list(val=NULL, type="med"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-8,8, length.out=50)))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma")), kernelXb=list(name="kern_rbf", pars=c(sigma="beta")), 
                      kernelYhk=list(name="kern_rbf", pars=c(sigma="kappa")), kernelRg=list(name="kern_rbf", pars=c(sigma="gamma")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, avgy=NULL)
optimizeParams <- list(losses=list(hsicLoss=list(func="hsicLoss"), 
                                   sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), numFolds=5, testTrain="test", max_iterations1=10, max_iterations2=100, num_init=30)
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
hsic <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.hsic, predict.hsic, resids=resids.add)

# GP.gptk - ML lambda, sigma, n_max = 50 (i.e. use fitc with 50 inducing points)


nonOptimizableParams <- list(modelInit=list(val=NULL))
optimizableParams <- list(model=list(val=NULL))
nonDataParams <- list(kernel=list("rbf", "white"))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list()
optimizeParams <- list(losses=list(sse=list(func="sse"), 
                                   rmse=list(func="rmse"), 
                                   corre=list(func="corre")), approx="fitc", numActive=100, fixInducing=FALSE, iters=2000)
heuristicSet <- "myGpCreate"
optimizeSet <- "myGpOptimise"
gptk1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.gptk, predict.gptk, resids=resids.add)

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

# quantile regression learner - cross validate lambda, fixed sigma val
nonOptimizableParams <- list(taus=list(val=list(c(0.1, 0.25, 0.5, 0.75, 0.9))), sigmaX=list(val=1))
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


# quantile regression learner - cross validate lambda, sigma data-driven seq
nonOptimizableParams <- list(taus=list(val=list(c(0.25, 0.5, 0.75))))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.1, 1,10, 100, 1000)), sigma.rbf.X=list(val=NULL, seq=NULL, length.out=5))
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigma.rbf.X")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alphas=NULL, bs=NULL)
optimizeParams <- list(losses=list(negLogLik=list(func="negLogLik"), 
                                   hsicLoss2=list(func="hsicLoss2"), 
                                   pinball=list(func="pinball")), numFolds=4, testTrain="test")
heuristicSet <- "getFixedParams_rbf2"
optimizeSet <- "optHP.CV"
kqr5 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.kqr, predict.kqr, resids=resids.kqr)

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


# logistic regression with vanilla feature - cross validate nout  
nonOptimizableParams <- list()
optimizableParams <- list()
nonDataParams <- list(featureX=list(name="feat_van", pars=NULL))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=10, testTrain="test")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
logReg1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logReg, predict.logReg, resids=resids.add, makeFeature=makePhi)

# logistic regression with rbf RFF - cross validate lambda internally - phi calculated externally 
nonOptimizableParams <- list(lambda=list(val=10^seq(-6,-1, length.out=6)))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=5, testTrain="test")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
logRegInt1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logRegInt, predict.logRegInt, resids=resids.add)

# logistic regression with lambda=0, for low dim regressions 
nonOptimizableParams <- list()
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")))
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
logRegInt2 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logRegInt0, predict.logRegInt0, resids=resids.add)

# logistic regression with rbf kernel - cross validate lambda internally - K calculated internally 
nonOptimizableParams <- list(lambda=list(val=10^seq(-6,-1, length.out=6)))
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=5, testTrain="test")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
logRegInt3 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logRegInt3, predict.logRegInt3, resids=resids.add)


# vanilla helper for conditional mean embedding based classification
# between fakes and trues
nonOptimizableParams <- list()
optimizableParams <- list()
nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")))
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
cmeClass1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cmeClass, predict.cmeClass)



# logistic regression with rbf RFF - cross validate lambda  
nonOptimizableParams <- list(sigma=list(val=10), num_f=list(val=1000), p_w=list(val="rnorm2"), seed=list(val=1234), map=list(val="cos"))
optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-4,1, length.out=6)))
nonDataParams <- list(featureX=list(name="rff_rbf", pars=c(sigma="sigma", p_w="p_w", seed="seed", map="map", num_f="num_f")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(model=NULL)
optimizeParams <- list(losses=list(negCE=list(func="negCE")), numFolds=10, testTrain="test")
heuristicSet <- "getFixedParams_rbf"
optimizeSet <- "optHP.CV"
logKReg1 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.logKReg, predict.logKReg, resids=resids.add, makeFeature=makePhi)


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

# laplacian KRR for semi super vised learning - W = K - no CV
nonOptimizableParams <- list(Gtype=list(val="adj_euc", pars=list(adjacency_k=7)), normL=list(val=TRUE),
                             kernelX=list(val="kern_rbf", seq=c("kern_rbf")), 
                             lambda1=list(val=54292.83), # laplacian - smoothness
                             lambda2=list(val=0.003377778), # smoothness on ambient space
                             sigma.rbf.X=list(val=1/(2*83.7678^2)))  
optimizableParams <- list()#10

nonDataParams <- list()
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL, K=NULL, W=NULL, L=NULL)
optimizeParams <- list(losses=list(myAUC=list(func="myAUC"), 
                                   rmse=list(func="rmse")), numFolds=5, testTrain="test", mainLoss="rmse")
heuristicSet <- NULL
optimizeSet <- "optHP.CV"
lapkrr_K0 <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSSLkrr1, predict.lapSSLkrr1, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)


# laplacian KRR for semi super vised learning - W = K
nonOptimizableParams <- list(Gtype=list(val="K", pars=list()), normL=list(val=TRUE)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
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
lapkrr_K <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSSLkrr1, predict.lapSSLkrr1, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)

# laplacian KRR for semi super vised learning - Adjacency matrix based on K for Laplacian
nonOptimizableParams <- list(Gtype=list(val="adj", pars=list(adjacency_k=6)), normL=list(val=TRUE)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
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
lapkrr_adj <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSSLkrr1, predict.lapSSLkrr1, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)

# laplacian KRR for semi super vised learning - Adjacency matrix based on isomap for Laplacian
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
lapkrr_iso <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSSLkrr1, predict.lapSSLkrr1, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)


# laplacian SVM for semi super vised learning - W=K
nonOptimizableParams <- list(Gtype=list(val="K", pars=list()), normL=list(val=TRUE), eps=list(val=1e-08)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
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
lapSVM_K <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSVM, predict.lapSVM, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)

# laplacian SVM for semi super vised learning - W = adjacency matrix based onK
nonOptimizableParams <- list(Gtype=list(val="adj", pars=list(adjacency_k=6)), normL=list(val=TRUE), eps=list(val=1e-08)) # Gtype = c("K","adj","isomap") if adj/isompa pars=list(adjacency_k) 
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
lapSVM_adj <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSVM, predict.lapSVM, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)

# laplacian SVM for semi super vised learning - W = adjacency matrix based onK
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
lapSVM_iso <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lapSVM, predict.lapSVM, optimizeSetPars=list(predFunc="pred.SSL.CV"), makeKernel=makeKernel, resids=resids.add)

# latent noise krr - lnkrr
nonOptimizableParams <- list(sigmaX=list(val=NULL),
                             L=list(val=1e6),
                             sigLevel=list(val=0.2),
                             gammas1N=list(val=10),
                             gammas2N=list(val=10),
                             logGammaRange=list(val=c(-14,0)),
                             ps=list(val=c(0.5)),
                             indxP=list(val=1),
                             krrRefModel=list(val=NULL,type="krr1"),
                             q=list(val=1))
optimizableParams <- list(lambda=list(val=NULL, seq=c(0.0001, 0.001, 0.01)),
                          eta=list(val=NULL, seq=c(0.0001, 0.001, 0.01)))
optLossFuncRmseHsicRX <- function(grid){ 
  indx <- which.max(grid[,"ksTestHsicRX"] > sigLevel)
  indxOpt <- which(grid[,"rmse2"]==min(grid[indx,"rmse2"]))
  return(indxOpt)
}
optLossFuncHsicRX <- function(grid) which.max(grid[,"ksTestHsicRX"])
optLossFuncHsicRXHsicResX <- function(grid) which.max(grid[,"ksTestHsicRX"]*grid[,"ksTestHsicResX"])
nonDataParams <- list(kernelXs=list(name="kern_rbf", pars=c(sigma="sigmaX")))
dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
hyperParams <- list(data=dataParams, non_data=nonDataParams)
learnParams <- list(alpha=NULL,gamma=NULL, rhat=NULL, sigma_r=NULL, 
                    ksPvalHsicRX=NULL, ksPvalHsicResX=NULL, ksPvalHsicResR=NULL, 
                    ksPvalHsicResZ=NULL, rmseOpt=NULL)
optimizeParams <- list(losses=list(rmse2=list(func="rmse2", aggFunc="mean"), 
                                   ksTestHsicRX=list(func="ksTestHsicRX", aggFunc="mean"),
                                   ksTestHsicResX=list(func="ksTestHsicResX", aggFunc="mean"),
                                   ksTestHsicResR=list(func="ksTestHsicResR", aggFunc="mean"),
                                   ksTestHsicResZ=list(func="ksTestHsicResZ", aggFunc="mean")), 
                       numFolds=1, numFoldsInt=10, mc_cores=3, numCoresFold=15, initOnly=FALSE, testTrain="test", seed=123, optLossFunc=optLossFuncHsicRXHsicResX)
heuristicSet <- "getFixedParams_rbf_ln"
optimizeSet <- "optHP.CV"
lnkrr <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.lnkrr, predict.lnkrr, resids=resids.add, optimizeSetPars=list(CVfunc="CV.boot"))