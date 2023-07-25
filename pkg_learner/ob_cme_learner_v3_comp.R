# Construct Conditional Mean Feature Measures

# remove(list=ls())
# repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
# setwd(repos)
# source("./pkg_learner/func_learners_v4.R")


# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("bin", "rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameters=c("num","sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                       paramFunc=c("getNumBins","getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

paramTab2 <- data.frame(kernel=c("rbf", "quad", "laplace1", "poly"), parameters=c("sigma", "offset", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_laplace1", "getFixedParams_poly"), stringsAsFactors=F)

#"rff_log"
#aux <- expand.grid(kernelX=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_poly"), 
#                   "kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_polyTnonlin"), stringsAsFactors=F) #, "feat_polyTnonlin"
aux <- expand.grid(kernelX=c("kern_bin","kern_rbf"), 
                   kernelY=c("kern_rbf"), stringsAsFactors=F) #, "feat_polyTnonlin"

aux$name <- with(aux, paste(kernelX, kernelY , sep="AND"))

betaLearns <- c("learnBlambda_L2", "learnBlambda_bin", "learnBlambda_KCMC", "learnBlambda_KCSC") 


kernelTab1 <- expand.grid(name=aux$name, betaLearns=betaLearns, centerLx=c(FALSE), centerKy=c(FALSE), optParms=c("none", "lambda", "lambda_kernParsX", "lambda_kernParsXY"),  stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$kernelY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")
kernelTab1$learnerName <- paste("cmem", kernelTab1$name, sapply(strsplit(kernelTab1$betaLearns,"_"), function(el) el[[2]]),centerCombo , kernelTab1$optParms, sep="_")

indx <- which(kernelTab1$kernelX == "kern_bin" & kernelTab1$betaLearns != "learnBlambda_bin")
kernelTab1 <- kernelTab1[-indx,]
indx <- which(kernelTab1$kernelX != "kern_bin" & kernelTab1$betaLearns == "learnBlambda_bin")
kernelTab1 <- kernelTab1[-indx,]
indx <- which(kernelTab1$kernelX != "kern_bin" & kernelTab1$optParms == "lambda")
kernelTab1 <- kernelTab1[-indx,]

dim(kernelTab1)


for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 1
  # lambda
  # i <- 120*1+7
  # lambda_kernPars
  # i <- 120*2+6
  #print("**********************")
  #print(kernelTab1[i,])
  kernelX <- kernelTab1$kernelX[i]
  aux <- strsplit(kernelX, "_")[[1]]
  # short name removes kern or feat part
  kernelX_short <- aux[length(aux)] 
  # short-short name removes part which specifies the name of parametrization ex poly.lin or poly.nonlin
  kernelX_short_short <- strsplit(kernelX_short, "T")[[1]][1]
  
  kernelY <- kernelTab1$kernelY[i]
  aux <- strsplit(kernelY, "_")[[1]]
  kernelY_short <- aux[length(aux)]
  kernelY_short_short <- strsplit(kernelY_short, "T")[[1]][1]
  
  indxX <- match(kernelX_short_short, paramTab1$kernel)
  indxY <- match(kernelY_short_short, paramTab2$kernel)
  # these are parameters needed to compute kernel: to feed to kernelMatrix or to phiy 
  # (not other such as center=T,or lambda,  which are not fed to kernelMatrix)
  parsX <- paramTab1$parameter[indxX]
  if(!is.na(indxX)){
    # if kernel or feature does have parameters then we give them the "last-name" of kernel to avoid ambiguity 
    # between different kernels with the same parameter name
    parsX <- paste( strsplit(parsX, "-")[[1]], kernelX_short_short, "X", sep=".") 
    names(parsX) <- strsplit(paramTab1$parameter[indxX], "-")[[1]]
  }
  parsY <- paramTab2$parameter[indxY]
  if(!is.na(indxY)){
    parsY <- c(paste(strsplit(parsY, "-")[[1]], kernelY_short_short, "Y", sep=".")) #
    names(parsY) <- strsplit(paramTab2$parameter[indxY], "-")[[1]]
  } 
  
  
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]
  
  lambdas_other <- list(val=NULL, seq=10^seq(-9,1, length.out=3))
  lambdasXY_other <- list(val=NULL, seq=10^seq(-9,1, length.out=3))
  lambdas_bin <- list(val=0, seq=0)
  numX <- 3#10
  numX_XY <- 3#10
  numY_XY <- 3#10
  
  lambdas <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=lambdas_bin, other=lambdas_other)
  lambdasXY <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=lambdas_bin, other=lambdasXY_other)
  
  # logistic regression comparing phiy(y_true) with 
  # phiy_hat(x_true) and phiy_hat(x_fake)
  classifier <- "logRegInt1"
  featFunc <- "makeLogRegFeats2"
  fakeDist_x <- "runif"
  fakeDist_y <- "runif"
  kappa <- 1
  fakeDistPars_x <- list(min=-0.5, max=0.5)
  fakeDistPars_y <- list(min=-0.5, max=0.5)
  NCE_learner <- list(classifier=classifier, featFunc=featFunc, fakeDist_x=fakeDist_x, fakeDistPars_x=fakeDistPars_x, fakeDist_y=fakeDist_y, fakeDistPars_y=fakeDistPars_y, kappa=kappa)
  
  # cme embedding classifier comparing phiy(y_true) with 
  # phiy_hat(x_true) and phiy_hat(x_fake)
  
  # classifier <- "cmeClass1"
  # featFunc <- "makeCME_cond_feats"
  # fakeDist_x <- "runif"
  # fakeDist_y <- "runif"
  # kappa <- 1
  # fakeDistPars_x <- list(min=-0.5, max=0.5)
  # fakeDistPars_y <- list(min=-0.5, max=0.5)
  # NCE_learner <- list(classifier=classifier, featFunc=featFunc, fakeDist_x=fakeDist_x, fakeDistPars_x=fakeDistPars_x, fakeDist_y=fakeDist_y, fakeDistPars_y=fakeDistPars_y, kappa=kappa)
  
  # cme embedding classifier comparing phiy(y_true) with 
  # phiy_hat(x_true) and phiy_hat(x_fake) but using logistic regression
  
  # classifier <- "logRegInt2"
  # featFunc <- "makeCME_cond_logRegfeats"
  # fakeDist_x <- "runif"
  # fakeDist_y <- "runif"
  # kappa <- 1
  # fakeDistPars_x <- list(min=-0.5, max=0.5)
  # fakeDistPars_y <- list(min=-0.5, max=0.5)
  # NCE_learner <- list(classifier=classifier, featFunc=featFunc, fakeDist_x=fakeDist_x, fakeDistPars_x=fakeDistPars_x, fakeDist_y=fakeDist_y, fakeDistPars_y=fakeDistPars_y, kappa=kappa)
  
  parsKernX <- switch(kernelX_short, bin=3, van=NA, rbf=1e2, quad=1e-3, log=0.98, poly=c(3,1,0) , laplace1=1e-3, laplace2=1e-3)
  parsKernY <- switch(kernelY_short, van=NA, rbf=1e2, quad=1e-3, log=0.98, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1e-3)
  
  
  if(kernelTab1$optParms[i]=="none"){
    
    
    nms <- c("lambda", "centerLx", "centerKy", "NCE_learner", parsX, parsY)
    nonOptimizableParams <- list(list(val=1e-4), list(val=centerLx), list(val=centerKy), list(val=NCE_learner))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j], seq=NULL)))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernY[j], seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list()
    
  } else if(kernelTab1$optParms[i]=="lambda"){
    nms <- c("centerLx", "centerKy", "NCE_learner",parsX, parsY)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy), list(val=NCE_learner))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j], seq=NULL)))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernY[j], seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list(lambda=lambdas)
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsX"){
    nms <- c("centerLx", "centerKy", "NCE_learner",parsY)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy), list(val=NCE_learner))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernY[j], seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda", parsX)
    optimizableParams <- list(lambdas)
    for(j in 1:length(parsX)) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsXY"){
    nms <- c("centerLx", "centerKy", "NCE_learner")
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy), list(val=NCE_learner))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda",parsX, parsY)
    optimizableParams <- list(lambdasXY)
    for(j in 1:(length(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX_XY)))
    for(j in 1:(length(parsY))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numY_XY)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  }
  
  nonDataParams <- list(kernelX=list(name=kernelX, pars=na.omit(parsX)), 
                        kernelY=list(name=kernelY, pars=na.omit(parsY)))
  
  heuristicSet <- na.omit(unique(c(paramTab1$paramFunc[c(indxX)], paramTab2$paramFunc[c(indxY)])))
  
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(Ky=NULL, Lx=NULL, Blambda=NULL, Cks=NULL, learnBlambda=kernelTab1$betaLearns[i])
  optimizeParams_bin <- list(losses=list(negCE=negCE, gauss_log_lik=gauss_log_lik, MisCR=MisCR, cmem_L2_f=cmem_L2_f, cmem_L2_k=cmem_L2_k, KCDC=KCDC, KCMC=KCMC, KCRDC=KCRDC), numFolds=5, testTrain="test", mainLoss="negCE")
  optimizeParams_other <- list(losses=list(negCE=negCE, gauss_log_lik=gauss_log_lik, MisCR=MisCR, cmem_L2_f=cmem_L2_f, cmem_L2_k=cmem_L2_k, KCDC=KCDC, KCMC=KCMC, KCRDC=KCRDC, KCSC=KCSC, KCNSC=KCNSC), numFolds=5, testTrain="test", mainLoss="negCE")
  
  optimizeParams <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=optimizeParams_bin, other=optimizeParams_other)
  optimizeSet <- "optHP.CV"
  measureParamsFull <- list(gll_tr=list(func="gll_tr", aggFunc="sum", pars=list()),
                            gll_te=list(func="gll_te", aggFunc="sum", pars=list()),
                            L2_f_tr=list(func="L2_f_tr", aggFunc="sum", pars=list()),
                            L2_f_te=list(func="L2_f_te", aggFunc="sum", pars=list()),
                            PCEtr=list(func="PCEtr", aggFunc="sum", pars=list()),
                            PCEte=list(func="PCEte", aggFunc="sum", pars=list()),
                            CCRtr=list(func="CCRtr", aggFunc="sum", pars=list()),
                            CCRte=list(func="CCRte", aggFunc="sum", pars=list()),
                            KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                            KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                            KCMC=list(func="KCMC", aggFunc="sum",pars=list()),
                            KCSC=list(func="KCSC", aggFunc="sum",pars=list()),
                            KCNSC=list(func="KCNSC", aggFunc="sum",pars=list()))

  measureParamsRed  <- list(gll_tr=list(func="gll_tr", aggFunc="sum", pars=list()),
                            gll_te=list(func="gll_te", aggFunc="sum", pars=list()),
                            L2_f_tr=list(func="L2_f_tr", aggFunc="sum", pars=list()),
                            L2_f_te=list(func="L2_f_te", aggFunc="sum", pars=list()),
                            PCEtr=list(func="PCEtr", aggFunc="sum", pars=list()),
                            PCEte=list(func="PCEte", aggFunc="sum", pars=list()),
                            CCRtr=list(func="CCRtr", aggFunc="sum", pars=list()),
                            CCRte=list(func="CCRte", aggFunc="sum", pars=list()),
                            KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                            KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                            KCMC=list(func="KCMC", aggFunc="sum",pars=list()))
  
  measureParams <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=measureParamsRed, other=measureParamsFull)
  
  cmem_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn=learn.cmem, predict=predict.cmem, msrs=measureParams, calcMsrs=calcMsrs,  makeKernel=makeKernel)
  assign(kernelTab1$learnerName[i], cmem_learner_aux)
}

cmem_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmem_learner_pack_lambda_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda")]
cmem_learner_pack_kernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX")]
cmem_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]

cmem_learner_pack_compCombos_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]
cmem_learner_pack_compCombos_1 <- sub("lambda_kernParsXY", "compCombo", cmem_learner_pack_compCombos_1)

# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

# do not add measures with underscore as this gets ignored... underscore are to identify variants of main measure eg. KCCC_ent vs KCCC_comp
ppTab <- data.frame(id=c("KCDC","KCDCrel","KCDCpval","KCRDC","KCMC", "KCMCpval","KCSC", "KCSCpval","KCCC", "KCCCpval", "KCNSC",
                         "MCX","EDML","EDMB","TRE","gll_tr","gll_te","L2_f_tr","L2_f_te","PCEtr","PCEte","CCRtr","CCRte","CVTeRNE", "CVTrRNE", "EAL", "DAL", "rankMajority","majority"))




# define rank2Func transf: how t
# KCDC              =                     in (0,Inf)   -> addQuot
# KCDCrel           =                     in (0,Inf)   -> addQuot
# KCDCpval          =                     in (0,Inf)   -> addQuot
# KCRDC             =                     in (0,Inf)   -> addQuot
# KCMC              =                     in (0,Inf)   -> addQuot
# KCMCpval          =                     in (0,Inf)   -> addQuot
# KCSC              =                     in (0,Inf)   -> addQuot
# KCSCpval          =                     in (0,Inf)   -> addQuot
# KCCC              =                     in (0,Inf)   -> addQuot
# KCCCpval          =                     in (0,Inf)   -> addQuot
# KCNSC             =                     in (0,Inf)   -> addQuot

numMsrs <- 29
ppTab$rank2Funcs <- c("quot","addQuot","differ")[rep(3,numMsrs)]
ppTab$rankFuncs <- rep("correctScoreToAdd", numMsrs)
ppTab$probFuncs <- rep("scoreToProb", numMsrs)
ppTab$argTypes <- c(rep("cmfs", numMsrs))


ppTabCMFM1 <- ppTab

