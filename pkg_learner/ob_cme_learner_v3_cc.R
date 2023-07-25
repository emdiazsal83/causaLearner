# Construct Conditional Mean Embedding Measures

# remove(list=ls())
# repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
# setwd(repos)
# source("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_learner/func_learners_v3.R")
# source("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_causaLearner/utilities/func_dagStuff.R")

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameter=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                       paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

paramTab2 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameter=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

aux <- expand.grid(kernelX=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_polyTnonlin"), 
                   kernelY=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_polyTnonlin"), stringsAsFactors=F)
aux$name <- with(aux, paste(kernelX, kernelY , sep="AND"))

kernelTab1 <- expand.grid(name=aux$name, centerLx=c(TRUE), centerKy=c(TRUE), optParms=c("none", "lambda", "lambda_kernParsX", "lambda_kernParsXY"), stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$kernelY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")
kernelTab1$learnerName <- paste("cmem", kernelTab1$name, centerCombo, "L2", kernelTab1$optParms, sep="_")

for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 1
  # lambda
  # i <- 10
  # lambda_kernPars
  # i <- 19
  #print("**********************")
  #print(kernelTab1[i,])
  kernelX <- kernelTab1$kernelX[i]
  kernelY <- kernelTab1$kernelY[i]
  
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
  
  
  lossFun <- cmem_L2
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]

  if(kernelTab1$optParms[i]=="none"){
    
    parsKernX <- switch(kernelX_short, van=NA, rbf=1e4, quad=1e-3, log=0.98, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1e-3)
    parsKernY <- switch(kernelY_short, van=NA, rbf=1e4, quad=1e-3, log=0.98, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1e-3)
    #parsKernY <- switch(kernelY_short, van=NA, rbf=1, quad=1, log=0.001, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1, laplace2=1)
    
    
    nms <- c("lambda", "centerLx", "centerKy", parsX, parsY)
    nonOptimizableParams <- list(list(val=1e-4), list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j])))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernY[j])))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list()
    
    
  } else if(kernelTab1$optParms[i]=="lambda"){
    nms <- c("centerLx", "centerKy",parsX, parsY)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-10,-6, length.out=9)))
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsX"){
    
    nms <- c("centerLx", "centerKy",parsY)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda", parsX)
    optimizableParams <- list(list(val=NULL, seq=10^seq(-10,-6, length.out=9)))
    for(j in 1:length(parsX)) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=10)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsXY"){
    
    nms <- c("centerLx", "centerKy")
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda",parsX, parsY)
    optimizableParams <- list(list(val=NULL, seq=10^seq(-10,-10, length.out=5)))
    for(j in 1:(length(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=4)))
    for(j in 1:(length(parsY))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=4)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
    
  }
  
  
  nonDataParams <- list(kernelX=list(name=kernelX, pars=na.omit(parsX)), 
                        kernelY=list(name=kernelY, pars=na.omit(parsY)))
  
  
  heuristicSet <- na.omit(unique(c(paramTab1$paramFunc[c(indxX)], paramTab2$paramFunc[c(indxY)])))
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(Ky=NULL, Lx=NULL, Blambda=NULL, Cks=NULL)
  optimizeParams <- list(losses=list(cmem_L2=lossFun), numFolds=5, testTrain="test")
  optimizeSet <- "optHP.CV"
  measureParams <- list(KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                        KCDCpval=list(func="KCDCpval", aggFunc="prod", pars=list(numPerms=2)),
                        KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                        KCMC=list(func="KCMC", aggFunc="sum",pars=list()),
                        KCSC=list(func="KCSC", aggFunc="sum",pars=list()),
                        KCNSC=list(func="KCNSC", aggFunc="sum",pars=list()))
  
  cmem_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cmem_L2, predict.cmem_L2, msrs=measureParams, calcMsrs=calcMsrs)
  assign(kernelTab1$learnerName[i], cmem_learner_aux)
}

cmem_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmem_learner_pack_lambda_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda")]
cmem_learner_pack_kernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX")]
cmem_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]


# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

ppTab <- data.frame(id=c("KCDC","KCDCrel","KCDCpval","KCRDC","KCMC", "KCMCpval","KCSC", "KCSCpval", "KCNSC", "majority", "rankMajority","rf","lasso"))




# define rank2Func transf: how t
# KCDC              =                     in (0,Inf)   -> addQuot
# KCDCrel           =                     in (0,Inf)   -> addQuot
# KCDCpval          =                     in (0,Inf)   -> addQuot
# KCRDC             =                     in (0,Inf)   -> addQuot
# KCMC              =                     in (0,Inf)   -> addQuot
# KCMCpval          =                     in (0,Inf)   -> addQuot
# KCSC              =                     in (0,Inf)   -> addQuot
# KCSCpval          =                     in (0,Inf)   -> addQuot
# KCNSC             =                     in (0,Inf)   -> addQuot
# score_rank_majority =                   in Naturals   -> differ
# score_majority      =                   in Naturals   -> differ       
# rf      =                                -> (0,1)    
# lasso      =                                -> (0,1)    

ppTab$rank2Funcs <- c("quot","addQuot","differ")[rep(3,13)]
ppTab$rankFuncs <- rep("correctScoreToAdd", 13)
ppTab$probFuncs <- rep("scoreToProb", 13)
ppTab$argTypes <- c(rep("cmes", 13))



ppTabCMEM1 <- ppTab