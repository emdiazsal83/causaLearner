# Construct Conditional Mean Feature Measures

# remove(list=ls())
# repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
# setwd(repos)
# source("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_learner/func_learners_v3.R")

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameters=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                       paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

paramTab2 <- data.frame(feature=c("rbf", "quad", "laplace1", "poly"), parameters=c("sigma", "offset", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_laplace1", "getFixedParams_poly"), stringsAsFactors=F)

#"rff_log"
aux <- expand.grid(kernelX=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_poly"), 
                   featureY=c("rff_rbf", "rff_quad","rff_laplace1","feat_toy"), stringsAsFactors=F) #, "feat_polyTnonlin"
aux$name <- with(aux, paste(kernelX, featureY , sep="AND"))

kernelTab1 <- expand.grid(name=aux$name, centerLx=c(TRUE), centerKy=c(TRUE), optParms=c("none", "lambda", "lambda_kernParsX", "lambda_kernParsXY"),  stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$featureY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")
kernelTab1$learnerName <- paste("cmfm", kernelTab1$name, centerCombo ,"L2", kernelTab1$optParms, sep="_")

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
  
  featureY <- kernelTab1$featureY[i]
  aux <- strsplit(featureY, "_")[[1]]
  featureY_short <- aux[length(aux)]
  featureY_short_short <- strsplit(featureY_short, "T")[[1]][1]
  
  indxX <- match(kernelX_short_short, paramTab1$kernel)
  indxY <- match(featureY_short_short, paramTab2$feature)
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
    parsY <- c(paste(strsplit(parsY, "-")[[1]], featureY_short_short, "Y", sep=".")) #
    names(parsY) <- strsplit(paramTab2$parameter[indxY], "-")[[1]]
  } 
  
  aux <- strsplit(featureY, "_")[[1]]
  if(aux[1]=="rff"){
    aux <- c("p_w","num_f", "seed", "map") #
    names(aux) <- aux
    parsYfeat <- c(aux)
  } else{
    parsYfeat <- c(NA, NA, NA, NA)
  }
  lossFun <- TNRE #cmem_L2
  num_f <- 100
  map <- "cos"
  
  dist_p_w <- switch(featureY_short_short, rbf="rnorm2", quad="rlaplace2", laplace1="rcauchy2", NA)
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]
  
  if(kernelTab1$optParms[i]=="none"){
    
    parsKernX <- switch(kernelX_short, van=NA, rbf=1e4, quad=1e-3, log=0.98, poly=c(3,1,0) , laplace1=1e-3, laplace2=1e-3)
      parsFeatY <- switch(featureY_short, van=NA, rbf=1e4, quad=1e-3, log=0.98, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1e-3)
    #parsFeatY <- switch(kernelY_short, van=NA, rbf=1, quad=1, log=0.001, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1, laplace2=1)
    
    nms <- c("lambda", "centerLx", "centerKy", parsX, parsY, parsYfeat)
    nonOptimizableParams <- list(list(val=1e-8), list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j])))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsFeatY[j])))
    nonOptimizableParams <- c(nonOptimizableParams, list(list(val=dist_p_w)),list(list(val=num_f)), list(list(val=1234)), list(list(val=map)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list()
    
  } else if(kernelTab1$optParms[i]=="lambda"){
    nms <- c("centerLx", "centerKy",parsX, parsY, parsYfeat)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    nonOptimizableParams <- c(nonOptimizableParams, list(list(val=dist_p_w)), list(list(val=num_f)), list(list(val=1234)), list(list(val=map)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list(lambda=list(val=NULL, seq=10^seq(-10,-6, length.out=9)))
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsX"){
    nms <- c("centerLx", "centerKy",parsY, parsYfeat)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy))
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL)))
    nonOptimizableParams <- c(nonOptimizableParams, list(list(val=dist_p_w)), list(list(val=num_f)), list(list(val=1234)), list(list(val=map)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda", parsX)
    optimizableParams <- list(list(val=NULL, seq=10^seq(-10,-6, length.out=9)))
    for(j in 1:length(parsX)) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=10)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsXY"){
    nms <- c("centerLx", "centerKy", parsYfeat)
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy),list(val=dist_p_w), list(val=num_f), list(val=1234), list(val=map))
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
                        featureY=list(name=featureY, pars=na.omit(c(parsY, parsYfeat))))
  
  heuristicSet <- na.omit(unique(c(paramTab1$paramFunc[c(indxX)], paramTab2$paramFunc[c(indxY)])))
  
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(phiy=NULL, Ky=NULL, Lx=NULL, Blambda=NULL, Cks=NULL)
  optimizeParams <- list(losses=list(cmem_L2=lossFun), numFolds=5, testTrain="test")
  optimizeSet <- "optHP.CV"
  measureParams <- list(KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                        KCDCpval=list(func="KCDCpval", aggFunc="prod", pars=list(numPerms=2)),
                        KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                        KCMC=list(func="KCMC", aggFunc="sum",pars=list()),
                        KCSC=list(func="KCSC", aggFunc="sum",pars=list()),
                        KCCC_ent=list(func="KCCC_ent", aggFunc="sum",pars=list()),
                        KCCC_comp=list(func="KCCC_comp", aggFunc="sum",pars=list()),
                        KCNSC=list(func="KCNSC", aggFunc="sum",pars=list()),
                        MCX=list(func="MCX", aggFunc="sum",pars=list()),
                        EDML=list(func="EDML", aggFunc="sum",pars=list()),
                        EDMB=list(func="EDMB", aggFunc="sum",pars=list()),
                        TRE=list(func="TRE", aggFunc="sum",pars=list()),
                        CVTeRNE=list(func="CVTeRNE", aggFunc="sum",pars=list()),
                        CVTrRNE=list(func="CVTrRNE", aggFunc="sum",pars=list()),
                        EAL=list(func="EAL", aggFunc="sum",pars=list()),
                        DAL=list(func="DAL", aggFunc="sum",pars=list()))

  cmfm_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cmfm_L2, predict.cmfm_L2, msrs=measureParams, calcMsrs=calcMsrs)
  assign(kernelTab1$learnerName[i], cmfm_learner_aux)
}

cmfm_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmfm_learner_pack_lambda_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda")]
cmfm_learner_pack_kernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX")]
cmfm_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]


# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

# do not add measures with underscore as this gets ignored... underscore are to identify variants of main measure eg. KCCC_ent vs KCCC_comp
ppTab <- data.frame(id=c("KCDC","KCDCrel","KCDCpval","KCRDC","KCMC", "KCMCpval","KCSC", "KCSCpval","KCCC", "KCCCpval", "KCNSC",
                         "MCX","EDML","EDMB","TRE","CVTeRNE", "CVTrRNE", "EAL", "DAL"))




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

ppTab$rank2Funcs <- c("quot","addQuot","differ")[rep(3,19)]
ppTab$rankFuncs <- rep("correctScoreToAdd", 19)
ppTab$probFuncs <- rep("scoreToProb", 19)
ppTab$argTypes <- c(rep("cmfs", 19))


ppTabCMFM1 <- ppTab

