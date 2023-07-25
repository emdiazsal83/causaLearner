# Construct Conditional Mean Feature Measures

# remove(list=ls())
# repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
# setwd(repos)
# source("./pkg_learner/func_learners_v3.R")

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel


paramTab <- data.frame(feature=c("rbf", "quad", "laplace1", "poly"), parameters=c("sigma", "offset", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_laplace1", "getFixedParams_poly"), stringsAsFactors=F)

#"rff_log"
aux <- expand.grid(featureY=c("rff_rbf"), stringsAsFactors=F) #, "feat_polyTnonlin"
aux$name <- with(aux, featureY)

kernelTab1 <- expand.grid(name=aux$name, centerKy=c(FALSE), optParms=c("none", "kernParsY"),  stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$featureY <- sapply(aux, function(el) el[[1]][1])
centerCombo <- c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))]
kernelTab1$learnerName <- paste("cmfm_bin", kernelTab1$name, centerCombo ,"L2", kernelTab1$optParms, sep="_")

for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 1
  #print("**********************")
  #print(kernelTab1[i,])
  
  featureY <- kernelTab1$featureY[i]
  aux <- strsplit(featureY, "_")[[1]]
  featureY_short <- aux[length(aux)]
  featureY_short_short <- strsplit(featureY_short, "T")[[1]][1]
  
  indxY <- match(featureY_short_short, paramTab$feature)
  # these are parameters needed to compute kernel: to feed to kernelMatrix or to phiy 
  # (not other such as center=T,or lambda,  which are not fed to kernelMatrix)
  
  parsY <- paramTab$parameter[indxY]
  if(!is.na(indxY)){
    parsY <- c(paste(strsplit(parsY, "-")[[1]], featureY_short_short, "Y", sep=".")) #
    names(parsY) <- strsplit(paramTab$parameter[indxY], "-")[[1]]
  } 
  
  aux <- strsplit(featureY, "_")[[1]]
  if(aux[1]=="rff"){
    aux <- c("p_w","num_f", "seed", "map") #
    names(aux) <- aux
    parsYfeat <- c(aux)
  } else{
    parsYfeat <- c(NA, NA, NA, NA)
  }
  lossFun <- negCE
  num_f <- 1000
  map <- "cos"
  
  dist_p_w <- switch(featureY_short_short, rbf="rnorm2", quad="rlaplace2", laplace1="rcauchy2", NA)
  numBins <- 3
  classifier <- "logRegInt1"
  fakeDist <- "runif"
  kappa <- 1
  fakeDistPars <- list()
  featFunc <- makePhi_f
  NCE_learner <- list(classifier=classifier, featFunc=makePhi_f, fakeDist=fakeDist, fakeDistPars=fakeDistPars, kappa=kappa)
  
  centerKy <- kernelTab1$centerKy[i]
  
  if(kernelTab1$optParms[i]=="none"){
    
    parsFeatY <- switch(featureY_short, van=NA, rbf=1e4, quad=1e-3, log=0.98, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1e-3)
    #parsFeatY <- switch(kernelY_short, van=NA, rbf=1, quad=1, log=0.001, polyTlin=c(1,1,0), polyTnonlin=c(3,1,0) , laplace1=1, laplace2=1)
    
    #nms <- c("numBins","centerKy", "classifier", "fakeDist","kappa", parsY, parsYfeat)
    #nonOptimizableParams <- list(list(val=numBins), list(val=centerKy), list(val=classifier), list(val=fakeDist), list(val=kappa))
    
    nms <- c("numBins","centerKy", "NCE_learner", parsY, parsYfeat)
    nonOptimizableParams <- list(list(val=numBins), list(val=centerKy), list(val=NCE_learner))
    
    
    for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsFeatY[j])))
    nonOptimizableParams <- c(nonOptimizableParams, list(list(val=dist_p_w)),list(list(val=num_f)), list(list(val=1234)), list(list(val=map)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list()
    
    
  } else if(kernelTab1$optParms[i]=="kernParsY"){
    nms <- c("numBins","centerKy", "NCE_learner", parsYfeat)
    nonOptimizableParams <- list(list(val=numBins), list(val=centerKy), list(val=NCE_learner), list(val=dist_p_w), list(val=num_f), list(val=1234), list(val=map))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c(parsY)
    optimizableParams <- list()
    for(j in 1:(length(parsY))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=12)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  }
  
  nonDataParams <- list(featureY=list(name=featureY, pars=na.omit(c(parsY, parsYfeat))))
  
  heuristicSet <- na.omit(unique(c(paramTab$paramFunc[c(indxY)])))
  
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(Ky=NULL, phiy=NULL, classifier=NULL, xqs=NULL)
  optimizeParams <- list(losses=list(negCE=negCE, KCDC=KCDC, KCRDC=KCRDC, KCMC=KCMC), numFolds=5, testTrain="test")
  optimizeSet <- "optHP.CV"
  measureParams <- list(KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                        KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                        KCMC=list(func="KCMC", aggFunc="sum",pars=list()))

  cmfm_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cmfm_bin, predict.cmfm_bin, msrs=measureParams, calcMsrs=calcMsrs, makeFeature=makePhi, makeLogRegFeature=makePhi)
  assign(kernelTab1$learnerName[i], cmfm_learner_aux)
  #learner <- cmfm_bin_rff_rbf_n_L2_none
}

cmfm_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmfm_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="kernParsY")]


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

