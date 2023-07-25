# Construct Conditional Mean Feature Measures
print("enters ob_cme_learner_v6g.R")

# remove(list=ls())
# repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
# setwd(repos)
# source("./pkg_learner/func_learners_v5.R")


# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("bin", "rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameters=c("num","sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                        paramFunc=c("getNumBins","getFixedParams_rbf4", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

paramTab2 <- data.frame(kernel=c("bin", "rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameters=c("num","sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                        paramFunc=c("getNumBins","getFixedParams_rbf4", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), stringsAsFactors=F)

# the H separates the diff kernels we want to try out as grid points
# S should separate different kernel we want to try simultaneously as a sum of kernels
aux <- expand.grid(kernelX=c("kern_rbf"), # kern_rbfHkern_quad 
                   kernelY=c("kern_rbf"), stringsAsFactors=F) #, "feat_polyTnonlin", kern_rbfHkern_quad
#aux <- aux[c(1,4),]

aux$name <- with(aux, paste(kernelX, kernelY , sep="AND"))

betaLearns <- c("learnBlambda_L2", "learnBlambda_bin", "learnBlambda_KCMC", "learnBlambda_KCSC")[c(1,3,4)] 
mainLosses <- c("negCE", "gauss_log_lik", "cmem_L2_f", "cmem_L2_f_rel")

kernelTab1 <- expand.grid(name=aux$name, mainLoss=mainLosses, betaLearns=betaLearns, 
                          centerLx=c(FALSE), centerKy=c(TRUE), 
                          optParms=c("none", "lambda", "lambda_kernParsX", "gamma_kernParsX","lambda_kernParsXY","gamma_kernParsXY"),  stringsAsFactors=F)
kernelTab1$feats <- "makeCME_cond_logRegfeatsK_x" #"makeLogRegFeats_x"

indxDel <- which(kernelTab1$optParms %in% c("gamma_kernParsX","gamma_kernParsXY") & !(kernelTab1$betaLearns %in% c("learnBlambda_KCMC", "learnBlambda_KCSC")) )
kernelTab1 <- kernelTab1[-indxDel,]

# for none, lambda, and lambda_kernParsX we only supervise with "cmem_L2_f" loss
#kernelTab1$mainLoss[which(! kernelTab1$optParms %in% c("lambda_kernParsX", "lambda_kernParsXY"))] <- "cmem_L2_f"

#, "makeCME_cond_logRegfeatsOdds"

indx <- which(kernelTab1$mainLoss=="negCE")
NCE_feats <- c("makeCME_cond_logRegfeatsK") #"makeLogRegFeats",  "makeCME_cond_logRegfeats", "makeCME_cond_logRegfeatsNorm", "makeCME_cond_feats"
NCE_feats <- c(paste(NCE_feats, "x", sep="_")) #,paste(NCE_feats, "y", sep="_")

kernelTabAux <- lapply(indx, function(i){
  #i <- 1
  expand.grid(c(as.list(kernelTab1[i,-which(colnames(kernelTab1)=="feats")]), 
                list(feats=NCE_feats)))
})
kernelTabAux <- do.call(rbind, kernelTabAux)
kernelTab1 <- rbind(kernelTab1[-indx,], kernelTabAux)

NCE_classifiers <- c("logRegInt3") #"logRegInt1", "logRegInt2","logRegInt2", "cmeClass1"
NCE_classifiers <- rep(NCE_classifiers, 1)

kernelTab1$classifier <- NCE_classifiers[match(kernelTab1$feats, NCE_feats)]

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$kernelY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")
kernelTab1$learnerName <- paste("cmem", kernelTab1$name, "DEL1", kernelTab1$mainLoss, kernelTab1$feats, kernelTab1$classifier, "DEL2",
                                centerCombo , kernelTab1$optParms, sapply(strsplit(kernelTab1$betaLearns,"_"), function(el) el[[2]]), sep="_")


# when we assigned L2 loss to all none lambda_kernParsXY learners we created duplicates
dups <- which(duplicated(kernelTab1$learnerName))
if(length(dups)>0) kernelTab1 <- kernelTab1[-dups,]

dim(kernelTab1)
table(kernelTab1$betaLearns, kernelTab1$mainLoss, kernelTab1$optParms)

for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 49
  # lambda
  # i <- 120*1+7
  # lambda_kernPars
  # i <- 120*2+6
  #print("**********************")
  #print(i)
  #print(kernelTab1[i,])
  kernelX <- kernelTab1$kernelX[i]
  aux <- strsplit(kernelX, "H")[[1]]
  aux <- lapply(aux, function(el) strsplit(el, "_")[[1]])
  kernelX <- sapply(aux, function(el) paste(el, collapse="_"))
  # short name removes kern or feat part
  kernelX_short <- sapply(aux, function(el) el[length(el)]) 
  # short-short name removes part which specifies the name of parametrization ex polyTlin or polyTnonlin
  kernelX_short_short <- sapply(kernelX_short, function(el) strsplit(el, "T")[[1]][1])
  
  kernelY <- kernelTab1$kernelY[i]
  aux <- strsplit(kernelY, "H")[[1]]
  aux <- lapply(aux, function(el) strsplit(el, "_")[[1]])
  kernelY <- sapply(aux, function(el) paste(el, collapse="_"))
  kernelY_short <- sapply(aux, function(el) el[length(el)])
  kernelY_short_short <- sapply(kernelY_short, function(el) strsplit(el, "T")[[1]][1])
  
  indxX <- match(kernelX_short_short, paramTab1$kernel)
  indxY <- match(kernelY_short_short, paramTab2$kernel)
  # these are parameters needed to compute kernel: to feed to kernelMatrix or to phiy 
  # (not other such as center=T,or lambda,  which are not fed to kernelMatrix)
  parsX <- paramTab1$parameter[indxX]
  parsX <- lapply(1:length(indxX), function(j){
    indx <- indxX[j]
    if(!is.na(indx)){
      # if kernel or feature does have parameters then we give them the "last-name" of kernel to avoid ambiguity 
      # between different kernels with the same parameter name
      parX <- paste( strsplit(parsX[j], "-")[[1]], kernelX_short_short[j], "X", sep=".") 
      names(parX) <- strsplit(paramTab1$parameter[indx], "-")[[1]]
    } else{
      parX <- NA
    }
    return(parX)
  })
  parsY <- paramTab2$parameter[indxY]
  parsY <- lapply(1:length(indxY), function(j){
    indx <- indxY[j]
    if(!is.na(indx)){
      parY <- c(paste(strsplit(parsY[j], "-")[[1]], kernelY_short_short[j], "Y", sep=".")) #
      names(parY) <- strsplit(paramTab2$parameter[indx], "-")[[1]]
    } else{
      parY <- NA
    }
    return(parY)
  }) 
  
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]
  
  gradX <- TRUE #FALSE
  
  # on my PC to try
  #gammas <- list(val=NULL, seq=10^seq(-9,1, length.out=2))
  #gammasBig <- list(val=NULL, seq=10^seq(-9,1, length.out=2))
  #lambdaNet <- list(val=1e-9)
  #lambdas_other <- list(val=NULL, seq=NULL, length.out=2) #seq=10^seq(-6,-1, length.out=2)
  #lambdasXY_other <- list(val=NULL, seq=NULL, length.out=2) # seq=10^seq(-6,-1, length.out=2)
  #lambdas_bin <- list(val=0, seq=0)
  #numX <- 2
  #numX_XY <- 2
  #numY_XY <- 2
  #on ERC serve run with
   gammas <- list(val=NULL, seq=c(0,10^seq(-9,1, length.out=9)))
   gammasBig <- list(val=NULL, seq=c(0,10^seq(-9,5, length.out=49)))
   lambdaNet <- list(val=1e-9)
   lambdas_other <- list(val=NULL, seq=NULL, length.out=15) #seq=10^seq(-6,-1, length.out=5)
   lambdasXY_other <- list(val=NULL, seq=NULL, length.out=5) #seq=10^seq(-6,-1, length.out=5)
   lambdas_bin <- list(val=0, seq=0)
   numX <- 15
   numX_XY <- 5
   numY_XY <- 10
  
  lambdas <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=lambdas_bin, other=lambdas_other)
  lambdasXY <- switch(c("other","bin")[(kernelTab1$betaLearns[i]=="learnBlambda_bin")*1+1], bin=lambdas_bin, other=lambdasXY_other)
  
  # logistic regression comparing phiy(y_true) with 
  # phiy_hat(x_true) and phiy_hat(x_fake)
  classifier <- NCE_classifiers #kernelTab1$classifier[i]
  featFunc <- NCE_feats #kernelTab1$feats[i]
  fakeDist_x <- "runif"
  fakeDist_y <- "runif"
  kappa <- 1
  fakeDistPars_x <- list(min=-0.5, max=0.5)
  fakeDistPars_y <- list(min=-0.5, max=0.5)
  NCE_learner <- list(classifier=classifier, featFunc=featFunc, fakeDist_x=fakeDist_x, fakeDistPars_x=fakeDistPars_x, fakeDist_y=fakeDist_y, fakeDistPars_y=fakeDistPars_y, kappa=kappa)
  densX <- list(val=NULL)
  
  parsKernX <- lapply(kernelX_short, function(el) switch(el, bin=3, van=NA, rbf=NULL, quad=NULL, log=0.98, poly=c(3,1,0) , laplace1=1e-3, laplace2=1e-3))
  parsKernY <- lapply(kernelX_short, function(el) switch(el, bin=3, van=NA, rbf=NULL, quad=NULL, log=0.98, poly=c(3,1,0) , laplace1=1e-3, laplace2=1e-3))
  
  
  if(kernelTab1$optParms[i]=="none"){
    
    nms <- c("lambda", "kernelX", "kernelY","centerLx", "centerKy", "gradX","NCE_learner", "densX",unlist(parsX), unlist(parsY))
    nonOptimizableParams <- list(list(val=0.01), list(val=kernelX[length(kernelX)]), list(val=kernelY[length(kernelY)]), list(val=centerLx), list(val=centerKy), list(val=gradX),list(val=NCE_learner), densX)
    if(kernelTab1$betaLearns[i] %in% c("learnBlambda_KCMC", "learnBlambda_KCSC")){
      nms <- c(nms, "gamma")
      nonOptimizableParams <- c(nonOptimizableParams, list(list(val=1e-4)))
    }
    for(j in 1:length(unlist(parsX))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    for(j in 1:length(unlist(parsY))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list()
    
  } else if(kernelTab1$optParms[i]=="lambda"){
    nms <- c("kernelX", "kernelY","centerLx", "centerKy", "gradX", "NCE_learner","densX", unlist(parsX), unlist(parsY))
    nonOptimizableParams <- list(list(val=kernelX[1]), list(val=kernelY[1]),list(val=centerLx), list(val=centerKy), list(val=gradX), list(val=NCE_learner), densX)
    for(j in 1:length(unlist(parsX))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    for(j in 1:length(unlist(parsY))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    optimizableParams <- list(lambda=lambdas)
    if(kernelTab1$betaLearns[i] %in% c("learnBlambda_KCMC", "learnBlambda_KCSC")){
      nms <- c(nms, "gamma")
      optimizableParams <- c(optimizableParams, list(gammas))
    }
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsX"){
    nms <- c("kernelY","centerLx", "centerKy", "gradX","NCE_learner","densX", unlist(parsY))
    nonOptimizableParams <- list(list(val=kernelY[1]), list(val=centerLx), list(val=centerKy), list(val=gradX), list(val=NCE_learner), densX)
    for(j in 1:length(unlist(parsY))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda", "kernelX")
    optimizableParams <- list(lambdas, list(val=NULL, seq=kernelX))
    if(kernelTab1$betaLearns[i] %in% c("learnBlambda_KCMC", "learnBlambda_KCSC")){
      nms <- c(nms, "gamma")
      optimizableParams <- c(optimizableParams, list(gammas))
    }
    nms <- c(nms, unlist(parsX))
    for(j in 1:length(unlist(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="gamma_kernParsX"){
    nms <- c("lambda","kernelY","centerLx", "centerKy", "gradX","NCE_learner", "densX",unlist(parsY))
    nonOptimizableParams <- list(lambdaNet,list(val=kernelY[1]), list(val=centerLx), list(val=centerKy), list(val=gradX), list(val=NCE_learner), densX)
    for(j in 1:length(unlist(parsY))) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=NULL, seq=NULL)))
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("gamma", "kernelX")
    optimizableParams <- list(gammasBig, list(val=NULL, seq=kernelX))
    nms <- c(nms, unlist(parsX))
    for(j in 1:length(unlist(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="lambda_kernParsXY"){
    nms <- c("centerLx", "centerKy", "gradX","NCE_learner","densX")
    nonOptimizableParams <- list(list(val=centerLx), list(val=centerKy), list(val=gradX),list(val=NCE_learner), densX)
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("lambda","kernelX","kernelY")
    optimizableParams <- list(lambdasXY, list(val=NULL, seq=kernelX), list(val=NULL, seq=kernelY))
    if(kernelTab1$betaLearns[i] %in% c("learnBlambda_KCMC", "learnBlambda_KCSC")){
      nms <- c(nms, "gamma")
      optimizableParams <- c(optimizableParams, list(gammas))
    }
    nms <- c(nms, unlist(parsX), unlist(parsY))
    for(j in 1:(length(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX_XY)))
    for(j in 1:(length(parsY))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numY_XY)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  } else if(kernelTab1$optParms[i]=="gamma_kernParsXY"){
    nms <- c("lambda","centerLx", "centerKy", "gradX","NCE_learner","densX")
    nonOptimizableParams <- list(lambdaNet,list(val=centerLx), list(val=centerKy), list(val=gradX),list(val=NCE_learner), densX)
    nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
    names(nonOptimizableParams) <- na.omit(nms)
    nms <- c("gamma","kernelX","kernelY")
    optimizableParams <- list(gammas, list(val=NULL, seq=kernelX), list(val=NULL, seq=kernelY))
    
    nms <- c(nms, unlist(parsX), unlist(parsY))
    for(j in 1:(length(parsX))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numX_XY)))
    for(j in 1:(length(parsY))) optimizableParams <- c(optimizableParams, list(list(val=NULL, seq=NULL, length.out=numY_XY)))
    optimizableParams <- optimizableParams[which(!is.na(nms))]
    names(optimizableParams) <- na.omit(nms)
    
  }
  
  nonDataParams <- list()
  
  heuristicSet <- na.omit(unique(c(paramTab1$paramFunc[c(indxX)], paramTab2$paramFunc[c(indxY)])))
  
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(phiy=NULL, Ky=NULL, Lx=NULL, Blambda=NULL, Cks=NULL, learnBlambda=kernelTab1$betaLearns[i])
  
  #MisCR=MisCR, , KCRDC=KCRDC, KCSC=KCSC, KCNSC=KCNSC
  losses <- list(negCE=list(func="negCE"),gauss_log_lik=list(func="gauss_log_lik"), cmem_L2_f_rel=list(func="cmem_L2_f_rel"),
                 cmem_L2_f=list(func="cmem_L2_f"), wcmem_L2_f=list(func="wcmem_L2_f"), 
                 cmem_L2_k=list(func="cmem_L2_k"), wcmem_L2_k=list(func="wcmem_L2_k"),
                 corrRKHS=list(func="corrRKHS"), wcorrRKHS=list(func="wcorrRKHS"),
                 KCDC=list(func="KCDC"), KCMC=list(func="KCMC"), 
                 WKCDC=list(func="WKCDC"), WKCMC=list(func="WKCMC"), KIIM=list(func="KIIM"), 
                 KCRDC=list(func="KCRDC"))# , 
                 #KCSC=list(func="KCSC"), KICSC=list(func="KICSC"), KCNSC=list(func="KCNSC"), 
                 #KICCC_pca_ent2b=list(func="KICCC_pca_ent2b"), KICCC_pca_entb=list(func="KICCC_pca_entb"))
  mainLossNm1 <-  kernelTab1$mainLoss[i]
  mainLossNm2 <- mainLossNm1
  if(mainLossNm1 == "negCE"){
    mainLossNm2 <- paste(mainLossNm1, paste(kernelTab1$feats[i], kernelTab1$classifier[i], sep="-") ,sep=".")
  }
  if(! mainLossNm1 %in% names(losses)){
    lossesAux <- list(list(func=mainLossNm1))
    names(lossesAux) <- mainLossNm1
    losses <- c(lossesAux, losses)
  }
  optLossFunc <- paste("function(grid) which.min(grid[,'",mainLossNm2,"'])", sep="")
  optLossFunc <- eval(parse(text=optLossFunc))
  optimizeParams <- list(losses=losses, numFolds=5, testTrain="test", optLossFunc=optLossFunc, mainLoss=mainLossNm2)
  optimizeSet <- "optHP.CV"
  measureParamsFull <- list(#gll_tr=list(func="gll_tr", aggFunc="sum", pars=list()),
                            gll_te=list(func="gll_te", aggFunc="sum", pars=list()),
                            #L2_f_tr=list(func="L2_f_tr", aggFunc="sum", pars=list()),
                            L2_f_te=list(func="L2_f_te", aggFunc="sum", pars=list()),
                            wL2_f_te=list(func="wL2_f_te", aggFunc="sum", pars=list()),
                            L2_k_te=list(func="L2_k_te", aggFunc="sum", pars=list()),
                            wL2_k_te=list(func="wL2_k_te", aggFunc="sum", pars=list()),
                            corrRKHS_k_te=list(func="corrRKHS_k_te", aggFunc="sum", pars=list()),
                            wcorrRKHS_k_te=list(func="wcorrRKHS_k_te", aggFunc="sum", pars=list()),
                            #PCEtr=list(func="PCEtr", aggFunc="sum", pars=list()),
                            PCEte=list(func="PCEte", aggFunc="sum", pars=list()),
                            #CCRtr=list(func="CCRtr", aggFunc="sum", pars=list()),
                            #CCRte=list(func="CCRte", aggFunc="sum", pars=list()),
                            KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                            WKCDC=list(func="WKCDC", aggFunc="sum",pars=list()),
                            #KCDCpval2=list(func="KCDCpval2", aggFunc="sum",pars=list(numPerms=100)),
                            KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                            KCMC=list(func="KCMC", aggFunc="sum",pars=list()),
                            WKCMC=list(func="WKCMC", aggFunc="sum",pars=list()),
                            KIIM=list(func="KIIM", aggFunc="sum",pars=list())) #,
                            #KCSC=list(func="KCSC", aggFunc="sum",pars=list()),
                            #KICSC=list(func="KICSC", aggFunc="sum",pars=list()),
                            #KCNSC=list(func="KCNSC", aggFunc="sum",pars=list()),
                            #KICCC_pca_ent2b=list(func="KICCC_pca_ent2b", aggFunc="sum",pars=list()), 
                            #KICCC_pca_entb=list(func="KICCC_pca_entb", aggFunc="sum",pars=list()))
                            
  
  measureParamsRed  <- list(#gll_tr=list(func="gll_tr", aggFunc="sum", pars=list()),
                            gll_te=list(func="gll_te", aggFunc="sum", pars=list()),
                            #L2_f_tr=list(func="L2_f_tr", aggFunc="sum", pars=list()),
                            L2_f_te=list(func="L2_f_te", aggFunc="sum", pars=list()),
                            #PCEtr=list(func="PCEtr", aggFunc="sum", pars=list()),
                            PCEte=list(func="PCEte", aggFunc="sum", pars=list()),
                            #CCRtr=list(func="CCRtr", aggFunc="sum", pars=list()),
                            #CCRte=list(func="CCRte", aggFunc="sum", pars=list()),
                            KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                            #KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                            KCMC=list(func="KCMC", aggFunc="sum",pars=list()))
  
  measureParamsRed2  <- list(KCDC=list(func="KCDC", aggFunc="sum",pars=list()),
                             WKCDC=list(func="WKCDC", aggFunc="sum",pars=list()),
                             #KCDCpval2=list(func="KCDCpval2", aggFunc="sum",pars=list(numPerms=100)),
                             KCRDC = list(func="KCRDC", aggFunc="sum", pars=list()),
                             KCMC=list(func="KCMC", aggFunc="sum",pars=list()),
                             KIIM=list(func="KIIM", aggFunc="sum",pars=list()),
                             KCSC=list(func="KCSC", aggFunc="sum",pars=list()),
                             KICSC=list(func="KICSC", aggFunc="sum",pars=list()),
                             KCNSC=list(func="KCNSC", aggFunc="sum",pars=list()),
                             KICCC_pca_ent2b=list(func="KICCC_pca_ent2b", aggFunc="sum",pars=list()), 
                             KICCC_pca_entb=list(func="KICCC_pca_entb", aggFunc="sum",pars=list()),
                            HSIC_cmem=list(func="HSIC_cmem", aggFunc="sum", pars=list()))
  
  
  measureParamsRed3  <- list(WKCDC_nonagg_loc=list(func="WKCDC_nonagg_loc", aggFunc="sum",pars=list()),
                             WKCDC_nonagg_loc=list(func="WKCDC_nonagg_superloc", aggFunc="sum",pars=list()),
                             WKCLC_nonagg_loc=list(func="WKCLC_nonagg_loc", aggFunc="sum",pars=list()),
                             WKCLC_nonagg_superloc=list(func="WKCLC_nonagg_superloc", aggFunc="sum",pars=list()),
                             WKCRDC_nonagg_loc=list(func="WKCRDC_nonagg_loc",aggFunc="sum",pars=list()),
                             WKCRDC_nonagg_loc=list(func="WKCRDC_nonagg_superloc",aggFunc="sum",pars=list()),
                             wL2_f_te_nonagg_loc=list(func="wL2_f_te_nonagg_loc", aggFunc="sum",pars=list()),
                             wL2_f_te_nonagg_superloc=list(func="wL2_f_te_nonagg_superloc", aggFunc="sum",pars=list()),
                             wL2_f_rel_nonagg_loc=list(func="wL2_f_rel_nonagg_loc", aggFunc="sum",pars=list()),
                             wL2_f_rel_nonagg_superloc=list(func="wL2_f_rel_nonagg_superloc", aggFunc="sum",pars=list())
)
  
  
  measureParams <- measureParamsRed3
  
  cmem_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn=learn.cmem, predict=predict.cmem, msrs=measureParams, calcMsrs=calcMsrs, makeKernel=makeKernel)
  assign(kernelTab1$learnerName[i], cmem_learner_aux)
}

cmem_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmem_learner_pack_lambda_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda")]
cmem_learner_pack_kernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX")]
cmem_learner_pack_gkernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="gamma_kernParsX")]
cmem_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]
cmem_learner_pack_gkernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="gamma_kernParsXY")]



# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

# do not add measures with underscore as this gets ignored... underscore are to identify variants of main measure eg. KCCC_ent vs KCCC_comp
ppTab <- data.frame(id=c("KCDC","WKCDC","WKCLC","WKCRDC","KCDCrel","KCDCpval2","KCRDC","KCMC","WKCMC", "KIIM","KCMCpval","KCSC","KICSC", "KCSCpval","KCCC", "KICCC","KCCCpval", "KCNSC",
                         "MCX","EDML","EDMB","TRE","gll_tr","gll_te","L2_f_tr","L2_f_te","wL2_f_te","L2_k_te","wL2_k_te","corrRKHS_k_te","wcorrRKHS_k_te","PCEtr","PCEte","CCRtr","CCRte","CVTeRNE", "CVTrRNE", "EAL", "DAL", "rankMajority","majority","HSIC"))




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

numMsrs <- nrow(ppTab)
ppTab$rank2Funcs <- c("quot","addQuot","differ")[rep(3,numMsrs)]
ppTab$rankFuncs <- rep("correctScoreToAdd", numMsrs)
ppTab$probFuncs <- rep("scoreToProb", numMsrs)
ppTab$argTypes <- c(rep("cmes", numMsrs))


ppTabCMEM1 <- ppTab

print("exits ob_cme_learner_v6g.R")

