# Construct Conditional Mean Embedding Measures

# remove(list=ls())
# repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
# setwd(repos)
# source("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_learner/func_learners_v3.R")

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameter=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                       paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), 
                       numPerPar=c("1","1","1","1","1","1-1-1"), stringsAsFactors=F)

paramTab2 <- data.frame(feature=c("rbf", "quad", "laplace1", "poly"), parameters=c("sigma", "offset", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_laplace1", "getFixedParams_poly"),
                        numPerPar=c("1","1","1","1-1-1"), stringsAsFactors=F)

aux <- expand.grid(kernelX=c("kern_rbf"), 
                   featureY=c("rff_rbf"), stringsAsFactors=F)
aux$name <- with(aux, paste(kernelX, featureY , sep="AND"))

lambdas <- c(10^seq(4,-9))
numKernParsX <- 3
numKernParsY <- 3

vanX <- NA
#rbfX <- list(sigma=c(0.0001, 0.001, 0.1, 1, 10, 100, 1000, 10000))
#rbfX <- list(sigma=c(0.1, 1, 10))
#rbfX <- list(sigma=c(1))
rbfX <- list(sigma=c(10, 100, 1000))

quadX <- list(offset=c(1e-3))
logX  <- list(degree=c(1))
polyTlinX <- list( degree=c(1,2,3), scale=c(1,2,3), offset=c(0,1,2))
polyTnonlinX <- list( degree=c(3,4,5), scale=c(1,2,3), offset=c(0,10,20))
laplace1X <- list(scale=c(1e-3))

toyY <- NA
vanY <- NA
#rbfY <- list(sigma=c(1))
rbfY <- list(sigma=c(1,100, 10000))
#rbfY <- list(sigma=c(0.1, 1, 10))
quadY <- list(offset=c(1))
logY  <- list(degree=c(0.001, 0.5, 0.98))
polyTlinY <- list( degree=c(1,2,3), scale=c(1,2,3), offset=c(0,1,2))
polyTnonlinY <- list( degree=c(3,4,5), scale=c(1,2,3), offset=c(0,10,20))
laplace1Y <- list(scale=c(1e-3))

kernelTab1 <- expand.grid(name=aux$name, centerLx=c( TRUE), centerKy=c(TRUE), optParms=c("none"), lambdaI=1:length(lambdas),  stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$featureY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")






aux <- strsplit(kernelTab1$kernelX, "_")
# short name removes kern or feat part
kernelTab1$kernelX_short <- sapply(aux, function(el) el[length(el)]) 
# short-short name removes part which specifies the name of parametrization ex poly.lin or poly.nonlin
kernelTab1$kernelX_short_short <- sapply(strsplit(kernelTab1$kernelX_short, "T"), function(el) el[1])


aux <- strsplit(kernelTab1$featureY, "_")
kernelTab1$featureY_short <- sapply(aux, function(el) el[length(el)])
kernelTab1$featureY_short_short <- sapply(strsplit(kernelTab1$featureY_short, "T"), function(el) el[1])

indxX <- match(kernelTab1$kernelX_short_short, paramTab1$kernel)
indxY <- match(kernelTab1$featureY_short_short, paramTab2$feature)

kernelTab1$kernParsX <- paramTab1$parameter[indxX]
kernelTab1$kernParsY <- paramTab2$parameter[indxY]
kernelTab1$paramFuncX <- paramTab1$paramFunc[indxX]
kernelTab1$paramFuncY <- paramTab2$paramFunc[indxY]

kernelTab1$numParsX <- sapply(strsplit(kernelTab1$kernParsX, "-"), function(el) sum(!is.na(el)))
kernelTab1$numParsY <- sapply(strsplit(kernelTab1$kernParsY, "-"), function(el) sum(!is.na(el)))

kernelTab1$numPerParX <- paramTab1$numPerPar[indxX]
kernelTab1$numPerParY <- paramTab2$numPerPar[indxY]


kernelTab1 <- lapply(1:nrow(kernelTab1), function(i){
  # i <- 16
  numParsX <- kernelTab1$numParsX[i]
  numParsY <- kernelTab1$numParsY[i]
  numPerParX <- kernelTab1$numPerParX[i]
  numPerParY <- kernelTab1$numPerParY[i]
  numPerParX <- as.numeric(strsplit(numPerParX, "-")[[1]])
  numPerParY <- as.numeric(strsplit(numPerParY, "-")[[1]])
  if(is.na(numPerParX)) numPerParX <- 1
  if(is.na(numPerParY)) numPerParY <- 1
  
  kernParsIndxs <- expand.grid(kernParsX_indx=1:prod(numPerParX), kernParsY_indx=1:prod(numPerParY))
  row <- kernelTab1[i,c("name","centerLx","centerKy","optParms","lambdaI", "kernelX", "featureY", "kernelX_short", 
                        "kernelX_short_short","featureY_short","featureY_short_short","kernParsX", 
                        "kernParsY", "paramFuncX", "paramFuncY","numParsX", "numParsY", "numPerParX", "numPerParY")] 
  
  row <- cbind(row, kernParsIndxs)
  
  return(row)
})
kernelTab1 <- do.call(rbind, kernelTab1)


kernelTab1$learnerName <- paste("cmfm", kernelTab1$name, centerCombo, "L2", kernelTab1$optParms, kernelTab1$lambdaI, 
                                kernelTab1$kernParsX_indx, kernelTab1$kernParsY_indx, sep="_")

dim(kernelTab1)

kernelTab1[, c("name","numParsX", "numParsY", "kernParsX_indx", "kernParsY_indx","lambdaI")]

for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 1
  print("********************")
  print(paste("i: ", i))
  print(kernelTab1[i,])
  
  # these are parameters needed to compute kernel: to feed to kernelMatrix or to phiy 
  # (not other such as center=T,or lambda,  which are not fed to kernelMatrix)
  parsX <- kernelTab1$kernParsX[i]
  if(!is.na(parsX)){
    # if kernel or feature does have parameters then we give them the "last-name" of kernel to avoid ambiguity 
    # between different kernels with the same parameter name
    parsX <- paste( strsplit(parsX, "-")[[1]], kernelTab1$kernelX_short_short[i], "X", sep=".") 
    names(parsX) <- strsplit(kernelTab1$kernParsX[i], "-")[[1]]
  }
  parsY <- kernelTab1$kernParsY[i]
  if(!is.na(parsY)){
    parsY <- c(paste(strsplit(parsY, "-")[[1]], kernelTab1$featureY_short_short[i], "Y", sep=".")) #
    names(parsY) <- strsplit(kernelTab1$kernParsY[i], "-")[[1]]
  } 
  
  aux <- strsplit(kernelTab1$featureY[i], "_")[[1]]
  if(aux[1]=="rff"){
    aux <- c("p_w","num_f", "seed", "map") #
    names(aux) <- aux
    parsYfeat <- c(aux)
  } else{
    parsYfeat <- c(NA, NA, NA, NA)
  }
  lossFun <-  "cmem_L2" #cmem_L2, TNRE
  lossFunList <- list(eval(parse(text=lossFun)))
  names(lossFunList) <- lossFun
  
  
  num_f <- 100
  map <- "cos"
  mySeed <- NULL # in toy example there seemed to be better results with diff rffs for forward and backward
  
  dist_p_w <- switch(kernelTab1$featureY_short_short[i], rbf="rnorm2", quad="rlaplace2", laplace1="rcauchy2", NA)
  
  
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]

  
  parsKernX <- switch(kernelTab1$kernelX_short[i], van=vanX, rbf=rbfX, quad=quadX, log=logX, polyTlin=polyTlinX, polyTnonlin=polyTnonlinX , laplace1=laplace1X)
  parsfeatY <- switch(kernelTab1$featureY_short[i], van=vanY, toy=toyY, rbf=rbfY, quad=quadY, log=logY, polyTlin=polyTlinY, polyTnonlin=polyTnonlinY , laplace1=laplace1Y)
  
  if(!is.na(parsKernX)){
    numPerParX <- as.numeric(strsplit(kernelTab1$numPerParX[i],"-")[[1]])
    indxArrX <- array(seq(prod(numPerParX)), dim=numPerParX)
    indxArrX <- indxArrX==kernelTab1$kernParsX_indx[i]
    #indxArrX
    #dim(indxArrX)
    #which(apply(indxArrX, 1, function(el) any(el)))
    indxX <- sapply(1:kernelTab1$numParsX[i], function(dm) which(apply(indxArrX, dm, function(mat) any(mat))))
    #indxX <- as.numeric(strsplit(dec2base(kernelTab1$kernParsX_indx[i]-1, base=max(2,kernelTab1$numParsX[i])),"")[[1]])+1
    parsKernX <- mapply(function(el, indx) el[indx], el=parsKernX, indx=indxX)
  }
  
  if(!is.na(parsfeatY)){ 
    numPerParY <- as.numeric(strsplit(kernelTab1$numPerParY[i],"-")[[1]])
    indxArrY <- array(seq(prod(numPerParY)), dim=numPerParY)
    indxArrY <- indxArrY==kernelTab1$kernParsY_indx[i]
    #indxArrY
    #dim(indxArrY)
    #which(apply(indxArrY, 1, function(el) any(el)))
    indxY <- sapply(1:kernelTab1$numParsY[i], function(dm) which(apply(indxArrY, dm, function(mat) any(mat))))
    #indxY <- as.numeric(strsplit(dec2base(kernelTab1$kernParsY_indx[i]-1, base=max(2,kernelTab1$numParsY[i])),"")[[1]])+1
    parsfeatY <- mapply(function(el, indx) el[indx], el=parsfeatY, indx=indxY)
  }
  
  

    
  nms <- c("lambda", "centerLx", "centerKy", parsX, parsY, parsYfeat)
  nonOptimizableParams <- list(list(val=lambdas[kernelTab1$lambdaI[i]]), list(val=centerLx), list(val=centerKy))
  for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j])))
  for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsfeatY[j])))
  nonOptimizableParams <- c(nonOptimizableParams, list(list(val=dist_p_w)),list(list(val=num_f)), list(list(val=mySeed)), list(list(val=map)))
  nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
  names(nonOptimizableParams) <- na.omit(nms)
  optimizableParams <- list()
    
  
  
  nonDataParams <- list(kernelX=list(name=kernelTab1$kernelX[i], pars=na.omit(parsX)), 
                        featureY=list(name=kernelTab1$featureY[i], pars=na.omit(c(parsY, parsYfeat))))
  
  
  heuristicSet <- na.omit(unique(c(kernelTab1$paramFuncX[i], kernelTab1$paramFuncY[i])))
  
  dataParams <- list(optimizable=optimizableParams, non_optimizable=nonOptimizableParams)
  hyperParams <- list(data=dataParams, non_data=nonDataParams)
  learnParams <- list(phiy=NULL, Ky=NULL, Lx=NULL, Blambda=NULL, alpha=NULL, Cks=NULL)
  optimizeParams <- list(losses=lossFunList, numFolds=5, testTrain="test")
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
                        CVTeRNSE=list(func="CVTeRNSE", aggFunc="sum",pars=list()),
                        #CVTrRNE=list(func="CVTrRNE", aggFunc="sum",pars=list()),
                        EAL=list(func="EAL", aggFunc="sum",pars=list()),
                        DAL=list(func="DAL", aggFunc="sum",pars=list())
                        #teErrNorm=list(func="teErrNorm", aggFunc="sum",pars=list()), 
                        #trErrNorm=list(func="trErrNorm", aggFunc="sum",pars=list()), 
                        #teNumUnd=list(func="teNumUnd", aggFunc="sum",pars=list()), 
                        #trNumUnd=list(func="trNumUnd", aggFunc="sum",pars=list())
                        )
  
  cmfm_learner_aux <- constructLearner(hyperParams, learnParams, optimizeParams, heuristicSet, optimizeSet, learn.cmfm_L2, predict.cmfm_L2, msrs=measureParams, calcMsrs=calcMsrs)
  assign(kernelTab1$learnerName[i], cmfm_learner_aux)
}


cmfm_learner_pack_none_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="none")]
cmfm_learner_pack_lambda_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda")]
cmfm_learner_pack_kernParsX_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX")]
cmfm_learner_pack_kernParsXY_1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY")]


# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

ppTab <- data.frame(id=c("KCDC","KCDCrel","KCDCpval","KCRDC","KCMC", "KCMCpval","KCSC", "KCSCpval","KCCC", "KCCCpval", "KCNSC",
                         "MCX","EDML","EDMB","TRE","CVTeRNE", "CVTeRNSE","CVTrRNE", "EAL", "DAL","bestTrainingError",
                         "teErrNorm", "trErrNorm", "teNumUnd", "trNumUnd"))



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

ppTab$rank2Funcs <- c("quot","addQuot","differ")[rep(3,25)]
ppTab$rankFuncs <- rep("correctScoreToAdd", 25)
ppTab$probFuncs <- rep("scoreToProb", 25)
ppTab$argTypes <- c(rep("cmfs", 25))



ppTabCMFM1 <- ppTab

