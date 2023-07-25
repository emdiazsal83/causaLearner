# Construct Conditional Mean Embedding Measures

# remove(list=ls())
# repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
# setwd(repos)
# source("/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_learner/func_learners_v3.R")
# library(oro.dicom)

# Notes
# 1) heuristic set is a function to set all parameters because one could perhaps
# want to set them in a joint fashion so you cant have individualized set parameter functions
# for each kernel

paramTab1 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameter=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                       paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), 
                       numPerPar=c("3","3","3","3","3","4-3-2"), stringsAsFactors=F)

paramTab2 <- data.frame(kernel=c("rbf", "quad", "log", "laplace1", "laplace2", "poly"), parameter=c("sigma", "offset", "degree", "scale", "scale", "degree-scale-offset"), 
                        paramFunc=c("getFixedParams_rbf2", "getFixedParams_quad", "getFixedParams_log", "getFixedParams_laplace1", "getFixedParams_laplace2", "getFixedParams_poly"), 
                        numPerPar=c("3","3","3","3","3","3-2-2"), stringsAsFactors=F)

aux <- expand.grid(kernelX=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_polyTnonlin"), 
                   kernelY=c("kern_van", "kern_rbf", "kern_quad", "kern_laplace1","kern_polyTnonlin"), stringsAsFactors=F)
aux$name <- with(aux, paste(kernelX, kernelY , sep="AND"))

lambdas <- c(0.001, 0.00001, 0.0000001)
#numKernParsX <- 4
#numKernParsY <- 5

vanX <- NA
rbfX <- list(sigma=c(0.0001, 1, 1e4))
quadX <- list(offset=c(1e-3, 1, 1e3))
logX  <- list(degree=c(0.001, 0.5, 0.98))
polyTlinX <- list( degree=c(1,2,3), scale=c(1,2,3), offset=c(0,1,2))
# 4-3-2
polyTnonlinX <- list( degree=c(2,3,4,5), scale=c(1,2,3), offset=c(0,10))
laplace1X <- list(scale=c(1e-3, 1, 1e3))

vanY <- NA
rbfY <- list(sigma=c(0.0001, 1, 1e4))
quadY <- list(offset=c(1e-3, 1, 1e3))
logY  <- list(degree=c(0.001, 0.5, 0.98))
polyTlinY <- list( degree=c(1,2,3), scale=c(1,2,3), offset=c(0,1,2))
# 3-2-2
polyTnonlinY <- list( degree=c(3,4,5), scale=c(1,2), offset=c(0,10))
laplace1Y <- list(scale=c(1e-3, 1, 1e3))

kernelTab1 <- expand.grid(name=aux$name, centerLx=c( TRUE), centerKy=c(TRUE), optParms=c("none"), lambdaI=1:length(lambdas),  stringsAsFactors=F)

aux <- strsplit(kernelTab1$name, split="AND")
kernelTab1$kernelX <- sapply(aux, function(el) el[[1]][1])
kernelTab1$kernelY <- sapply(aux, function(el) el[[2]][1])
centerCombo <- paste(c("n","c")[match(kernelTab1$centerLx, c(FALSE, TRUE))], c("n","c")[match(kernelTab1$centerKy, c(FALSE, TRUE))] , sep="")






aux <- strsplit(kernelTab1$kernelX, "_")
# short name removes kern or feat part
kernelTab1$kernelX_short <- sapply(aux, function(el) el[length(el)]) 
# short-short name removes part which specifies the name of parametrization ex poly.lin or poly.nonlin
kernelTab1$kernelX_short_short <- sapply(strsplit(kernelTab1$kernelX_short, "T"), function(el) el[1])


aux <- strsplit(kernelTab1$kernelY, "_")
kernelTab1$kernelY_short <- sapply(aux, function(el) el[length(el)])
kernelTab1$kernelY_short_short <- sapply(strsplit(kernelTab1$kernelY_short, "T"), function(el) el[1])

indxX <- match(kernelTab1$kernelX_short_short, paramTab1$kernel)
indxY <- match(kernelTab1$kernelY_short_short, paramTab2$kernel)

kernelTab1$kernParsX <- paramTab1$parameter[indxX]
kernelTab1$kernParsY <- paramTab2$parameter[indxY]
kernelTab1$paramFuncX <- paramTab1$paramFunc[indxX]
kernelTab1$paramFuncY <- paramTab2$paramFunc[indxY]

kernelTab1$numParsX <- sapply(strsplit(kernelTab1$kernParsX, "-"), function(el) sum(!is.na(el)))
kernelTab1$numParsY <- sapply(strsplit(kernelTab1$kernParsY, "-"), function(el) sum(!is.na(el)))

kernelTab1$numPerParX <- paramTab1$numPerPar[indxX]
kernelTab1$numPerParY <- paramTab2$numPerPar[indxY]


kernelTab1[,c("name","kernelX","kernelY")]

kernelTab1 <- lapply(1:nrow(kernelTab1), function(i){
  # i <- 25
  #print(paste("i: ", i))
  numParsX <- kernelTab1$numParsX[i]
  numParsY <- kernelTab1$numParsY[i]
  numPerParX <- kernelTab1$numPerParX[i]
  numPerParY <- kernelTab1$numPerParY[i]
  numPerParX <- as.numeric(strsplit(numPerParX, "-")[[1]])
  numPerParY <- as.numeric(strsplit(numPerParY, "-")[[1]])
  if(is.na(numPerParX)) numPerParX <- 1
  if(is.na(numPerParY)) numPerParY <- 1
  kernParsIndxs <- expand.grid(kernParsX_indx=1:prod(numPerParX), kernParsY_indx=1:prod(numPerParY))
  row <- kernelTab1[i,c("name","centerLx","centerKy","optParms","lambdaI", "kernelX", "kernelY", "kernelX_short", 
                        "kernelX_short_short","kernelY_short","kernelY_short_short","kernParsX", 
                        "kernParsY", "paramFuncX", "paramFuncY","numParsX", "numParsY", "numPerParX", "numPerParY")] 
  row <- cbind(row, kernParsIndxs)
  return(row)
})
kernelTab1 <- do.call(rbind, kernelTab1)


kernelTab1$learnerName <- paste("cmem", kernelTab1$name, centerCombo, "L2", kernelTab1$optParms, kernelTab1$lambdaI, 
                                kernelTab1$kernParsX_indx, kernelTab1$kernParsY_indx, sep="_")

dim(kernelTab1)

# A: 1 X and Y kernels with 0 parameter and  0 parameter levels per parameter: 
# B: 3 X and Y kernels with 1 parameter and  3 parameter levels per parameter: 
# Cx: 1 X  kernel with 3 parameters and 4-3-2 parameter levels for each one
# Cy: 1 Y  kernel with 3 parameters and 3-2-2 parameter levels for each one

# Ax vs Ay: 1*1 *3 =3
# Ax vs By: 1*(3*3) *3 = 27
# Ax vs Cy: 1*(3*2*2) *3 = 36
# Bx vs Ay: (3*3)*1 *3= 27
# Bx vs By: (3*3)*(3*3) *3 = 243
# Bx vs Cy: (3*3)*(3*2*2) *3= 324
# Cx vs Ay: (4*3*2)*1 *3 = 72
# Cx vs By: (4*3*2)*(3*3) *3= 648
# Cx vs Cy: (4*3*2)*(3*2*2) *3= 864
# total: 3 +27+36+27+243+324+72+648+864=2244




kernelTab1[, c("name","numParsX", "numParsY", "lambdaI", "kernParsX_indx", "kernParsY_indx")]
kernelTab1[which(kernelTab1$numParsX==3 & kernelTab1$numParsY==3 & kernelTab1$lambdaI>1), c("name","numParsX", "numParsY", "lambdaI","kernParsX_indx", "kernParsY_indx")]


for(i in 1:nrow(kernelTab1)){
  # none
  # i <- 1218
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
    parsY <- c(paste(strsplit(parsY, "-")[[1]], kernelTab1$kernelY_short_short[i], "Y", sep=".")) #
    names(parsY) <- strsplit(kernelTab1$kernParsY[i], "-")[[1]]
  } 
  
  
  lossFun <- cmem_L2
  centerLx <- kernelTab1$centerLx[i]
  centerKy <- kernelTab1$centerKy[i]

  
  parsKernX <- switch(kernelTab1$kernelX_short[i], van=vanX, rbf=rbfX, quad=quadX, log=logX, polyTlin=polyTlinX, polyTnonlin=polyTnonlinX , laplace1=laplace1X)
  parsKernY <- switch(kernelTab1$kernelY_short[i], van=vanY, rbf=rbfY, quad=quadY, log=logY, polyTlin=polyTlinY, polyTnonlin=polyTnonlinY , laplace1=laplace1Y)
  
  # formatC(as.numeric(sapply(0:26, dec2base, base=3)), width=3, flag="0")  
  
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
   
  if(!is.na(parsKernY)){ 
    numPerParY <- as.numeric(strsplit(kernelTab1$numPerParY[i],"-")[[1]])
    indxArrY <- array(seq(prod(numPerParY)), dim=numPerParY)
    indxArrY <- indxArrY==kernelTab1$kernParsY_indx[i]
    #indxArrY
    #dim(indxArrY)
    #which(apply(indxArrY, 1, function(el) any(el)))
    indxY <- sapply(1:kernelTab1$numParsY[i], function(dm) which(apply(indxArrY, dm, function(mat) any(mat))))
    #indxY <- as.numeric(strsplit(dec2base(kernelTab1$kernParsY_indx[i]-1, base=max(2,kernelTab1$numParsY[i])),"")[[1]])+1
    parsKernY <- mapply(function(el, indx) el[indx], el=parsKernY, indx=indxY)
  }
  

    
  nms <- c("lambda", "centerLx", "centerKy", parsX, parsY)
  nonOptimizableParams <- list(list(val=lambdas[kernelTab1$lambdaI[i]]), list(val=centerLx), list(val=centerKy))
  for(j in 1:length(parsX)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernX[j])))
  for(j in 1:length(parsY)) nonOptimizableParams <- c(nonOptimizableParams, list(list(val=parsKernY[j])))
  nonOptimizableParams <- nonOptimizableParams[which(!is.na(nms))]
  names(nonOptimizableParams) <- na.omit(nms)
  optimizableParams <- list()
    
    
  
  nonDataParams <- list(kernelX=list(name=kernelTab1$kernelX[i], pars=na.omit(parsX)), 
                        kernelY=list(name=kernelTab1$kernelY[i], pars=na.omit(parsY)))
  
  
  heuristicSet <- na.omit(unique(c(kernelTab1$paramFuncX[i], kernelTab1$paramFuncY[i])))
  
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


# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

ppTab <- data.frame(id=c("KCDC","KCDCrel","KCDCpval","KCRDC","KCMC", "KCMCpval","KCSC", "KCSCpval", "KCNSC"))




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

ppTab$rank2Funcs <- c("quot","addQuot","differ")[c(2,2,2,2,2,2,2,2,2)]
ppTab$rankFuncs <- rep("correctScoreToAdd", 9)
ppTab$probFuncs <- rep("scoreToProb", 9)
ppTab$argTypes <- c(rep("cmes", 9))



ppTabCMEM1 <- ppTab

