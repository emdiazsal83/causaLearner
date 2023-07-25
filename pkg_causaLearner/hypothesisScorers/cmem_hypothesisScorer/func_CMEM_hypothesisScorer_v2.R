# Conditional Mean Embeddding Measure (CMEM) Hypothesis Scorer functions
print("in func_CMEM_hypothesisScorer_v2.R")
library(psych) #geometric.mean
library(EnvStats) #rtri
library(rjson) # fromJSON
library(jsonlite)
library(dimRed)
library(pdfCluster) #kepdf
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_learner/func_learners_v5.R")

########################################################################################*
# Hypothesis scorer functions - a hypothesis score consists of a function that
# takes in data and parameters and gives back a score
#
########################################################################################*

# Conditional Mean embedding Measure hypothesis scorer
# approximate p(y|x_i) for all data with conditional mean embedding,  evaluate complexity
#
# normal hypothesis scorer params data -x, and approx hypothesis set -hypArray, plus
# 1. kernel packages for dependent variables (X) and independent variables (Y)
# 2. measures packages to apply to conditional mean embedding

# this one includes regressions between non parents, but for the Markov Equivalence Class
# approach (only compare scores for dags in same ME class) this is not necessary
cmem_hypScorer_noParents_deprecated <- function(x, hypArray, cmemLearner, ppTab=NULL, plot=FALSE){
  
  
  uniqueNoParentsList <- getUniqueNoParentsList(dags=hypArray)
  # we don't need parents lists of size just one - no hsic to perform
  indx <- which(sapply(uniqueNoParentsList, length)<=1)
  uniqueNoParentsList <- uniqueNoParentsList[-indx] 
  # we need to make sure that the regressions  implied between every uniqueNoParents are included in uniqueRegsList
  # ie. for a given unique set of nonParents the regressions implied in the fully connected directed graph are
  # i nuniqueRegs list. We'll use these as a kind of 2 way indpendence test but keeping same units as measure
  # as opposed to taking the hsic
  # a) First for each element in uniqueNoParentsList make a DG (not DAG) where nodes in the element are connected 
  # both ways
  p <- dim(hypArray)[2]
  if(length(uniqueNoParentsList)>0){
    dgArray <- sapply(uniqueNoParentsList, function(nodes){
      # i <- 1; nodes <- uniqueNoParentsList[[i]]
      dg <- matrix(0, p, p)
      indxCon <- combn(nodes, 2)
      indxCon <- cbind(indxCon, indxCon[2:1,])
      indxCon <- t(indxCon)
      dg[indxCon] <- 1
      return(dg)
    }, simplify="array")
  } else{
    dgArray <- array(0, dim=c(p,p,0))
    
  }
  
  dimnames(dgArray)[1:2] <- dimnames(hypArray)[1:2]
  names(dimnames(dgArray))[1:2] <- names(dimnames(hypArray))[1:2]
  
  # we can now eget unique regression list with both the original dag-array of hypothesis and
  # the dg-array we need to make sure all no-parent regressions are included
  dagDgArray <- abind(hypArray, dgArray, along=3)
  uniqueRegsList <- getUniqueRegsList(dags=dagDgArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # calculate CME measures needed for each unique regression
  print(paste("calculating cmes needed to evaluate all hypotheses"))
  cmemSet <- getMeasuresGivenDAGSet(uniqueRegsList, x, cmemLearner)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  print(paste("scoring all hypotheses"))
  
  
  
  scores <- assembleDagScores(dags=hypArray, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE)
  
  #scores <- melt(scores)
  #scores <- cast(scores, dag~measures+kernels)
  #dagNms <- scores$dag
  #scores <- scores[,-1]
  #scoreNms <- colnames(scores)
  #scores <- as.matrix(scores)
  #dimnames(scores) <- list(dag=dagNms, score=scoreNms)
  
  return(scores)
}

cmem_hypScorer <- function(x, hypArray, cmemLearner, ppTab=NULL, plot=FALSE, dataNm, folderSave,...){
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # calculate CME measures needed for each unique regression
  #print(paste("calculating cmes needed to evaluate all hypotheses"))
  cmemSet <- getMeasuresGivenDAGSet(uniqueRegsList, x, cmemLearner, dataNm, folderSave)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  #print(paste("scoring all hypotheses"))
  
  
  
  scores <- assembleDagScores(dags=hypArray, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE)
  
 
  
  return(scores)
}

boot_cmem_hypScorer <- function(x, hypArray, ppTab=NULL, plot=FALSE, dataNm, folderSave){
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # calculate CME measures needed for each unique regression
  #print(paste("calculating cmes needed to evaluate all hypotheses"))
  cmemSet <- getMeasuresGivenDAGSet_boot(uniqueRegsList, x)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  #print(paste("scoring all hypotheses"))
  
  
  
  scores <- assembleDagScores_boot(dags=hypArray, uniqueRegsList, cmemSet, prnt=FALSE)
  dmnms <- dimnames(scores)
  scores <- as.matrix(scores)
  dimnames(scores) <- dmnms
  
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}

boot_cmem_hypScorer_eqSig <- function(x, hypArray, jointFeats=FALSE, smoothFeats=FALSE, ppTab=NULL, plot=FALSE, dataNm, folderSave, hypSc_char= NULL){
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # calculate CME measures needed for each unique regression
  #print(paste("calculating cmes needed to evaluate all hypotheses"))
  fileSave <- paste("cmemBootEqSig",dataNm,".RData",sep="")
  fileSave2 <- paste("modMat",dataNm,".RData",sep="")
  if(paste(fileSave) %in% dir(folderSave)){
   load(file=paste(folderSave, fileSave, sep="")) 
    
    if(fileSave2 %in% dir(folderSave)){
      load(file=paste(folderSave, fileSave2, sep=""))
    } else{
      
      cmemSet2 <- cmemSet$measureList
      print("reformat")
      res <- reformat_msrs(cmemSet=cmemSet2, hypArray)
      modMatNew <- res$modMatNew
      testPts <- res$testPts
      
      save(list=c("modMatNew","testPts"), file=paste(folderSave, fileSave2, sep=""))
    }
  # cmemSetRead <- cmemSet  
  } else{
    pm <- proc.time()
    cmemSet <- getMeasuresGivenDAGSet_boot_eqSig(uniqueRegsList, x, jointFeats, smoothFeats)
    proc.time()-pm
    save("cmemSet", file=paste(folderSave, fileSave, sep=""))
    
    print("reformat")
    cmemSet2 <- cmemSet$measureList
    res <- reformat_msrs(cmemSet=cmemSet2, hypArray)
    modMatNew <- res$modMatNew
    testPts <- res$testPts
    
    save(list=c("modMatNew","testPts"), file=paste(folderSave, fileSave2, sep=""))  
  }
  
  #plot(cmemSet$measureList$y$`1`$wwkcrdc, wwKCRDCs_xy); abline(a=0, b=1, col="red")
  #plot(cmemSet$measureList$x$`1`$wwkcrdc, wwKCRDCs_yx); abline(a=0, b=1, col="red")
  #sum(cmemSet$measureList$y$`1`$wwkcrdc<cmemSet$measureList$x$`1`$wwkcrdc)/length(cmemSet$measureList$y$`1`$wwkcrdc)
  #sum(wwKCRDCs_xy<wwKCRDCs_yx)/length(wwKCRDCs_xy)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  #print(paste("scoring all hypotheses"))
  
  #colnames(cmemSet$x$`1`)
  #plot(cmemSet$y$`1`[,"wwkcrdc"],cmemSet$x$`1`[,"wwkcrdc"]); abline(a=0, b=1, col="red")
  #sum(cmemSet$y$`1`[,"wwkcrdc"]<cmemSet$x$`1`[,"wwkcrdc"])/nrow(cmemSet$x$`1`)
  
  cmemSet2 <- cmemSet$measureList
  print("assemble scores")
  #scores <- assembleDagScores_boot_eqSig(dags=hypArray, uniqueRegsList, cmemSet=cmemSet2, prnt=FALSE)
  
  scores <- rnd_hypScorer(x, hypArray, ppTab=NULL, plot=FALSE, dataNm=dataNm,  folderSave=folderSave)
  
  dmnms <- dimnames(scores)
  scores <- as.matrix(scores)
  dimnames(scores) <- dmnms
  
  
  names(dimnames(scores)) <- c("dag", "score")
  
  
  return(scores)
}

cmem_hypScorer_confounder_isomap <- function(x, hypArray, cmemLearner, ppTab=NULL, plot=FALSE, dataNm, folderSave){
  
  
  if(!"T" %in% colnames(x)) x <- cbind(x, T=NA)
  
  
  # need 33 (x->y, T),129 (x<-y, T), 73 (x<-T->y), 67 (x <- T<- y) and 13 (x->T->y) dags on hypArray
  hypArray <- sapply(c(33, 129, 73,67,13), function(id){
    mat <- dagIDtoMatrix(id,3)
    colnames(mat) <- rownames(mat) <-  colnames(x)
    return(mat)
  }, simplify="array")
  
  dimnames(hypArray)[[3]] <- as.character(c(33, 129, 73,67,13))
  names(dimnames(hypArray)) <- c("from","to","dag")
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # replace T with ISOMAP estimate
  
  
  
  #ggpairs(as.data.frame(x))
  num <- 5
  valve <- T
  emb <- 1
  while(class(emb)=="try-error" | valve){
    #print(num)
    #print(class(emb))
    valve <- FALSE
    emb <- try(embed(x[,c("x","y")], "Isomap", .mute = c("message"), knn = num))
    num <- num + 1
  }
  
  #dats <- cbind(x[,"T"],emb@data@data)
  dats <- cbind(x,emb@data@data)
  colnames(dats) <- c(colnames(dats)[1:3],"f1","f2")
  # dats <- as.data.frame(dats)
  # df <- as.data.frame(cbind(f1=dats$f1))
  # 
  # o <- order(dats$f1)
  # Tord <- dats$T[o]
  # f1ord <- dats$f1[o]
  # #plot(f1ord, Tord)
  # f1brk <- f1ord[which.max(abs(diff(Tord)))+1]
  # df$fac <- (dats$f1<f1brk)*1
  # df$T <- dats$T
  # mod <- lm(T~f1+fac, data=df)
  # summary(mod)
  # #plot(predict(mod), dats[,"T"])
  # 
  # dats$That <- predict(mod)
  
  #par(mfrow=c(2,2))
  #plot(emb, type = "2vars")
  #plot(x[,"T"], emb@data@data[,1])
  #plot(emb@data@data[,2], x[,"T"])
  #plot(dats$That, dats$T)
  #p <- ggpairs(dats)
  #print(p)
  #plot(dats$T, dats$f1)
  x[,"T"] <- dats[,"f1"] #dats$That
  
  # calculate CME measures needed for each unique regression
  #print(paste("calculating cmes needed to evaluate all hypotheses"))
  cmemSet <- getMeasuresGivenDAGSet(uniqueRegsList, x, cmemLearner, dataNm, folderSave)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  #print(paste("scoring all hypotheses"))
  
  #return kcdc_T_x, kcdc_x_T, kcdc_T_y, kcdc_y_T 
  
  names(cmemSet)
  
  cmemSet2 <- mapply(function(scrsNod, regsNod){
    # nod <- "x"; scrsNod <- cmemSet[[nod]]; regsNod <- uniqueRegsList[[nod]]
    nms <- rownames(regsNod)[apply(regsNod,2, function(col) which(col==1))]
    names(scrsNod) <- nms
    return(scrsNod)
  }, scrsNod=cmemSet, regsNod=uniqueRegsList, SIMPLIFY=FALSE)
  
  scores <- unlist(cmemSet2)
  
  
  
  #names(scores) <- gsub("x.1", "T_2_x", names(scores))
  #names(scores) <- gsub("y.1", "T_2_y", names(scores))
  #names(scores) <- gsub("T.1", "y_2_T", names(scores))
  #names(scores) <- gsub("T.2.KCDC", "x_2_T.KCDC", names(scores))
  #names(scores) <- gsub("T.2.KCMC", "x_2_T.KCMC", names(scores))
  
  #scores <- assembleDagScores(dags=hypArray, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE)
  
  #print(cor(dats$T, dats$That, method="spearman"))
  #scores <- c(scores, cor=cor(dats$T, dats$That, method="spearman"))
  
  return(scores)
}


cmemJoint_hypScorer <- function(x, hypArray, cmemLearner, ppTab=NULL, plot=FALSE){
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  # calculate CME measures needed for each unique regression
  #print(paste("calculating cmes needed to evaluate all hypotheses"))
  cmemSet <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner)
  
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  #print(paste("scoring all hypotheses"))
  
  
  
  scores <- assembleDagCMEM(x, dags=hypArray, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE)
  #scores2 <- assembleDagCMEM(x, dags=hypArray, uniqueRegsList, cmemSet2, cmemLearner, prnt=FALSE)
  
  #dev.new()
  #plot(as.numeric(scores), as.numeric(scores2)); abline(a=0, b=1, col="red")
  
  return(scores)
}

# Note for this learner all measures in measurePack have to be loss functions also
# but there can't be any measures in measure Pack not in loss functions
cmem_hypScorer_comp_deprecated <- function(x, hypArray, cmemLearner, ppTab=NULL, plot=FALSE, dataNm, folderSave){
  
  cmemLearner1 <- sub("compCombo", "lambda_kernParsXY", cmemLearner)
  cmemLearner2 <- sub("compCombo", "lambda_kernParsX", cmemLearner)
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner=cmemLearner1, dataNm=paste(dataNm,"1",sep="_"), folderSave)
  
  print("calculating opt kernParsXY measures")
  cmemSet1 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet)
  
  print("assembling opt kernParsXY measures into dag scores")
  scores1 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet1, cmemLearner1, prnt=FALSE)
  
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  loss <- "cmem_L2_f"
  reglrs <- names(eval(parse(text=cmemLearner2))$optimizeParams$losses) 
  aux <- c("KCDC", "KCMC", "KCRDC", "KCSC", "KCNSC", "KCCC_ent")
  reglrs <- intersect(reglrs, aux)
  

  print("calculating opt kernParsX measures per size")
  msrs_perSize <- lapply(unique(uniqueSizeRegs), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    indx <- which(sizeRegs==size)
    
    # step 2 - choose Hy : only regressions with same number of 
    # regressors will be compared so we choose sum of Hys for
    # same size regressions. Actually best do it by node?
    # one Hy per node
    cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner2, cmeSet, nodeIndx[indx], regIndx[indx])
    # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
    curves_size <- getCurves_comp(cmemLearnerSum, uniqueRegsList, x, loss, reglrs, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx], dataNm=paste(dataNm,"2",sep="_"), folderSave)
  
    msrs_size <- getMeasures_comp(curves=curves_size, reglrs)
    
    return(msrs_size)
  })
  
  msrs <- do.call(rbind, msrs_perSize)
  aux <- strsplit(rownames(msrs), "\\.")
  nod  <- sapply(aux, function(el) el[1]) 
  regr <- sapply(aux, function(el) el[2])
  indxNodReg <- order(nod,regr) 
  msrs <- msrs[indxNodReg,]
  nod <- nod[indxNodReg]
  regr <- regr[indxNodReg]
  
  print("arranging opt kernParsX measures per node")
  # step 5: assign scores back to the node-regression to prepare for assemblage.
  msrs_perNodePerReg <- lapply(nodes, function(nodeTo){
    msr_nod <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indx_nod_reg <- which(nod == nodeTo & regr==numReg)
      msr_nod_reg <- msrs[indx_nod_reg,] 
      return(msr_nod_reg)
    })
    names(msr_nod) <- numRegs[[nodeTo]]
    return(msr_nod)
  })
  names(msrs_perNodePerReg) <- nodes
  
  print(class(msrs_perNodePerReg))
  print(length(msrs_perNodePerReg))
  print(class(msrs_perNodePerReg[[1]]))
  print(msrs_perNodePerReg)
  
  print("assembling opt kernParsX measures into dags")
  # step 6 - assemble dags with scores as before
  scores2 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                              msrs_perNodePerReg, cmemLearner2, prnt=FALSE)
  
  colnames(scores2) <- paste(colnames(scores2), "comp", sep="_")
  scores <- cbind(scores1, scores2)
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}

cmem_hypScorer_comp <- function(x, hypArray, cmemLearner1, cmemLearner2, ppTab=NULL, plot=FALSE, dataNm,  folderSave){
  
  
  aux <- strsplit(cmemLearner1, "_")[[1]]
  indxDEL1 <- which(aux=="DEL1")
  indxDEL2 <- which(aux=="DEL2")
  cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
  #ini <- gregexpr("compCombo",cmemLearner_save)
  #cmemLearner1_save <- paste(substr(cmemLearner_save, 1, as.numeric(ini)-2), aux1, sep="_")
  #cmemLearner2_save <- paste(substr(cmemLearner_save, 1, as.numeric(ini)-2), aux2, sep="_")
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner=cmemLearner1, nm=paste(dataNm, cmemLearner1_save, "1",sep="_"), folderSave, round="first")
  
  print("calculating opt kernParsXY measures")
  cmemSet1 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet)
  
  print("assembling opt kernParsXY measures into dag scores")
  scores1 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet=cmemSet1, cmemLearner=cmemLearner1, prnt=FALSE)
  colnames(scores1) <- paste(colnames(scores1), "XY", sep="_")
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  loss <- "cmem_L2_f"
  reglrs <- names(eval(parse(text=cmemLearner2))$optimizeParams$losses) 
  aux <- c("KCDC", "KCMC", "KCRDC", "KCSC", "KCNSC", "KCCC_ent","KCCC_pca_ent")
  reglrs <- intersect(reglrs, aux)
  
  
  print("building the appropriate y-hilbert space measures per size")
  cmemLearnersSum_preHyperParamSelec_perSize <- lapply(sort(unique(uniqueSizeRegs)), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    indx <- which(sizeRegs==size)
    
    # step 2 - choose Hy : only regressions with same number of 
    # regressors will be compared so we choose sum of Hys for
    # same size regressions. Actually best do it by node?
    # one Hy per node
    cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
    cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
    
    return(cmemLearnerSum)
  })
  names(cmemLearnersSum_preHyperParamSelec_perSize) <- sort(unique(uniqueSizeRegs))
  
  nodes <- names(uniqueRegsList)
  
  print("putting back into nodeTo-numReg structure")
  cmeSet2  <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      whatSize <- sizeRegs[match(paste(nodeTo, numReg, sep="."), names(sizeRegs))]
      cmemLearnerSumAux <- cmemLearnersSum_preHyperParamSelec_perSize[[as.character(whatSize)]]
      
      return(cmemLearnerSumAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  names(cmeSet2) <- nodes
  
  # step 3 - regression: get regressions for each node in dags
  # with CV on L2-loss for lambda and parsKernX
  print("calculating kernParsX learners")
  cmeSet2 <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner=cmeSet2, nm=paste(dataNm, cmemLearner2, "2",sep="_"), folderSave)
  
  print("calculating opt kernParsX measures")
  cmemSet2 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet=cmeSet2)
  
  print("assembling opt kernParsX measures into dag scores")
  scores2 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet2, cmemLearner2, prnt=FALSE)
  colnames(scores2) <- paste(colnames(scores2), "X", sep="_")
  
  print("obtain maximum efficiency curves for kernParsX learners")
  
  curves <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  
  
  print("comparing curves of sames size to get measures per size")
  msrs_perSize <- lapply(unique(uniqueSizeRegs), function(size){
    # size <- 2
    #print("*******")
    #print(paste("size: ", size))
    
    indx <- which(sizeRegs==size)
    curves_size <- mapply(function(node, reg) curves[[node]][[reg]], node=nodeIndx[indx], reg=regIndx[indx], SIMPLIFY = F)
    names(curves_size) <- names(indx)
    msrs_size <- getMeasures_comp(curves=curves_size, reglrs)
    #rownames(msrs_size) <- names(indx)
    
    return(msrs_size)
  })
  
  msrs <- do.call(rbind, msrs_perSize)
  aux <- strsplit(rownames(msrs), "\\.")
  nod  <- sapply(aux, function(el) el[1]) 
  regr <- sapply(aux, function(el) el[2])
  indxNodReg <- order(nod,regr) 
  msrs <- msrs[indxNodReg,]
  nod <- nod[indxNodReg]
  regr <- regr[indxNodReg]
  
  msrs_nms <- colnames(msrs)
  colnames(msrs) <- sapply(strsplit(msrs_nms, "_"), function(el) paste(el[1:(length(el)-1)], collapse="_"))
  
  print("arranging opt kernParsX measures per node")
  # step 5: assign scores back to the node-regression to prepare for assemblage.
  msrs_perNodePerReg <- lapply(nodes, function(nodeTo){
    msr_nod <- lapply(numRegs[[nodeTo]], function(numReg){
      # nodeTo <- "1"; numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indx_nod_reg <- which(nod == nodeTo & regr==numReg)
      msr_nod_reg <- msrs[indx_nod_reg,] 
      return(msr_nod_reg)
    })
    names(msr_nod) <- numRegs[[nodeTo]]
    return(msr_nod)
  })
  names(msrs_perNodePerReg) <- nodes
  
  
  print("assembling opt kernParsX measures into dags")
  # step 6 - assemble dags with scores as before
  scores3 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet=msrs_perNodePerReg, cmemLearner=cmemLearner2, prnt=FALSE)
  
  dimnames(scores3)$score <- msrs_nms
  colnames(scores3) <- paste(colnames(scores3), "comp", sep="_")
  scores <- cbind(scores1, scores2, scores3)
  #scores <- scores3
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}

# PARETO REGION
cmem_hypScorer_comp <- function(x, hypArray, cmemLearner1, cmemLearner2, noiseLearner=NULL, augmentData=FALSE, ppTab=NULL, plot=FALSE, dataNm,  folderSave){
  
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[3]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  
  # this is to avoid duplicating calculations... if 2 cmem_comp hyp scorers share a first
  # learner up to kernel params, alpha (kcmc, l2 or kcsc) AND loss function then we would
  # recycle.. now we loss function doesnt have to be the same as we calculate all loss
  # functions the first time round and then just change the optimizing function and main loss
  # and pass through setParams which doesnt re-calculate the grid just the optimal params
  aux <- strsplit(cmemLearner1, "_")[[1]]
  indxDEL1 <- which(aux=="DEL1")
  indxDEL2 <- which(aux=="DEL2")
  cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
  nm1 <- paste(dataNm, cmemLearner1_save, "1",sep="_")
  #cmemLearner=cmemLearner1; nm=nm1; round="first"
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                               cmemLearner=cmemLearner1, 
                               noiseLearner=noiseLearner,
                               augmentData=augmentData,
                               dataNm=dataNm,
                               nm=nm1, 
                               folderSave, round="first")
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
  preds <- lapply(cmeSet, function(cme){
    # cme <- cmeSet[[2]]
    auxLearn <- cme[[1]]$learn(cme[[1]])
    pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
    plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
  })
  }
  
  
  print("calculating opt kernParsXY measures")
  cmemSet1 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet)
  
  print("assembling opt kernParsXY measures into dag scores")
  scores1 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet=cmemSet1, cmemLearner=cmemLearner1, prnt=FALSE)
  colnames(scores1) <- paste(colnames(scores1), "XY", sep="_")
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  
  
  print("building the appropriate y-hilbert space measures per size")
  cmemLearnersSum_preHyperParamSelec_perSize <- lapply(sort(unique(uniqueSizeRegs)), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    if(size==0){
      cmemLearnerSum <- NA
    } else{
      indx <- which(sizeRegs==size)
      indxNotZero <- which(sizeRegs>0)
      
      # step 2 - choose Hy : only regressions with same number of 
      # regressors will be compared so we choose sum of Hys for
      # same size regressions. Actually best do it by node?
      # one Hy per node
    
      #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    
      # modifying to do a direct sum where the best parameter from all regressions
      # from the point of view that it doesn't make sense to add kcmcs or kcds from different
      # Hy spaces
      #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
      # nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
      cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
                                                nodeIndx=nodeIndx[indxNotZero], 
                                                regIndx=regIndx[indxNotZero])
    
      # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
      cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
    }
      return(cmemLearnerSum)
  })
  names(cmemLearnersSum_preHyperParamSelec_perSize) <- sort(unique(uniqueSizeRegs))
  
  nodes <- names(uniqueRegsList)
  
  print("putting back into nodeTo-numReg structure")
  cmeSet2  <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      whatSize <- sizeRegs[match(paste(nodeTo, numReg, sep="."), names(sizeRegs))]
      cmemLearnerSumAux <- cmemLearnersSum_preHyperParamSelec_perSize[[as.character(whatSize)]]
      
      
      return(cmemLearnerSumAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  names(cmeSet2) <- nodes
  
  
  # step 3 - regression: get regressions for each node in dags
  # with CV on L2-loss for lambda and parsKernX
  print("calculating kernParsX learners")
  nm1 <- paste(dataNm, cmemLearner1, "1",sep="_")
  nm2 <- paste(cmemLearner2, "2",sep="_")
  #cmemLearner=cmeSet2; nm=paste("nm1",nm1,"nm2",nm2, sep="_"); round="last"
  cmeSet2 <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                                cmemLearner=cmeSet2, 
                                noiseLearner=noiseLearner,
                                augmentData=augmentData,
                                dataNm=dataNm,
                                nm=paste("nm1",nm1,"nm2",nm2, sep="_"), 
                                folderSave)
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
    preds <- lapply(cmeSet2, function(cme){
     # cme <- cmeSet[[1]]
      auxLearn <- cme[[1]]$learn(cme[[1]])
      pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
      plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
    })
  }

  print("calculating opt kernParsX measures")
  cmemSet2 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet=cmeSet2)
  
  print("assembling opt kernParsX measures into dag scores")
  scores2 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet2, cmemLearner2, prnt=FALSE)
  colnames(scores2) <- paste(colnames(scores2), "X", sep="_")
  
  
  loss <- "cmem_L2_f"
  #reglrs <- names(eval(parse(text=cmemLearner2))$optimizeParams$losses) 
  reglrs <- names(cmeSet2[[1]][[1]]$optimizeParams$losses)
  aux <- c("KCDC", "WKCDC","KCMC","WKCMC", "KIIM","KCRDC", "KCSC","KICSC", "KCNSC", "KCCC_ent","KCCC_pca_ent","KCCC_pca_ent2","KICCC_pca_ent2","KICCC_pca_ent2b","KICCC_pca_entb")
  reglrs <- intersect(reglrs, aux)
  print("obtain maximum efficiency curves for kernParsX learners")
  curves_nodeReg <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  print("assembling opt kernParsX measures into dags")
  # step 6 - assemble dags with scores as before
  scores3 <- assembleDagCurvesIntoDagCurvesThenScores(dags=hypArray, uniqueRegsList, 
                               curves=curves_nodeReg, reglrs, prnt=FALSE, plot=NULL)
  colnames(scores3) <- paste(colnames(scores3), "comp", sep="_")
  scores3 <- scores3[match(rownames(scores1), rownames(scores3)),]
  loss <- "wcmem_L2_f"
  print("obtain maximum efficiency curves for kernParsX learners 2")
  curves_nodeReg <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  print("assembling opt kernParsX measures into dags 2")
  # step 6 - assemble dags with scores as before
  scores4 <- assembleDagCurvesIntoDagCurvesThenScores(dags=hypArray, uniqueRegsList, 
                                                      curves=curves_nodeReg, reglrs, prnt=FALSE, plot=NULL)
  colnames(scores4) <- paste(colnames(scores4), "wcomp", sep="_")
  scores4 <- scores4[match(rownames(scores1), rownames(scores4)),]
  
  
  scores <- cbind(scores1, scores2, scores3, scores4)
  #scores <- scores3
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}

# GRANULAR COMPARISONS
cmem_hypScorer_comp <- function(x, hypArray, cmemLearner1, cmemLearner2, noiseLearner=NULL, augmentData=FALSE, ppTab=NULL, plot=FALSE, dataNm,  folderSave, hypSc_char= NULL){
  
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[3]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  
  # this is to avoid duplicating calculations... if 2 cmem_comp hyp scorers share a first
  # learner up to kernel params, alpha (kcmc, l2 or kcsc) AND loss function then we would
  # recycle.. now we loss function doesnt have to be the same as we calculate all loss
  # functions the first time round and then just change the optimizing function and main loss
  # and pass through setParams which doesnt re-calculate the grid just the optimal params
  aux <- strsplit(cmemLearner1, "_")[[1]]
  indxDEL1 <- which(aux=="DEL1")
  indxDEL2 <- which(aux=="DEL2")
  cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
  nm1 <- paste(dataNm, cmemLearner1_save, "1",sep="_")
  #cmemLearner=cmemLearner1; nm=nm1; round="first"
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                               cmemLearner=cmemLearner1, 
                               noiseLearner=noiseLearner,
                               augmentData=augmentData,
                               dataNm=dataNm,
                               nm=nm1, 
                               folderSave, round="first")
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
    preds <- lapply(cmeSet, function(cme){
      # cme <- cmeSet[[2]]
      auxLearn <- cme[[1]]$learn(cme[[1]])
      pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
      plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
    })
  }
  
  
  #print("calculating opt kernParsXY measures")
  #cmemSet1 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet, cmemLearner=cmemLearner1)
  
  #print("assembling opt kernParsXY measures into dag scores")
  #scores1 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
  #                             cmemSet=cmemSet1, cmemLearner=cmemLearner1, prnt=FALSE)
  #colnames(scores1) <- paste(colnames(scores1), "XY", sep="_")
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  
  
  print("building the appropriate y-hilbert space measures per size")
  cmemLearnersSum_preHyperParamSelec_perSize <- lapply(sort(unique(uniqueSizeRegs)), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    if(size==0){
      cmemLearnerSum <- NA
    } else{
      indx <- which(sizeRegs==size)
      indxNotZero <- which(sizeRegs>0)
      
      # step 2 - choose Hy : only regressions with same number of 
      # regressors will be compared so we choose sum of Hys for
      # same size regressions. Actually best do it by node?
      # one Hy per node
      
      #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
      
      # modifying to do a direct sum where the best parameter from all regressions
      # from the point of view that it doesn't make sense to add kcmcs or kcds from different
      # Hy spaces
      #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
      # nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
      cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
                                                  nodeIndx=nodeIndx[indxNotZero], 
                                                  regIndx=regIndx[indxNotZero])
      
      # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
      cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
    }
    return(cmemLearnerSum)
  })
  names(cmemLearnersSum_preHyperParamSelec_perSize) <- sort(unique(uniqueSizeRegs))
  
  nodes <- names(uniqueRegsList)
  
  print("putting back into nodeTo-numReg structure")
  cmeSet2  <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      whatSize <- sizeRegs[match(paste(nodeTo, numReg, sep="."), names(sizeRegs))]
      cmemLearnerSumAux <- cmemLearnersSum_preHyperParamSelec_perSize[[as.character(whatSize)]]
      
      
      return(cmemLearnerSumAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  names(cmeSet2) <- nodes
  
  
  # step 3 - regression: get regressions for each node in dags
  # with CV on L2-loss for lambda and parsKernX
  print("calculating kernParsX learners")
  nm1 <- paste(dataNm, cmemLearner1, "1",sep="_")
  nm2 <- paste(cmemLearner2, "2",sep="_")
  #cmemLearner=cmeSet2; nm=paste("nm1",nm1,"nm2",nm2, sep="_"); round="last"
  cmeSet2 <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                                cmemLearner=cmeSet2, 
                                noiseLearner=noiseLearner,
                                augmentData=augmentData,
                                dataNm=dataNm,
                                nm=paste("nm1",nm1,"nm2",nm2, sep="_"), 
                                folderSave)
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
    preds <- lapply(cmeSet2, function(cme){
      # cme <- cmeSet[[1]]
      auxLearn <- cme[[1]]$learn(cme[[1]])
      pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
      plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
    })
  }
  
  print("calculating opt kernParsX measures")
  #cmeSet=cmeSet2; cmemLearner=cmemLearner2
  cmemSet2 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet=cmeSet2, cmemLearner=cmemLearner2)
  #print("cmemSet2")
  #print(cmemSet2)
  print("assembling opt kernParsX measures into dag scores")
  scores2 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet2, cmemLearner2, prnt=FALSE)
  #colnames(scores2) <- paste(colnames(scores2), "X", sep="_")
  
  
  scrs <- sapply(strsplit(colnames(scores2), "_"), function(el) paste(el[1:(length(el)-1)], collapse="_"))
  scores3 <- sapply(unique(scrs), function(uni_scr){
    indx <- which(scrs==uni_scr)
  
    #col <- apply(scores2[,indx],2,function(col) col[1]<col[2])+1
    #plot(dat$x[,1], dat$x[,2], col=col)
  
    res <- length(indx)-apply(apply(scores2[,indx], 2, function(col) col==min(col)),1,sum)
  })
    
  scores <- cbind(scores3)
  #scores <- scores3
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}


cmem_hypScorer_comp_boot <- function(x, hypArray, cmemLearner1, cmemLearner2, numBoots, numPerBoot, noiseLearner=NULL, augmentData=FALSE, ppTab=NULL, plot=FALSE, dataNm,  folderSave){
  
  
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  
  aux <- strsplit(cmemLearner1, "_")[[1]]
  indxDEL1 <- which(aux=="DEL1")
  indxDEL2 <- which(aux=="DEL2")
  cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
  
  set.seed(12)
  xs <- lapply(1:numBoots, function(i){
    smpl <- sample(1:nrow(x), numPerBoot)
    res <- x[smpl,]
    return(res)
  })
  
  counts <- 1:numBoots
  
  scores <- mcmapply(FUN=function(x, cnt){
    # i <- 3; x <- xs[[i]]; cnt <- counts[i]
    print("-------------------------------------")
    print(paste("boot number: ", cnt, " of ", numBoots))
    print(paste("nrow(x): ", nrow(x)))
  # this is to avoid duplicating calculations... if 2 cmem_comp hyp scorers share a first
  # learner up to kernel params, alpha (kcmc, l2 or kcsc) AND loss function then we would
  # recycle.. now we loss function doesnt have to be the same as we calculate all loss
  # functions the first time round and then just change the optimizing function and main loss
  # and pass through setParams which doesnt re-calculate the grid just the optimal params
  
  scores <- try({
  nm1 <- paste(dataNm, cmemLearner1_save,numPerBoot,cnt,"1",sep="_") #, "bt"
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                               cmemLearner=cmemLearner1, 
                               noiseLearner=noiseLearner,
                               augmentData=augmentData,
                               dataNm=dataNm,
                               nm=nm1, 
                               folderSave, round="first")
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
    preds <- lapply(cmeSet, function(cme){
      # cme <- cmeSet[[2]]
      auxLearn <- cme[[1]]$learn(cme[[1]])
      pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
      plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
    })
  }
  
  
  print("calculating opt kernParsXY measures")
  cmemSet1 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet)
  
  print("assembling opt kernParsXY measures into dag scores")
  scores1 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet=cmemSet1, cmemLearner=cmemLearner1, prnt=FALSE)
  colnames(scores1) <- paste(colnames(scores1), "XY", sep="_")
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  
  
  print("building the appropriate y-hilbert space measures per size")
  cmemLearnersSum_preHyperParamSelec_perSize <- lapply(sort(unique(uniqueSizeRegs)), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    indx <- which(sizeRegs==size)
    
    # step 2 - choose Hy : only regressions with same number of 
    # regressors will be compared so we choose sum of Hys for
    # same size regressions. Actually best do it by node?
    # one Hy per node
    
    #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    
    # modifying to do a direct sum where the best parameter from all regressions
    # from the point of view that it doesn't make sense to add kcmcs or kcds from different
    # Hy spaces
    #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
    # nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
                                                nodeIndx=nodeIndx, regIndx=regIndx)
    
    # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
    cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
    
    return(cmemLearnerSum)
  })
  names(cmemLearnersSum_preHyperParamSelec_perSize) <- sort(unique(uniqueSizeRegs))
  
  nodes <- names(uniqueRegsList)
  
  print("putting back into nodeTo-numReg structure")
  cmeSet2  <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      whatSize <- sizeRegs[match(paste(nodeTo, numReg, sep="."), names(sizeRegs))]
      cmemLearnerSumAux <- cmemLearnersSum_preHyperParamSelec_perSize[[as.character(whatSize)]]
      
      
      return(cmemLearnerSumAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  names(cmeSet2) <- nodes
  
  
  # step 3 - regression: get regressions for each node in dags
  # with CV on L2-loss for lambda and parsKernX
  print("calculating kernParsX learners")
  nm1 <- paste(dataNm, cmemLearner1, "1",sep="_")
  nm2 <- paste(cmemLearner2, "2",sep="_")
  #cmemLearner=cmeSet2; nm=paste("nm1",nm1,"nm2",nm2, sep="_"); round="last"
  cmeSet2 <- getCMEsGivenDAGSet(uniqueRegsList, x, 
                                cmemLearner=cmeSet2, 
                                noiseLearner=noiseLearner,
                                augmentData=augmentData,
                                dataNm=dataNm,
                                nm=paste(nm1,nm2,numPerBoot,cnt, sep="_"), #, "bt","nm2","nm1"
                                folderSave)
  
  # for plotting when its x->y vs x<-y
  if(FALSE){
    preds <- lapply(cmeSet2, function(cme){
      # cme <- cmeSet[[1]]
      auxLearn <- cme[[1]]$learn(cme[[1]])
      pred <- auxLearn$predict(learner=auxLearn, data=auxLearn$hyperParams$trainData, forLoss=F)
      plotPredCMEM(auxLearn, pred=pred, var="k",indx=1:9)
    })
  }
  
  print("calculating opt kernParsX measures")
  cmemSet2 <- getMeasuresGivenCMEsList(uniqueRegsList, x, cmeSet=cmeSet2)
  
  print("assembling opt kernParsX measures into dag scores")
  scores2 <- assembleDagScores(dags=hypArray, uniqueRegsList, 
                               cmemSet2, cmemLearner2, prnt=FALSE)
  colnames(scores2) <- paste(colnames(scores2), "X", sep="_")
  
  
  loss <- "cmem_L2_f"
  #reglrs <- names(eval(parse(text=cmemLearner2))$optimizeParams$losses) 
  reglrs <- names(cmeSet2[[1]][[1]]$optimizeParams$losses)
  aux <- c("KCDC", "WKCDC","KCMC","WKCMC", "KIIM","KCRDC", "KCSC","KICSC", "KCNSC", "KCCC_ent","KCCC_pca_ent","KCCC_pca_ent2","KICCC_pca_ent2","KICCC_pca_ent2b","KICCC_pca_entb")
  reglrs <- intersect(reglrs, aux)
  print("obtain maximum efficiency curves for kernParsX learners")
  curves_nodeReg <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  print("assembling opt kernParsX measures into dags")
  # step 6 - assemble dags with scores as before
  scores3 <- assembleDagCurvesIntoDagCurvesThenScores(dags=hypArray, uniqueRegsList, 
                                                      curves=curves_nodeReg, reglrs, prnt=FALSE, plot=NULL)
  colnames(scores3) <- paste(colnames(scores3), "comp", sep="_")
  scores3 <- scores3[match(rownames(scores1), rownames(scores3)),]
  loss <- "wcmem_L2_f"
  print("obtain maximum efficiency curves for kernParsX learners 2")
  curves_nodeReg <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  print("assembling opt kernParsX measures into dags 2")
  # step 6 - assemble dags with scores as before
  scores4 <- assembleDagCurvesIntoDagCurvesThenScores(dags=hypArray, uniqueRegsList, 
                                                      curves=curves_nodeReg, reglrs, prnt=FALSE, plot=NULL)
  colnames(scores4) <- paste(colnames(scores4), "wcomp", sep="_")
  scores4 <- scores4[match(rownames(scores1), rownames(scores4)),]
  
  
  scores <- cbind(scores1, scores2, scores3, scores4)
  #scores <- scores3
  names(dimnames(scores)) <- c("dag", "score")
  scores
  })

  return(scores)
  
  }, x=xs, cnt=counts, mc.cores=1, SIMPLIFY="array")

  dimnames(scores)[[3]] <- 1:numBoots
  names(dimnames(scores))[3] <- "boot"
  scores2 <- melt(scores)
  scores3 <- cast(scores2, dag~score+boot, value="value")
  scores3 <- scores3[match(scores3$dag, dimnames(scores)$dag),]
  scores4 <- as.matrix(scores3[,2:ncol(scores3)])
  dimnames(scores4) <- list(dag=scores3$dag, score=colnames(scores3)[2:ncol(scores3)])
  return(scores4)
}


cmem_hypScorer_comp_nn <- function(x, hypArray, ppTab=NULL, plot=FALSE, dataNm,  folderSave){
  

  print("read in data jason file")
  #print(folderSave)
  #print(dir(folderSave))
  #print(dir(paste(folderSave, "problem_",dataNm, "/", sep="")))
  file <- c("10_12.json", "40_30.json", "60_50.json")[1]
  
  #dir(paste(folderSave, dataNm, "/",sep=""))
  char <- paste(folderSave, "problem_", dataNm, "/", file, sep="")
  aux <- readLines(con=char)
  aux <- gsub("NaN", "null", aux)
  curves <- jsonlite:::fromJSON(txt=aux)
  
  
  curves <- lapply(curves$metricskdcd, function(el) na.omit(as.data.frame(el)))
  ids <- getHypID(hypArray)
  ids$dir <- c("Y->X","X->Y")
  
  names(curves) <- ids$id[match(names(curves), ids$dir)]
  
  scores2 <- sapply(curves, function(el){
    # el <- curves[[2]]
    indx <- which.min(el$RMSE)
    res <- as.numeric(el[indx,c("KCDC","KCMC")])
    return(res)
  })
  scores2 <- t(scores2)
  colnames(scores2) <- c("KCDC","KCMC")
  colnames(scores2) <- paste(colnames(scores2), "X", sep="_")
  
  reglrs <- c("KCDC","KCMC")
  loss  <- "cmem_L2_f"
  curves_dags <- lapply(curves, function(curve){
    # curve <- curves[[1]]
    curves_dags_reglr <- lapply(reglrs, function(reglr){
      # reglr <- reglrs[1]
      envPts <- extractLeftEnvPts(log(curve[,reglr],10), curve[,"RMSE"])
      if(length(envPts$x)==1){
      envPts$x <- c(envPts$x, envPts$x+0.00001)
      envPts$y <- c(envPts$y, envPts$y-0.00001)
    }
      #plot(log(curve[,reglr],10), curve[,"RMSE"])
      #lines(envPts$x, envPts$y, type="p", cex=2, col="red")
      splf_yx <- splinefun(envPts$y, envPts$x, method="monoH.FC")
      splf_xy <- splinefun(envPts$x, envPts$y, method="monoH.FC")
      #yy_1 <- seq(envPts$y[1], envPts$y[length(envPts$y)], length.out=100)
      #xx_1 <- splf_yx(yy_1)
      #xx_2 <- seq(envPts$x[1], envPts$x[length(envPts$x)], length.out=100)
      #yy_2 <- splf_xy(xx_2)
      #plot(log(curve[,reglr],10), curve[,"RMSE"])
      #lines(xx_1, yy_1, col="red", lty=2)
      #lines(xx_2, yy_2, col="blue", lty=3)
      curve <- list(reglrRng=range(envPts$x), lossRng=range(envPts$y), fun_xy=splf_xy, fun_yx=splf_yx, loss=loss, reglr=reglr, logReglr=T)
      return(curve)
    })
    names(curves_dags_reglr) <- reglrs
    return(curves_dags_reglr)
  })
  
  print("assembling opt kernParsX measures into dags")
  
  scores <- sapply(reglrs, function(reglr){
    # reglr <- reglrs[1]
    curves_dags_reglr <- lapply(curves_dags, function(el) el[[reglr]])
    res <- getComparableReglrVal(curveList=curves_dags_reglr)
    res <- cbind(reglr=res$logReglrs, loss=res$losses)
    rownames(res) <- names(curves_dags_reglr)
    return(res)
  }, simplify="array")
  
  names(dimnames(scores)) <- c("dag","regLoss","msr")
  scores <- melt(scores)
  scores <- cast(scores, dag~msr+regLoss, value="value")
  rownames(scores) <- scores$dag
  scores <- as.matrix(as.data.frame(scores[,2:ncol(scores)]))
  c.nms <- colnames(scores)
  indxChng <- grep("loss",c.nms)
  c.nms[indxChng] <- paste("L2", sapply(strsplit(c.nms[indxChng],"_"), function(el) el[1]), sep="_")
  colnames(scores) <- c.nms
  names(dimnames(scores)) <- c("dag","score") #"measures"
  
  
  
  scores[,c("KCDC_reglr", "KCMC_reglr")]
  colnames(scores) <- paste(colnames(scores), "comp", sep="_")
  scores <- scores[match(rownames(scores2), rownames(scores)),]
  
  
  scores <- cbind(scores2, scores)
  #scores <- scores3
  names(dimnames(scores)) <- c("dag", "score")
  
  return(scores)
}



plotCurveUAI_master <- function(x, hypArray, cmemLearner1, cmemLearner2, ppTab=NULL, plot=FALSE, dataNm,  folderSave,x0=NULL, xN=NULL, y0=NULL, yN=NULL){
  
  aux <- strsplit(cmemLearner1, "_")[[1]]
  indxDEL1 <- which(aux=="DEL1")
  indxDEL2 <- which(aux=="DEL2")
  cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
  #ini <- gregexpr("compCombo",cmemLearner_save)
  #cmemLearner1_save <- paste(substr(cmemLearner_save, 1, as.numeric(ini)-2), aux1, sep="_")
  #cmemLearner2_save <- paste(substr(cmemLearner_save, 1, as.numeric(ini)-2), aux2, sep="_")
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  
  print("calculating kernParsXY learners")
  # step 1 - regressions : get regression for each node in dags
  # as usual with CV on NCE for lambda, parsKernX, parsKernY
  
  cmeSet <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner=cmemLearner1, nm=paste(dataNm, cmemLearner1_save, "1",sep="_"), folderSave, round="first")
  
  sizeRegs <- lapply(uniqueRegsList, function(nodeReg) apply(nodeReg,2,sum))
  sizeRegs <- unlist(sizeRegs)
  aux <- strsplit(names(sizeRegs), "\\.")
  nodeIndx <- sapply(aux, function(el) el[1])
  regIndx <- sapply(aux, function(el) el[2])
  uniqueSizeRegs <- unlist(sizeRegs)
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  loss <- "cmem_L2_f"
  reglrs <- names(eval(parse(text=cmemLearner2))$optimizeParams$losses) 
  aux <- c("KCDC", "KCMC", "KCRDC", "KCSC", "KCNSC", "KCCC_ent","KCCC_pca_ent")
  reglrs <- intersect(reglrs, aux)
  
  
  print("building the appropriate y-hilbert space measures per size")
  cmemLearnersSum_preHyperParamSelec_perSize <- lapply(sort(unique(uniqueSizeRegs)), function(size){
    # size <- 1
    #print("*******")
    #print(paste("size: ", size))
    indx <- which(sizeRegs==size)
    
    # step 2 - choose Hy : only regressions with same number of 
    # regressors will be compared so we choose sum of Hys for
    # same size regressions. Actually best do it by node?
    # one Hy per node
    
    #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    
    # modifying to do a direct sum where the best parameter from all regressions
    # from the point of view that it doesn't make sense to add kcmcs or kcds from different
    # Hy spaces
    #cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
    # nodeIndx=nodeIndx[indx], regIndx=regIndx[indx])
    cmemLearnerSum <- directSumCMEMlearner_comp(cmemLearner=cmemLearner2, cmeSet, 
                                                nodeIndx=nodeIndx, regIndx=regIndx)
    
    # cmemLearnerSum$hyperParams$data$non_optimizable$sigma.rbf.Y
    cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc
    
    return(cmemLearnerSum)
  })
  names(cmemLearnersSum_preHyperParamSelec_perSize) <- sort(unique(uniqueSizeRegs))
  
  nodes <- names(uniqueRegsList)
  
  print("putting back into nodeTo-numReg structure")
  cmeSet2  <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      whatSize <- sizeRegs[match(paste(nodeTo, numReg, sep="."), names(sizeRegs))]
      cmemLearnerSumAux <- cmemLearnersSum_preHyperParamSelec_perSize[[as.character(whatSize)]]
      
      return(cmemLearnerSumAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  names(cmeSet2) <- nodes
  
  # step 3 - regression: get regressions for each node in dags
  # with CV on L2-loss for lambda and parsKernX
  print("calculating kernParsX learners")
  cmeSet2 <- getCMEsGivenDAGSet(uniqueRegsList, x, cmemLearner=cmeSet2, nm=paste(dataNm, cmemLearner2, "2",sep="_"), folderSave)
  
  
  
  print("obtain maximum efficiency curves for kernParsX learners")
  
  curves <- getCurves_comp(cmeLearnerSumList=cmeSet2, uniqueRegsList, loss, reglrs)
  
  dags <- hypArray
  numDags <- dim(dags)[3]
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag

  curves_dags <- lapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 2; dag <- dagsList[[i]]
    
    # obtain curves that will make up generla curve of dag hypothesis
    curves_components <- lapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[5])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        curve <- curves[[nodeTo]][[indxReg]]
        
      } else{
        curve <- NA
      }
      
      return(curve)
    })
    
    curves_components <- curves_components[which(!sapply(curves_components, function(el) is.na(el[1])))]
    curve_dag <- estimateOverallCurve(curves=curves_components, reglrs)
  
    return(curve_dag)
  })
  curvesPlot <- lapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 2; dag <- dagsList[[i]]
    
    
    # obtain curves that will make up generla curve of dag hypothesis
    curves_components <- lapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[5])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        curve <- curves[[nodeTo]][[indxReg]]
        
      } else{
        curve <- NA
      }
      
      return(curve)
    })
    
    curves_components <- curves_components[which(!sapply(curves_components, function(el) is.na(el[1])))]
    
    
    return(curves_components)
  })
  curvesPlot <- lapply(curvesPlot, function(el) el[[1]])
  reglr <- "KCMC"
  
  # step 4b - compare regressions: compare regressions of same size
  # (there should be at least two otherwise all dags in hypothesis
  # set have the same regression and we don't have to assess it)
  # get envelope curves for each set of regressions, compare to get 
  # measures 
  j <- which(reglr == reglrs)
  curves_reglr <- lapply(1:length(curvesPlot), function(i) curvesPlot[[i]][[j]])
  names(curves_reglr) <- names(curvesPlot)
  compReglrVals <- getComparableReglrVal(curveList=curves_reglr)
  plotCurvesUAI(compReglrVals, curveList=curves_reglr,x0=x0, xN=xN, y0=y0, yN=yN)
  
  
  
  return(NULL)
}

getMeasuresGivenDAGSet <- function(uniqueRegsList, x, cmemLearner, dataNm, folderSave){
  
  cmemLearner <- eval(parse(text=cmemLearner))
  
  
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  #count <- 0
  #pm0 <- proc.time()
  measureList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      if(length(indxPreds)==0){
        cmems <- rep(NA, length(cmemLearner$msrs))
        names(cmems) <- names(cmemLearner$msrs)
      } else{
      
      #count <- count + 1
      regressors <- nodes[indxPreds] 
      dataX <- x[,regressors, drop=F]
      dataY <- x[,nodeTo]
      
      trainDataO <- constructData(x=dataX, y=dataY)
      
      regressorsChar <- paste(regressors, collapse="-")
      regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
      fileSave <- paste(dataNm, "_", regressionChar, ".RData", sep="")
    
      if(!is.null(folderSave)){
        if(fileSave %in% dir(folderSave) ){
          load(file=paste(folderSave, fileSave, sep=""))
        } else{
          cmemLearnerAux <- setParams(learner=cmemLearner, trainData=trainDataO)
          save(cmemLearnerAux, file=paste(folderSave, fileSave, sep=""))  
        }
        }
      else{
        
        cmemLearnerAux <- setParams(learner=cmemLearner, trainData=trainDataO)
      }
      
      
      cmemLearnerAux <- cmemLearnerAux$learn(learner=cmemLearnerAux, forLoss=F)
      
      cmems <- cmemLearnerAux$calcMsrs(cmemLearner=cmemLearnerAux)
      }
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(cmems)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(measureList) <- nodes
  
  return(measureList)
}

getMeasuresGivenDAGSet_boot <- function(uniqueRegsList, x){
  
  
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  
  #count <- 0
  #pm0 <- proc.time()
  measureList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      if(length(indxPreds)==0){
        cnms <- c("corrRKHS","kcrdc","wkcrdc","l2rel","wl2rel","wwkcrdc","wVarCorrRKHS","wVarL2rels")
        cmems <- matrix(NA, nrow(x), cnms)
        names(cmems) <- cnms
      } else{
        
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        xTe <- X #[51:100,,drop=F]
        xTr <- X #[1:50,,drop=F]
        yTe <- as.matrix(Y)#[51:100,])
        yTr <- as.matrix(Y)#[1:50,])
        nTr <- nrow(xTr)
        nTe <- nrow(xTe)
        n <- nrow(X) 
        Itr <- diag(nTr)  
        I <- diag(n)  
        
        nboots <- 100
        set.seed(123)
        smpl_tr <- sapply(1:nboots, function(i) sample(nTr, replace=T), simplify="matrix")
        smpl_te <- sapply(1:nboots, function(i) sample(setdiff(1:n, unique(smpl_tr[,i])), nTe+1, replace=T), simplify="matrix")
        
        dists2 <- as.numeric(dist(X)^2)
        dists2 <- dists2[which(dists2>0)]
        
        sigma1 <- 1/quantile(dists2, 0.99)  
        sigma2 <- 1/quantile(dists2, 0.2)
        
        sigmaxSeq <- 10^seq(log(sigma1,10), log(sigma2,10), length.out=10)
        lambdaSeq <- 10^seq(-5,-1,1)
        
        sigmasy <- 10^seq(-2,8,1)
        
        #i <- 0
        #j <- 0
        print("get corrsRKHS_boot")
        pm <- proc.time()
        corrsRKHS_boot <- sapply(lambdaSeq, function(lam){ 
          #i <<- i + 1
          #j <<- 0
          sapply(sigmaxSeq , function(sigmax){
            # i <- 1; j <- 1; lam <- unique(lambda)[i]; sigmax <- unique(sigma.rbf.X)[j]
            # lam <- 0.1; sigmax <- 1
            #j <<- j + 1
            #print(paste("i-j: ", i, "-",j))
            
            
            
            difYs <- sapply(sigmasy, function(sigmay){
              K <- kern_rbf(yTr, sigma=sigmay)
              min(apply(K,2, function(col) length(unique(col))))
            })/nTr
            
            pm <- proc.time()
            corrsRKHS_te_boot <- sapply(1:nboots, function(i){
              #i <- 1
              
              xTe_b <- X[smpl_te[,i],,drop=F]
              xTr_b <- X[smpl_tr[,i],,drop=F]
              yTe_b <- as.matrix(Y[smpl_te[,i],])
              yTr_b <- as.matrix(Y[smpl_tr[,i],])
              nTr_b <- nrow(xTr_b)
              nTe_b <- nrow(xTe_b)
              Itr_b <- diag(nTr_b)  
              
              Ltr <- kern_rbf(xTr_b, sigma=sigmax)
              Lte <- kern_rbf(xTe_b, sigma=sigmax)
              Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
              Blambda_tr <- solve(Ltr+nTr*lam*Itr)
              if(FALSE){
                LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
                plot(c(Lte), c(LBK)) 
                plot(c(Lte[21,]), c(LBK[21,])) 
                cor(c(Lte[21,]), c(LBK[21,]))
              }
              
              corrsRKHS_te <- sapply(sigmasy, function(sigmay) wcorrRKHS_col(sigmay, B=Blambda_tr, Lte_tr=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b ), simplify="array")
              
              return(corrsRKHS_te)
            }, simplify="array")
            proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
            dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nTe+1), sigmay=sigmasy, boot=1:nboots)
            
            corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
            
            return(corrsRKHS_te_boot)
            
          }, simplify="array")}, simplify="array")
        print(proc.time() - pm) #21 mins
        
        
        print("get opt df")
        pm <- proc.time()
        dimnames(corrsRKHS_boot)[4:5] <- list(sigmax=sigmaxSeq, lambda=lambdaSeq)
        names(dimnames(corrsRKHS_boot))[4:5] <- c("sigmax", "lambda")
        bestSigmay <- apply(corrsRKHS_boot, c("testPt","boot","lambda","sigmax"), function(vec) which.max(vec))
        corrBestSigmay <- apply(corrsRKHS_boot, c("testPt","boot","lambda","sigmax"), function(vec) max(vec,na.rm=T))
        dfBestSigmay <- melt(bestSigmay)
        dfcorrBestSigmay <- melt(corrBestSigmay)
        dfBestSigmay <- cbind(dfBestSigmay, corr=dfcorrBestSigmay$value)
        colnames(dfBestSigmay)[5] <- "indxSigmay"
        dfBestSigmay$indxSigmax <- match(round(dfBestSigmay$sigmax,6), round(sigmaxSeq,6))
        corrsRKHS_boot_sigmay <- cast(dfBestSigmay,testPt~boot~lambda~sigmax, value="corr")
        bestSigmax <- apply(corrsRKHS_boot_sigmay, c("testPt","boot","lambda"), function(vec) which.max(vec))
        corrBestSigmax <- apply(corrsRKHS_boot_sigmay, c("testPt","boot","lambda"), function(vec) max(vec,na.rm=T))
        dfBestSigmax <- melt(bestSigmax)
        dfcorrBestSigmax <- melt(corrBestSigmax)
        dfBestSigmax <- cbind(dfBestSigmax, corr=dfcorrBestSigmax$value)
        head(dfBestSigmax)
        colnames(dfBestSigmax)[4] <- "indxSigmax"
        dfBestSigmax <- merge(dfBestSigmax, dfBestSigmay,by=c("testPt","boot","lambda","indxSigmax"),all.x=T)
        dfBestSigmax$indxLambda <- match(dfBestSigmax$lambda, lambdaSeq)
        corrsRKHS_boot_lambda <- cast(dfBestSigmax,testPt~boot~lambda, value="corr.x")
        bestLambda <- apply(corrsRKHS_boot_lambda, c("testPt","boot"), function(vec) which.max(vec))
        corrBestLambda <- apply(corrsRKHS_boot_lambda, c("testPt","boot"), function(vec) max(vec,na.rm=T))
        dfBestLambda <- melt(bestLambda)
        dfcorrBestLambda <- melt(corrBestLambda)
        dfBestLambda <- cbind(dfBestLambda, corr=dfcorrBestLambda$value)
        head(dfBestLambda)
        colnames(dfBestLambda)[3] <- "indxLambda"
        dfBestLambda <- merge(dfBestLambda, dfBestSigmax,by=c("testPt","boot","indxLambda"),all.x=T)
        dfBestLambda$sigmay <- sigmasy[dfBestLambda$indxSigmay]
        indxMat <- matrix(c(dfBestLambda$testPt,dfBestLambda$boot),nrow(dfBestLambda),2)
        dfBestLambda$y <- Y[smpl_te[indxMat]]
        dfBestLambda$x <- X[smpl_te[indxMat]]
        dim(unique(dfBestLambda[,c("testPt","boot","indxLambda","indxSigmax")]))
        print(proc.time()-pm)
        
        print("KCRDC, wKCRDC, l2rel, wL2rel")
        pm <- proc.time()
        
        msrs1 <- apply(dfBestLambda, 1, function(row){
          # row <- as.numeric(dfBestLambda[1,]); names(row) <- colnames(dfBestLambda) 
          lam <- row["lambda"]
          sigmax <- row["sigmax"]
          sigmay <- row["sigmay"]
          boot <- row["boot"]
          testPt <- row["testPt"]
          xTe_b <- X[smpl_te[,boot],,drop=F]
          xTr_b <- X[smpl_tr[,boot],,drop=F]
          yTr_b <- as.matrix(Y[smpl_tr[,boot],])
          yTe_b <- as.matrix(Y[smpl_te[,boot],])
          nTr_b <- nrow(xTr_b)
          Itr_b <- diag(nTr_b)  
          Ltr <- kern_rbf(xTr_b, sigma=sigmax)
          Lte <- kern_rbf(xTe_b, sigma=sigmax)
          Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
          Blambda_tr <- solve(Ltr+nTr*lam*Itr)
          Ktr <- kern_rbf(yTr_b, sigma=sigmay)
          Kte <- kern_rbf(yTe_b, sigma=sigmay)
          Kte_tr <- kern_rbf(yTe_b, yTr_b, sigma=sigmay)
          Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
          res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
          res <- res^0.5
          weights <- Lte[testPt,]
          weights <- weights/sum(weights)
          kcrdc <- res[testPt]/mean(res)
          wkcrdc <- res[testPt]/weighted.mean(res, w=weights)
          resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
          l2rel <- sqrt(resids[testPt]/mean(diag(Kte)^2))
          wl2rel <- sqrt(resids[testPt]/weighted.mean(diag(Kte)^2,w=weights))
          return(c(kcrdc=kcrdc, wkcrdc=wkcrdc, l2rel=l2rel, wl2rel=wl2rel))
        })
        msrs1 <- t(msrs1)
        
        print("wwKCRDC, wVarCorrRKHS, wVarl2rel")
        msrs2 <- sapply(1:nrow(dfBestLambda), function(i){
          # i <- 47
          row <- as.numeric(dfBestLambda[i,]); names(row) <- colnames(dfBestLambda)
          lam <- row["lambda"]
          sigmax <- row["sigmax"]
          sigmay <- row["sigmay"]
          boot <- row["boot"]
          testPt <- row["testPt"]
          xPt <- matrix(row["x"],1,1)
          xTr <- X
          Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
          weights <- Lte_tr[,match(dfBestLambda$x, xTr)]
          weights <- weights/sum(weights)
          # hist(weights)
          wwkcrdc <- weighted.var(msrs1[,"wkcrdc"],w=weights)
          mu <- weighted.mean(dfBestLambda$corr,w=weights)
          wVarCorrRKHS <- weighted.var(dfBestLambda$corr,w=weights)/mu
          wVarL2rels <- weighted.var(msrs1[,"wl2rel"], w=weights)
          res <- c(wwkcrdc=wwkcrdc, wVarCorrRKHS=wVarCorrRKHS, wVarL2rels=wVarL2rels)
          return(res)
        })
        msrs2 <- t(msrs2)
        
        res <- cbind(corrRKHS=dfBestLambda$corr, msrs1, msrs2)
        
        print(proc.time()-pm)
        
        
      
      }
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(res)
    })
    
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(measureList) <- nodes
  
  return(measureList)
}

featsCalc <- function(res, testPt, weights1, weights2, weights){
  # local value - 1
  normCME <- res[testPt]
  # global mean - 2
  normCME_smth <- mean(res)
  # local weighted mean (smooth the mean) - 3 weights - 3-5
  normCME_smth_loc <- weighted.mean(res, w=weights1)
  normCME_smth_ad <- weighted.mean(res, w=weights2)
  normCME_smth_loc_ad <- weighted.mean(res, w=weights)
  
  # relative (to global mean) local weighted mean, 3 weights- 6-8
  normCME_smth_rel_loc <- normCME_smth_loc/normCME_smth
  normCME_smth_rel_ad <- normCME_smth_ad/normCME_smth
  normCME_smth_rel_loc_ad <- normCME_smth_loc_ad/normCME_smth
  
  # relative local value (to global mean, and to local weightedmeans) 
  # , 4 different ways of relativizing - 9-12
  normCME_rel <- normCME/normCME_smth
  normCME_rel_loc <- normCME/normCME_smth_loc
  normCME_rel_ad <- normCME/normCME_smth_ad
  normCME_rel_loc_ad <- normCME/normCME_smth_loc_ad
  
  # global variance - 13
  kcdc <- var(res)
  
  # local weighted variance (smooth the var) - 3 weights - 14-16 
  kcdc_loc <- weighted.var(res, w=weights1)
  kcdc_ad <- weighted.var(res, w=weights2)
  kcdc_loc_ad <- weighted.var(res, w=weights)
  
  
  kcrdc <- kcdc/normCME_smth
  kcrdc_loc <- kcdc_loc/normCME_smth_loc
  kcrdc_ad <- kcdc_ad/normCME_smth_ad
  kcrdc_loc_ad <- kcdc_loc_ad/normCME_smth_loc_ad
  
  res <- c(normCME,normCME_smth,normCME_smth_loc,normCME_smth_ad,normCME_smth_loc_ad,
           normCME_smth_rel_loc,normCME_smth_rel_ad,normCME_smth_rel_loc_ad,
           normCME_rel,normCME_rel_loc,normCME_rel_ad,normCME_rel_loc_ad,kcdc,
           kcdc_loc,kcdc_ad,kcdc_loc_ad,kcrdc,kcrdc_loc,kcrdc_ad,kcrdc_loc_ad)
  
  names(res) <-  c("v_or","v_smth","v_smth_loc" ,"v_smth_dist"   ,"v_smth_loc_dist",
                                 "v_excSmth_loc","v_excSmth_dist","v_excSmth_loc_dist",
                   "v_exc",        "v_exc_loc","v_exc_dist","v_exc_loc_dist",
                   "v_var",        "v_var_loc","v_var_dist","v_var_loc_dist",
                   "v_var_rel", "v_var_rel_loc","v_var_rel_dist","v_var_rel_loc_dist")
  
  return(res)
}

reformat_msrs <- function(cmemSet, hypArray){
  print("enters reformat_msrs")
  
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  # we don't need dummy "regressions" with no dep variables
  uniqueRegsList <- lapply(uniqueRegsList, function(mat){
    # mat <- uniqueRegsList[[1]]
    #print(mat)
    # if its only regression we leave it so that we get a vector of NAs for all measures
    indx <- which(apply(mat, 2, function(col) all(col==0)))
    if(ncol(mat)>1 & length(indx)>0) mat <- mat[,-indx, drop=FALSE]
    colnames(mat) <- 1:ncol(mat)
    return(mat)
  })
  dags <- hypArray
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  testPts_x <- cmemSet$x$`1`[,c("testPt","boot","indxPt")]
  testPts_y <- cmemSet$y$`1`[,c("testPt","boot","indxPt")]
  head(testPts_x)
  head(testPts_y)
  all(testPts_x==testPts_y)
  testPts <- testPts_x[,c("testPt","boot","indxPt")]
  
  print("scores local")
  prnt <- F
  pm <- proc.time()
  scoresLocal <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 2; dag <- dagsList[[i]]
    
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[2])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      } else{
        indxReg <- 1
        correct <- NA
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      
      
      indxColsOwnSig <- grep(paste(".",nodeTo, ".",1, sep=""), colnames(cmemSet2))
      # length(indxColsOwnSig)
      cmemSet2best <- cmemSet2[,indxColsOwnSig]
      
      colnames(cmemSet2best) <- gsub(paste(".",nodeTo, ".",1, sep=""), "_best", colnames(cmemSet2best))
      
      cmemSet2 <- cbind(cmemSet2, cmemSet2best)
      indxCol <- which(colnames(cmemSet2) %in% c("testPt","boot","indxPt"))
      
      cmemSet3 <- cmemSet2[,-indxCol]
      
      fac <- as.matrix(cmemSet3)*correct
      
      
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    #aux <- melt(facs)
    
    
    pm <- proc.time()
    facs <- apply(facs, c(1,2), sum, na.rm=T)
    proc.time() - pm # 3 seconds for 1 core
    
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(facs)
  }, simplify="array")
  names(dimnames(scoresLocal))[3] <- "dag"
  print((proc.time()-pm)[3])
  
  
  #difference between the score of a dag-pt and the min of the rest of the dags
  print("scores local 2")
  pm <- proc.time()
  scoresLocal2  <- apply(scoresLocal, c(1,2), function(vec) vec[2]-vec[1])
  
  
  # here i need to take aggregate
  head(colnames(scoresLocal2))
  colnames(scoresLocal2)[grep("sigmay", colnames(scoresLocal2))]
  
  indx_x <- grep(".x.1", colnames(scoresLocal2))
  indx_y <- grep(".y.1", colnames(scoresLocal2))
  indx_best <- grep("best", colnames(scoresLocal2))
  indx_joint <- grep("joint", colnames(scoresLocal2))
  head(colnames(scoresLocal2)[indx_x])
  head(colnames(scoresLocal2)[indx_y])
  head(colnames(scoresLocal2)[indx_best])
  head(colnames(scoresLocal2)[indx_joint])
  length(indx_x); length(indx_y);length(indx_best); length(indx_joint)
  scoresLocal2_best <- scoresLocal2[,indx_best]
  scoresLocal2_joint <- scoresLocal2[,indx_joint]
  
  print((proc.time()-pm)[3])
  
  print("scores local 3")
  pm <- proc.time()
  cnms <- colnames(scoresLocal2)[indx_x]
  cnms <- gsub(".x.1","",cnms)
  scoresLocal3 <- abind(scoresLocal2[,indx_x],scoresLocal2[,indx_y], along=3)
  dimnames(scoresLocal3)[[2]] <- cnms
  dimnames(scoresLocal3)[[3]] <- c("x","y")
  names(dimnames(scoresLocal3)) <- c(names(dimnames(scoresLocal2)),"sigyType")
  print((proc.time()-pm)[3])
  
  print("scores local 4")
  pm <- proc.time()
  scoresLocal4 <- apply(scoresLocal3, c("pts", "measures"), function(msrs){
    msrs2 <- c(-Inf, msrs)
    msrs3 <- c(Inf, msrs)
    maxM <- sign(msrs2[which.max(msrs2)])*max(abs(msrs2),na.rm=T)
    minM <- sign(msrs3[which.min(msrs3)])*min(abs(msrs3),na.rm=T)
    res <- c(mean(msrs), maxM, minM, max(msrs)-min(msrs))
    names(res) <- c("mean", "max", "min", "range")
    return(res)
  })
  names(dimnames(scoresLocal4))[1] <- c("stat")
  print((proc.time()-pm)[3])
  
  if(FALSE){
    print("scores local 5")
    scoresLocal5 <- melt(scoresLocal4)
    print("scores local 6")
    scoresLocal6 <- cast(scoresLocal5, pts~measures+stat, value="value")
    scoresLocal6 <- scoresLocal6[,2:ncol(scoresLocal6)]
  } else{
    
    pm <- proc.time()
    scoresLocal5 <- aperm(scoresLocal4, c(2, setdiff(order(dim(scoresLocal4)),2)) )
    aux <- expand.grid(dimnames(scoresLocal4)[[setdiff(order(dim(scoresLocal4)),2)[1]]], 
                       dimnames(scoresLocal4)[[setdiff(order(dim(scoresLocal4)),2)[2]]])
    colnames(aux) <- names(dimnames(scoresLocal4))[setdiff(order(dim(scoresLocal4)),2)]
    scoresLocal6 <- scoresLocal5
    dim(scoresLocal6) <- c(dim(scoresLocal5)[1], dim(scoresLocal5)[2]*dim(scoresLocal5)[3])
    colnames(scoresLocal6) <- paste(aux$measures, aux$stat, sep="_")
    #all(sort(colnames(scoresLocal6b))==colnames(scoresLocal6))
    scoresLocal6 <- scoresLocal6[,order(colnames(scoresLocal6))]
    #all(colnames(scoresLocal6b)==colnames(scoresLocal6))
    #scoresLocal6[1:3,1:10]
    #scoresLocal6b[1:3,1:10]
    #summary(c(as.matrix(scoresLocal6)-scoresLocal6b))
    print((proc.time()-pm)[3])
    
  }
  
  pm <- proc.time()
  ncol(scoresLocal6)
  ncol(scoresLocal2_best)*4
  print("scores local 7")
  if(length(indx_joint)>0){
    scoresLocal7 <- cbind(scoresLocal2_best, scoresLocal2_joint, scoresLocal6)
  } else{
    scoresLocal7 <- cbind(scoresLocal2_best, scoresLocal6)
  }
  print((proc.time()-pm)[3])
  #modMatNew <- scoresLocal2#*respNew
  modMatNew <- scoresLocal7#*respNew
  print("exits reformat_msrs")
  return(list(modMatNew=modMatNew, testPts=testPts))
}

getMeasuresGivenDAGSet_boot_eqSig <- function(uniqueRegsList, x, jointFeats, smoothFeats){
  
  
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  # 100 on ERC
  # 50 on my lap
  nPerBoot <- 50
  
  repl <- FALSE
  if(n < 100) repl <- TRUE 
  # 100 or more on ERC
  # 10 on my lap
  #nboots <- max(100, floor(n/nPerBoot)*10)
  nboots <- 10
  RNGversion("3.5.0")
  set.seed(1234)
  smpl_tr <- sapply(1:nboots, function(i) sample(1:n,nPerBoot, replace=repl), simplify="matrix")
  smpl_te <- sapply(1:nboots, function(i) sample(setdiff(1:n, unique(smpl_tr[,i])), nPerBoot, replace=repl), simplify="matrix")
  #all(smpl_tr2==smpl_tr)
  #all(smpl_te2==smpl_te)
  
  lambdaSeq <- 10^seq(-5,-1,1)
  #sigmasy <- 10^seq(-2,8,1)
  # length.out 10 on ERC
  #            5 on lap
  sigmasy <- 10^seq(-1,3, length.out=5) 
  
  #plot(cmemSetRead$x[,"x"], x[,"x"])
  #plot(cmemSetRead$x[,"y"], x[,"y"])
  # o1 <- order(cmemSetRead$x[,"x"]); o2 <- order(x[,"x"])
  #plot(cmemSetRead$x[o1,"x"], x[o2,"x"])
  #plot(cmemSetRead$x[o1,"y"], x[o2,"y"])
  #plot(cmemSetRead$x[,"x"],cmemSetRead$x[,"y"])
  #plot(x[,"x"], x[,"y"])
  # x <- x[match(cmemSetRead$x[,"x"],x[,"x"]),]
  #plot(cmemSetRead$x[,"x"], x[,"x"])
  #plot(cmemSetRead$x[,"y"], x[,"y"])
  
  # 10 on ERC
  # 5 on lap
  numSigx <- 5
  sigmasx_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      
      #count <- count + 1
      regressors <- nodes[indxPreds] 
      dataX <- x[,regressors, drop=F]
      dataY <- x[,nodeTo]
      trainDataO <- constructData(x=dataX, y=dataY)
      X <- as.matrix(trainDataO$x)
      Y <- as.matrix(trainDataO$y)
      n <- nrow(X)
      dists2 <- as.numeric(dist(X)^2)
      dists2 <- dists2[which(dists2>0)]
      sigma1 <- 1/quantile(dists2, 0.99)  
      sigma2 <- 1/quantile(dists2, 0.2)
      # 10 on ERC
      # 5 on lap
      sigmaxSeq <- 10^seq(log(sigma1,10), log(sigma2,10), length.out=numSigx)
      
      return(sigmaxSeq)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  names(sigmasx_by_node) <- nodes
  Ls_by_node <- mcmapply(function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #j <- 0
        
        #mc_cores_loc <- detectCores()/2
        mc_cores_loc <- 1
        pm <- proc.time()
        Ls <- mcmapply(function(sigmax){
          #  j <- 1;  sigmax <- sigmaxSeq[j]
          #j <<- j + 1
          #print(paste("j: ","-",j))
          L <- kern_rbf(X, sigma=sigmax)
          return(L)
            
          }, sigmax=sigmaxSeq, SIMPLIFY=FALSE, mc.cores=mc_cores_loc)
        print(proc.time() - pm) #0.254 1 core
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(Ls)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)}, nodeTo=nodes, mc.cores=1, SIMPLIFY=FALSE)
  names(Ls_by_node) <- nodes
  
  pm <- proc.time()
  Bs_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #j <- 0
        pm <- proc.time()
        Bs <- lapply(lambdaSeq, function(lam) lapply(1:length(sigmaxSeq) , function(sigmax){
          #  i <- 1; lam <- lambdaSeq[i]; sigmax <- 1;  
          #j <<- j + 1
          #print(paste("j: ","-",j))
          L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
          Bs <- lapply(1:nboots, function(k){
            # k <- 1
            indxTr_b <- smpl_tr[,k]
            Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
            Itr_b <- diag(nPerBoot)
            Blambda_tr <- solve(Ltr+nPerBoot*lam*Itr_b)  
            return(Blambda_tr)
          })
          return(Bs)
          
        }))
        print(proc.time() - pm) # 0.3 secs
    
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(Bs)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  names(Bs_by_node) <- nodes
  proc.time() - pm # 1.321
  
  Ks_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        I <- diag(n)
        H <- I-matrix(1/n,n,n)
        #sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #j <- 0
        
        #pm <- proc.time()
        Ks <- lapply(sigmasy , function(sigmay){
          #  j <- 1;  sigmax <- sigmaxSeq[j]
          #j <<- j + 1
          #print(paste("j: ","-",j))
          K <- kern_rbf(Y, sigma=sigmay)
          K <- H %*% K %*% H
          return(K)
          
        })
        #print(proc.time() - pm) #21 mins
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(Ks)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  names(Ks_by_node) <- nodes
  # standardize y locally - weighted mean 0, weighted sd = 1
  # weights defined according to L and sigmx
  yNorm_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #j <- 0
        
        #pm <- proc.time()
        ysNorm <- lapply(1:length(sigmaxSeq) , function(indxSigmax){
          #  indxSigmax <- 1
          #j <<- j + 1
          #print(paste("j: ","-",j))
          L <- Ls_by_node[[nodeTo]][[numReg]][[indxSigmax]]
          weights <- apply(L, 2, function(col) col/sum(col))
          # table(apply(weights, 2, sum))
          ysNorm <- apply(weights, 2, function(col) (Y-weighted.mean(Y, w=col))/sqrt(weighted.var(Y, w=col)))
          # j <- 2; weighted.mean(ysNorm[,j], w=weights[,j]); weighted.var(ysNorm[,j], w=weights[,j])
          return(ysNorm)
          
        })
        #print(proc.time() - pm) #21 mins
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(ysNorm)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  names(yNorm_by_node) <- nodes
  KsNorm_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        I <- diag(n)
        H <- I-matrix(1/n,n,n)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #j <- 0
        #pm <- proc.time()
        KsNorm <- lapply(1:length(sigmaxSeq) , function(indxSigmax){
          #  indxSigmax <- 3
          #j <<- j + 1
          #print(paste("indxSigmax: ","-",indxSigmax))
          ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[indxSigmax]]
          
          Ks <- lapply(sigmasy , function(sigmay){
            #  j <- 1;  sigmax <- sigmaxSeq[j]
            #j <<- j + 1
            #print(paste("j: ","-",j))
            K <- kern_rbf(ysNorm, sigma=sigmay)
            K <- H %*% K %*% H
            return(K)
            
          })
          
          return(Ks)
          
        })
        #print(proc.time() - pm) #21 mins
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(KsNorm)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  names(KsNorm_by_node) <- nodes
  
  #count <- 0 
  #pm0 <- proc.time()
  if(FALSE){
  dfs_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #i <- 0
        #j <- 0
        print("get corrsRKHS_boot")
        pm <- proc.time()
        corrsRKHS_boot <- sapply(lambdaSeq, function(lam){ 
          #i <<- i + 1
          #j <<- 0
          sapply(sigmaxSeq , function(sigmax){
            # i <- 1; j <- 1; lam <- lambdaSeq[i]; sigmax <- sigmaxSeq[j]
            # lam <- 0.1; sigmax <- 1
            #j <<- j + 1
            #print(paste("i-j: ", i, "-",j))
            
            
            
            difYs <- sapply(sigmasy, function(sigmay){
              K <- kern_rbf(Y, sigma=sigmay)
              min(apply(K,2, function(col) length(unique(col))))
            })/n
            
            pm <- proc.time()
            corrsRKHS_te_boot <- sapply(1:nboots, function(i){
              #i <- 1
              
              xTe_b <- X[smpl_te[,i],,drop=F]
              xTr_b <- X[smpl_tr[,i],,drop=F]
              yTe_b <- as.matrix(Y[smpl_te[,i],])
              yTr_b <- as.matrix(Y[smpl_tr[,i],])
              nTr_b <- nrow(xTr_b)
              nTe_b <- nrow(xTe_b)
              Itr_b <- diag(nTr_b)  
              
              Ltr <- kern_rbf(xTr_b, sigma=sigmax)
              Lte <- kern_rbf(xTe_b, sigma=sigmax)
              Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
              Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
              if(FALSE){
                LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
                plot(c(Lte), c(LBK)) 
                plot(c(Lte[21,]), c(LBK[21,])) 
                cor(c(Lte[21,]), c(LBK[21,]))
              }
              
              corrsRKHS_te <- sapply(sigmasy, function(sigmay) wcorrRKHS_col(sigmay, B=Blambda_tr, Lte_tr=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b ), simplify="array")
              
              return(corrsRKHS_te)
            }, simplify="array")
            proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
            dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nPerBoot), sigmay=sigmasy, boot=1:nboots)
            
            #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
            
            return(corrsRKHS_te_boot)
            
          }, simplify="array")}, simplify="array")
        print(proc.time() - pm) #21 mins
        #all(corrsRKHS_boot==corrsRKHS_boot_xy)
        #all(corrsRKHS_boot==corrsRKHS_boot_yx)
        # smplPlt <- sample(prod(dim(corrsRKHS_boot)),1000)
        #plot(c(corrsRKHS_boot)[smplPlt], c(corrsRKHS_boot_xy)[smplPlt])
        #plot(c(corrsRKHS_boot)[smplPlt], c(corrsRKHS_boot_yx)[smplPlt])
        
        print("get opt df")
        pm <- proc.time()
        dimnames(corrsRKHS_boot)[4:5] <- list(sigmax=sigmaxSeq, lambda=lambdaSeq)
        names(dimnames(corrsRKHS_boot))[4:5] <- c("sigmax", "lambda")
        
        indxOpt <- apply(corrsRKHS_boot, c("testPt","boot"), function(arr){
          # arr <- corrsRKHS_boot[1,,1,,]
          maxArr <- arr==max(arr,na.rm=T)
          indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
          #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
          return(indxs)
        })
        dimnames(indxOpt)[[1]] <- c("sigmay","sigmax","lambda")
        names(dimnames(indxOpt))[1] <- "indx"
        # all(indxOpt==indxOpt_xy)
        # all(indxOpt==indxOpt_yx)
        df <- melt(indxOpt)
        df <- cast(df, testPt+boot~indx, value="value")
        df$corr <- corrsRKHS_boot[as.matrix(df[,names(dimnames(corrsRKHS_boot))])]
        colnames(df) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
        df$node <- nodeTo
        df$numReg <- numReg
        df$lambda <- lambdaSeq[df$indxLambda]
        df$sigmax <- sigmaxSeq[df$indxSigmax]
        df$sigmay <- sigmasy[df$indxSigmay]
        df <- as.data.frame(df)
        #all(df==df_xy[which(df_xy$node=="y"),colnames(df)])
        #all(df==df_yx[which(df_yx$node=="x"),colnames(df)])
        
        print(proc.time()-pm)
        
        res <- list(df=df, corrsRKHS_boot=corrsRKHS_boot, indxOpt=indxOpt, sigmaxSeq=sigmaxSeq)
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(res)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  #print(proc.time()-pm0) #  
  names(dfs_by_node) <- nodes
  }
  
  print("dfs_by_node")
  # obtain a  dataframe  of size nboots x nPtsPerBoot x 2 x 2 where for each boot_pt we have 
  # the optimal sigma_x, sigma_y and lambda according to locally (according to pt and diff
  # sigma_x considered) weighted correlation of similarities. We take two diff sigma_ys
  # the raw arg_max for correlation and a loess smoothed (in sigmay~x) one. We also use
  # two yTypes to calculate Ks, the raw 0-1 normalized one, and the locally standardized
  # one meaning for any given m x n boot_pts we have 4 optimal arg_max and max values. 
  pm <- proc.time()
  dfs_by_node <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        #i <- 0
        #j <- 0
        print("get corrsRKHS_boot")
        # For every bootstrap, lambda, sigmax, sigmay, calculate locally weighted corr
        # (one for every pt in the bootstrap) using both raw y and locally standardized y 
        pm <- proc.time()
        corrsRKHS_boot <- sapply(1:length(lambdaSeq), function(lam){ 
          #i <<- i + 1
          #j <<- 0
          sapply(1:length(sigmaxSeq) , function(sigmax){
            # i <- 1; j <- 1; lam <- i; sigmax <- j
            #j <<- j + 1
            #print(paste("i-j: ", i, "-",j))
            
            L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
            
            mc_cores_loc <- detectCores()/2
            pm <- proc.time()
            corrsRKHS_te_boot <- mcmapply(function(i){
              #i <- 1
              
              indxTr_b <- smpl_tr[,i]
              indxTe_b <- smpl_te[,i]
              Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
              Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
              Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
              Itr_b <- diag(nPerBoot)
              Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[i]]
              
              if(FALSE){
                LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
                plot(c(Lte), c(LBK)) 
                plot(c(Lte[21,]), c(LBK[21,])) 
                cor(c(Lte[21,]), c(LBK[21,]))
              }
              
              LB <- Lte_tr%*%Blambda_tr
              
              # calculate corrs and corrsNorm for all sigmays and sigmax, lambda and bootstrap from above
              # corrs is the weighted correlation using Ks_by_node- based on y
              # corrsNorm is the weigthed correlation using  KsNorm_by_node - basedo on locally standardized y 
              
              corrsRKHS_te <- sapply(1:length(sigmasy), function(sigmay){
                # sigmay <- 4
                K <- Ks_by_node[[nodeTo]][[numReg]][[sigmay]]
                Ktr_te <- K[indxTr_b,]; Ktr_te <- Ktr_te[,indxTe_b]
                LBK <- LB%*%Ktr_te
                
                mc_cores_loc2 <- 1#detectCores()/2
                corrs <- mcmapply(function(Lte_col,LBK_col){
                    res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                    return(res)
                  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=mc_cores_loc2)
                  
                K <- KsNorm_by_node[[nodeTo]][[numReg]][[sigmax]][[sigmay]]
                Ktr_te <- K[indxTr_b,]; Ktr_te <- Ktr_te[,indxTe_b]
                LBK <- LB%*%Ktr_te
                corrsNorm <- mcmapply(function(Lte_col,LBK_col){
                    res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                    return(res)
                  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=mc_cores_loc2)
                  
                corrs <- cbind(corr=corrs, corrNorm=corrsNorm)
                  
                  return(corrs)
                }, simplify="array")
                
               
              
              return(corrsRKHS_te)
            }, i=1:nboots, mc.cores=mc_cores_loc, SIMPLIFY="array")
            proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
            dim(corrsRKHS_te_boot)
            dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nPerBoot), yType=c("raw","nrm"),sigmay=sigmasy, boot=1:nboots)
            
            #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
            
            return(corrsRKHS_te_boot)
            
          }, simplify="array")}, simplify="array")
        print(proc.time() - pm) #21 mins
        #all(corrsRKHS_boot==corrsRKHS_boot_xy)
        #all(corrsRKHS_boot==corrsRKHS_boot_yx)
        # smplPlt <- sample(prod(dim(corrsRKHS_boot)),1000)
        #plot(c(corrsRKHS_boot)[smplPlt], c(corrsRKHS_boot_xy)[smplPlt])
        #plot(c(corrsRKHS_boot)[smplPlt], c(corrsRKHS_boot_yx)[smplPlt])
        
        
        pm <- proc.time()
        dimnames(corrsRKHS_boot)[5:6] <- list(sigmax=sigmaxSeq, lambda=lambdaSeq)
        names(dimnames(corrsRKHS_boot))[5:6] <- c("sigmax", "lambda")
        
        # alternative way of optimizing sigmax, sigmay and lambda
        # based not only on maximizing corr but looking at its variance
        if(FALSE){
        dfTest <- cmemSet$dfs_by_node$y[[1]]$corrsRKHS_boot
        dfTest <- melt(dfTest)
        dfTest <- dfTest[which(dfTest$yType=="raw"),]
        indxMat <- dfTest[,c("testPt","boot")]
        dfTest$indxPt <- cmemSet$smpl_te[as.matrix(indxMat)]
        dfTest_mean <- cast(dfTest, indxPt+sigmax+lambda+sigmay~., value="value", fun.aggregate="mean")
        dfTest_sd <- cast(dfTest, indxPt+sigmax+lambda+sigmay~., value="value", fun.aggregate="sd")
        colnames(dfTest_mean)[5] <- "meanCorr"
        colnames(dfTest_sd)[5] <- "sdCorr"
        dfTest2 <- cbind(dfTest_mean, dfTest_sd[,"sdCorr"])
        colnames(dfTest2)[6] <- "sdCorr"
        
        pt <- 993
        p <- ggplot(dfTest2[which(dfTest2$indxPt==pt),])
        p <- p+  geom_point(aes(x=meanCorr, y=sdCorr, colour=factor(sigmax)))
        p
        
        dfTestPt <- dfTest2[which(dfTest2$indxPt==pt),]
        
        # angle between two pts and horizontal cos(theta)= a/h
        # where a is diff in x (mean_corr) and h is dist
        h <- sqrt(matNorm2(as.matrix(dfTestPt[,c("meanCorr","sdCorr")])))
        a <- sqrt(matNorm2(as.matrix(dfTestPt[,c("meanCorr")])))
        
        a <- sapply(dfTestPt[,c("meanCorr")], function(i) sapply(dfTestPt[,c("meanCorr")], function(j) i-j))
        o <- sapply(dfTestPt[,c("sdCorr")], function(i) sapply(dfTestPt[,c("sdCorr")], function(j) i-j))
        
        ang <- sign(o)*acos(a/h)
        summary(c(ang))
        
        
        
        ptRef <- 100
        dfTestPt$ang <- ang[,ptRef]#<0&ang[,ptRef]>-pi/2
        dfTestPt$ang2 <- findInterval(dfTestPt$ang, c(-pi, -pi/2, 0, pi/2,  pi))
        dfTestPt$ang3 <- ang[,ptRef]<pi&ang[,ptRef]>pi/2
        p <- ggplot()
        p <- p+  geom_point(aes(x=meanCorr, y=sdCorr, colour=factor(ang3)), data=dfTestPt)
        p <- p+  geom_point(aes(x=meanCorr, y=sdCorr), colour="blue", size=2, data=dfTestPt[ptRef,])
        p
        
        indxKeep <- apply(ang, 2, function(col) sum(col<=pi & col>=pi/2, na.rm=T))
        dfTestPt$inferior <- TRUE
        dfTestPt$inferior[indxKeep==0] <- FALSE
        
        p <- ggplot(dfTestPt)
        p <- p+  geom_point(aes(x=meanCorr, y=sdCorr, colour=inferior))
        p
        
        dfTestPt <- dfTestPt[-which(dfTestPt$inferior),]
        p <- ggplot(dfTestPt)
        p <- p+  geom_point(aes(x=meanCorr, y=sdCorr, colour=factor(sigmax)))
        p
        
        
        
        # lets see how many efficient par-combos this would generate 
        # for all pts
        
        numNotInf <- sapply(unique(dfTest2$indxPt), function(pt){ 
          # pt <- unique(dfTest2$indxPt)[1]
          print(paste("pt: ", pt))
          dfTestPt <- dfTest2[which(dfTest2$indxPt==pt),]
          h <- sqrt(matNorm2(as.matrix(dfTestPt[,c("meanCorr","sdCorr")])))
          a <- sapply(dfTestPt[,c("meanCorr")], function(i) sapply(dfTestPt[,c("meanCorr")], function(j) i-j))
          o <- sapply(dfTestPt[,c("sdCorr")], function(i) sapply(dfTestPt[,c("sdCorr")], function(j) i-j))
          ang <- sign(o)*acos(a/h)
          indxKeep <- apply(ang, 2, function(col) sum(col<=pi & col>=pi/2, na.rm=T))
          dfTestPt$inferior <- TRUE
          dfTestPt$inferior[indxKeep==0] <- FALSE
          indxDel <- which(dfTestPt$inferior)
          if(length(indxDel)>0) dfTestPt <- dfTestPt[-indxDel,]
          return(nrow(dfTestPt))
        })
        range(numNotInf)
        hist(numNotInf)
        table(numNotInf)
        sum(numNotInf==35)
        which(numNotInf==35)[1]
        sum(numNotInf); mean(numNotInf)
        
        dfTest3 <- lapply(unique(dfTest2$indxPt), function(pt){ 
          # pt <- unique(dfTest2$indxPt)[1]
          print(paste("pt: ", pt))
          dfTestPt <- dfTest2[which(dfTest2$indxPt==pt),]
          h <- sqrt(matNorm2(as.matrix(dfTestPt[,c("meanCorr","sdCorr")])))
          a <- sapply(dfTestPt[,c("meanCorr")], function(i) sapply(dfTestPt[,c("meanCorr")], function(j) i-j))
          o <- sapply(dfTestPt[,c("sdCorr")], function(i) sapply(dfTestPt[,c("sdCorr")], function(j) i-j))
          ang <- sign(o)*acos(a/h)
          indxKeep <- apply(ang, 2, function(col) sum(col<=pi & col>=pi/2, na.rm=T))
          dfTestPt$inferior <- TRUE
          dfTestPt$inferior[indxKeep==0] <- FALSE
          indxDel <- which(dfTestPt$inferior)
          if(length(indxDel)>0) dfTestPt <- dfTestPt[-indxDel,]
          return(dfTestPt)
        })
        dfTest3 <- do.call(rbind, dfTest3)
        
        dfTest4 <- cbind(dfTest3, cmemSet$x[dfTest3$indxPt,])
        
        p <- ggplot(dfTest4)
        p <- p + geom_point(aes(x=x, y=y,colour=factor(sigmax)), size=0.5, alpha=0.3)
        p
        p <- ggplot(dfTest4)
        p <- p + geom_point(aes(x=x, y=y,colour=factor(sigmay)), size=0.5, alpha=0.3)
        p
        
        dfTest5 <- cmemSet$dfs_by_node$y[[1]]$df
        dfTest5 <- dfTest5[which(dfTest5$yType=="raw"),]
        indxMat <- dfTest5[,c("testPt","boot")]
        dfTest5$indxPt <- cmemSet$smpl_te[as.matrix(indxMat)] 
        dfTest6 <- cbind(dfTest5, cmemSet$x[dfTest5$indxPt,])
        
        p <- ggplot(dfTest6)
        p <- p + geom_point(aes(x=x, y=y,colour=factor(sigmax)), size=0.5, alpha=0.3)
        p
        p <- ggplot(dfTest6)
        p <- p + geom_point(aes(x=x, y=y,colour=factor(sigmay)), size=0.5, alpha=0.3)
        p
        
        }
        
        print("get opt df")
        
        # for each bootstrap, testPtPerBootstrap and yType (raw, locally standardized)
        # obtaind the corrRKHS optimal index of sigma_x, sigma_y and lambda
        indxOpt <- apply(corrsRKHS_boot, c("testPt","boot","yType"), function(arr){
          # arr <- corrsRKHS_boot[1,,1,,]
          maxArr <- arr==max(arr,na.rm=T)
          indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
          #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
          return(indxs)
        })
        dimnames(indxOpt)[[1]] <- c("sigmay","sigmax","lambda")
        names(dimnames(indxOpt))[1] <- "indx"
        # all(indxOpt==indxOpt_xy)
        # all(indxOpt==indxOpt_yx)
        df <- melt(indxOpt)
        df <- cast(df, testPt+boot+yType~indx, value="value")
        indxMat <- df[,names(dimnames(corrsRKHS_boot))]
        indxMat$yType <- c(1,2)[match(indxMat$yType, c("raw","nrm"))]
        df$corr <- corrsRKHS_boot[as.matrix(indxMat)]
        colnames(df) <- c("testPt","boot","yType","indxLambda","indxSigmax","indxSigmay","corr")
        df$node <- nodeTo
        df$numReg <- numReg
        df$principal <- TRUE
        df$lambda <- lambdaSeq[df$indxLambda]
        df$sigmax <- sigmaxSeq[df$indxSigmax]
        df$sigmay <- sigmasy[df$indxSigmay]
        df <- as.data.frame(df)
        df$sigyType <- "rough"
        #all(df==df_xy[which(df_xy$node=="y"),colnames(df)])
        #all(df==df_yx[which(df_yx$node=="x"),colnames(df)])
        
        res <- list(df=df,  corrsRKHS_boot=corrsRKHS_boot,
                    indxOpt=indxOpt, sigmaxSeq=sigmaxSeq)
        
        # smooth - for non-normalized and normalized
        if(smoothFeats){ 
        print("smooth sig y")
        # For every pt of the original N pts obtain a mean sig_y assigned for
        # both types of yType used (raw and locally standardized)
        optSigmay <- sapply(c("raw","nrm"), function(type){
          # type <- unique(df$yType)[2]
          indxType <- which(df$yType==type)
          indxMat <- matrix(c(df$testPt[indxType],df$boot[indxType]),length(indxType),2)
          indxPt <- smpl_te[indxMat]
          #plot(X[indxPt], log(df$sigmay[indxType],10))
          #o <- order(X)
          mod <- sapply(1:length(X), function(pt) mean(log(df$sigmay[indxType][which(indxPt==pt)],10)))
          mod[is.nan(mod)] <- mean(mod, na.rm=T)
          #lines(X[o,], mod[o], col="red")
          mod_smth <- loess(mod~X, span=0.25)
          optSigmay <- predict(mod_smth, X)
          #lines(X[o,], optSigmay[o], col="green")
          #plot(X,Y)
          #lines(X[o,], norml(mod[o]), col="red")
          #lines(X[o,], norml(optSigmay[o]), col="green")
          optSigmay <- 10^optSigmay
          #plot(X, log(optSigmay,10))
          return(optSigmay)
        }, simplify="array")
        
        
        #i <- 0
        #j <- 0
        print("get corrsRKHS_boot_smooth")
        pm <- proc.time()
        corrsRKHS_boot_smooth <- sapply(1:length(lambdaSeq), function(lam){ 
          #i <<- i + 1
          #j <<- 0
          sapply(1:length(sigmaxSeq) , function(sigmax){
            # i <- 1; j <- 2; lam <- i; sigmax <- j
            
            
            L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
            ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[sigmax]]
            
            #weighted.mean(ysNorm[,3], w=L[,3]/sum(L[,3])
            #weighted.var(ysNorm[,5], w=L[,5]/sum(L[,5])
            
            pm <- proc.time()
            corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
              #k <- 1
              
              optSigmay_b <- optSigmay[smpl_te[,k],1]
              optSigmay_norm_b <- optSigmay[smpl_te[,k],2]
              
              indxTr_b <- smpl_tr[,k]
              indxTe_b <- smpl_te[,k]
              Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
              Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
              Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
              Itr_b <- diag(nPerBoot)
              Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[k]]
              
              
              LB <- Lte_tr%*%Blambda_tr
              I <- diag(nPerBoot)
              H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
              
              yTr_b <- as.matrix(Y[smpl_tr[,k]])
              yTe_b <- as.matrix(Y[smpl_te[,k]])
              
              ysTrNorm_b <- ysNorm[smpl_tr[,k],]; ysTrNorm_b <- ysTrNorm_b[,smpl_tr[,k]] 
              ysTeNorm_b <- ysNorm[smpl_te[,k],]; ysTeNorm_b <- ysTeNorm_b[,smpl_te[,k]] 
              # l <- 4; weighted.mean(ysNorm[,l], w=L[,l]/sum(L[,l]))
              
              corrs <- mcmapply(function(sigmay, pt){ 
                # l <- 1; sigmay <- optSigmay_b[l]; pt <- l
                
                
                Ktr_te <- kern_rbf(yTr_b, yTe_b, sigma=sigmay)
                Ktr_te <- H %*% Ktr_te %*% H
                LBK <- LB%*%Ktr_te
                Lte_col <- Lte[pt,]
                LBK_col <- LBK[pt,]
                res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                return(res)
              }, sigmay=optSigmay_b, pt=1:length(optSigmay_b), SIMPLIFY="array")
              
              corrsNorm <- mcmapply(function(sigmay, pt){ 
                # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
                
                
                Ktr_te <- kern_rbf(ysTrNorm_b[,pt,drop=F], ysTeNorm_b[,pt,drop=F], sigma=sigmay)
                Ktr_te <- H %*% Ktr_te %*% H
                LBK <- LB%*%Ktr_te
                Lte_col <- Lte[pt,]
                LBK_col <- LBK[pt,]
                res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                return(res)
              }, sigmay=optSigmay_norm_b, pt=1:length(optSigmay_norm_b), SIMPLIFY="array")
              
              corrs <- cbind(corr=corrs, corrNorm=corrsNorm)
              
              return(corrs)
            }, simplify="array")
            proc.time() - pm # 17 secs
            
            #dim(corrsRKHS_te_boot)
            dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nPerBoot),yType=c("raw","nrm") ,boot=1:nboots)
            
            
            return(corrsRKHS_te_boot)
            
          }, simplify="array")}, simplify="array")
        proc.time() - pm # 15 mins mins
        dim(corrsRKHS_boot_smooth)
        dimnames(corrsRKHS_boot_smooth)[4:5] <- list(sigmaxSeq, lambdaSeq)
        names(dimnames(corrsRKHS_boot_smooth))[4:5] <- c("sigmax", "lambda")
        print(proc.time()-pm)
        
        print("get opt df smooth")
        
        indxOptSmooth <- apply(corrsRKHS_boot_smooth, c("testPt","boot","yType"), function(arr){
          # arr <- corrsRKHS_boot_xy[1,1,,]
          maxArr <- arr==max(arr,na.rm=T)
          indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
          #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
          return(indxs)
        })
        dimnames(indxOptSmooth)[[1]] <- c("sigmax","lambda")
        names(dimnames(indxOptSmooth))[1] <- "indx"
        
        
        df_smooth <- melt(indxOptSmooth)
        df_smooth <- cast(df_smooth, testPt+boot+yType~indx, value="value")
        indxMat <- df_smooth[,names(dimnames(corrsRKHS_boot_smooth))]
        indxMat$yType <- c(1,2)[match(indxMat$yType, c("raw","nrm"))]
        df_smooth$corr <- corrsRKHS_boot_smooth[as.matrix(indxMat)]
        colnames(df_smooth) <- c("testPt","boot","yType","indxLambda","indxSigmax","corr")
        df_smooth$node <- nodeTo
        df_smooth$numReg <- numReg
        df_smooth$principal <- TRUE
        df_smooth$lambda <- lambdaSeq[df_smooth$indxLambda]
        df_smooth$sigmax <- sigmaxSeq[df_smooth$indxSigmax]
        indxMat <- matrix(c(df_smooth$testPt,df_smooth$yType),nrow(df_smooth),2)
        df_smooth$sigmay <- optSigmay[indxMat]
        df_smooth <- as.data.frame(df_smooth)
        df_smooth$sigyType <- "smooth"
        
        res <- c(res, list(df_smooth=df_smooth, corrsRKHS_boot_smooth=corrsRKHS_boot_smooth,
                           indxOptSmooth=indxOptSmooth, optSigmay=optSigmay))
        }
        
        
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(res)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  print((proc.time()-pm)[3]) #  
  names(dfs_by_node) <- nodes
  
  # all(dfs_by_node$y$`1`$df ==df_xy[which(df_xy$node=="y"),colnames(dfs_by_node$y$`1`$df)])
  # all(dfs_by_node$x$`1`$df ==df_yx[which(df_xy$node=="x"),colnames(dfs_by_node$x$`1`$df)])
  # all(dfs_by_node$y$`1`$corrsRKHS_boot==corrsRKHS_boot_xy)
  # all(dfs_by_node$x$`1`$corrsRKHS_boot==corrsRKHS_boot_yx)
  # all(dfs_by_node$y$`1`$indxOpt==indxOpt_xy)
  # all(dfs_by_node$x$`1`$indxOpt==indxOpt_yx)
  # all(dfs_by_node$y$`1`$sigmaxSeq==sigmaxSeq_xy)
  # all(dfs_by_node$x$`1`$sigmaxSeq==sigmaxSeq_yx)
  
  #all(dfs_by_node$x$`1`$sigmaxSeq==cmemSetRead$dfs_by_node$x$`1`$sigmaxSeq)
  #max(abs(dfs_by_node$x$`1`$corrsRKHS_boot-cmemSetRead$dfs_by_node$x$`1`$corrsRKHS_boot))
  #plot(dfs_by_node$x$`1`$corrsRKHS_boot,cmemSetRead$dfs_by_node$x$`1`$corrsRKHS_boot)
  
  if(jointFeats){
  print("dfs_by_node joint- stack RKHS (sum kerns) and stack feat space (product kerns)")
  # sum -> stack non-linear rbf features
  # prod -> stack feautres first then calculate non-linear rbf features
  indxOptJoint <- unlist(dfs_by_node, recursive=FALSE)
  indxOptJoint <- lapply(indxOptJoint, function(el) el$indxOpt)
  indxOptJoint <- do.call(abind, c(indxOptJoint, along=5))
  dimnames(indxOptJoint)[[5]] <- seq(dim(indxOptJoint)[5])
  names(dimnames(indxOptJoint)) <- c(names(dimnames(dfs_by_node[[1]][[1]]$indxOpt)), "numSigy")
  
  print("df_joints")
  dfs_by_node_joint <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      {
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]
        
        trainDataO <- constructData(x=dataX, y=dataY)
        
        X <- as.matrix(trainDataO$x)
        Y <- as.matrix(trainDataO$y)
        n <- nrow(X)
        
        sigmaxSeq <- sigmasx_by_node[[nodeTo]][[numReg]]
        #sigmaxSeq; sigmaxSeq_xy
        #sigmaxSeq; sigmaxSeq_yx
        
        
        
        #i <- 0
        #j <- 0
        print("get corrsRKHS_boot_joint")
        pm <- proc.time()
        corrsRKHS_boot_joint <- sapply(1:length(lambdaSeq), function(lam){ 
          #i <<- i + 1
          #j <<- 0
          sapply(1:length(sigmaxSeq) , function(sigmax){
            # i <- 1; j <- 2; lam <- i; sigmax <- j
            
            
            L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
            ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[sigmax]]
            
            #weighted.mean(ysNorm[,3], w=L[,3]/sum(L[,3])
            #weighted.var(ysNorm[,5], w=L[,5]/sum(L[,5])
            
            pm <- proc.time()
            corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
              #k <- 1
              
              # indxOptJoint dims: 
              #c("indx","testPt","boot","yType","numSigy")
              indxSigy_b <- indxOptJoint["sigmay",,k,"raw",]
              indxSigy_norm_b <- indxOptJoint["sigmay",,k,"nrm",]
              
              
              indxTr_b <- smpl_tr[,k]
              indxTe_b <- smpl_te[,k]
              Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
              Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
              Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
              Itr_b <- diag(nPerBoot)
              Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[k]]
              
              
              LB <- Lte_tr%*%Blambda_tr
              I <- diag(nPerBoot)
              H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
              
              yTr_b <- as.matrix(Y[smpl_tr[,k]])
              yTe_b <- as.matrix(Y[smpl_te[,k]])
              
              ysTrNorm_b <- ysNorm[smpl_tr[,k],]; ysTrNorm_b <- ysTrNorm_b[,smpl_tr[,k]] 
              ysTeNorm_b <- ysNorm[smpl_te[,k],]; ysTeNorm_b <- ysTeNorm_b[,smpl_te[,k]] 
              # l <- 4; weighted.mean(ysNorm[,l], w=L[,l]/sum(L[,l]))
              
              corrs <- mcmapply(function(indxy, pt){ 
                # l <- 1; indxy <- indxSigy_b[l,]; pt <- l
                
                sigsy <- unique(sigmasy[indxy])
                Ktr_te <- sapply(sigsy, function(sigy){
                  Ktr_te <- kern_rbf(yTr_b, yTe_b, sigma=sigy)
                  }, simplify="array")
                Ktr_te_sum <- apply(Ktr_te, c(1,2), sum)
                Ktr_te_prod <- apply(Ktr_te, c(1,2), prod)
                Ktr_te_sum <- H %*% Ktr_te_sum  %*% H
                Ktr_te_prod <- H %*% Ktr_te_prod  %*% H
               
                
                LBK_sum <- LB%*%Ktr_te_sum
                LBK_prod <- LB%*%Ktr_te_prod
                Lte_col <- Lte[pt,]
                LBK_sum_col <- LBK_sum[pt,]
                LBK_prod_col <- LBK_prod[pt,]
                res_sum <- weightedCorr(Lte_col,LBK_sum_col, method="pearson", weights=Lte_col/sum(Lte_col))
                res_prod <- weightedCorr(Lte_col,LBK_prod_col, method="pearson", weights=Lte_col/sum(Lte_col))
                res <- c(sum=res_sum, prod=res_prod)
                return(res)
              }, indxy=as.list(as.data.frame(t(indxSigy_b))), pt=1:nrow(indxSigy_b), SIMPLIFY="array")
              
              corrsNorm <- mcmapply(function(indxy, pt){ 
                # l <- 1; indxy <- indxSigy_norm_b[l,]; pt <- l
                
                sigsy <- unique(sigmasy[indxy])
                Ktr_te <- sapply(sigsy, function(sigy){
                  Ktr_te <- kern_rbf(ysTrNorm_b[,pt,drop=F], ysTeNorm_b[,pt,drop=F], sigma=sigy)
                  
                }, simplify="array")
                Ktr_te_sum <- apply(Ktr_te, c(1,2), sum)
                Ktr_te_prod <- apply(Ktr_te, c(1,2), prod)
                Ktr_te_sum <- H %*% Ktr_te_sum %*% H
                Ktr_te_prod <- H %*% Ktr_te_prod %*% H
                
                
                LBK_sum <- LB%*%Ktr_te_sum
                LBK_prod <- LB%*%Ktr_te_prod
                Lte_col <- Lte[pt,]
                LBK_sum_col <- LBK_sum[pt,]
                LBK_prod_col <- LBK_prod[pt,]
                res_sum <- weightedCorr(Lte_col,LBK_sum_col, method="pearson", weights=Lte_col/sum(Lte_col))
                res_prod <- weightedCorr(Lte_col,LBK_prod_col, method="pearson", weights=Lte_col/sum(Lte_col))
                res <- c(sum=res_sum, prod=res_prod)
                
                return(res)
              }, indxy=as.list(as.data.frame(t(indxSigy_norm_b))), pt=1:nrow(indxSigy_norm_b), SIMPLIFY="array")
              
              
              
              corrs <- abind(corr=corrs, corrNorm=corrsNorm, along=3)
              
              return(corrs)
            }, simplify="array")
            proc.time() - pm # 17 secs
            
            #dim(corrsRKHS_te_boot)
            dimnames(corrsRKHS_te_boot) <- list(jointType=c("sum","prod"), testPt=1:(nPerBoot), yType=c("raw","nrm") ,boot=1:nboots)
            
            
            return(corrsRKHS_te_boot)
            
          }, simplify="array")}, simplify="array")
        proc.time() - pm # 15 mins mins
        dim(corrsRKHS_boot_joint)
        dimnames(corrsRKHS_boot_joint)[5:6] <- list(sigmaxSeq, lambdaSeq)
        names(dimnames(corrsRKHS_boot_joint))[5:6] <- c("sigmax", "lambda")
        print(proc.time()-pm)
        
        print("get opt df joint")
        
        indxOptJoint2 <- apply(corrsRKHS_boot_joint, c("testPt","boot","yType","jointType"), function(arr){
          # arr <- corrsRKHS_boot_xy[1,1,,]
          maxArr <- arr==max(arr,na.rm=T)
          indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
          #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
          return(indxs)
        })
        dimnames(indxOptJoint2)[[1]] <- c("sigmax","lambda")
        names(dimnames(indxOptJoint2))[1] <- "indx"
        
        
        df_joint <- melt(indxOptJoint2)
        df_joint <- cast(df_joint, testPt+boot+yType+jointType~indx, value="value")
        indxMat <- df_joint[,names(dimnames(corrsRKHS_boot_joint))]
        indxMat$yType <- c(1,2)[match(indxMat$yType, c("raw","nrm"))]
        indxMat$jointType <- c(1,2)[match(indxMat$jointType, c("sum","prod"))]
        df_joint$corr <- corrsRKHS_boot_joint[as.matrix(indxMat)]
        colnames(df_joint) <- c("testPt","boot","yType","jointType","indxLambda","indxSigmax","corr")
        df_joint$node <- nodeTo
        df_joint$numReg <- numReg
        df_joint$principal <- NA
        df_joint$lambda <- lambdaSeq[df_joint$indxLambda]
        df_joint$sigmax <- sigmaxSeq[df_joint$indxSigmax]
        
        sigsy <- sapply(dimnames(indxOptJoint)$numSigy, function(sigyNm){
          # sigyNm <- dimnames(indxOptJoint)$numSigy[1]
          indxMat <- df_joint[,c("testPt", "boot", "yType")]
          indxMat$yType <- c(1,2)[match(indxMat$yType, c("raw","nrm"))]
          indxMat$numSigy <- as.numeric(sigyNm)
          indxMat$indx <- 1
          indxMat <- indxMat[, names(dimnames(indxOptJoint))]
          indxOpt <- indxOptJoint[as.matrix(indxMat)]
          sigsy <- sigmasy[indxOpt]
          return(sigsy)
        })
        
        colnames(sigsy) <- paste("sigmay", dimnames(indxOptJoint)$numSigy, sep="_")
        df_joint <- cbind(df_joint, sigsy)
        df_joint <- as.data.frame(df_joint)
        df_joint$sigyType <- "rough vector"
        indxMat <- matrix(c(df_joint$testPt,df_joint$boot),nrow(df_joint),2)
        df_joint$indxPt <- smpl_te[indxMat]
        df_joint <- df_joint[order(df_joint$yType, df_joint$node, df_joint$numReg,df_joint$testPt, df_joint$boot),]
        
        res <- list(df_joint=df_joint, corrsRKHS_boot_joint=corrsRKHS_boot_joint, 
                    indxOptJoint=indxOptJoint,  sigmaxSeq=sigmaxSeq)
        
      }  
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(res)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)})
  #print(proc.time()-pm0) #  
  names(dfs_by_node_joint) <- nodes
  }
  
  print("dfs_by_node othr - use sigy of other direction")
  pm <- proc.time()
  dfs_by_node2 <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[2]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      #count <- count + 1
      regressors <- nodes[indxPreds] 
      dataX <- x[,regressors, drop=F]
      dataY <- x[,nodeTo]
        
      trainDataO <- constructData(x=dataX, y=dataY)
        
      X <- as.matrix(trainDataO$x)
      Y <- as.matrix(trainDataO$y)
      
      df <- dfs_by_node[[nodeTo]][[numReg]]$df
      if(smoothFeats) df_smooth <- dfs_by_node[[nodeTo]][[numReg]]$df_smooth
      corrsRKHS_boot  <- dfs_by_node[[nodeTo]][[numReg]]$corrsRKHS_boot 
      
      sigmaxSeq <- dfs_by_node[[nodeTo]][[numReg]]$sigmaxSeq
      otherNodes <- setdiff(nodes, nodeTo)
      df_othr <- lapply(otherNodes, function(othrNode){
        # othrNode <- otherNodes[1]
        regs <- names(dfs_by_node[[othrNode]])
        df_othr <- lapply(regs, function(othrReg){
          #othrReg <- regs[1]
          indxOpt_othr <- dfs_by_node[[othrNode]][[othrReg]]$indxOpt
          indxOther <- sapply(dimnames(corrsRKHS_boot)[["testPt"]], function(testPt) sapply(dimnames(corrsRKHS_boot)[["boot"]], function(boot) sapply(dimnames(corrsRKHS_boot)[["yType"]], function(type){
            # i<-1; j<-1; k<-1; testPt <- dimnames(corrsRKHS_boot)[["testPt"]][i]; boot <- dimnames(corrsRKHS_boot)[["boot"]][j]; type <- dimnames(corrsRKHS_boot)[["yType"]][k]
            
            indxSigmay <- indxOpt_othr["sigmay",testPt,boot,type]
            mat <- corrsRKHS_boot[testPt,type,indxSigmay,boot,,]
            maxMat <- mat==max(mat,na.rm=T)
            indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
            indxs <- c(indxSigmay, indxs)
            return(indxs)
          }, simplify="array"),simplify="array"),simplify="array")
          indxOther <- aperm(indxOther, c(1,2,4,3))
          dimnames(indxOther)[[1]] <- c("sigmay","sigmax","lambda")
          names(dimnames(indxOther)) <- c("indx","yType", "testPt","boot")
          df_othr <- melt(indxOther)
          df_othr <- cast(df_othr, testPt+boot+yType~indx, value="value")
          indxMat <- df_othr[,names(dimnames(corrsRKHS_boot))]
          indxMat$yType <- c(1,2)[match(indxMat$yType,c("raw","nrm"))]
          df_othr$corr <- corrsRKHS_boot[as.matrix(indxMat)]
          colnames(df_othr) <- c("testPt","boot","yType","indxLambda","indxSigmax","indxSigmay","corr")
          df_othr$node <- othrNode
          df_othr$numReg <- othrReg
          df_othr$principal <- FALSE
          df_othr$lambda <- lambdaSeq[df_othr$indxLambda]
          df_othr$sigmax <- sigmaxSeq[df_othr$indxSigmax]
          df_othr$sigmay <- sigmasy[df_othr$indxSigmay]
          df_othr$sigyType <- "rough"
          df_othr <- as.data.frame(df_othr)
          
          return(df_othr)
        })
        return(df_othr)
      })
      df_othr <- unlist(df_othr, recursive=F)
      df_othr <- do.call(rbind, df_othr)
      df <- rbind(df, df_othr)
      
      if(smoothFeats){
      corrsRKHS_boot_smoothL <- lapply(otherNodes, function(othrNode){
        # othrNode <- otherNodes[1]
        regs <- names(dfs_by_node[[othrNode]])
        corrsRKHS_boot_smooth <- lapply(regs, function(othrReg){
          #othrReg <- regs[1]
          
          
          optSigmay <- dfs_by_node[[othrNode]][[othrReg]]$optSigmay
           
          pm <- proc.time()
          corrsRKHS_boot_smooth <- sapply(1:length(lambdaSeq), function(lam){ 
            #i <<- i + 1
            #j <<- 0
            sapply(1:length(sigmaxSeq) , function(sigmax){
              # i <- 1; j <- 2; lam <- i; sigmax <- j
              
              
              L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
              ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[sigmax]]
              
              #weighted.mean(ysNorm[,3], w=L[,3]/sum(L[,3])
              #weighted.var(ysNorm[,5], w=L[,5]/sum(L[,5])
              
              pm <- proc.time()
              corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
                #k <- 1
                
                optSigmay_b <- optSigmay[smpl_te[,k],1]
                optSigmay_norm_b <- optSigmay[smpl_te[,k],2]
                
                indxTr_b <- smpl_tr[,k]
                indxTe_b <- smpl_te[,k]
                Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
                Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
                Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
                Itr_b <- diag(nPerBoot)
                Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[k]]
                
                
                LB <- Lte_tr%*%Blambda_tr
                I <- diag(nPerBoot)
                H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
                
                yTr_b <- as.matrix(Y[smpl_tr[,k]])
                yTe_b <- as.matrix(Y[smpl_te[,k]])
                
                ysTrNorm_b <- ysNorm[smpl_tr[,k],]; ysTrNorm_b <- ysTrNorm_b[,smpl_tr[,k]] 
                ysTeNorm_b <- ysNorm[smpl_te[,k],]; ysTeNorm_b <- ysTeNorm_b[,smpl_te[,k]] 
                # l <- 4; weighted.mean(ysNorm[,l], w=L[,l]/sum(L[,l]))
                
                corrs <- mcmapply(function(sigmay, pt){ 
                  # l <- 1; sigmay <- optSigmay_b[l]; pt <- l
                  
                  
                  Ktr_te <- kern_rbf(yTr_b, yTe_b, sigma=sigmay)
                  Ktr_te <- H %*% Ktr_te %*% H
                  LBK <- LB%*%Ktr_te
                  Lte_col <- Lte[pt,]
                  LBK_col <- LBK[pt,]
                  res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                  return(res)
                }, sigmay=optSigmay_b, pt=1:length(optSigmay_b), SIMPLIFY="array")
                
                corrsNorm <- mcmapply(function(sigmay, pt){ 
                  # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
                  
                  
                  Ktr_te <- kern_rbf(ysTrNorm_b[,pt,drop=F], ysTeNorm_b[,pt,drop=F], sigma=sigmay)
                  Ktr_te <- H %*% Ktr_te %*% H
                  LBK <- LB%*%Ktr_te
                  Lte_col <- Lte[pt,]
                  LBK_col <- LBK[pt,]
                  res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
                  return(res)
                }, sigmay=optSigmay_norm_b, pt=1:length(optSigmay_norm_b), SIMPLIFY="array")
                
                corrs <- cbind(corr=corrs, corrNorm=corrsNorm)
                
                return(corrs)
              }, simplify="array")
              proc.time() - pm # 17 secs
              
              #dim(corrsRKHS_te_boot)
              dimnames(corrsRKHS_te_boot) <- list(testPt=1:(nPerBoot),yType=c("raw","nrm") ,boot=1:nboots)
              
              
              return(corrsRKHS_te_boot)
              
            }, simplify="array")}, simplify="array")
          proc.time() - pm # 15 mins mins
          dim(corrsRKHS_boot_smooth)
          dimnames(corrsRKHS_boot_smooth)[4:5] <- list(sigmaxSeq, lambdaSeq)
          names(dimnames(corrsRKHS_boot_smooth))[4:5] <- c("sigmax", "lambda")
            
      
          
          return(corrsRKHS_boot_smooth)
        })
        names(corrsRKHS_boot_smooth) <- regs
        return(corrsRKHS_boot_smooth)
      })
      names(corrsRKHS_boot_smoothL) <- otherNodes
      
      indxOtherSmoothL <- lapply(otherNodes, function(othrNode){
        # othrNode <- otherNodes[1]
        regs <- names(dfs_by_node[[othrNode]])
        indxOther <- lapply(regs, function(othrReg){
          #othrReg <- regs[1]
          
          
          optSigmay <- dfs_by_node[[othrNode]][[othrReg]]$optSigmay
          corrsRKHS_boot_smooth <- corrsRKHS_boot_smoothL[[othrNode]][[othrReg]]
          
          
          indxOther <- sapply(dimnames(corrsRKHS_boot_smooth)[["testPt"]], function(testPt) sapply(dimnames(corrsRKHS_boot_smooth)[["boot"]], function(boot) sapply(dimnames(corrsRKHS_boot_smooth)[["yType"]], function(type){
            # i<-1; j<-1; k <- 1; testPt <- dimnames(corrsRKHS_boot)[["testPt"]][i]; boot <- dimnames(corrsRKHS_boot)[["boot"]][j]; type <- dimnames(corrsRKHS_boot)[["yType"]][k]
            
            mat <- corrsRKHS_boot_smooth[testPt,type,boot,,]
            maxMat <- mat==max(mat,na.rm=T)
            indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
            return(indxs)
          }, simplify="array"),simplify="array"),simplify="array")
          indxOther <- aperm(indxOther, c(1,2,4,3))
          dimnames(indxOther)[[1]] <- c("sigmax","lambda")
          names(dimnames(indxOther)) <- c("indx", "yType","testPt","boot")
          
          return(indxOther)
        })
        names(indxOther) <- regs
        return(indxOther)
      })
      names(indxOtherSmoothL) <- otherNodes
      
      df_othr_smooth <- lapply(otherNodes, function(othrNode){
        # othrNode <- otherNodes[1]
        regs <- names(dfs_by_node[[othrNode]])
        df_othr_smooth <- lapply(regs, function(othrReg){
          #othrReg <- regs[1]
          
          
          optSigmay <- dfs_by_node[[othrNode]][[othrReg]]$optSigmay
          corrsRKHS_boot_smooth <- corrsRKHS_boot_smoothL[[othrNode]][[othrReg]]
          indxOther <- indxOtherSmoothL[[othrNode]][[othrReg]]
          
          
          df_othr_smooth <- melt(indxOther)
          df_othr_smooth <- cast(df_othr_smooth, testPt+boot+yType~indx, value="value")
          indxMat <- df_othr_smooth[,names(dimnames(corrsRKHS_boot_smooth))]
          indxMat$yType <- c(1,2)[match(indxMat$yType, c("raw","nrm"))]
          df_othr_smooth$corr <- corrsRKHS_boot_smooth[as.matrix(indxMat)]
          colnames(df_othr_smooth) <- c("testPt","boot","yType","indxLambda","indxSigmax","corr")
          df_othr_smooth$node <- othrNode
          df_othr_smooth$numReg <- othrReg
          df_othr_smooth$principal <- FALSE
          df_othr_smooth$lambda <- lambdaSeq[df_othr_smooth$indxLambda]
          df_othr_smooth$sigmax <- sigmaxSeq[df_othr_smooth$indxSigmax]
          indxMat <- matrix(c(df_othr_smooth$testPt,df_othr_smooth$yType),nrow(df_othr_smooth),2)
          df_othr_smooth$sigmay <- optSigmay[indxMat]
          df_othr_smooth$sigyType <- "smooth"
          df_othr_smooth <- as.data.frame(df_othr_smooth)
          
          
          
          return(df_othr_smooth)
        })
        names(df_othr_smooth) <- regs
        return(df_othr_smooth)
      })
      names(df_othr_smooth) <- otherNodes
      df_othr_smooth <- unlist(df_othr_smooth, recursive=F)
      df_othr_smooth <- do.call(rbind, df_othr_smooth)
      df_smooth <- rbind(df_smooth, df_othr_smooth)
      }
            
      
      
      #all(df_othr==df_xy[which(df_xy$node=="x"),colnames(df_othr)])
      #all(df_othr==df_yx[which(df_yx$node=="y"),colnames(df_othr)])
      
      
      
      otherRegs <- setdiff(numRegs[[nodeTo]], numReg)
      if(length(otherRegs)>0){
        
        df_othr_reg <- lapply(otherRegs, function(othrReg){
        #othrReg <- regs[1]
        indxOpt_othr <- dfs_by_node[[nodeTo]][[othrReg]]$indxOpt
        indxOther <- sapply(dimnames(corrsRKHS_boot)[["testPt"]], function(testPt) sapply(dimnames(corrsRKHS_boot)[["boot"]], function(boot){
          # i<-1; j<-1; testPt <- dimnames(corrsRKHS_boot)[["testPt"]][i]; boot <- dimnames(corrsRKHS_boot)[["boot"]][j]
          
          indxSigmay <- indxOpt_othr["sigmay",testPt,boot]
          mat <- corrsRKHS_boot[testPt,indxSigmay,boot,,]
          maxMat <- mat==max(mat,na.rm=T)
          indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
          indxs <- c(indxSigmay, indxs)
          return(indxs)
        }, simplify="array"),simplify="array")
        
        indxOther <- aperm(indxOther, c(1,3,2))
        dimnames(indxOther)[[1]] <- c("sigmay","sigmax","lambda")
        names(dimnames(indxOther)) <- c("indx", "testPt","boot")
        df_othr <- melt(indxOther)
        df_othr <- cast(df_othr, testPt+boot~indx, value="value")
        df_othr$corr <- corrsRKHS_boot[as.matrix(df_othr[,names(dimnames(corrsRKHS_boot))])]
        colnames(df_othr) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
        df_othr$node <- othrNode
        df_othr$numReg <- othrReg
        df_othr$lambda <- lambdaSeq[df_othr$indxLambda]
        df_othr$sigmax <- sigmaxSeq[df_othr$indxSigmax]
        df_othr$sigmay <- sigmasy[df_othr$indxSigmay]
        df_othr <- as.data.frame(df_othr)
        return(df_othr)
      })
        df_other_reg <- do.call(rbind, df_othr_reg)
        df <- rbind(df, df_other_reg)
        
        df_othr_reg_smooth <- lapply(otherRegs, function(othrReg){
          #othrReg <- regs[1]
          
          indxOpt_othr <- dfs_by_node[[nodeTo]][[othrReg]]$indxOptSmooth
          optSigmay <- dfs_by_node[[nodeTo]][[othrReg]]$optSigmay
          indxOther <- sapply(dimnames(corrsRKHS_boot_smooth)[["testPt"]], function(testPt) sapply(dimnames(corrsRKHS_boot_smooth)[["boot"]], function(boot){
            # i<-1; j<-1; testPt <- dimnames(corrsRKHS_boot)[["testPt"]][i]; boot <- dimnames(corrsRKHS_boot)[["boot"]][j]
            
            indxSigmay <- indxOpt_othr["sigmay",testPt,boot]
            mat <- corrsRKHS_boot_smooth[testPt,indxSigmay,boot,,]
            maxMat <- mat==max(mat,na.rm=T)
            indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
            indxs <- c(indxSigmay, indxs)
            return(indxs)
          }, simplify="array"),simplify="array")
          indxOther <- aperm(indxOther, c(1,3,2))
          dimnames(indxOther)[[1]] <- c("sigmay","sigmax","lambda")
          names(dimnames(indxOther)) <- c("indx", "testPt","boot")
          df_othr_smooth <- melt(indxOther)
          df_othr_smooth <- cast(df_othr_smooth, testPt+boot~indx, value="value")
          df_othr_smooth$corr <- corrsRKHS_boot_smooth[as.matrix(df_othr_smooth[,names(dimnames(corrsRKHS_boot_smooth))])]
          colnames(df_othr_smooth) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
          df_othr_smooth$node <- othrNode
          df_othr_smooth$numReg <- othrReg
          df_othr_smooth$principal <- FALSE
          df_othr_smooth$lambda <- lambdaSeq[df_othr_smooth$indxLambda]
          df_othr_smooth$sigmax <- sigmaxSeq[df_othr_smooth$indxSigmax]
          df_othr_smooth$sigmay <- optSigmay[df_othr_smooth$indxSigmay]
          df_othr_smooth <- as.data.frame(df_othr_smooth)
          return(df_othr_smooth)
        })
      }
        
      indxMat <- matrix(c(df$testPt,df$boot),nrow(df),2)
      df$indxPt <- smpl_te[indxMat]
      df <- df[order(df$yType, df$node, df$numReg,df$testPt, df$boot),]
      res <- list(df=df)
      if(smoothFeats){
        indxMat <- matrix(c(df_smooth$testPt,df_smooth$boot),nrow(df_smooth),2)
        df_smooth$indxPt <- smpl_te[indxMat]
        df_smooth <- df_smooth[order(df_smooth$yType,df_smooth$node, df_smooth$numReg,df_smooth$testPt, df_smooth$boot),]
        print(df_smooth$sigmay[1:10])
        res <- c(res, list(df_smooth=df_smooth,indxOtherSmoothL=indxOtherSmoothL))
      }
      #all(df[,colnames(df)[which(colnames(df)%in%colnames(df_xy))]]==df_xy[,colnames(df)[which(colnames(df)%in%colnames(df_xy))]])
      #all(df[,colnames(df)[which(colnames(df)%in%colnames(df_yx))]]==df_yx[,colnames(df)[which(colnames(df)%in%colnames(df_yx))]])
      
      #print(table(df$sigmay))
      
      
        
        
      print(proc.time()-pm)
        
      return(res)
    })  
      
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  print((proc.time()-pm)[3])
  names(dfs_by_node2) <- nodes
  
  print("measureList")
  pm <- proc.time()
  measureList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      #count <- count + 1
      regressors <- nodes[indxPreds] 
      dataX <- x[,regressors, drop=F]
      dataY <- x[,nodeTo]
      
      trainDataO <- constructData(x=dataX, y=dataY)
      
      X <- as.matrix(trainDataO$x)
      Y <- as.matrix(trainDataO$y)
      
      
      
      df <- dfs_by_node2[[nodeTo]][[numReg]]$df
      df_rnd <- df
      
      set.seed(123456)
      smpl1 <- sample(nrow(df))
      smpl2 <- sample(nrow(df))
      smpl3 <- sample(nrow(df))
      
      df_rnd$sigmay <- df_rnd$sigmay[smpl1]
      df_rnd$indxSigmay <- df_rnd$indxSigmay[smpl1]
      df_rnd$sigmax <- df_rnd$sigmax[smpl2]
      df_rnd$indxSigmax <- df_rnd$indxSigmax[smpl2]
      df_rnd$lambda <- df_rnd$lambda[smpl3]
      df_rnd$lambda <- df_rnd$indxLambda[smpl3]
      
      df_raw <- df[which(df$yType=="raw"),]
      df_nrm <- df[which(df$yType=="nrm"),]
      df_rnd_raw <- df_rnd[which(df_rnd$yType=="raw"),]
      df_rnd_nrm <- df_rnd[which(df_rnd$yType=="nrm"),]
      dim(df_raw)  
      dim(df_nrm)
      all(df_raw[,c("testPt","boot","node","numReg","indxPt")]==df_nrm[,c("testPt","boot","node","numReg","indxPt")])
      
      
      if(smoothFeats){
      df_smooth <- dfs_by_node2[[nodeTo]][[numReg]]$df_smooth
      df_smooth_raw <- df_smooth[which(df_smooth$yType=="raw"),]
      df_smooth_nrm <- df_smooth[which(df_smooth$yType=="nrm"),]
      dim(df_smooth_raw)  
      dim(df_smooth_nrm)
      all(df_raw[,c("testPt","boot","node","numReg","indxPt")]==df_smooth_raw[,c("testPt","boot","node","numReg","indxPt")])
      all(df_raw[,c("testPt","boot","node","numReg","indxPt")]==df_smooth_nrm[,c("testPt","boot","node","numReg","indxPt")])
      
      }
      
      
      
      if(jointFeats){
        df_joint <- dfs_by_node_joint[[nodeTo]][[numReg]]$df
        df_joint_sum_raw <- df_joint[which(df_joint$yType=="raw" & df_joint$jointType=="sum"),]
        df_joint_sum_nrm <- df_joint[which(df$yType=="nrm" & df_joint$jointType=="sum"),]
        df_joint_prod_raw <- df_joint[which(df_joint$yType=="raw" & df_joint$jointType=="prod"),]
        df_joint_prod_nrm <- df_joint[which(df_joint$yType=="nrm" & df_joint$jointType=="prod"),]
        dim(df_joint_sum_raw)  
        dim(df_joint_sum_nrm) 
        dim(df_joint_prod_raw)  
        dim(df_joint_prod_nrm) 
        
        all(df_raw[(nrow(df_joint_sum_raw)+1):nrow(df_raw),c("testPt","boot","node","numReg","indxPt")]==df_joint_sum_raw[,c("testPt","boot","node","numReg","indxPt")])
        all(df_raw[(nrow(df_joint_sum_raw)+1):nrow(df_raw),c("testPt","boot","node","numReg","indxPt")]==df_joint_sum_nrm[,c("testPt","boot","node","numReg","indxPt")])
        all(df_raw[(nrow(df_joint_sum_raw)+1):nrow(df_raw),c("testPt","boot","node","numReg","indxPt")]==df_joint_prod_raw[,c("testPt","boot","node","numReg","indxPt")])
        all(df_raw[(nrow(df_joint_sum_raw)+1):nrow(df_raw),c("testPt","boot","node","numReg","indxPt")]==df_joint_prod_nrm[,c("testPt","boot","node","numReg","indxPt")])
        
        head(df_raw[1:nrow(df_joint_sum_raw),c("testPt","boot","node","numReg","indxPt")])
        head(df_joint_sum_raw[,c("testPt","boot","node","numReg","indxPt")])
        
      }
      
      
      
      
      
      
      # non-smooth - can read off Ks from list Ks_by_node (raw) or KsNorm_by_node (norml)
      # smooth - sigmay is not so discrete -> Ks must be calculated
      
      # how to reference:
      # Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
      # Bs_by_node[[nodeTo]][[numReg]][[lambda]][[sigmax]][[boot]] 
      # Ks_by_node[[nodeTo]][[numReg]][[sigmay]]
      # yNorm_by_node[[nodeTo]][[numReg]][[sigmax]][,testPt] 
      # KsNorm_by_node[[nodeTo]][[numReg]][[sigmax]][[sigmay]]
      
      modDens <- kepdf(X, eval.points = X, kernel = "gaussian", bwtype = "adaptive")
      
      dfList <- list(df_raw, df_nrm, df_rnd_raw, df_rnd_nrm) 
      
      nms <- c("raw","nrm","raw_rnd","nrm_rnd")
      
      if(TRUE){
        mc_cores_loc <- 1
        mc_cores_loc2 <- detectCores()/2
        print("msrs1")
        pm <- proc.time()
        msrs1List <- mcmapply(function(df, nm){
          # df <- dfList[[1]]; nm <- nms[[1]]
          if(strsplit(nm, "_")[[1]][1]=="raw"){
            KList <- Ks_by_node[[nodeTo]][[numReg]]  
            KList <-  lapply(1:numSigx, function(i) KList)
          } else{
            KList <- KsNorm_by_node[[nodeTo]][[numReg]]    
          }
        
          pm <- proc.time()
          msrs1 <- t(mcmapply(function(row){
            # row <- df[1,]; names(row) <- colnames(df) 
            lam <- as.numeric(df[row,"indxLambda"])
            sigmax <- as.numeric(df[row,"indxSigmax"])
            sigmay <- as.numeric(df[row,"indxSigmay"])
            boot <- as.numeric(df[row, "boot"])
            testPt <- as.numeric(df[row, "testPt"])
            indxTe_b <- smpl_te[,boot]
            indxTr_b <- smpl_tr[,boot]
            L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
            Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
            Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
            Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
            Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
            K <- KList[[sigmax]][[sigmay]]
            Ktr <- K[indxTr_b,]; Ktr <- Ktr[,indxTr_b]
            Kte <- K[indxTe_b,]; Kte <- Kte[,indxTe_b]
            Kte_tr <- K[indxTe_b,]; Kte_tr <- Kte_tr[,indxTr_b]
            Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
            normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
            normsCME <- normsCME^0.5
            resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
            # local weights
            weights1 <- Lte[testPt,]
            weights1 <- weights1/sum(weights1)
            # distribution of hypothesis cause
            weights2 <- 1/modDens@estimate[indxTe_b]
            weights2 <- weights2/sum(weights2)
            # both
            weights <- weights1*weights2
            weights <- weights/sum(weights)
        
        
        
            ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
            ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
            names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
            names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
            lam <- as.numeric(df[row, "lambda"])
            sigmax <- as.numeric(df[row, "sigmax"])
            sigmay <- as.numeric(df[row, "sigmay"])
            pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(df[row, "corr"]))
        
            return(c(pars, ftsNorm, ftsResid))
      }, row=1:nrow(df), mc.cores=mc_cores_loc2, SIMPLIFY="array"))
          proc.time() - pm # 5secs
        
          return(msrs1)
      }, df=dfList, nm=nms, mc.cores=mc_cores_loc, SIMPLIFY=FALSE)
        
        print((proc.time()-pm)[3])
        #length(msrs1List); length(msrs1List2)
        #i <- 4; summary(c(msrs1List[[i]]-msrs1List2[[i]]))
      } else{
      
        print("msrs1")
        pm <- proc.time()
        print("raw")
        pm <- proc.time()
        msrs1_raw <- t(apply(df_raw, 1, function(row){
        # row <- df_raw[1,]; names(row) <- colnames(df_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["indxSigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        K <- Ks_by_node[[nodeTo]][[numReg]][[sigmay]]
        Ktr <- K[indxTr_b,]; Ktr <- Ktr[,indxTr_b]
        Kte <- K[indxTe_b,]; Kte <- Kte[,indxTe_b]
        Kte_tr <- K[indxTe_b,]; Kte_tr <- Kte_tr[,indxTr_b]
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        # local weights
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        # distribution of hypothesis cause
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        # both
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        proc.time() - pm # 5secs
        print("nrm")
        msrs1_nrm <- t(apply(df_nrm, 1, function(row){
        # row <- df_raw[1,]; names(row) <- colnames(df_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["indxSigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        K <- KsNorm_by_node[[nodeTo]][[numReg]][[sigmax]][[sigmay]]
        Ktr <- K[indxTr_b,]; Ktr <- Ktr[,indxTr_b]
        Kte <- K[indxTe_b,]; Kte <- Kte[,indxTe_b]
        Kte_tr <- K[indxTe_b,]; Kte_tr <- Kte_tr[,indxTr_b]
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        print("raw rnd")
        msrs1_rnd_raw <- t(apply(df_rnd_raw, 1, function(row){
        # row <- df_raw[1,]; names(row) <- colnames(df_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["indxSigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        K <- Ks_by_node[[nodeTo]][[numReg]][[sigmay]]
        Ktr <- K[indxTr_b,]; Ktr <- Ktr[,indxTr_b]
        Kte <- K[indxTe_b,]; Kte <- Kte[,indxTe_b]
        Kte_tr <- K[indxTe_b,]; Kte_tr <- Kte_tr[,indxTr_b]
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        # local weights
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        # distribution of hypothesis cause
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        # both
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        print("nrm rnd")
        msrs1_rnd_nrm <- t(apply(df_rnd_nrm, 1, function(row){
        # row <- df_raw[1,]; names(row) <- colnames(df_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["indxSigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        K <- KsNorm_by_node[[nodeTo]][[numReg]][[sigmax]][[sigmay]]
        Ktr <- K[indxTr_b,]; Ktr <- Ktr[,indxTr_b]
        Kte <- K[indxTe_b,]; Kte <- Kte[,indxTe_b]
        Kte_tr <- K[indxTe_b,]; Kte_tr <- Kte_tr[,indxTr_b]
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
      
        msrs1List <- list(msrs1_raw, msrs1_nrm, msrs1_rnd_raw, msrs1_rnd_nrm)
      
      }
      
      if(smoothFeats){
      print("smth_raw")
      msrs1_smooth_raw <- t(apply(df_smooth_raw, 1, function(row){
        # row <- df_smooth_raw[1,]; names(row) <- colnames(df_smooth_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["sigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        Ktr <- kern_rbf(yTr_b, sigma=sigmay)
        Kte <- kern_rbf(yTe_b, sigma=sigmay)
        Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmay)
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
      print("smth_nrm")
      msrs1_smooth_nrm <- t(apply(df_smooth_nrm, 1, function(row){
        # row <- df_smooth_nrm[1,]; names(row) <- colnames(df_smooth_nrm) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        sigmay <- as.numeric(row["sigmay"])
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[sigmax]]
        yTr_b <- ysNorm[indxTr_b,]; yTr_b <- yTr_b[,indxTr_b[testPt], drop=F] 
        yTe_b <- ysNorm[indxTe_b,]; yTe_b <- yTe_b[,indxTe_b[testPt], drop=F]
        Ktr <- kern_rbf(yTr_b, sigma=sigmay)
        Kte <- kern_rbf(yTe_b, sigma=sigmay)
        Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmay)
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- as.numeric(row["sigmay"])
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=log(sigmay,10), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
      dfList <- c(dfList, list(df_smooth_raw, df_smooth_nrm))
      msrs1List <- c(msrs1List, list(msrs1_smooth_raw, msrs1_smooth_nrm))
      nms <- c(nms,"smth_raw","smth_nrm")
      }
      
      if(jointFeats){
        print("joint sum raw")
        msrs1_joint_sum_raw <- t(apply(df_joint_sum_raw, 1, function(row){
        # row <- df_joint_sum_raw[1,]; names(row) <- colnames(df_joint_sum_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        
        indxSigyVars <- grep("sigmay", names(row))
        sigsy <- unique(as.numeric(row[indxSigyVars]))
        
        
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        I <- diag(nPerBoot)
        H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
        
        Ktr <- sapply(sigsy, function(sigy){
          # sigy <- sigsy[1]
          Ktr <- kern_rbf(yTr_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Ktr <- apply(Ktr, c(1,2), sum)
        Ktr <- H %*% Ktr %*% H
        
        Kte <- sapply(sigsy, function(sigy){
          Kte <- kern_rbf(yTe_b, yTe_b, sigma=sigy)
        }, simplify="array")
        Kte <- apply(Kte, c(1,2), sum)
        Kte <- H %*% Kte %*% H
        
        Kte_tr <- sapply(sigsy, function(sigy){
          Kte_tr <- kern_rbf(yTe_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Kte_tr <- apply(Kte_tr, c(1,2), sum)
        Kte_tr <- H %*% Kte_tr %*% H
        
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        # local weights
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        # distribution of hypothesis cause
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        # both
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- sigsy
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=mean(log(sigmay,10)), corr=-as.numeric(row["corr"]))
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        print("joint sum nrm")
        msrs1_joint_sum_nrm <- t(apply(df_joint_sum_nrm, 1, function(row){
        # row <- df_joint_sum_nrm[1,]; names(row) <- colnames(df_joint_sum_nrm)
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        
        indxSigyVars <- grep("sigmay", names(row))
        sigsy <- unique(as.numeric(row[indxSigyVars]))
        
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        I <- diag(nPerBoot)
        H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
        
        Ktr <- sapply(sigsy, function(sigy){
          Ktr <- kern_rbf(yTr_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Ktr <- apply(Ktr, c(1,2), sum)
        Ktr <- H %*% Ktr %*% H
        
        Kte <- sapply(sigsy, function(sigy){
          Kte <- kern_rbf(yTe_b, yTe_b, sigma=sigy)
        }, simplify="array")
        Kte <- apply(Kte, c(1,2), sum)
        Kte <- H %*% Kte %*% H
        
        Kte_tr <- sapply(sigsy, function(sigy){
          Kte_tr <- kern_rbf(yTe_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Kte_tr <- apply(Kte_tr, c(1,2), sum)
        Kte_tr <- H %*% Kte_tr %*% H
        
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- sigsy
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=mean(log(sigmay,10)), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        print("joint prod raw")
        msrs1_joint_prod_raw <- t(apply(df_joint_prod_raw, 1, function(row){
        # row <- df_smooth_raw[1,]; names(row) <- colnames(df_smooth_raw) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        
        indxSigyVars <- grep("sigmay", names(row))
        sigsy <- unique(as.numeric(row[indxSigyVars]))
        
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        I <- diag(nPerBoot)
        H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
        
        Ktr <- sapply(sigsy, function(sigy){
          Ktr <- kern_rbf(yTr_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Ktr <- apply(Ktr, c(1,2), prod)
        Ktr <- H %*% Ktr %*% H
        
        Kte <- sapply(sigsy, function(sigy){
          Kte <- kern_rbf(yTe_b, yTe_b, sigma=sigy)
        }, simplify="array")
        Kte <- apply(Kte, c(1,2), prod)
        Kte <- H %*% Kte %*% H
        
        Kte_tr <- sapply(sigsy, function(sigy){
          Kte_tr <- kern_rbf(yTe_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Kte_tr <- apply(Kte_tr, c(1,2), prod)
        Kte_tr <- H %*% Kte_tr %*% H
        
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- sigsy
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=mean(log(sigmay,10)), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
        print("joint prod nrm")
        msrs1_joint_prod_nrm <- t(apply(df_joint_prod_nrm, 1, function(row){
        # row <- df_smooth_nrm[1,]; names(row) <- colnames(df_smooth_nrm) 
        lam <- as.numeric(row["indxLambda"])
        sigmax <- as.numeric(row["indxSigmax"])
        
        indxSigyVars <- grep("sigmay", names(row))
        sigsy <- unique(as.numeric(row[indxSigyVars]))
        
        boot <- as.numeric(row["boot"])
        testPt <- as.numeric(row["testPt"])
        indxTe_b <- smpl_te[,boot]
        indxTr_b <- smpl_tr[,boot]
        yTr_b <- as.matrix(Y[indxTr_b])
        yTe_b <- as.matrix(Y[indxTe_b])
        L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
        Ltr <- L[indxTr_b,]; Ltr <- Ltr[,indxTr_b]
        Lte <- L[indxTe_b,]; Lte <- Lte[,indxTe_b]
        Lte_tr <- L[indxTe_b,]; Lte_tr <- Lte_tr[,indxTr_b]
        Blambda_tr <- Bs_by_node[[nodeTo]][[numReg]][[lam]][[sigmax]][[boot]]
        ysNorm <- yNorm_by_node[[nodeTo]][[numReg]][[sigmax]]
        yTr_b <- ysNorm[indxTr_b,]; yTr_b <- yTr_b[,indxTr_b[testPt], drop=F] 
        yTe_b <- ysNorm[indxTe_b,]; yTe_b <- yTe_b[,indxTe_b[testPt], drop=F]
        I <- diag(nPerBoot)
        H <- I-matrix(1/nPerBoot,nPerBoot,nPerBoot)
        
        Ktr <- sapply(sigsy, function(sigy){
          Ktr <- kern_rbf(yTr_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Ktr <- apply(Ktr, c(1,2), prod)
        Ktr <- H %*% Ktr %*% H
        
        Kte <- sapply(sigsy, function(sigy){
          Kte <- kern_rbf(yTe_b, yTe_b, sigma=sigy)
        }, simplify="array")
        Kte <- apply(Kte, c(1,2), prod)
        Kte <- H %*% Kte %*% H
        
        Kte_tr <- sapply(sigsy, function(sigy){
          Kte_tr <- kern_rbf(yTe_b, yTr_b, sigma=sigy)
        }, simplify="array")
        Kte_tr <- apply(Kte_tr, c(1,2), prod)
        Kte_tr <- H %*% Kte_tr %*% H
        
        Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
        normsCME <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
        normsCME <- normsCME^0.5
        resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
        weights1 <- Lte[testPt,]
        weights1 <- weights1/sum(weights1)
        weights2 <- 1/modDens@estimate[indxTe_b]
        weights2 <- weights2/sum(weights2)
        weights <- weights1*weights2
        weights <- weights/sum(weights)
        
        
        
        ftsNorm <- featsCalc(normsCME, testPt, weights1, weights2, weights)
        ftsResid <- featsCalc(resids, testPt, weights1, weights2, weights)
        names(ftsNorm) <- gsub("v_","normCME_",names(ftsNorm))
        names(ftsResid) <- gsub("v_","rsds_",names(ftsResid))
        
        lam <- as.numeric(row["lambda"])
        sigmax <- as.numeric(row["sigmax"])
        sigmay <- sigsy
        pars <- c(lambda=log(lam,10), sigmax=log(sigmax,10), sigmay=mean(log(sigmay,10)), corr=-as.numeric(row["corr"]))
        
        
        return(c(pars, ftsNorm, ftsResid))
      }))
      }
      
      
      # plot(msrs1[,"kcrdc"], KCRDCs_xy); abline(a=0, b=1, col="red")
      # plot(msrs1[,"wkcrdc"], wKCRDCs_xy); abline(a=0, b=1, col="red")
      # plot(msrs1[,"wl2rel"], wL2_f_xy); abline(a=0, b=1, col="red")
      # 
      # plot(msrs1[,"kcrdc"], KCRDCs_yx); abline(a=0, b=1, col="red")
      # plot(msrs1[,"wkcrdc"], wKCRDCs_yx); abline(a=0, b=1, col="red")
      # plot(msrs1[,"wl2rel"], wL2_f_yx); abline(a=0, b=1, col="red")
      
      print("rearrange 1")
      df_msrs1_List <- mcmapply(FUN=function(df, msrs1, nm){
        # i <- 1; df <- dfList[[i]]; msrs1 <- msrs1List[[i]]; nm <- nms[i]
        print(nm)
        # separate by sigmay view
        df_msrs1 <- lapply(nodes, function(nd){
          # nd <- nodes[2]
          regs <- names(dfs_by_node2[[nd]])
          df_msrs1 <- lapply(regs, function(rg){
            # rg <- regs[[1]]
            indxNR <- which(df$node==nd & df$numReg==rg)
            
            res <- list(df=df[indxNR,], msrs1=msrs1[indxNR,])
            
          })
          names(df_msrs1) <- regs
          return(df_msrs1)
        })
        names(df_msrs1) <- nodes
        
        df_msrs1 <- unlist(df_msrs1, recursive=F)
        
        dfs_el <- lapply(df_msrs1, function(el) el$df)
        msrs1_el <- lapply(df_msrs1, function(el) el$msrs1)
      
        res <- list(dfs=dfs_el, msrs1=msrs1_el)
        return(res)
        
      }, df=dfList, msrs1=msrs1List, nm=nms, mc.cores=length(dfList), SIMPLIFY=FALSE)
      names(df_msrs1_List) <- nms
      
      print("rearrange 2a")
      dfList <- lapply(df_msrs1_List, function(el) el$dfs)
      
      print("rearrange 2b")
      msrs1List <- lapply(df_msrs1_List, function(el) el$msrs1)
      
      print("rearrange 2c")
      dfList <- unlist(dfList, recursive=F)
      msrs1List <- unlist(msrs1List, recursive=F)
      
      names(dfList)
      names(msrs1List)
      sapply(dfList, dim)
      sapply(msrs1List, dim)
      
      if(jointFeats){
        nms2 <- c("joint_sum_raw","joint_sum_nrm", "joint_prod_raw","joint_prod_nrm")
        dfList2 <- list(df_joint_sum_raw, df_joint_sum_nrm, df_joint_prod_raw, df_joint_prod_nrm)
        msrs1List2 <- list(msrs1_joint_sum_raw, msrs1_joint_sum_nrm, msrs1_joint_prod_raw, msrs1_joint_prod_nrm)
        names(dfList2) <- nms2
        names(msrs1List2) <- nms2
      
        dfList <- c(dfList, dfList2)
        msrs1List <- c(msrs1List, msrs1List2)
      }
      
      
      
      print("msrs2")
      mc_cores_loc <- 1
      mc_cores_loc2 <- detectCores()/2
      msrs2List <- mcmapply(FUN=function(df, msrs1, nm){
        # i <- 1; df <- dfList[[i]]; msrs1 <- msrs1List[[i]]; nm <- names(dfList)[i]
        print(nm)
        msrs2 <- t(mcmapply(function(i){
            # i <- 1
            row <- df[i,]; names(row) <- colnames(df)
            sigmax <- as.numeric(row["indxSigmax"])
            testPt <- as.numeric(row["testPt"])
            xPt <- as.numeric(row["indxPt"])
            
            L <- Ls_by_node[[nodeTo]][[numReg]][[sigmax]]
            Lte_tr <- L[xPt,]; 
            # get the weights for the all bootstrap x pointsPerBoot
            weights1 <- Lte_tr[df$indxPt]
            weights1 <- weights1/sum(weights1)
            # hist(weights1)
            # plot(X[df$indxPt[indxNR]], weights1); abline(v=X[xPt], col="red")
            weights2 <- 1/modDens@estimate[df$indxPt]
            weights2 <- weights2/sum(weights2)
            # plot(X[df$indxPt[indxNR]], weights2); abline(v=X[xPt], col="red")
            weights <- weights1*weights2
            weights <- weights/sum(weights)
            # plot(X[$indxPt[indxNR]], weights); abline(v=X[xPt], col="red")
            
            locIndx <- grep("loc",colnames(msrs1))
            distIndx <- grep("dist",colnames(msrs1))
            locDistIndx <- intersect(locIndx,distIndx)
            locIndx <- setdiff(locIndx, locDistIndx)
            distIndx <- setdiff(distIndx, locDistIndx)
            allIndx <- setdiff(1:ncol(msrs1), c(locIndx, distIndx, locDistIndx))
            # all(sort((c(locIndx, distIndx,locDistIndx, allIndx)))==1:ncol(msrs1))
            indxList <- list(locIndx, distIndx,locDistIndx, allIndx)
            weightsList <- list(cbind(loc=weights1), cbind(dist=weights2), cbind(loc_dist=weights), cbind(loc=weights1, dist=weights2, loc_dist=weights))
            
            msrs2 <- mcmapply(FUN=function(indx, wMat){
              # i <- 2; indx <- indxList[[i]]; wMat <- weightsList[[i]]
              msrs2 <- sapply(indx, function(j){
                # j <- indx[1]
                msrs2 <- apply(wMat, 2, function(col) weighted.var(msrs1[,j], w=col))
                
                return(msrs2)
              }, simplify="array")
              
              if(is.null(dim(msrs2))){
                msrs2 <- matrix(msrs2, 1, length(msrs2))
                rownames(msrs2) <- colnames(wMat)
              }
              colnames(msrs2) <- colnames(msrs1)[indx]
              names(dimnames(msrs2)) <- c("wType","var")
              msrs2 <- melt(msrs2)
              #nms <- paste(msrs2$var, "wVar", msrs2$wType, sep="_") # this is how it shd be i think
              nmAux <- strsplit(nm,"\\.")[[1]][1]
              
              nms <- paste(msrs2$var, nmAux,"wVar", msrs2$wType, sep="_") # doubled upon the ytype-sigytype label to be backwards compatible
              msrs2 <- msrs2$value
              names(msrs2) <- nms
              return(msrs2)
            }, indx=indxList, wMat=weightsList, mc.cores=1)
            
            
            msrs2 <- unlist(msrs2)
            indxs <- mcmapply(function(indxEl, wMat) rep(indxEl, rep(ncol(wMat),length(indxEl))), indxEl=indxList, wMat=weightsList, mc.cores=1)
            indxs <- unlist(indxs)
            msrs2 <- msrs2[order(indxs)]
            
            return(msrs2)
          }, i=1:nrow(df), mc.cores=mc_cores_loc2, SIMPLIFY="array"))
        colnames(msrs2) <- paste(colnames(msrs2), nm, sep="_")  
          
          # plot(msrs2[,"wwkcrdc"], wwKCRDCs_xy1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_xy1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarL2rels"], wwL2_f_xy1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wwkcrdc"], wwKCRDCs_xy0); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_xy0); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarL2rels"], wwL2_f_xy0); abline(a=0, b=1, col="red")
          # 
          # plot(msrs2[,"wwkcrdc"], wwKCRDCs_yx1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_yx1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarL2rels"], wwL2_f_yx1); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wwkcrdc"], wwKCRDCs_yx0); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_yx0); abline(a=0, b=1, col="red")
          # plot(msrs2[,"wVarL2rels"], wwL2_f_yx0); abline(a=0, b=1, col="red")
          
       return(msrs2)
        
      }, df=dfList, msrs1=msrs1List, 
      nm=names(dfList), mc.cores=mc_cores_loc, SIMPLIFY=FALSE)
      
      msrs1List <- mcmapply(FUN=function(msrs1,nm){
        # i <- 1; msrs1 <- msrs1List[[i]]; nm <- names(msrs1List)[i]
        res <- msrs1
        colnames(res) <- paste(colnames(res), nm, sep="_")
        return(res)
      }, msrs1=msrs1List, nm=names(msrs1List), mc.cores=1, SIMPLIFY=FALSE)
      
      
      sapply(msrs1List, dim)
      sapply(msrs2List, dim)
      
      sapply(msrs1List, colnames)
      sapply(msrs2List, colnames)
      
      msrs1 <- do.call(cbind, msrs1List)
      msrs2 <- do.call(cbind, msrs2List)
      
      dim(msrs1)
      dim(msrs2)
      colnames(msrs1)
      colnames(msrs2)
      
      # plot(msrs2[,"wwkcrdc"], wwKCRDCs_xy); abline(a=0, b=1, col="red")
      # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_xy); abline(a=0, b=1, col="red")
      # plot(msrs2[,"wVarL2rels"], wwL2_f_xy); abline(a=0, b=1, col="red")
      # 
      # plot(msrs2[,"wwkcrdc"], wwKCRDCs_yx); abline(a=0, b=1, col="red")
      # plot(msrs2[,"wVarCorrRKHS"], wVarCorrs_yx); abline(a=0, b=1, col="red")
      # plot(msrs2[,"wVarL2rels"], wwL2_f_yx); abline(a=0, b=1, col="red")
      
      res <- cbind(df_raw[c("testPt","boot","indxPt")], as.data.frame(msrs1), as.data.frame(msrs2))
      
      
      print(proc.time()-pm)
      
      return(res)
    })  
    
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  print((proc.time()-pm)[3])
  names(measureList) <- nodes
  
  #plot(measureList$y$`1`$wwkcrdc, wwKCRDCs_xy); abline(a=0, b=1, col="red")
  #plot(measureList$x$`1`$wwkcrdc, wwKCRDCs_yx); abline(a=0, b=1, col="red")
  #sum(measureList$y$`1`$wwkcrdc<measureList$x$`1`$wwkcrdc)/length(measureList$x$`1`$wwkcrdc)
  res <- list(x=x,smpl_tr=smpl_tr, smpl_te=smpl_te, dfs_by_node=dfs_by_node,  dfs_by_node2=dfs_by_node2, measureList=measureList)
  if(jointFeats){
    res <- c(res, dfs_by_node_joint=dfs_by_node_joint)
  }
  return(res)
}

getCMEsGivenDAGSet <- function(uniqueRegsList, x, cmemLearner, noiseLearner=NULL, augmentData=FALSE, dataNm, nm, folderSave, round="last"){
  
  
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  mc_cores <- 8#min(40, detectCores()-4)
  numCoresFold <- 5
  #mc_cores <- 1; numCoresFold <- 1
  
  #count <- 0
  pm0 <- proc.time()
  cmemList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      if(length(indxPreds)==0){
        cmemLearnerAux <- NA
      } else{
        
        if(!is.list(cmemLearner)){
          cmemLearnerAux <- eval(parse(text=cmemLearner))
        } else{
          cmemLearnerAux <- cmemLearner[[nodeTo]][[numReg]]
        }
        #count <- count + 1
        regressors <- nodes[indxPreds] 
        dataX <- x[,regressors, drop=F]
        dataY <- x[,nodeTo]

        # hack to try transorming (cause, effect) -> (copula(cause), effect)
        if(FALSE){
          #print("enters copula cause hack")
          #print(paste("dim(dataX): ", dim(dataX)))
          xnew <- dataX
          xnew <- apply(xnew,2, function(col){
            o <- order(col)
            aux <- (1:length(o))/length(o)
            aux <- aux[order(o)]
            # so that its between -0.5 and 0.5
            res <- aux-0.5
            return(res)
          })
          #print(paste("dim(xnew): ", dim(xnew)))
          #print(summary(xnew))
          
          colnames(xnew) <- colnames(dataX)
          # sample 100 pts
          set.seed(12345)
          smpl <- sample(1:nrow(xnew), 100, replace=F)
          dataX <- xnew[smpl,,drop=F]
          dataY <- dataY[smpl]
          #print(paste("dim(dataX): ", dim(dataX)))
        }
        
        # hack to try transorming (cause, effect) -> (resample_unif(cause), effect)
        if(FALSE){
          xnew <- dataX
          ynew <- dataY
          #plot(dataX[,1], dataY[,2])
          #hist(xnew[,1])
          for(i in 1:10){
            
            #print(paste("i: ", i))
            #print(paste("nrow(xnew): ", nrow(xnew)))
            #print(paste("length(ynew): ", length(ynew)))
            
            
            mod <- kepdf(xnew, eval.points = xnew, kernel = "gaussian", bwtype = "adaptive")@estimate
            
            probs <- 1/mod  
            probs <- probs/sum(probs)
            
            set.seed(12345)
            smpl <- sample(1:nrow(xnew), size=10000, prob=probs, replace=T)
            xnew <- xnew[smpl,,drop=F]
            ynew <- ynew[smpl]
            
            #plot(x[,1],x[,2])
            #plot(xnew[,1], xnew[,2])
            #hist(x[,1])
            #hist(xnew[,1])
            #pval <- ks.test(xnew[,1], y="punif", min=min(xnew[,1]), max=max(xnew[,1]))$p.value
            #print(paste("pval: ", pval))
          }
          
          
          set.seed(1234)
          smpl <- sample(1:nrow(xnew), size=100)
          dataX <- xnew[smpl,,drop=F]
          dataY <- ynew[smpl]
          #print(paste("nrow(dataX): ", nrow(dataX)))
          #print(paste("length(dataY): ", length(dataY)))
        }
        
        regressorsChar <- paste(regressors, collapse="-")
        regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
        fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
        fileSaveAugData <- paste(dataNm, "_","augData","_",regressionChar, ".RData", sep="")
        fileSaveNoiseLearner <- paste(dataNm, "_", "noisLearner", "_", regressionChar, ".RData", sep="")
        trainDataO <- constructData(x=dataX, y=dataY)
        
        if(augmentData){
          if(!fileSaveAugData %in% dir(folderSave)){
            print("imputes aug data")
            # augment data using noiseLearner, save trainData to cmemLearnerAux, write noiseLearner
            noiseLearner1 <- eval(parse(text=noiseLearner))
            print("set noiseLearner params")
            pm <- proc.time()
            noiseLearner1 <- setParams(learner=noiseLearner1, trainData=trainDataO, mc_cores=1)
            pm <- proc.time()-pm
            print(paste("set params took: ", round(pm["elapsed"]/60/60,2), " hours", sep=""))
            print("find best initialization")
            pm <- proc.time()
            noiseLearner1 <- noiseLearner1$learn(learner=noiseLearner1, inits=TRUE)
            pm <- proc.time()-pm
            print(paste("finding best ini took: ", round(pm["elapsed"]/60/60,2), " hours", sep=""))
            rhat <- noiseLearner1$hyperParams$trainData$r
            # we dont need to normalize per the multi output kernel ridge reg model
            # only the y var must be centered
            #rhat <- apply(rhat, 2, norml)
            dataX <- cbind(dataX, rhat)
            trainDataO <- constructData(x=dataX, y=dataY)
            save(trainDataO, file=paste(folderSave, fileSaveAugData, sep=""))
            save(noiseLearner1, file=paste(folderSave, fileSaveNoiseLearner, sep=""))
          } else{
            print("reads aug data")
            load(file=paste(folderSave, fileSaveAugData, sep=""))
          }
        } 
            
   
        
        if(fileSaveCmemLearner %in% dir(folderSave)){
          print("reads")
          cmemLearnerAuxOrig <- cmemLearnerAux 
          load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
          
          # in the first round if a learner with the same alpha (L2, Kcmc, kcsc) and
          # learner params has already been run (irrespective of the sigmay criterion)
          # we can use it and just change the criterion (if all relevant critera are
          # calculated the first time -> check learner)
          if(round=="first"){
            cmemLearnerAux$optimizeParams$mainLoss <- cmemLearnerAuxOrig$optimizeParams$mainLoss
            cmemLearnerAux$optimizeParams$optLossFunc <- cmemLearnerAuxOrig$optimizeParams$optLossFunc
            print("re-sets optimal params based on ready made grid and right mainLoss")
            pm <- proc.time()
            cmemLearnerAux <- setParams(learner=cmemLearnerAux, trainData=trainDataO, mc_cores=mc_cores, numCoresFold=numCoresFold)
            pm <- proc.time() - pm
            print(paste("shd be quick: ", round(pm["elapsed"],3), " secs"))
          }
          
        } else{
          print("calculates")
          if(round=="first"){
            # if we're running this function in a two-time scheme dont bother calculating kernel grads
            # calculating regularization losses or measures
            losses <- c("RMSE2","negCE", "gauss_log_lik", "cmem_L2_f", "cmem_L2_f_rel", "cmem_L2_k","corrRKHS","wcorrRKHS")
            msrs <- c("gll_tr", "gll_te", "L2_f_tr", "L2_f_te", "RMSEte","PCEtr", "PCEte", "CCRtr", "CCRte","corrRKHS_k_te", "wcorrRKHS_k_te")
            lossesOr <- names(cmemLearnerAux$optimizeParams$losses) 
            msrsOr <- names(cmemLearnerAux$msrs)
            cmemLearnerAux$optimizeParams$losses <- cmemLearnerAux$optimizeParams$losses[which(lossesOr %in% losses)]
            cmemLearnerAux$msrs <- cmemLearnerAux$msrs[which(msrsOr %in% msrs)] 
            
            cmemLearnerAux$hyperParams$data$non_optimizable$gradX$val <- FALSE
          } else{
            lossesOr <- names(cmemLearnerAux$optimizeParams$losses)
            lossesDont <- c("negCE", "gauss_log_lik")
            indxDont <- which(lossesOr %in% lossesDont)
            if(length(indxDont)>0) cmemLearnerAux$optimizeParams$losses <- cmemLearnerAux$optimizeParams$losses[-indxDont]
            
          }
          #learner=cmemLearnerAux;trainData=trainDataO;
          cmemLearnerAux <- setParams(learner=cmemLearnerAux, trainData=trainDataO, mc_cores=mc_cores, numCoresFold=numCoresFold)
          
          save(cmemLearnerAux, file=paste(folderSave, fileSaveCmemLearner, sep=""))  
        }
        
        
        #print(dimnames(cmemLearnerAux$hyperParams$data$grid)$var)
        
        cmemLearnerAux <- cmemLearnerAux$learn(learner=cmemLearnerAux, forLoss=F)
        
        
      }
      
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(cmemLearnerAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  print(proc.time()-pm0) #  
  
  names(cmemList) <- nodes
  
  return(cmemList)
}

getMeasuresGivenCMEsList <- function(uniqueRegsList, x, cmeSet, cmemLearner){
  
  n <- dim(x)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  cmemLearner <- eval(parse(text=cmemLearner))
  
  #count <- 0
  #pm0 <- proc.time()
  measureList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))

      cmemLearnerAux <- cmeSet[[nodeTo]][[numReg]] 
      
      if(is.na(cmemLearnerAux)){
        indxAux2 <- which(sapply(unlist(cmeSet, recursive=F), function(el) !is.na(el[1])))[1]
        cmemLearnerAux2 <- unlist(cmeSet, recursive=F)[[indxAux2]]
        cmems <- rep(NA, length(cmemLearnerAux2$msrs))
        names(cmems) <- names(cmemLearnerAux2$msrs)
      } else{
      
        cmemLearnerAux <- cmemLearnerAux$learn(learner=cmemLearnerAux, forLoss=F)
        cmemLearnerAux$msrs <- cmemLearner$msrs
        cmems <- cmemLearner$calcMsrs(cmemLearner=cmemLearnerAux) # learner <- cmemLearnerAux
      }
      #print("estimated time to completion:")
      #numRegsLeft <- numRegsTot - count
      #numRegsDone <- numRegsTot -numRegsLeft
      #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      
      return(cmems)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(measureList) <- nodes
  
  return(measureList)
}

logSum <- function(x, na.rm) sum(log(x), na.rm=T)

assembleDagScores <- function(dags, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE){
  
  
  cmemLearner <- eval(parse(text=cmemLearner))
  measurePack <- cmemLearner$msrs
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  scores <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        
      } else{
        indxReg <- 1
        correct <- NA
        
      }
      
      fac <- cmemSet[[nodeTo]][[indxReg]]*correct
      return(fac)
    }, simplify="array")
    
    names(dimnames(facs)) <- c("measures","nodeFactors")
    
    
    indx <- match(sapply(strsplit(dimnames(facs)$measures,"_"), function(el) el[1]), sapply(strsplit(names(measurePack),"_"), function(el) el[1]))
    measurePack <- measurePack[indx]
    
    aggFuncs <- sapply(measurePack, function(el) el$aggFunc)
    
    
    score <- sapply(1:length(aggFuncs), function(i){
      # i <- 1
      apply(facs[i,,drop=F], 1, aggFuncs[i], na.rm=T)
    })
    
    
    #colnames(score) <- dimnames(facs)$measures
    #names(dimnames(score)) <- names(dimnames(facs))[1]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  
  
  names(dimnames(scores)) <- c("score","dag") #"measures"
  
  
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

assembleDagScores_boot <- function(dags, uniqueRegsList, cmemSet, prnt=FALSE){
  
  
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  scoresNodeWise <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        
      } else{
        indxReg <- 1
        correct <- NA
        
      }
      
      fac <- cmemSet[[nodeTo]][[indxReg]]*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    q90 <- function(x, ...){
      res <- quantile(x, probs=0.9, ...)
      names(res) <- NULL
      return(res)
    }
    q100 <- function(x, ...){
      res <- quantile(x, probs=1, ...)
      names(res) <- NULL
      return(res)
    }
    aggFuncs <- list(corrRKHS=c("mean"), kcrdc=c("var"), wkcrdc=c("var"), 
                     l2rel=c("mean","var"), wl2rel=c("mean","var"),
                     wwkcrdc=c("median","q90","q100"), wVarCorrRKHS=c("mean"), wVarL2rels=c("mean"))
    
    
    scoreNodeWise <- lapply(1:length(aggFuncs), function(i){
      # i <- 4
      res <- sapply(1:length(aggFuncs[[i]]), function(j) apply(facs[,i,,drop=F], c(2,3), aggFuncs[[i]][j], na.rm=T))
      colnames(res) <- paste(names(aggFuncs)[i], aggFuncs[[i]], sep="_")
      rownames(res) <- dimnames(facs)$nodeFactors
      return(res)
    })
    scoreNodeWise <- do.call(cbind,scoreNodeWise)
    scoreNodeWise[is.nan(scoreNodeWise)] <- NA
    scoreNodeWise <- apply(scoreNodeWise, 2, sum, na.rm=T)
    
    
    
    
    
    #colnames(score) <- dimnames(facs)$measures
    #names(dimnames(score)) <- names(dimnames(facs))[1]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(scoreNodeWise)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  names(dimnames(scoresNodeWise)) <- c("score","dag") #"measures"
  scoresNodeWise <- aperm(scoresNodeWise, c(2,1))
  
  scoresLocal <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        
      } else{
        indxReg <- 1
        correct <- NA
        
      }
      
      fac <- cmemSet[[nodeTo]][[indxReg]]*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    facs <- apply(facs, c(1,2), sum, na.rm=T)
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(facs)
  }, simplify="array")
  names(dimnames(scoresLocal))[3] <- "dag"
  
  #difference between the score of a dag-pt and the min of the rest of the dags
  scoresLocal2 <- apply(scoresLocal, 2, function(mat){
    # mat <- scoresLocal[,1,]
    difs <- sapply(1:nrow(mat), function(i) sapply(1:ncol(mat), function(j) mat[i,j]-min(mat[i,setdiff(1:ncol(mat),j)])))
    difs <- t(difs)
    res <- apply(difs, 2, function(col) c(sum(col<0)/length(col), sum(col)))
    return(res)
  })
  dim(scoresLocal2) <- c(2,2,dim(scoresLocal2)[2])
  dimnames(scoresLocal2) <- list(aggFunc=c("pctPos","sum"), dag=dimnames(scoresLocal)$dag ,measures=dimnames(scoresLocal)$measures)
  scoresLocal3 <- melt(scoresLocal2)
  scoresLocal4 <- cast(scoresLocal3, dag~measures+aggFunc, value="value")
  scoresLocal5 <- scoresLocal4[,2:ncol(scoresLocal4)]
  rownames(scoresLocal5) <- scoresLocal4$dag
  names(dimnames(scoresLocal5)) <- c("dag","score")
  scores <- cbind(scoresNodeWise, scoresLocal5[match(rownames(scoresLocal5), rownames(scoresNodeWise)),])
  return(scores)
}

assembleDagScores_boot_eqSig_old_vertical <- function(dags, uniqueRegsList, cmemSet, prnt=FALSE){
  
  #table(cmemSet$x$`1`$lambda_nrm,cmemSet$y$`1`$lambda_nrm)
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  nodeRegs <- apply(cmemSet[[1]][[1]][,c("node","numReg")],1,  paste, collapse=".")
  
  print("scoresNodeWise")
  scoresNodeWise <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        
      } else{
        indxReg <- 1
        correct <- NA
        
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      cmemSet2 <- cmemSet2[,-which(colnames(cmemSet2)%in%c("node","numReg","testPt","boot","indxPt","principal"))]
      
      fac <- as.matrix(cmemSet2)*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    q90 <- function(x, ...){
      res <- quantile(x, probs=0.9, ...)
      names(res) <- NULL
      return(res)
    }
    q100 <- function(x, ...){
      res <- quantile(x, probs=1, ...)
      names(res) <- NULL
      return(res)
    }
    aggFuncs <- c("median","var","q90")
    
    # obtain the agg measure  - separateley for each node/sigmay pair
    scoreNodeWise <- sapply(unique(nodeRegs), function(nd){
      # nd <- unique(nodeRegs)[1]
      ks <- which(nodeRegs==nd)
      # to each factor apply all factor-specific functions 
      res <- apply(facs[ks,,,drop=F], c(2,3), function(vec) sapply(aggFuncs, function(func) do.call(func, list(x=vec, na.rm=T))))
      
      names(dimnames(res))[1] <- "stat"
      return(res)
    }, simplify="array")
    
    scoreNodeWise <- apply(scoreNodeWise, c(1,2,3), mean, na.rm=T)
    scoreNodeWise[is.nan(scoreNodeWise)] <- NA
    
    
    scoreNodeWise <- apply(scoreNodeWise, c(1,2), sum, na.rm=T)
    
    scoreNodeWise <- melt(scoreNodeWise)
    scoreNodeWise <- cast(scoreNodeWise, ~measures+stat, value="value")
    scoreNodeWise <- as.matrix(scoreNodeWise)[1,]
    
    
    #colnames(score) <- dimnames(facs)$measures
    #names(dimnames(score)) <- names(dimnames(facs))[1]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(scoreNodeWise)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  names(dimnames(scoresNodeWise)) <- c("score","dag") #"measures"
  scoresNodeWise <- aperm(scoresNodeWise, c(2,1))
  
  print("scores local")
  scoresLocal <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      } else{
        indxReg <- 1
        correct <- NA
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      cmemSet2 <- cmemSet2[,-which(colnames(cmemSet2)%in%c("node","numReg","testPt","boot","indxPt","principal"))]
      
      fac <- as.matrix(cmemSet2)*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    facs <- apply(facs, c(1,2), sum, na.rm=T)
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(facs)
  }, simplify="array")
  names(dimnames(scoresLocal))[3] <- "dag"
  
  #difference between the score of a dag-pt and the min of the rest of the dags
  scoresLocal2 <- apply(scoresLocal, 2, function(mat){
    # mat <- scoresLocal[,3,]
    difs <- sapply(1:nrow(mat), function(i) sapply(1:ncol(mat), function(j) mat[i,j]-min(mat[i,setdiff(1:ncol(mat),j)])))
    difs <- t(difs)
    res <- apply(difs, 2, function(col) c(-sum(col<0)/length(col), sum(col)))
    return(res)
  })
  dim(scoresLocal2) <- c(2,2,dim(scoresLocal2)[2])
  dimnames(scoresLocal2) <- list(aggFunc=c("pctPos","sum"), dag=dimnames(scoresLocal)$dag ,measures=dimnames(scoresLocal)$measures)
  scoresLocal3 <- melt(scoresLocal2)
  scoresLocal4 <- cast(scoresLocal3, dag~measures+aggFunc, value="value")
  scoresLocal5 <- scoresLocal4[,2:ncol(scoresLocal4)]
  rownames(scoresLocal5) <- scoresLocal4$dag
  names(dimnames(scoresLocal5)) <- c("dag","score")
  scores <- cbind(scoresNodeWise, scoresLocal5[match(rownames(scoresLocal5), rownames(scoresNodeWise)),])
  
  #Now re ordering according to best node
  print("scores local- best")
  scoresLocal <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      } else{
        indxReg <- 1
        correct <- NA
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      cmemSet2 <- cmemSet2[which(cmemSet2$principal),]
      cmemSet2 <- cmemSet2[order(cmemSet2$testPt, cmemSet2$boot),]
      cmemSet2 <- cmemSet2[,-which(colnames(cmemSet2)%in%c("node","numReg","testPt","boot","indxPt","principal"))]
      
      fac <- as.matrix(cmemSet2)*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    facs <- apply(facs, c(1,2), sum, na.rm=T)
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(facs)
  }, simplify="array")
  names(dimnames(scoresLocal))[3] <- "dag"
  
  #difference between the score of a dag-pt and the min of the rest of the dags
  scoresLocal2 <- apply(scoresLocal, 2, function(mat){
    # mat <- scoresLocal[,3,]
    difs <- sapply(1:nrow(mat), function(i) sapply(1:ncol(mat), function(j) mat[i,j]-min(mat[i,setdiff(1:ncol(mat),j)])))
    difs <- t(difs)
    res <- apply(difs, 2, function(col) c(-sum(col<0)/length(col), sum(col)))
    return(res)
  })
  dim(scoresLocal2) <- c(2,2,dim(scoresLocal2)[2])
  dimnames(scoresLocal2) <- list(aggFunc=c("pctPos","sum"), dag=dimnames(scoresLocal)$dag ,measures=dimnames(scoresLocal)$measures)
  scoresLocal3 <- melt(scoresLocal2)
  scoresLocal4 <- cast(scoresLocal3, dag~measures+aggFunc, value="value")
  scoresLocal5 <- scoresLocal4[,2:ncol(scoresLocal4)]
  rownames(scoresLocal5) <- scoresLocal4$dag
  colnames(scoresLocal5) <- paste(colnames(scoresLocal5),"best",sep="_")
  names(dimnames(scoresLocal5)) <- c("dag","score")
  scores <- cbind(scores, scoresLocal5[match(rownames(scoresLocal5), rownames(scores)),])
  
  
  
  return(scores)
}


assembleDagScores_boot_eqSig <- function(dags, uniqueRegsList, cmemSet, prnt=FALSE){
  print("enters assembleDagScores_boot_eqSig")
  #table(cmemSet$x$`1`$lambda_nrm,cmemSet$y$`1`$lambda_nrm)
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  print("scoresNodeWise")
  scoresNodeWise <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        
      } else{
        indxReg <- 1
        correct <- NA
        
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      cmemSet2 <- cmemSet2[,-which(colnames(cmemSet2)%in%c("testPt","boot","indxPt"))]
      
      fac <- as.matrix(cmemSet2)*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    q90 <- function(x, ...){
      res <- quantile(x, probs=0.9, ...)
      names(res) <- NULL
      return(res)
    }
    q100 <- function(x, ...){
      res <- quantile(x, probs=1, ...)
      names(res) <- NULL
      return(res)
    }
    aggFuncs <- c("median","var","q90")
    
    # obtain the agg measure  - separateley for each node/sigmay pair
    scoreNodeWise <- apply(facs, c(2,3), function(vec) sapply(aggFuncs, function(func) do.call(func, list(x=vec, na.rm=T))))
    
    names(dimnames(scoreNodeWise))[1] <- "stat"
    
    scoreNodeWise[is.nan(scoreNodeWise)] <- NA
    
    
    scoreNodeWise <- apply(scoreNodeWise, c(1,2), sum, na.rm=T)
    
    scoreNodeWise <- melt(scoreNodeWise)
    scoreNodeWise <- cast(scoreNodeWise, ~measures+stat, value="value")
    scoreNodeWise <- as.matrix(scoreNodeWise)[1,]
    
    
    #colnames(score) <- dimnames(facs)$measures
    #names(dimnames(score)) <- names(dimnames(facs))[1]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(scoreNodeWise)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  names(dimnames(scoresNodeWise)) <- c("score","dag") #"measures"
  scoresNodeWise <- aperm(scoresNodeWise, c(2,1))
  
  print("scores local")
  scoresLocal <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      correct <- 1
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      } else{
        indxReg <- 1
        correct <- NA
      }
      
      cmemSet2 <- cmemSet[[nodeTo]][[indxReg]]
      cmemSet2 <- cmemSet2[,-which(colnames(cmemSet2)%in%c("testPt","boot","indxPt"))]
      
      fac <- as.matrix(cmemSet2)*correct
      return(fac)
    }, simplify="array")
    
    dimnames(facs)[[1]] <- 1:(dim(facs)[1])
    names(dimnames(facs)) <- c("pts","measures","nodeFactors")
    
    facs <- apply(facs, c(1,2), sum, na.rm=T)
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(facs)
  }, simplify="array")
  names(dimnames(scoresLocal))[3] <- "dag"
  
  #difference between the score of a dag-pt and the min of the rest of the dags
  print("scores local 2")
  scoresLocal2 <- apply(scoresLocal, 2, function(mat){
    # mat <- scoresLocal[,3,]
    difs <- sapply(1:nrow(mat), function(i) sapply(1:ncol(mat), function(j) mat[i,j]-min(mat[i,setdiff(1:ncol(mat),j)])))
    difs <- t(difs)
    res <- apply(difs, 2, function(col) c(-sum(col<0)/length(col), sum(col)))
    return(res)
  })
  dim(scoresLocal2) <- c(2,2,dim(scoresLocal2)[2])
  dimnames(scoresLocal2) <- list(aggFunc=c("pctPos","sum"), dag=dimnames(scoresLocal)$dag ,measures=dimnames(scoresLocal)$measures)
  print("scores local 3")
  scoresLocal3 <- melt(scoresLocal2)
  print("scores local 4")
  scoresLocal4 <- cast(scoresLocal3, dag~measures+aggFunc, value="value")
  scoresLocal5 <- scoresLocal4[,2:ncol(scoresLocal4)]
  rownames(scoresLocal5) <- scoresLocal4$dag
  names(dimnames(scoresLocal5)) <- c("dag","score")
  scores <- cbind(scoresNodeWise, scoresLocal5[match(rownames(scoresLocal5), rownames(scoresNodeWise)),])
  
  
  print("exits assembleDagScores_boot_eqSig")
  return(scores)
}


assembleDagCurvesIntoDagCurvesThenScores <- function(dags, uniqueRegsList, curves, reglrs,  prnt=FALSE, plot=NULL){
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  
  
  curves_dags <- lapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain curves that will make up generla curve of dag hypothesis
    # one component per node
    curves_components <- lapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[5])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        curve <- curves[[nodeTo]][[indxReg]]
        
      } else{
        curve <- NA
      }
      
      return(curve)
    })
    #names(curves_components[[1]][["KCDC"]])
    #curves_components[[1]][["KCDC"]]$reglrRng
    #curves_components[[1]][["KCDC"]]$lossRng
    
    curves_components <- curves_components[which(!sapply(curves_components, function(el) is.na(el[1])))]
    
    
    
    curve_dag <- estimateOverallCurve(curves=curves_components, reglrs)
    names(curve_dag[["KCDC"]])
    #curve_dag[["KCDC"]]$reglrRng
    #curve_dag[["KCDC"]]$lossRng
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(curve_dag)
  })
  # length(curves_dags)
  # names(curves_dags[[1]][["KCDC"]])
  # curves_dags[[1]][["KCDC"]]$reglrRng
  # curves_dags[[1]][["KCDC"]]$lossRng
  
  
  if(!is.null(plot)){
    curves_dags_plot <- lapply(curves_dags, function(curve) curve[[plot]])  
    plotCurvesSimple(curves_dags_plot)
  } 
  
  # strictly for plotting
  
  # curvesPlot <- lapply(dagsList, function(dag){
  #   #get residuals for the four regressions implied in each column
  #   # i <- 2; dag <- dagsList[[i]]
  #   
  #   count <<- count + 1
  #   if(prnt){
  #     print("*******************")
  #     print(paste("dag # ",count, sep=""))
  #   }
  #   
  #   # obtain curves that will make up generla curve of dag hypothesis
  #   curves_components <- lapply(nodes, function(nodeTo){
  #     # (nodeTo <- nodes[5])
  #     # print(paste("nodeTo: ", nodeTo))
  #     parents <- dag[,nodeTo]
  #     
  #     if(sum(parents)>0){
  #       indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
  #       curve <- curves[[nodeTo]][[indxReg]]
  #       
  #     } else{
  #       curve <- NA
  #     }
  #     
  #     return(curve)
  #   })
  #   
  #   curves_components <- curves_components[which(!sapply(curves_components, function(el) is.na(el[1])))]
  #   
  #   
  #   return(curves_components)
  # })
  # curvesPlot <- lapply(curvesPlot, function(el) el[[1]])
  # getMeasures_comp(curves=curvesPlot, reglrs)
  
  
  scores <- sapply(reglrs, function(reglr){
    # reglr <- reglrs[1]
    curves_dags_reglr <- lapply(curves_dags, function(el) el[[reglr]])
    # length(curves_dags_reglr)
    # names(curves_dags_reglr[[1]])
    # curves_dags_reglr[[1]]$reglrRng
    # curves_dags_reglr[[1]]$lossRng
    # plotCurvesSimple(curveList=curves_dags_reglr)
    
    res <- getComparableReglrVal(curveList=curves_dags_reglr)
    
    
    res <- cbind(reglr=res$logReglrs, loss=res$losses)
    rownames(res) <- names(curves_dags_reglr)
    return(res)
  }, simplify="array")
  
  names(dimnames(scores)) <- c("dag","regLoss","msr")
  scores <- melt(scores)
  scores <- cast(scores, dag~msr+regLoss, value="value")
  rownames(scores) <- scores$dag
  scores <- as.matrix(as.data.frame(scores[,2:ncol(scores)]))
  c.nms <- colnames(scores)
  indxChng <- grep("loss",c.nms)
  c.nms[indxChng] <- paste("L2", sapply(strsplit(c.nms[indxChng],"_"), function(el) el[1]), sep="_")
  colnames(scores) <- c.nms
  
  names(dimnames(scores)) <- c("dag","score") #"measures"
  
  
  if(prnt) print(proc.time()-pm0) #
  
  
  
  return(scores)
}



estimateOverallCurve <- function(curves, reglrs){
  # for each regularizer
  
  curve_overall <- lapply(reglrs, function(reglr){
    # reglr <- reglrs[1]
    # 1. obtain initial (complexity, loss) point: sum min complexities and max losses.
    # 2. obtain final (complexity, loss) point: sum max complexities and min losses
    #print(reglr)
    lossRngs <- sapply(curves, function(el) el[[reglr]]$lossRng)
    reglrRngs <- sapply(curves, function(el) el[[reglr]]$reglrRng)
    #lossIniFin  <- apply(lossRngs, 1, sum) 
    #reglrIniFin <- apply(reglrRngs, 1, sum)
    
    
    #3.  sample n x p 
    n <- 10000
    p <- length(curves)
    
    reglr_smpl <- matrix(rtri(n*p, min=0, max=1, mode=0.99999), n, p)
    # make sure max and min of each curve will be in sample
    reglr_smpl <- rbind(reglr_smpl, rep(0, p), rep(1, p))
    
    #3. Convert along columns into reglr in Support_j j=1,...,p
    reglr_smpl <- sapply(1:p, function(j){
        # j <- 1
        res <- reglr_smpl[,j]
        res <- res*(reglrRngs[2,j]-reglrRngs[1,j])+reglrRngs[1,j]
        # range(res); reglrRngs[,j]
        return(res)
    }, simplify="array")
    
    loss_smpl <- sapply(1:p, function(j){
      # j <- 1
      res <- curves[[j]][[reglr]]$fun_xy(reglr_smpl[,j])
      return(res)
    }, simplify="array")
    
    reglr_smpl <- apply(reglr_smpl, 1, function(row) log(sum(10^row)))
    loss_smpl <- apply(loss_smpl, 1, sum)
    
    # plot(reglr_smpl, loss_smpl)
    
    envPts <- extractLeftEnvPts(reglr_smpl, loss_smpl)
    if(length(envPts$x)==1){
      envPts$x <- c(envPts$x, envPts$x+0.00001)
      envPts$y <- c(envPts$y, envPts$y-0.00001)
    }
    #plot(reglr_smpl, loss_smpl)
    #lines(envPts$x, envPts$y, type="p", cex=2, col="red")
    splf_yx <- splinefun(envPts$y, envPts$x, method="monoH.FC")
    splf_xy <- splinefun(envPts$x, envPts$y, method="monoH.FC")
    #yy_1 <- seq(envPts$y[1], envPts$y[length(envPts$y)], length.out=100)
    #xx_1 <- splf_yx(yy_1)
    #xx_2 <- seq(envPts$x[1], envPts$x[length(envPts$x)], length.out=100)
    #yy_2 <- splf_xy(xx_2)
    #plot(reglr_smpl, loss_smpl)
    #lines(xx_1, yy_1, col="red", lty=2)
    #lines(xx_2, yy_2, col="blue", lty=3)
    curve <- list(reglrRng=range(envPts$x), lossRng=range(envPts$y), fun_xy=splf_xy, fun_yx=splf_yx, loss=curves[[1]][[reglr]]$loss, reglr=reglr, logReglr=curves[[1]][[reglr]]$logReglr)
      
    
    return(curve)
  
  })
  
  names(curve_overall) <- reglrs
  return(curve_overall)
}


getLBK <- function(cmem){
  nms <- names(cmem$learnParams)
  Lx <- cmem$learnParams$Lx
  Ky <- cmem$learnParams$Ky
  Blambda <- cmem$learnParams$Blambda
  #Cks <- cmem$learnParams$Cks
  if("phiy" %in% nms){
    phiy <- cmem$learnParams$phiy
    res <- list(Lx=Lx, phiy=phiy, Ky=Ky, Blambda=Blambda) #, Cks=Cks
  } else{
    res <- list(Lx=Lx, Ky=Ky, Blambda=Blambda) # , Cks=Cks
  }
  return(res)
} 

buildSuperCMEM <- function(LBKs, cmemLearner, x, xAug){
  nms <- names(cmemLearner$learnParams) #names(LBKs[[1]])
  indx <- which(!sapply(LBKs, function(el) is.na(el[1])))
  LBKs <- LBKs[indx] 
  Lxs <- lapply(LBKs, function(el) el$Lx)
  Kys <- sapply(LBKs, function(el) el$Ky, simplify="array")
  Blambdas <- lapply(LBKs, function(el) el$Blambda)
  
  # Cks must be Rn x Rn x n dim where R is the number of nodes with at least one parent in the graph
  #Ckss <- lapply(LBKs, function(el) el$Cks)
  
  if("phiy" %in% nms){
    phiys <- lapply(LBKs, function(el) el$phiy)
    Lx <- do.call(cbind, Lxs)
    Ky <- apply(Kys, c(1,2), sum)
    phiy <- do.call(cbind, phiys)
    #Ky2 <- phiy %*% t(phiy)
    #dim(Ky); dim(Ky)
    #plot(Ky, Ky2); abline(a=0, b=1, col="red")  
    Blambda <-  do.call(rbind, Blambdas)
    #Cks <- abind(Ckss, along=2)
    cmemLearner$learnParams$phiy <- phiy
    cmemLearner$learnParams$Lx <- Lx
    cmemLearner$learnParams$Ky <- Ky
    cmemLearner$learnParams$Blambda <- Blambda
    #cmemLearner$learnParams$Cks <- Cks
  } else{
    Lx <- do.call(cbind, Lxs)
    Ky <- apply(Kys, c(1,2), sum)
    Blambda <-  do.call(rbind, Blambdas)
    #Cks <- abind(Ckss, along=2)
    cmemLearner$learnParams$Lx <- Lx
    cmemLearner$learnParams$Ky <- Ky
    cmemLearner$learnParams$Blambda <- Blambda
    #cmemLearner$learnParams$Cks <- Cks
  }
  
  parsX <- getKernelPars(cmemLearner, kernelName="kernelX")
  kernelX_char <- cmemLearner$hyperParams$non_data$kernelX$name
  kernelX_char <- strsplit(kernelX_char, split="T")[[1]][1]
  Cks <- kernelGradNormMatrix(kernel=kernelX_char, x=xAug, y=xAug, z=x,  pars=parsX)
  # any NA was generated in previous step by calculating dot products between variables where 
  # one was not in one of the nodes making up the graph so the gradient is zero so we force it here
  # (forcing it before is not possible as for some kernels k(0,0) != 0)
  Cks[is.na(Cks)] <- 0
  cmemLearner$learnParams$Cks <- Cks
  
  return(cmemLearner)
}

assembleDagCMEM <- function(x, dags, uniqueRegsList, cmemSet, cmemLearner, prnt=FALSE){
  
  
  cmemLearner <- eval(parse(text=cmemLearner))
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  
  scores <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # i <- 1; dag <- dagsList[[i]]
    
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    LBKs <- lapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      # print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      
      if(sum(parents)>0){
        indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
        # names(cmemSet[[nodeTo]])
        
        cmem <- cmemSet[[nodeTo]][[indxReg]]
        LBK <- getLBK(cmem)
        
      } else{
        LBK <- NA
      }
      return(LBK)
    })
    names(LBKs) <- nodes
    
    
    xAug <- lapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[4])
       #print(paste("nodeTo: ", nodeTo))
      parents <- dag[,nodeTo]
      if(sum(parents)>0){
        xDag <- x#t(t(x)*parents)
        xDag[,which(parents==0)] <- NA 
      } else{
        xDag <- matrix(0, 0, ncol(x))
        colnames(xDag) <- colnames(x)
      }
      return(xDag)
    })
    
    names(xAug) <- nodes
    xAug <- do.call(rbind, xAug)
    
    cmemLearnerAux <- cmemLearner
  
    superCmemLearner <- buildSuperCMEM(LBKs, cmemLearner=cmemLearnerAux, x, xAug)
    
    
    score <- superCmemLearner$calcMsrs(cmemLearner=superCmemLearner)
    
    
    #colnames(score) <- dimnames(facs)$measures
    #names(dimnames(score)) <- names(dimnames(facs))[1]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  
  names(dimnames(scores)) <- c("score","dag") #"measures"
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

# Get unique (index) list of nodes with no parents in a set of dags
getUniqueNoParentsList <- function(dags){
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  noParents <- lapply(dagsList, function(mat) which(apply(mat, 2, function(col) all(col==0))))
  uniqueNoParents <- unique(noParents)
  return(uniqueNoParents)
}

directSumCMEMlearner_comp_deprecated <- function(cmemLearner, cmeSet, nodeIndx, regIndx){
  
  cmemLearnerSum <- eval(parse(text=cmemLearner))
  parsY <- mapply(function(nodeTo, numReg){
    #cme <- cmeSet[["1"]][["1"]]
    cme <- cmeSet[[nodeTo]][[numReg]]
    indxPars <- which(names(cme$hyperParams$data$optimizable) %in% names(cmemLearnerSum$hyperParams$data$non_optimizable))
    nmPars <- names(cme$hyperParams$data$optimizable)[indxPars]
    # nm <- nmPars[1]
    #res <- lapply(nmPars, function(nm) getHyperPar(cmemAux, nm))
    # ok just gonna do it for one par! for now assuming rbf!!!
    res <- getHyperPar(cme, nmPars)
    return(res)
  }, nodeTo=nodeIndx, numReg=regIndx, SIMPLIFY=TRUE)
  
  parsY <- sort(unique(parsY))
  indxPars <- which(names(cmeSet[[1]][[1]]$hyperParams$data$optimizable) %in% names(cmemLearnerSum$hyperParams$data$non_optimizable))
  nmPars <- names(cmeSet[[1]][[1]]$hyperParams$data$optimizable)[indxPars]
  cmemLearnerSum$hyperParams$data$non_optimizable[[nmPars]]$val <- parsY
  
  return(cmemLearnerSum)
}
directSumCMEMlearner_comp <- function(cmemLearner, cmeSet, nodeIndx, regIndx){
  
  cmemLearnerSum <- eval(parse(text=cmemLearner))
  
  # here we cant apply more than one classifier coz we rely on the y parameters
  # from previous step
  
  aux <- cmemLearnerSum$optimizeParams$mainLoss
  mainLoss <- strsplit(aux, "\\.")[[1]][1]
  if(mainLoss=="negCE"){
    aux <- strsplit(aux,"\\.")[[1]][2]
    aux <- strsplit(aux, "-")[[1]]
    featFunc <- aux[1]
    classifier <- aux[2]
    cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$featFunc <- featFunc
    cmemLearnerSum$hyperParams$data$non_optimizable$NCE_learner$val$classifier <- classifier
  }
  # get the kernels and the parameters
  
  parsY <- mapply(function(nodeTo, numReg){
    # i <- 1; nodeTo <- nodeIndx[i]; numReg <- regIndx[i]; cme <- cmeSet[[nodeTo]][[numReg]]
    cme <- cmeSet[[nodeTo]][[numReg]]
    
    aux1 <- (!is.null(getHyperPar(cme, "kernelY")))*1
    aux2 <- (!is.null(getHyperPar(cme, "featureY")))*2
    kernFeatY <- c("kernelY","featureY")[aux1+aux2]
    krnY <- getHyperPar(cme, kernFeatY)
    krnY_regexpr <- strsplit(krnY,"_")[[1]][2]
    krnY_regexpr <- paste(krnY_regexpr,"Y",sep=".")
    aux1 <- (kernFeatY %in% names(cme$hyperParams$data$optimizable))*1
    aux2 <- (kernFeatY %in% names(cme$hyperParams$data$non_optimizable))*2
    optimNonoptim <- c("optimizable","non_optimizable")[aux1+aux2]
    krnY_pars_indx <- grep(krnY_regexpr, names(cme$hyperParams$data[[optimNonoptim]]))
    parY <- names(cme$hyperParams$data[[optimNonoptim]])[krnY_pars_indx]
    valY <- lapply(parY, function(par) cme$hyperParams$data[[optimNonoptim]][[par]]$val)
    names(valY) <- parY
    
    res <- list(krn=krnY, par=parY, val=valY)
    
    return(res)
  }, nodeTo=nodeIndx, numReg=regIndx, SIMPLIFY=FALSE)
  
  krnsY <- sapply(parsY, function(el) el$krn)
  valsY <- lapply(parsY, function(el) el$val)
  parsY <- lapply(parsY, function(el) el$par)
  
  aux1 <- (!is.null(getHyperPar(cmemLearnerSum, "kernelY")))*1
  aux2 <- (!is.null(getHyperPar(cmemLearnerSum, "featureY")))*2
  kernFeatY <- c("kernelY","featureY")[aux1+aux2]
  
  cmemLearnerSum$hyperParams$data$non_optimizable[[kernFeatY]]$val <- unique(krnsY)
  
  parms <- unique(unlist(parsY))
  for(parm in parms) cmemLearnerSum$hyperParams$data$non_optimizable[[parm]]$val <- NULL
  parms <- unlist(parsY)
  valsY <- unlist(valsY, recursive=F)
  names(valsY) <- parms
  for(parm in unique(parms)){
    # parm <- parms[1]
    # print(parm)
    indx <- which(parm==names(valsY))
    
    cmemLearnerSum$hyperParams$data$non_optimizable[[parm]]$val <- unique(unlist(valsY[indx]))
  }
  
  
  return(cmemLearnerSum)
}

getCurves_comp_deprecated <- function(cmemLearnerSum, uniqueRegsList, x, loss, reglrs, nodeIndx, regIndx, dataNm, folderSave){
  
  nodes <- names(uniqueRegsList)
  curves <- mapply(function(nodeTo, numReg){
    # nodeTo <- "2"; numReg <- "1"
    print("**********")
    print(paste(nodeTo, numReg, sep="-"))
    
    
    # step 3 - regression: get regressions for each node in dags
    # with CV on L2-loss for lambda and parsKernX
    indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
    regressors <- nodes[indxPreds]
    size <- length(regressors)
    dataX <- x[,regressors, drop=F]
    dataY <- x[,nodeTo]
    trainDataO <- constructData(x=dataX, y=dataY)
    
    regressorsChar <- paste(regressors, collapse="-")
    regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
    fileSave <- paste(dataNm, "_", size, "_",regressionChar, ".RData", sep="")
    
    if(fileSave %in% dir(folderSave)){
      load(file=paste(folderSave, fileSave, sep=""))
    } else{
      cmemLearnerSumAux <- setParams(learner=cmemLearnerSum, trainData=trainDataO)
      save(cmemLearnerAux, file=paste(folderSave, fileSave, sep=""))  
    }
    
    
    # step 4a - compare regressions: compare regressions of same size
    # (there should be at least two otherwise all dags in hypothesis
    # set have the same regression and we don't have to assess it)
    # get envelope curves for each set of regressions, compare to get 
    # measures.
    curves <- lapply(reglrs, function(reglr){
      # reglr <- reglrs[6]
      print(reglr)
      # temporary this should be assigned from measure pack
      logReglr <- switch(reglr, KCCC_ent=F, T)
      curve <- getLossReglrCurve(learner=cmemLearnerSumAux, loss=loss, reglr=reglr, logReglr=logReglr)
      return(curve)
    })
    
    names(curves) <- names(reglrs)
    return(curves)
  }, 
  nodeTo=nodeIndx, numReg=regIndx, SIMPLIFY=FALSE)
  
  names(curves) <- paste(nodeIndx, regIndx, sep=".")
  return(curves)
}
getCurves_comp <- function(cmeLearnerSumList, uniqueRegsList, loss, reglrs){
  
  nodes <- names(uniqueRegsList)
  
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims <- dims[c(2,1)]
    return(dims)
  })
  
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  curves <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      
      cmemLearnerSumAux <- cmeLearnerSumList[[nodeTo]][[numReg]]
      
      # step 4a - compare regressions: compare regressions of same size
      # (there should be at least two otherwise all dags in hypothesis
      # set have the same regression and we don't have to assess it)
      # get envelope curves for each set of regressions, compare to get 
      # measures.
      if(!is.na(cmemLearnerSumAux)){
        curves <- lapply(reglrs, function(reglr){
        # reglr <- reglrs[1]
        #print(reglr)
        # temporary this should be assigned from measure pack
        logReglr <- switch(reglr, KCCC_ent=F, KCCC_pca_ent=F,T)
        curve <- getLossReglrCurve(learner=cmemLearnerSumAux, loss=loss, reglr=reglr, logReglr=logReglr)
        return(curve)
      })
        names(curves) <- reglrs
      } else{
        curves <- NA
      }
      return(curves)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  
  
  names(curves) <- nodes 
  return(curves)
}



getMeasures_comp <- function(curves, reglrs){
  
  msrs <- sapply(reglrs, function(reglr){
    # reglr <- reglrs[2]
    
    # step 4b - compare regressions: compare regressions of same size
    # (there should be at least two otherwise all dags in hypothesis
    # set have the same regression and we don't have to assess it)
    # get envelope curves for each set of regressions, compare to get 
    # measures 
    j <- which(reglr == reglrs)
    curves_reglr <- lapply(1:length(curves), function(i) curves[[i]][[j]])
    names(curves_reglr) <- names(curves)
    compReglrVals <- getComparableReglrVal(curveList=curves_reglr)
    #plotCurvesUAI(compReglrVals, curveList=curves_reglr)
    #plotCurves(compReglrVals, curveList=curves_reglr)
    #plotCurves2(compReglrVals, curveList=curves_reglr, xlim=c(-5,0), ylim=c(0,2))
    #res <- cbind(reglr=compReglrVals$reglrs, loss=compReglrVals$losses)
    res <- cbind(reglr=compReglrVals$logReglrs, loss=compReglrVals$losses)
    rownames(res) <- names(curves)
    return(res)  
  }, simplify="array")
  names(dimnames(msrs)) <- c("curves","regLoss","msr")
  msrs <- melt(msrs)
  msrs <- cast(msrs, curves~msr+regLoss, value="value")
  rownames(msrs) <- msrs$curves
  msrs <- as.matrix(as.data.frame(msrs[,2:ncol(msrs)]))
  c.nms <- colnames(msrs)
  indxChng <- grep("loss",c.nms)
  c.nms[indxChng] <- paste("L2", sapply(strsplit(c.nms[indxChng],"_"), function(el) el[1]), sep="_")
  colnames(msrs) <- c.nms
  
  return(msrs)
}



extractLeftEnvPts <- function(x, y){
  o <- order(x)
  x2 <- x[o]
  y2 <- y[o]
  indxMin <- which.min(y2)
  #indxMin2 <- order(o)[which.min(y)]
  o2 <- o[1:indxMin]
  yEnv <- y2[o2]
  xEnv <- x2[o2]
  #plot(x2[order(o)], x)
  #plot(y2[order(o)], y)
  #plot(xEnv, x2[o[1:indxMin]])
  #plot(xEnv, x[o[o2]])
  indxKeep <- o[o2]
  
  i <- 1
  #plot(x2,y2)
  while(i < indxMin){
    #print(i)
    if(yEnv[i]<yEnv[i+1]){
      yEnv <- yEnv[-(i+1)]
      xEnv <- xEnv[-(i+1)]
      indxKeep <- indxKeep[-(i+1)]
      indxMin <- indxMin - 1
    } else{
      i <- i + 1
    }
    #lines(xEnv, yEnv, type="p",col="red", cex=2)
  }
  #plot(xEnv, x[indxKeep])
  #plot(yEnv, y[indxKeep])
  
  return(list(x=xEnv, y=yEnv, indx=indxKeep))
}
getLossReglrCurve <- function(learner, loss, reglr, logReglr=T){
  grid <- learner$hyperParams$data$grid
  dimnames(grid)$params <- 1:length(dimnames(grid)$params)
  df <- melt(grid)
  form <- as.formula(paste(paste(setdiff(colnames(df), c("var","value")), collapse="+"), "~var", sep=""))
  chk <- cast(data=df, formula=form, value="value", fun.aggregate="length")
  if(any(chk[,3:ncol(chk)]!=1)) stop("Error shd be one value per cell")
  df2 <- cast(df, form, value="value")
  df2$obs <- 1:nrow(df2)
  indx <- which(df2$trainTest=="test"  & !is.nan(df2[,reglr]) )
  x <- df2[indx, reglr]
  if(logReglr) x <- log(x, 10)
  y <- df2[indx, loss]
  indxValid <- which(!is.nan(x) & !is.na(x) & !is.nan(y) & !is.na(y) &x<Inf &x>-Inf &y<Inf &y>-Inf)
  x <- x[indxValid]; y <- y[indxValid]
  #o <- order(x)
  #x <- x[o]; y <- y[o]
  envPts <- extractLeftEnvPts(x, y)
  if(length(unique(envPts$x))==1){
    envPts$x <- c(envPts$x, envPts$x+0.00001)
    envPts$y <- c(envPts$y, envPts$y-0.00001)
  }
  #plot(x,y, ylim=c(0.2,1))
  #lines(envPts$x, envPts$y, type="p", cex=2, col="red")
  splf_yx <- splinefun(envPts$y, envPts$x, method="monoH.FC")
  splf_xy <- splinefun(envPts$x, envPts$y, method="monoH.FC")
  #yy_1 <- seq(envPts$y[1], envPts$y[length(envPts$y)], length.out=100)
  #xx_1 <- splf_yx(yy_1)
  #xx_2 <- seq(envPts$x[1], envPts$x[length(envPts$x)], length.out=100)
  #yy_2 <- splf_xy(xx_2)
  #plot(x,y, ylim=c(0.5,1))
  #lines(xx_1, yy_1, col="red", lty=2)
  #lines(xx_2, yy_2, col="blue", lty=3)
  return(list(df=df2, reglrRng=range(envPts$x), lossRng=range(envPts$y), fun_xy=splf_xy, fun_yx=splf_yx, loss=loss, reglr=reglr, logReglr=logReglr))
}
getComparableReglrVal <- function(curveList){
  # el <- curveList[[1]]
  minRngLosses <- sapply(curveList, function(el) el$lossRng[1])
  maxRngLosses <- sapply(curveList, function(el) el$lossRng[2])
  minRngReglrs <- sapply(curveList, function(el) el$reglrRng[1])
  maxRngReglrs <- sapply(curveList, function(el) el$reglrRng[2])
  refCurve <- which.min(minRngLosses)
  validCurves <- which(minRngLosses <= maxRngLosses[refCurve])
  invalidCurves <- setdiff(1:length(curveList), validCurves)
  maxMinLoss <- max(minRngLosses[validCurves])
  minMaxReglr <- min(maxRngReglrs[validCurves])
  
  
  # initial default value
  #logReglrs <- rep(curveList[[refCurve]]$reglrRng[2], length(curveList))
  #losses <- rep(curveList[[refCurve]]$lossRng[2], length(curveList))
  logReglrs <- rep(Inf, length(curveList))
  losses <- rep(Inf, length(curveList))
  logReglrs[validCurves] <- sapply(curveList[validCurves], function(el) el$fun_yx(maxMinLoss))
  losses[validCurves] <- sapply(curveList[validCurves], function(el) el$fun_xy(minMaxReglr))
  
  #logReglrs <- sapply(curveList, function(el) el$fun_yx(maxMinLoss))
  #losses <- sapply(curveList, function(el) el$fun_xy(minMaxReglr))
  
  #data.frame(logReglrsVal=logReglrs2, logReglrs=logReglrs, lossesVal=losses2, losses=losses)
  #data.frame(valid=minRngLosses <= maxRngLosses[refCurve], minReglr=minRngReglrs, maxReglr=maxRngReglrs, minLoss=minRngLosses, maxLoss=maxRngLosses)
  
  reglrs <- mapply(function(logReglr, curve){
    # i <- 1; curve <- curveList[[i]]; logReglr <- logReglrs[i]
    res <- logReglr
    if(curve$logReglr) res <- 10^logReglr
    return(res)
  }, logReglr=logReglrs, curve=curveList)
  res <- list(reglrs=reglrs, logReglrs=logReglrs, losses=losses, minMaxReglr=minMaxReglr, maxMinLoss=maxMinLoss, validCurves=validCurves, refCurve=refCurve)
  return(res)
}

plotCurvesSimple <- function(curveList){
  
  dfs <- lapply(1:length(curveList), function(i){
    # i <- 1
    loss_yy <- seq(curveList[[i]]$lossRng[1], curveList[[i]]$lossRng[2], length.out=100)
    reglr_xx <- curveList[[i]]$fun_yx(loss_yy)
    #plot(reglr_xx, loss_yy, type="l")
    df <- data.frame(reglr=reglr_xx, loss=loss_yy, curve=names(curveList)[[i]])
    return(df)
  })
  df <- do.call(rbind, dfs)
  p <- ggplot(df)
  p <- p + geom_line(aes(x=reglr, y=loss, colour=curve))
  print(p)

}

plotCurves <- function(compReglrVals, curveList){
  loss <- unique(sapply(curveList, function(el) el$loss))
  reglr <- unique(sapply(curveList, function(el) el$reglr))
  logReglr <- unique(sapply(curveList, function(el) el$logReglr))
  if(logReglr) reglr2 <- paste("log(",reglr,",10)") else reglr2 <- reglr
  dfs <- mapply(function(curve, nm) cbind(curve$df, node=nm), curve=curveList, nm=names(curveList), SIMPLIFY=FALSE)
  df <- do.call(rbind, dfs)
  funs_yx <- lapply(curveList, function(el) el$fun_yx)
  rngsLoss <- sapply(curveList, function(el) el$lossRng)
  rngLoss <- range(rngsLoss)
  indx <- which(df$trainTest=="test" & !is.na(df[,reglr])& !is.na(df[,loss]) )
  
  minReglr <- curveList[[compReglrVals$refCurve]]$reglrRng[1] #min(df[indx,reglr],na.rm=T)
  maxReglr <- max(compReglrVals$logReglrs)# curveList[[compReglrVals$refCurve]]$reglrRng[2] #max(df[indx,reglr], na.rm=T)
  minLoss <- curveList[[compReglrVals$refCurve]]$lossRng[1] 
  maxLoss <- curveList[[compReglrVals$refCurve]]$lossRng[2] 
  
  
  validCurves <- compReglrVals$validCurves
  
  dfAux <- data.frame(loss1=rep(min(df[indx,loss],na.rm=T), length(validCurves)),
                      loss2=rep(compReglrVals$maxMinLoss, length(validCurves)),
                      reglr1=compReglrVals$logReglrs[validCurves],
                      reglr2=compReglrVals$logReglrs[validCurves])
  
  
  p <- ggplot()
  p <- p + geom_polygon( aes_string(x = loss, y=reglr2 , group="obs"), colour = "grey", alpha=0.5, data=df[indx,])
  p <- p + geom_segment(aes(x = compReglrVals$maxMin, y =min(compReglrVals$reglrs) , 
                            xend =compReglrVals$maxMin , yend = max(compReglrVals$reglrs)), colour = "light yellow", alpha=0.8, size=2)
  p <- p + geom_segment(aes(x = loss1, y =reglr1 , xend =loss2 , yend = reglr2), 
                        colour = "light yellow", alpha=0.8, data=dfAux, size=2)
  
  
  p <- p + geom_point(aes_string(x=loss, y=reglr2, colour="node"), size=0.6, data=df[indx,])
  p <- p + xlim(minLoss, maxLoss) 
  p <- p + ylim(minReglr, maxReglr)
  for(i in 1:length(funs_yx)) p <- p + stat_function(fun = funs_yx[[i]], data=df[indx,],
                                                  colour="black", linetype=5  ) 
  p <- p + coord_flip()
  print(p)
  
  
}
plotCurves2 <- function(compReglrVals, curveList, xlim=NULL, ylim=NULL){
  loss <- unique(sapply(curveList, function(el) el$loss))
  reglr <- unique(sapply(curveList, function(el) el$reglr))
  logReglr <- unique(sapply(curveList, function(el) el$logReglr))
  if(logReglr) reglr2 <- paste("log(",reglr,",10)") else reglr2 <- reglr
  dfs <- mapply(function(curve, nm) cbind(curve$df, node=nm), curve=curveList, nm=names(curveList), SIMPLIFY=FALSE)
  df <- do.call(rbind, dfs)
  funs_yx <- lapply(curveList, function(el) el$fun_yx)
  rngsLoss <- sapply(curveList, function(el) el$lossRng)
  rngLoss <- range(rngsLoss)
  indx <- which(df$trainTest=="test" & !is.na(df[,reglr])& !is.na(df[,loss]) )
  
  minReglr <- curveList[[compReglrVals$refCurve]]$reglrRng[1] #min(df[indx,reglr],na.rm=T)
  maxReglr <- max(compReglrVals$logReglrs)# curveList[[compReglrVals$refCurve]]$reglrRng[2] #max(df[indx,reglr], na.rm=T)
  minLoss <- curveList[[compReglrVals$refCurve]]$lossRng[1] 
  maxLoss <- curveList[[compReglrVals$refCurve]]$lossRng[2] 
  
  
  validCurves <- compReglrVals$validCurves
  
  dfAux <- data.frame(loss1=rep(min(df[indx,loss],na.rm=T), length(validCurves)),
                      loss2=rep(compReglrVals$maxMinLoss, length(validCurves)),
                      reglr1=compReglrVals$logReglrs[validCurves],
                      reglr2=compReglrVals$logReglrs[validCurves])
  
  if(is.null(xlim)) xlim <- range(df[indx,loss], na.rm=T)
  if(is.null(ylim)) ylim <- range(log(df[indx,reglr],10), na.rm=T)
  
  p <- ggplot()
  p <- p + geom_polygon( aes_string(x = loss, y=reglr2 , group="obs"), colour = "grey", alpha=0.5, data=df[indx,])
  p <- p + geom_point(aes_string(x=loss, y=reglr2, colour="node"), size=0.6, data=df[indx,])
 # for(i in 1:length(funs_yx)) p <- p + stat_function(fun = funs_yx[[i]], data=df[indx,],
  #                                                   colour="black", linetype=5  ) 
  p <- p + xlim(xlim[1], xlim[2]) 
  p <- p + ylim(ylim[1], ylim[2])
  
  p <- p + coord_flip()
  print(p)
  
  
}

plotCurvesUAI <- function(compReglrVals, curveList, x0=NULL, xN=NULL, y0=NULL, yN=NULL){
  loss <- unique(sapply(curveList, function(el) el$loss))
  reglr <- unique(sapply(curveList, function(el) el$reglr))
  logReglr <- unique(sapply(curveList, function(el) el$logReglr))
  if(logReglr) reglr2 <- paste("log(",reglr,",10)") else reglr2 <- reglr
  loss2 <- paste("log(",loss,",10)")
  dfs <- mapply(function(curve, nm) cbind(curve$df, node=nm), curve=curveList, nm=names(curveList), SIMPLIFY=FALSE)
  df <- do.call(rbind, dfs)
  
  summary(df[,reglr])
  min(df[,reglr])
  df[which(df[,reglr]==0),reglr]
  
  funs_yx <- lapply(curveList, function(el) el$fun_yx)
  rngsLoss <- sapply(curveList, function(el) el$lossRng)
  rngLoss <- range(rngsLoss)
  indx <- which(df$trainTest=="test" & !is.na(df[,reglr])& !is.na(df[,loss]) )
  df <- df[indx,]
  minReglr <- curveList[[compReglrVals$refCurve]]$reglrRng[1] #min(df[indx,reglr],na.rm=T)
  maxReglr <- max(compReglrVals$logReglrs)# curveList[[compReglrVals$refCurve]]$reglrRng[2] #max(df[indx,reglr], na.rm=T)
  minLoss <- curveList[[compReglrVals$refCurve]]$lossRng[1] 
  maxLoss <- curveList[[compReglrVals$refCurve]]$lossRng[2] 
  
  #x0 <- xN <- y0 <- yN <- NULL
  if(!is.null(x0)) minReglr <- x0
  if(!is.null(xN)) maxReglr <- xN
  if(!is.null(y0)) minLoss <- y0
  if(!is.null(yN)) maxLoss <- yN
  
  # direction
  df$g.truth <- c("x->y","y->x")[match(df$node, c(3,5))]
  df$dir <- df$g.truth
  
  df$outOfPlot <- 1
  indx <- which(log(df[,reglr],10)<minReglr)
  if(length(indx)>0){
    df[indx,reglr] <- 10^minReglr
    df$outOfPlot[indx] <- 0
    df$g.truth[indx] <- "out"
  }
  
  indx <- which(log(df[,reglr],10)>maxReglr)
  if(length(indx)>0){
    df[indx,reglr] <- 10^maxReglr
    df$outOfPlot[indx] <- 0
    df$g.truth[indx] <- "out"
  }
  
  indx <- which(df[,loss]<minLoss)
  if(length(indx)>0){
    df[indx,loss] <- minLoss
    df$outOfPlot[indx] <- 0
    df$g.truth[indx] <- "out"
  } 
  
  indx <- which(df[,loss]>maxLoss)
  if(length(indx)>0){
    df[indx,loss] <- maxLoss
    df$outOfPlot[indx] <- 0
    df$g.truth[indx] <- "out"
  }
  
  
  
  df$logKCMC <- log(df$KCMC, 10)
  
  # obtain gamma value
  indxStart <- regexpr("gamma=",df$params)+6
  indxFinish <- regexpr("sigma.rbf.X=",df$params)-2
  df$gamma <- as.numeric(substr(df$params, indxStart, indxFinish))
  
  # see if corresponding (sigma_x, lambda, gamma) KCMCs are "right"
  # to add colour to segments
  
  cast(df, obs~dir, value="KCMC", fun.aggregate = length)
  indxObsKCMC <- cast(df, obs~dir, value="KCMC")
  indxObsKCMC$dir <- c(2,1)[(indxObsKCMC[,2]<indxObsKCMC[,3])*1+1]
  
  df$pred <- as.factor(c("x->y","x<-y")[indxObsKCMC$dir[match(df$obs, indxObsKCMC$obs)]])
  
  validCurves <- compReglrVals$validCurves
  
  # min(df[,loss],na.rm=T)
  dfAux <- data.frame(loss1=rep(minLoss, length(validCurves)),
                      loss2=rep(compReglrVals$maxMinLoss, length(validCurves)),
                      reglr1=compReglrVals$logReglrs[validCurves],
                      reglr2=compReglrVals$logReglrs[validCurves])
  
  softCol1 <- "#c6dbef"
  softCol2 <- "#c7e9c0"    
  col1 <- "#31a354"
  col2 <- "#3182bd"
  colOutOfPlot <- "light grey"
  
  p <- ggplot()
  p <- p + geom_polygon( aes_string(x = loss, y=reglr2 , group="obs", colour="pred"), show.legend = F ,alpha=0.3, data=df,  fill="black")
  p <- p + theme(legend.title=element_blank())
  #p <- p  + scale_colour(guide = 'none')
  p <- p +  theme(legend.position="none") 
  p <- p + guides(fill=FALSE,size = FALSE, shape=FALSE)
  p <- p + geom_segment(aes(x = compReglrVals$maxMin, y =min(dfAux$reglr1) , 
                            xend =compReglrVals$maxMin , yend = max(dfAux$reglr1)), 
                        colour = "light yellow", alpha=0.8, size=2)
  p <- p + geom_segment(aes(x = loss1, y =reglr1 , xend =loss2 , yend = reglr2), 
                        colour = "light yellow", alpha=0.8,data=dfAux, size=2)
  p <- p + geom_segment(aes(x = compReglrVals$maxMin, y =min(dfAux$reglr1) , 
                            xend =compReglrVals$maxMin , yend = max(dfAux$reglr1)), 
                        colour = "black", alpha=0.8, size=0.5, linetype=5)
  p <- p + geom_segment(aes(x = loss1, y =reglr1 , xend =loss2 , yend = reglr2), 
                        colour = "black", alpha=0.8,data=dfAux, size=0.5, linetype=5)
  p <- p + scale_color_manual(values=c(softCol2,softCol1))
  p <- p + new_scale_color()
  p <- p + scale_color_manual(values=c(colOutOfPlot,col1,col2))
  p <- p + geom_point(aes(x=cmem_L2_f, y=logKCMC,  colour=g.truth, alpha=outOfPlot), show.legend=T, data=df, size=0.7)
  #p <- p + scale_color_hue(l=40, c=35)
  #p <- p + scale_color_manual(values=c(col1, col2)) #, #56B4E9, #999999
  #p <- p + theme(legend.title="ground truth")
  #p <- p + new_scale_color()
  #p <- p + geom_point(aes(x=cmem_L2_f, y=logKCMC, colour=log(gamma,10)), size=4, data=df, alpha=0.1)
  p <- p + xlim(minLoss, maxLoss) 
  p <- p + geom_blank()
  
  p <- p + ylim(minReglr, maxReglr)
  for(i in 1:length(funs_yx)) p <- p + stat_function(fun = funs_yx[[i]], data=df,
                                                     colour=c(col2, col1)[i], linetype=5  ,n=10000) 
  p <- p + coord_flip()
  p <- p + xlab("RMSE")
  print(p)
  
  
}
