# Additive Noise model hypothesis scorer functions

########################################################################################*
# Hypothesis scorer functions - a hypothesis score consists of a function that
# takes in data and parameters and gives back a score
#
########################################################################################*

# Additive noise model hypothesis scorer
# fit y = f(x) + e to dags, obtain residuals, evaluate complexity
#
# takes:
# normal hypothesis scorer params data -x, and approx hypothesis set -hypArray, plus
# 1. data regiment (split/recycle) - dataReg
# 2. f(x) estimation method - learner

print("loading learner functions")
source("./pkg_learner/func_learners_v5.R")


# 3. complexity function(s) - complexityPack

print("loading complexity functions")
source("./pkg_cmplxtyScores/func_complexity_v1.R")



anm_hypScorer <- function(x, hypArray, dataReg, learner, complexityPack, ppTab=NULL, plot=FALSE, dataNm, folderSave, hypSc_char){
  
  cmplxScorePack <- eval(parse(text=complexityPack))
  
  # apply learner-data setting to each hypothesis to get residuals
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  
  # get data
  dataRg <- dataRegime(x, type=dataReg)
  xTrain <- dataRg$train
  xTest <- dataRg$test
  # get learner
  file <- paste(dataNm, "_",hypSc_char,".RData", sep="")
  learner <- eval(parse(text=learner))  
  
  
  # learn
  #print(paste("fitting regressions needed to evaluate all hypotheses"))
  if(! file %in% dir(folderSave)){
    print("learns")
    #trainData=xTrain
    modsClass <- fitSEMSetGivenDAGSet(uniqueRegsList, trainData=xTrain, learner=learner) 
    save(list=c("modsClass"), file=paste(folderSave, file, sep=""))
  } else{
    print("reads")
    load(paste(folderSave, file, sep=""))
    #mcqrnn.predict(as.matrix(rnorm(10)), modsClass$x[[1]]$learnParams$model)
    #mcqrnn.predict(as.matrix(rnorm(10)), modsClass$y[[2]]$learnParams$model)
  }
  # predict
  #print(paste("obtaining residuals needed to evaluate all hypotheses"))
  
  # data=xTest; learnerList=modsClass
  predsClass <- predictSEMSet(uniqueRegsList, data=xTest, learnerList=modsClass, plot=plot) #data=xTest
  
  #print(paste("scoring all hypotheses"))
  #predsSet <- lapply(predsClass, function(el) adrop(el[,"resid", ,drop=FALSE], drop=2))
  # dags=hypArray; vars=predsClass; cmplxScores=cmplxScorePack
  scores <- complexityScoreList(dags=hypArray, uniqueRegsList=uniqueRegsList, vars=predsClass, cmplxScores=cmplxScorePack, prnt=FALSE)
  dmnms <- dimnames(scores)
  timeLearns <- timeScoreList(dags=hypArray, uniqueRegsList, learnersList=modsClass)
  
  scores <- cbind(scores, timeLearn=timeLearns[match(rownames(scores), names(timeLearns))])
  dmnms$score <- c(dmnms$score, "timeLearn")
  dimnames(scores) <- dmnms
  
  return(scores)
}

anm_hypScorer_boot <- function(x, hypArray, dataReg, learner, complexityPack, numBoots, numPerBoot, ppTab=NULL, plot=FALSE, dataNm, folderSave){
  print("enters anm_hypScorer_boot")
  print(paste("learner:",learner))
  cmplxScorePack <- eval(parse(text=complexityPack))
  
  # apply learner-data setting to each hypothesis to get residuals
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  
  set.seed(12)
  xs <- lapply(1:numBoots, function(i){
    smpl <- sample(1:nrow(x), numPerBoot)
    res <- x[smpl,]
    
    # get data
    dataRg <- dataRegime(res, type=dataReg)
    xTrain <- dataRg$train
    xTest <- dataRg$test
    res <- list(xTrain=xTrain, xTest=xTest)
    return(res)
  })
  
  xTrains <- lapply(xs, function(el) el$xTrain)
  xTests <- lapply(xs, function(el) el$xTest)
  
  # get learner
  learner <- eval(parse(text=learner))  
  
  counts <- 1:numBoots
  
  scores <- mcmapply(FUN=function(xTrain, xTest, cnt){
    # i <- 1; xTrain <- xTrains[[i]]; xTest <- xTests[[i]]; cnt <- counts[i]
    print("-------------------------------------")
    print(paste("boot number: ", cnt, " of ", numBoots))
    
    # learn
    #print(paste("fitting regressions needed to evaluate all hypotheses"))
    print("model")
    modsClass <- fitSEMSetGivenDAGSet(uniqueRegsList, trainData=xTrain, learner=learner) #trainData=xTrain
    # predict
    #print(paste("obtaining residuals needed to evaluate all hypotheses"))
    print("predictions")
    predsClass <- predictSEMSet(uniqueRegsList, data=xTest, learnerList=modsClass, dataReg, plot=plot) #data=xTest
  
    #print(paste("scoring all hypotheses"))
    #predsSet <- lapply(predsClass, function(el) adrop(el[,"resid", ,drop=FALSE], drop=2))
    print("scores")
    scores <- complexityScoreList(dags=hypArray, uniqueRegsList=uniqueRegsList, vars=predsClass, cmplxScores=cmplxScorePack, prnt=FALSE)
  
    return(scores)
  }, xTrain=xTrains, xTest=xTests, cnt=counts, mc.cores=1, SIMPLIFY="array")
  
  dimnames(scores)[[3]] <- 1:numBoots
  names(dimnames(scores))[3] <- "boot"
  scores2 <- melt(scores)
  scores3 <- cast(scores2, dag~score+boot, value="value")
  scores3 <- scores3[match(scores3$dag, dimnames(scores)$dag),]
  scores4 <- as.matrix(scores3[,2:ncol(scores3)])
  dimnames(scores4) <- list(dag=scores3$dag, score=colnames(scores3)[2:ncol(scores3)])
  
  print("exits anm_hypScorer_boot")
  return(scores4)
}


#################*
# data regiment
#################*
dataRegime <- function(data, type=c("recycle", "holdout")){
  dataReg <- switch(type, 
                    recycle={
                      xTrain <- data
                      xTest <- data
                      return(list(train=xTrain, test=xTest))
                    },
                    holdout={
                      midpoint <- ceiling(nrow(data)/2)
                      xTrain <- data[1:midpoint,]
                      xTest <- data[(midpoint+1):nrow(data),]
                      return(list(train=xTrain, test=xTest))
                    })
  return(dataReg)
}


###############################################*
# fit learner to regressions 
# implied by a set of graphs
#
# sets of dags share certain regressions so its 
# best to apply this method directly to groups 
# of dags instead of to 1 dag to avoid fitting 
# the same regression more than once
################################################*


# obtains all possible unique regressions for each node (ie for all possible m-node dags)
getUniqueRegs <- function(allDags){
  # slice <- "w"
  uniqueRegs <- sapply(dimnames(allDags)[[2]], function(slice) unique(allDags[,slice,], MARGIN=2), simplify="array")
  dimnames(uniqueRegs) <- list(nodeFrom=dimnames(uniqueRegs)[[1]], numReg=1:(dim(uniqueRegs)[2]), nodeTo=dimnames(uniqueRegs)[[3]])
  return(uniqueRegs)
}

# obtains the unique regressions for a certain set of candidate dags
getUniqueRegsList <- function(dags){
  # slice <- "w"
  uniqueRegsList <- lapply(dimnames(dags)[[2]], function(slice){ 
    # slice <- "1"
    res <- unique( adrop(dags[,slice,,drop=FALSE], drop=2), MARGIN=2)
    dimnames(res) <- list(nodeFrom=dimnames(res)[[1]], numReg=1:(dim(res)[2]))
    return(res)
  })
  names(uniqueRegsList) <-c(dimnames(dags)[[2]])
  return(uniqueRegsList)
}

# obtains the models for unique regressions implied in a set of dags (output of getUniqueRegsList)
fitSEMSetGivenDAGSet <- function(uniqueRegsList, trainData, learner, learnerMarg="qr_marg"){
  
  
  n <- dim(trainData)[1]
  dims <- lapply(uniqueRegsList, function(el){
    dims <- dimnames(el)[2]
    dims$sim <- 1:n
    dims <- dims[c(2,1,3)]
    return(dims)
  })
  
  nodes <- names(uniqueRegsList)
  numRegs <- lapply(dims, function(el) el$numReg)
  numRegsTot <- sum(sapply(uniqueRegsList, function(el) sum(apply(el,2, function(col) any(col!=0)*1))))
  
  mc_cores <- 1
  numCoresFold <- detectCores()/2
  #mc_cores <- numCoresFold <- 1
  
  #count <- 0
  #pm0 <- proc.time()
  learnerList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      # numReg <- "2"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>=0){
        #count <- count + 1
        if(length(indxPreds)==0) lrnr <- eval(parse(text=learnerMarg)) else lrnr <- learner
        trainDataO <- constructData(x=as.matrix(trainData[,nodes[indxPreds]]), y=trainData[,nodeTo])
        pm <- proc.time()
        learnerAux <- setParams(lrnr, trainData=trainDataO, mc_cores=mc_cores, numCoresFold=numCoresFold)
        learnerAux <- learnerAux$learn(learner=learnerAux)
        pm <- (proc.time()-pm)[3]
        learnerAux$timeLearn <- pm
        print("test learner")
        print(try(head(mcqrnn.predict(as.matrix(rnorm(10)), learnerAux$learnParams$model))))
        
        #print("estimated time to completion:")
        #numRegsLeft <- numRegsTot - count
        #numRegsDone <- numRegsTot -numRegsLeft
        #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        learnerAux <- list(NULL)
      }
      return(learnerAux)
    })
    names(res) <- numRegs[[nodeTo]]
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(learnerList) <- nodes
  
  return(learnerList)
}


# uses models to obtain predictions and residuals for set of candidate dags
predictSEMSet <- function(uniqueRegsList, data, learnerList, dataReg, plot=TRUE){
  
  
  n <- dim(data)[1]
  dims <- lapply(uniqueRegsList, function(el){
    # el <- uniqueRegsList[[1]]
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
  predResids <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    # nodeTo <- nodes[2]
    #print("*********************")
    #print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    
    res <- sapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      # numReg <- "2"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>=0){
        #count <- count + 1
        
        learnerAux <- learnerList[[nodeTo]][[numReg]]
        
        #if(dataReg=="recycle"){
        #  dataO <- learnerAux$hyperParams$trainData
        #} else{
          dataO <- constructData(x=as.matrix(data[,nodes[indxPreds]]), y=data[,nodeTo])  
        #}
        
        
        pred <- learnerAux$predict(learner=learnerAux, data=dataO)
        
        # if(dataReg=="recycle"){
        #   pred <- learnerAux$predict(learnerAux, data=dataO)
        # } 
        # if(dataReg=="holdout"){
        #   pred <- predict.CV(learnerAux, data=dataO)
        #   pred <- pred$test
        # }
        
        if(plot){
          predList <- list(pred)
          names(predList) <- nodeTo
          plot.emeley.1D(predList=predList)
        } 
        #preds  <-  pred$gyh
        #resids <- pred$gyh - pred$gy
        res <- learnerAux$resids(learner=learnerAux, pred)
        #hist(res[,"resid"])
        #plot(dataO$x, res[,"resid"])
        #dhsic.test(dataO$x, res[,"resid"])
        #qs <- quantile(dataO$x, c(0,learnerAux$hyperParams$data$non_optimizable$taus$val[[1]],1))
        #grp <- findInterval(dataO$x, qs, rightmost.closed = T)
        #table(grp)
        #ksTests <- aggregate(res[,"resid"], by=list(grp), FUN=function(x) ks.test(x, y="punif")$p.value)
        # hist(res[which(grp==3),"resid"])
        #hiscTests <- sapply(ksTests$Group.1, function(gr){
          ## gr <- ksTests$Group.1[1]
        #  indx <- which(grp==gr)
        #  plot(dataO$x[indx,],res[indx,"resid"], main=gr)
        #  dhsic.test(res[indx,"resid"], dataO$x[indx,])$p.value
        #})
        
        #grpy <- sapply(1:nrow(pred$gyh), function(i){ findInterval(dataO$y[i], sort(pred$gyh[i,]), rightmost.closed=T)})
        #table(grpy)
        
        #hiscTestsY <- sapply(1:6, function(gr){
        #  # gr <- ksTests$Group.1[1]
        #  indx <- which(res[,"grp"]==gr)
        #  plot(dataO$x[indx,],res[indx,"resid"], main=gr)
        #  dhsic.test(res[indx,"resid"], dataO$x[indx,])$p.value
        #})
        
        #mean(ksTests$x)
        #mean(hiscTests)
        #mean(hiscTestsY)
        #hist(res[which(grp==1),"resid"])
        #hist(res[which(grp==2),"resid"])
        #hist(res[which(grp==3),"resid"])
        #hist(res[which(grp==4),"resid"])
        
        #res <- cbind(pred=preds, resid=resids)
      
  
        #print("estimated time to completion:")
        #numRegsLeft <- numRegsTot - count
        #numRegsDone <- numRegsTot -numRegsLeft
        #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        
        
        
        res <- cbind(pred=rep(0,n), resid=data[,nodeTo], grp=rep(1,n))
      }
      return(res)
    }, simplify="array")
    
    
    dimnames(res) <- list(numObs=1:n, predResid=c("pred","resid","grp"), numReg=numRegs[[nodeTo]])
    return(res)
  })
  #print(proc.time()-pm0) #  
  
  names(predResids) <- nodes
  
  return(predResids)
}


##############################################################################################*
# Functions for applying corresponding unique regressions implied by a class/set of dags
# arranging residuals per dag and then evaluating complexity
##############################################################################################*

# NOTA!: CREO QUE VALE LA PENA PONER CMPLX FUNCS, Y PARS EN UNA LISTA Y SI SON VARIAS EN UNA LISTA DE LISTAS

complexityScore <- function(allDags, uniqueRegs, vars, cmplxScores){
  
  #indxDag <- 543
  #allDags[,,indxDag]
  
  matlabFunctions <- c("Shannon_vME")
  # Functions which are in matlab but which I have now programmed in R
  # "Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", 
  # "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", "Shannon_spacing_LL", "Shannon_PSD_SzegoT", 
  #"Shannon_Edgeworth", "Shannon_MaxEnt1", "Shannon_MaxEnt2", "Shannon_expF", "Shannon_KDP"
  
  matlab <- rep(FALSE, length(cmplxScores))
  indxEntropy <- which(sapply(cmplxScores, function(el) el$func) == "score_sumMarginalEntropies")
  
  type <- matlab
  
  entropyEstimators <- sapply( cmplxScores[indxEntropy], function(el) el$pars$type)
  
  if(length(indxEntropy)>0) type[indxEntropy] <- entropyEstimators
  
  indxMatlab <- which(entropyEstimators %in% matlabFunctions)
  matlab[indxEntropy][indxMatlab] <- TRUE
  
  data.frame(sapply(cmplxScores, function(el) el$func), matlab, type)
  
  if(any(matlab)){
    Matlab$startServer(port=9998)
    ## Connect to Matlab session
    matlabSession <- Matlab(port=9998)
    setVerbose(matlabSession, threshold=200000)
    open(matlabSession)
    # Load ITE library into Matlab session
    evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/matlab'))")
    #evaluate(matlab, "help HShannon_kNN_k_initialization")
    
  } else{
    matlabSession <- NULL
  }
  
  cmplxFuncs <- sapply(cmplxScores, function(el) el$func)
  pars <- lapply(cmplxScores, function(el) el$pars)
  
  nodes <- dimnames(allDags)[[1]]
  numDags <- dim(allDags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(allDags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- allDags[,,1]
    count <<- count + 1
    #print("*******************")
    #print(paste("dag # ",count, sep=""))
    
    # obtain residuals/vars to evaluate for each dag
    varsDag <- sapply(nodes, function(nodeTo){
      children <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegs[,,nodeTo],2, function(col) all(col==children)))
      rsNode <- vars[,indxReg, nodeTo]
      return(rsNode)
    }, simplify="array")
    
    
    #score <- dhsic.test(residsDag, matrix.input=TRUE, pairwise=FALSE)[[type]]
    
    #PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
    # ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
    # ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
    # COMBINAR COMO CORRESPONDE
    
    
    
    score <- mapply(FUN=function(func, p, m, matlabSession){
      # i<-4; func <- cmplxFuncs[i]; p <- pars[[i]]; m <- matlab[i]
      #print(paste("function: ", func))
      #print("parameters:")
      #print(p)
      if(m){
        ps <- c(p, list(vars=varsDag), matlabSession=list(matlabSession))
      } else{
        ps <- c(p, list(vars=varsDag))
      }
      do.call(func, ps)
    }, func=cmplxFuncs, p=pars, m=matlab, MoreArgs=list(matlabSession=matlabSession)) 
    
    
    if(length(indxEntropy)>0) names(score)[indxEntropy] <- entropyEstimators
    
    #print("estimated time until completion")
    #print((proc.time()-pm0)*(numDags-count)/count)
    return(score)
  })
  #print(proc.time()-pm0) #
  
  if(any(matlab)){
    close(matlabSession)
  }
  
  
  dims <- dimnames(scores)
  dims <- list(score=dims[[1]], dag=dims$dag)
  dimnames(scores) <- dims
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

complexityScoreList <- function(dags, uniqueRegsList, vars, cmplxScores, prnt=FALSE){
  
  ids <- sapply(cmplxScores, function(el) el$id)
  cmplxFuncs <- sapply(cmplxScores, function(el) el$func)
  pars <- lapply(cmplxScores, function(el) el$pars)
  args <- sapply(cmplxScores, function(el) el$arg)
  
   
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,2]
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    varsDag <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      parents <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      #rsNode <- vars[[nodeTo]][,indxReg]
      rsNode <- vars[[nodeTo]][,,indxReg]
      return(rsNode)
    }, simplify="array")
    
    #PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
    # ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
    # ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
    # COMBINAR COMO CORRESPONDE
    
    score <- mapply(FUN=function(func, p, a){
      # i<-4; func <- cmplxFuncs[i]; p <- pars[[i]]; a <- args[[i]]
      if(prnt){
        print(paste("function: ", func))
        print("parameters:")
        print(p)
      }
 
      ps <- c(p, list(vars=varsDag[,a, ]))
      
      # score_entropy(vars=varsDag[,a,], type=ps$type)
      # score_HSIC_fix(vars=ps$vars, method=ps$method)
      # score_pdHSICgMin(vars=ps$vars, method=ps$method)
      res <- do.call(func, ps)
      return(res)
    }, func=cmplxFuncs, p=pars, a=args) 
    
    names(score) <- ids
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  })
  if(prnt) print(proc.time()-pm0) #
  
 
  
  
  dims <- dimnames(scores)
  dims <- list(score=dims[[1]], dag=dims$dag)
  dimnames(scores) <- dims
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

timeScoreList <- function(dags, uniqueRegsList, learnersList, prnt=FALSE){
  
  
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,1]
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    time <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      parents <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==parents)))
      res <- learnersList[[nodeTo]][[indxReg]]$timeLearn
      return(res)
    }, simplify="array")
    time <- sum(unlist(time))
    
    return(time)
  })
  
  return(scores)
}
