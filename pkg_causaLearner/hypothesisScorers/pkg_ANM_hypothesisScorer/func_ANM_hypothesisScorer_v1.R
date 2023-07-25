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
source("/func_learners_v2.R")


# 3. complexity function(s) - complexityPack

print("loading complexity functions")
source("./hypothesisScorers/pkg_ANM_hypothesisScorer/func_complexity_v1.R")



anm_hypScorer <- function(x, hypArray, dataReg, learner, complexityPack, ppTab=NULL, plot=FALSE){
  
  cmplxScorePack <- eval(parse(text=complexityPack))
  
  # apply learner-data setting to each hypothesis to get residuals
  uniqueRegsList <- getUniqueRegsList(dags=hypArray)
  
  # get data
  dataRg <- dataRegime(x, type=dataReg)
  xTrain <- dataRg$train
  xTest <- dataRg$test 
  # get learner
  learner <- eval(parse(text=learner))  
  
  # learn
  print(paste("fitting regressions needed to evaluate all hypotheses"))
  modsClass <- fitSEMSetGivenDAGSet(uniqueRegsList, trainData=xTrain, learner=learner)
  # predict
  print(paste("obtaining residuals needed to evaluate all hypotheses"))
  predsClass <- predictSEMSet(uniqueRegsList, data=xTest, learnerList=modsClass, plot=plot)
  
  print(paste("scoring all hypotheses"))
  predsSet <- lapply(predsClass, function(el) adrop(el[,"resid", ,drop=FALSE], drop=2))
  scores <- complexityScoreList(dags=hypArray, uniqueRegsList=uniqueRegsList, vars=predsSet, cmplxScores=cmplxScorePack, prnt=FALSE)
  
  
  
  return(scores)
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

# given train data and a graph, estimate corresponding non-linear additive SEM returning a  named list of models
fitSEMgivenDAG <- function(G, trainData, learner){
  
  
  if(class(G)=="graphNEL"){
    G <- amat(as.bn(G))
    Nodes <- colnames(G)
  } else if(class(G)=="matrix"){
    Nodes <- colnames(G) 
  } else{
    stop("G must be a graph or a matrix")
  }
  
  
  if(all(colnames(trainData)!=Nodes)) stop("names (G, trainData) dont match")
  
  
  
  topoNodes <- topoSort(G)
  p <- length(Nodes)
  nTrain <- nrow(trainData)
  
  learnerList <- list()
  learnerList[Nodes] <- list(NULL)
  
  
  indxNoParents <- which(apply(G,2, function(col) all(col==0)))
  
  
  print(paste("no parent nodes: ", Nodes[indxNoParents], sep=""))
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # nd <- "x"
    print("*********************************")
    print(paste("node: ", nd,sep=""))
    indxParents <- which(G[,nd]==1)
    Parents <- Nodes[indxParents]
    
    trainDataO <- constructData(as.matrix(trainData[,Parents]), trainData[,nd])
    
    
    learnerAux <- setParams(learner, trainDataO)
    learnerAux <- learnerAux$learn(learnerAux)
    
    learnerList[[nd]] <- learnerAux
    
  }
  
  return(learnerList)
  
}

# given data a graph and a list of models for each node in the graph, predict node variables
# obtain residuals and return predictions and residuals

predictSEM <- function(G, data, learnerList, plot=TRUE){
  
  
  if(class(G)=="graphNEL"){
    G <- amat(as.bn(G))
    Nodes <- colnames(G)
  } else if(class(G)=="matrix"){
    Nodes <- colnames(G) 
  } else{
    stop("G must be a graph or a matrix")
  }
  
  
  if(all(colnames(data)!=Nodes)) stop("names (G, data) dont match")
  if(all(names(learnerList)!=Nodes)) stop("names (learnerList, G) dont match")
  
  topoNodes <- topoSort(G)
  p <- length(Nodes)
  n <- nrow(data)
  
  
  
  indxNoParents <- which(apply(G,2, function(col) all(col==0)))
  
  
  preds <- matrix(NA, nrow(data), ncol(data))
  resids <- matrix(NA, nrow(data), ncol(data))
  colnames(preds) <- colnames(data)
  colnames(resids) <- colnames(data)
  
  
  preds[,indxNoParents] <- 0
  resids[,indxNoParents] <- data[,indxNoParents]
  
  print(paste("no parent nodes: ", Nodes[indxNoParents], sep=""))
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # nd <- "x"
    print("*********************************")
    print(paste("node: ", nd,sep=""))
    indxParents <- which(G[,nd]==1)
    Parents <- Nodes[indxParents]
    
    learner <- learnerList[[nd]]
    
    dataOb <- constructData(as.matrix(data[,Parents]), data[,nd])
    
    pred <- learner$predict(learner, data=dataOb)
    
    if(plot){
      predList <- list(pred)
      names(predList) <- nd  
      plot.emeley.1D(predList=predList)
    } 
    
    preds[,nd] <- pred$gyh
    resids[,nd] <- pred$gyh - pred$gy
    
  }
  
  
  return(list(preds=preds, resids=resids))
  
}


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


# obtains the models for all possible unique regressions implied in all m-node dags (output of getUniqueRegs)
fitSEMClassGivenDAGClass <- function(uniqueRegs, trainData, learner){
  
  
  nodes <- dimnames(uniqueRegs)[[1]]
  numRegs <- dimnames(uniqueRegs)[[2]]
  numRegsTot <- (length(numRegs)-1)*length(nodes)
  
  pm0 <- proc.time()
  learnerList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- nodes[1]
    print("*********************")
    print(paste(nodeTo, " regressions:", sep=""))
    res <- lapply(numRegs, function(numReg){
      # numReg <- numRegs[2]
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegs[,numReg,nodeTo]==1)
      if(length(indxPreds)>0){
        
        trainDataO <- constructData(as.matrix(trainData[,nodes[indxPreds]]), trainData[,nodeTo])
        learnerAux <- setParams(learner, trainDataO)
        learnerAux <- learnerAux$learn(learnerAux)
        
        print("estimated time to completion:")
        numRegsLeft <- ((length(nodes)-match(nodeTo, nodes))*(length(numRegs)-1)+(length(numRegs)-match(numReg,numRegs)))
        numRegsDone <- numRegsTot -numRegsLeft
        print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        learnerAux <- list(NULL)
      }
      return(learnerAux)
    })
    names(res) <- numRegs
    return(res)
  })
  print(proc.time()-pm0) # 5.2 mins for krr1, 250 train points 
  
  names(learnerList) <- nodes
  
  
  return(learnerList)
}

# uses models to obtain predictions and residuals for all possible unique regression implied in all m-node dags
predictSEMClass <- function(uniqueRegs, data, learnerList, plot=TRUE){
  
  dims <- dimnames(uniqueRegs)[-1]
  n <- dim(data)[1]
  dims$sim <- 1:n
  dims$predResid <- c("pred","resid")
  dims <- dims[c(3,4,1,2)]
  
  nodes <- dimnames(uniqueRegs)[[1]]
  numRegs <- dimnames(uniqueRegs)[[2]]
  numRegsTot <- (length(numRegs)-1)*length(nodes)
  
  pm0 <- proc.time()
  predResid <- sapply(nodes,  function(nodeTo){
    print("*********************")
    print(paste(nodeTo, " regressions:", sep=""))
    sapply(numRegs, function(numReg){
      # nodeTo <- "w"; numReg <- "2"
      print("******")
      print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegs[,numReg,nodeTo]==1)
      if(length(indxPreds)>0){
        
        dataOb <- constructData(as.matrix(data[,nodes[indxPreds]]), data[,nodeTo])
        learnerAux <- learnerList[[nodeTo]][[numReg]]
        pred <- learnerAux$predict(learnerAux, data=dataOb)
        if(plot){
          predList <- list(pred)
          names(predList) <- nodeTo
          plot.emeley.1D(predList=predList)
        } 
        preds  <- pred$gyh
        resids <- pred$gyh - pred$gy
        res <- cbind(preds, resids)
        
        
        print("estimated time to completion:")
        numRegsLeft <- ((length(nodes)-match(nodeTo, nodes))*(length(numRegs)-1)+(length(numRegs)-match(numReg,numRegs)))
        numRegsDone <- numRegsTot -numRegsLeft
        print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        res <- cbind(rep(0,n), data[,nodeTo])
      }
      return(res)
    }, simplify="array")
  }, simplify="array")
  print(proc.time()-pm0) # 30 mins for krr, 
  
  dimnames(predResid) <- dims 
  
  return(predResid)
}

# obtains the models for unique regressions implied in a set of dags (output of getUniqueRegsList)
fitSEMSetGivenDAGSet <- function(uniqueRegsList, trainData, learner){
  
  
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
  
  
  #count <- 0
  #pm0 <- proc.time()
  learnerList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "2"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>0){
        #count <- count + 1
        
        trainDataO <- constructData(as.matrix(trainData[,nodes[indxPreds]]), trainData[,nodeTo])
        learnerAux <- setParams(learner, trainData=trainDataO)
        learnerAux <- learnerAux$learn(learnerAux)
        
        
        
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
predictSEMSet <- function(uniqueRegsList, data, learnerList, plot=TRUE){
  
  
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
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    
    res <- sapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "2"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      if(length(indxPreds)>0){
        #count <- count + 1
        
        dataO <- constructData(as.matrix(data[,nodes[indxPreds]]), data[,nodeTo])
        learnerAux <- learnerList[[nodeTo]][[numReg]]
        pred <- learnerAux$predict(learnerAux, data=dataO)
        if(plot){
          predList <- list(pred)
          names(predList) <- nodeTo
          plot.emeley.1D(predList=predList)
        } 
        preds  <- pred$gyh
        resids <- pred$gyh - pred$gy
        res <- cbind(pred=preds, resid=resids)
        
        #print("estimated time to completion:")
        #numRegsLeft <- numRegsTot - count
        #numRegsDone <- numRegsTot -numRegsLeft
        #print(numRegsLeft*(proc.time()-pm0)/numRegsDone)
      } else{
        res <- cbind(pred=rep(0,n), resid=data[,nodeTo])
      }
      return(res)
    }, simplify="array")
    
    
    dimnames(res) <- list(numObs=1:n, predResid=c("pred","resid"), numReg=numRegs[[nodeTo]])
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
    print("*******************")
    print(paste("dag # ",count, sep=""))
    
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
      print(paste("function: ", func))
      print("parameters:")
      print(p)
      if(m){
        ps <- c(p, list(vars=varsDag), matlabSession=list(matlabSession))
      } else{
        ps <- c(p, list(vars=varsDag))
      }
      do.call(func, ps)
    }, func=cmplxFuncs, p=pars, m=matlab, MoreArgs=list(matlabSession=matlabSession)) 
    
    
    if(length(indxEntropy)>0) names(score)[indxEntropy] <- entropyEstimators
    
    print("estimated time until completion")
    print((proc.time()-pm0)*(numDags-count)/count)
    return(score)
  })
  print(proc.time()-pm0) #
  
  if(any(matlab)){
    close(matlabSession)
  }
  
  
  dims <- dimnames(scores)
  dims <- list(score=dims[[1]], dag=dims$dag)
  dimnames(scores) <- dims
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

complexityScoreListMatlabDeprecated <- function(dags, uniqueRegsList, vars, cmplxScores, prnt=FALSE){
  
  
  matlabFunctions <- c("Shannon_vME")
  # Functions which are in matlab but which I have now programmed in R
  # "Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", 
  # "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", "Shannon_spacing_LL", "Shannon_PSD_SzegoT", 
  #"Shannon_Edgeworth", "Shannon_MaxEnt1", "Shannon_MaxEnt2", "Shannon_expF","Shannon_KDP"
  
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
  
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,31]
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    varsDag <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      children <- dag[,nodeTo]
      indxReg <- which(apply(uniqueRegsList[[nodeTo]], 2, function(col) all(col==children)))
      rsNode <- vars[[nodeTo]][,indxReg]
      return(rsNode)
    }, simplify="array")
    
    #PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
    # ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
    # ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
    # COMBINAR COMO CORRESPONDE
    
    score <- mapply(FUN=function(func, p, m, matlabSession){
      # i<-7; func <- cmplxFuncs[i]; p <- pars[[i]]; m <- matlab[i]
      if(prnt){
        print(paste("function: ", func))
        print("parameters:")
        print(p)
      }
      if(m){
        ps <- c(p, list(vars=varsDag), matlabSession=list(matlabSession))
      } else{
        ps <- c(p, list(vars=varsDag))
      }
      return(do.call(func, ps))
    }, func=cmplxFuncs, p=pars, m=matlab, MoreArgs=list(matlabSession=matlabSession)) 
    
    if(length(indxEntropy)>0) names(score)[indxEntropy] <- entropyEstimators
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  })
  if(prnt) print(proc.time()-pm0) #
  
  if(any(matlab)){
    close(matlabSession)
  }
  
  
  dims <- dimnames(scores)
  dims <- list(cmplxFunc=dims[[1]], dag=dims$dag)
  dimnames(scores) <- dims
  scores <- aperm(scores, c(2,1))
  
  return(scores)
}

complexityScoreList <- function(dags, uniqueRegsList, vars, cmplxScores, prnt=FALSE){
  
  ids <- sapply(cmplxScores, function(el) el$id)
  cmplxFuncs <- sapply(cmplxScores, function(el) el$func)
  pars <- lapply(cmplxScores, function(el) el$pars)
  
  
   
  nodes <- dimnames(dags)[[1]]
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  scores <- apply(dags, "dag", function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,31]
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
      rsNode <- vars[[nodeTo]][,indxReg]
      return(rsNode)
    }, simplify="array")
    
    #PARA LOS SCORES DE ENTROPIA ESTA FORMA NO ES EFICIENTE YA QUE SE PUEDE OBTENER LA 
    # ENTROPIA DE CADA VARIABLE/RESIDUAL DE FORMA INDEPENDIENTE Y LUEGO COMBINAR LAS 
    # ENTROPIAS COMO CORRESPONDA, ES DECIR UNA ENTROPIA POR UNIQUE-REGRESSION Y LUEGO
    # COMBINAR COMO CORRESPONDE
    
    score <- mapply(FUN=function(func, p){
      # i<-7; func <- cmplxFuncs[i]; p <- pars[[i]]
      if(prnt){
        print(paste("function: ", func))
        print("parameters:")
        print(p)
      }
 
      ps <- c(p, list(vars=varsDag))
      
      return(do.call(func, ps))
    }, func=cmplxFuncs, p=pars) 
    
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

