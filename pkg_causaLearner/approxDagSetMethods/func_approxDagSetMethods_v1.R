# Approximation dag set methods

library(pcalg) # pc, skeleton, dag2cpdag, 
library(kpcalg) # kernelCItest
library(igraph) #V(), graph in get_n_cycles_directed_E
source("./pkg_causaLearner/utilities/func_dagStuff.R")

# Mooij 2009 - Regression by dependence minimization algorithm
causalOrdering <- function(Xtr, Xte, learner,  method, alpha=0){
  #print("enters causalOrdering")
  d <- ncol(Xtr)
  S <- colnames(Xtr) # 1:d
  ord <- as.numeric()
  
  totRegs <- d*(d+1)/2 -1
  #count <- 0
  # we'll choose one node at a time to be next in back to front causal ordering
  #pm0 <- proc.time()
  for (iter in d:2){
    # iter <- 3
    #print("**********************")
    #print(paste("iter: ",iter))
    
    # for each  nodes not chosen that are left we use the other nodes not chosen 
    # and see if we get independent residuals. We choose "most" independent
    # regression and take the dependent variable to be next in back to front
    # causal ordering
    pvals <- sapply(S, function(effectNodeTry){
      # effectNodeTry <- S[1]
      # print(paste("iter: ", iter,"effectNodeTry: ", effectNodeTry))
      #count <<- count  + 1
      causeNodesTry <- setdiff(S, effectNodeTry)
      
      trainDataO <- constructData(as.matrix(Xtr[,causeNodesTry]), Xtr[,effectNodeTry])
      
      learnerAux <- setParams(learner, trainData=trainDataO)
      learnerAux <- learnerAux$learn(learnerAux)
      testDataO <- constructData(as.matrix(Xte[,causeNodesTry]), Xte[,effectNodeTry])
      pred <- learnerAux$predict(learnerAux, data=testDataO)
      rs <- pred$gyh - pred$gy
      
      
      
      pval <- dhsic.test(Xte[, causeNodesTry], rs, method=method)$p.value
      #print(paste("p-value: ", pval))
      
      #timePast <- proc.time()-pm0
      #avgTime <- timePast/count
      #regsLeft <- totRegs - count
      #timeLeft <- regsLeft*avgTime
      
      #print(paste("Estimated time left causal ordering: ", round(timeLeft[3]/60,2), " mins."))
      
      return(pval)
    })
    
    mostEffectNode <- S[which.max(pvals)]
    if(pvals[mostEffectNode] < alpha){
      print("no consitent DAGs")
      return(NULL)
    }
    ord <- c(mostEffectNode, ord)
    S <- setdiff(S, mostEffectNode)
    #print(paste("causal order so far: ", paste(ord, collapse="-> ")))
  }
  ord <- c(S, ord)
  
  print(paste("final causal order: ", paste(ord, collapse="-> ")))
  
  #print("exits causalOrdering")
  return(ord)
}

minimalDAG <- function(Xtr, Xte, learner, method, co, alpha=0.05){
  print("enters minimalDAG")
  d <- ncol(Xtr)
  print("calculating causal ordering")
  
  ord <- causalOrdering(Xtr, Xte, learner, method=method, alpha=0)
  
  if(is.null(ord)){
    print("no consitent DAGs")
    return(NULL)
  }  
  
  if(co){
    parents <- lapply(1:d, function(j){
      node <- ord[j]
      if(j == 1){
        prnts_node <- as.numeric()
      } else{
        prnts_node <- ord[1:(j-1)]
      }
    })
  } else{
  
  totRegs <- d*(d-1)/2
  #count <- 0
  
  print("calculating parents")
  # for each node starts with all the parents and removes certain parents if we keep
  # independence of residuals with reduced inputs
  #pm0 <- proc.time()
  parents <- lapply(1:d, function(j){
    node <- ord[j]
    #print(paste("node: ", node))
    if(j  == 1){
      prnts_node <- as.numeric()
      #print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    } else{
      
      prnts_node <- ord[1:(j-1)]
      # we will try and take off as many inputs as we can maintaining independence
      for(k in (j-1):1){
        #print(paste("for node ", node, "testing taking off input: ", ord[k]))
        #count <<- count + 1
        notParent_node_try <- ord[k]
        prnts_node_try <- setdiff(prnts_node, notParent_node_try)
        
        if(length(prnts_node_try)==0){
          # if theres only one potential parent then to test if we should eliminate it we
          # check if the response variable "node" is independent of the explanatory variable
          # "prnts_node"
          pval <- dhsic.test(Xte[, prnts_node], Xte[,node])$p.value
        } else{
          
          trainDataO <- constructData(as.matrix(Xtr[,prnts_node_try]), Xtr[,node])
          learnerAux <- setParams(learner, trainDataO)
          learnerAux <- learnerAux$learn(learnerAux)
          testDataO <- constructData(as.matrix(Xte[,prnts_node_try]), Xte[,node])
          pred <- learnerAux$predict(learnerAux, data=testDataO)
          rs <- pred$gyh - pred$gy
          
          pval <- dhsic.test(Xte[,prnts_node_try], rs, method=method)$p.value
          
          
          #timePast <- proc.time()-pm0
          #avgTime <- timePast/count
          #regsLeft <- totRegs - count
          #timeLeft <- regsLeft*avgTime
          
          #print(paste("Estimated time left parents: ", round(timeLeft[3]/60,2), " mins."))
        }
        
        # if we maintain independence of residuals and inputs (p-value big) while not
        # using input k we take it off the parents list
        #print(paste("pval: ", pval))
        if(pval >= alpha){
          #print(paste("for node ", node,"taking off input: ", notParent_node_try))
          prnts_node <- setdiff(prnts_node, notParent_node_try)
        }
      }
      
    }
    print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    return(prnts_node)
  })
  
  }
  
  # make into a matrix and plot
  V <- ord
  parents2 <- sapply(parents, function(el) (V %in% el)*1)
  
  rownames(parents2) <- V
  colnames(parents2) <- V
  
  parents2 <- parents2[,colnames(Xtr)]
  parents2 <- parents2[colnames(Xtr),]
  
  parents3 <- getGraph(parents2)
  #plot(parents3)
  
  print("exits minimalDAG")
  return(parents2)
}

causalOrdering_cmem <- function(x, cmemLearner, msr){
  #print("enters causalOrdering_cmem")
  d <- ncol(x)
  S <- colnames(x) # 1:d
  ord <- as.numeric()
  
  cmemLearner <- eval(parse(text=cmemLearner))
  
  totRegs <- d*(d+1)/2 -1
  #count <- 0
  # we'll choose one node at a time to be next in back to front causal ordering
  #pm0 <- proc.time()
  for (iter in d:2){
    # iter <- 4
    print("**********************")
    print(paste("iter: ",iter))
    
    # for each  nodes not chosen that are left we use the other nodes not chosen 
    # and see if we get independent residuals. We choose "most" independent
    # regression and take the dependent variable to be next in back to front
    # causal ordering
    cmems <- sapply(S, function(effectNodeTry){
      # effectNodeTry <- S[1]
      #print(paste("iter: ", iter,"effectNodeTry: ", effectNodeTry))
      #count <<- count  + 1
      causeNodesTry <- setdiff(S, effectNodeTry)
      
      causes <- as.matrix(x[,causeNodesTry,drop=F])
      effects <- x[,effectNodeTry, drop=F]
      
      
      trainDataO <- constructData(causes, effects)
      cmemLearnerAux <- setParams(learner=cmemLearner, trainData=trainDataO)
      cmemLearnerAux <- cmemLearnerAux$learn(cmemLearnerAux)
      
      mesr <- do.call(msr, list(learner=cmemLearnerAux))
  
  
      
      #timePast <- proc.time()-pm0
      #avgTime <- timePast/count
      #regsLeft <- totRegs - count
      #timeLeft <- regsLeft*avgTime
      
      #print(paste("Estimated time left causal ordering: ", round(timeLeft[3]/60,2), " mins."))
      
      return(mesr)
    })
    
    mostEffectNode <- S[which.min(cmems)]
    # if(pvals[mostEffectNode] < alpha){
    #   print("no consitent DAGs")
    #   return(NULL)
    # }
    ord <- c(mostEffectNode, ord)
    S <- setdiff(S, mostEffectNode)
    #print(paste("causal order so far: ", paste(ord, collapse="-> ")))
  }
  ord <- c(S, ord)
  
  print(paste("final causal order: ", paste(ord, collapse="-> ")))
  
  #print("exits causalOrdering_cmem")
  return(ord)
}

permTestVarDep <- function(cmemLearner, effects, node, causes, causesNotPerm, causesPerm, msr, numPerms){
  print("enters permTestVarDep")
  cmemLearner <- eval(parse(text=cmemLearner))
  trainDataO <- constructData(causes, effects)
  cmemLearnerAux <- setParams(learner=cmemLearner, trainData=trainDataO)
  cmemLearnerAux <- cmemLearnerAux$learn(cmemLearnerAux)
  mesr <- do.call(msr, list(learner=cmemLearnerAux))
  
  n <- nrow(causes)
  
  set.seed(12)
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  print("calculating under null: the tested parent not true parent")
  msrPerm <- apply(rperms, 2, function(col){
    # i<-1; col <- rperms[,i]
    causesPerms <- causes[col,causesPerm,drop=FALSE]
    causesNotPerms <- causes[,causesNotPerm, drop=FALSE]
    trainDataO <- constructData(cbind(causesNotPerms, causesPerms), effects)
    cmemLearnerAux <- setParams(learner=cmemLearner, trainData=trainDataO)
    cmemLearnerAux <- cmemLearnerAux$learn(cmemLearnerAux)
    mesr <- do.call(msr, list(learner=cmemLearnerAux))
    return(mesr)
  })
  
  #hist(msrPerm, 30)
  #abline(v=mesr, col="red")
  
  #hist(msrPerm, xlim=range(c(msrPerm, mesr)), main=paste("node: ", node, " prnt test take-off: ", causesPerm))
  #abline(v=mesr, lwd=2, col="red")
  
  pval <- sum(mesr<msrPerm)/numPerms
  return(pval)
  print("exits permTestVarDep")
}

permTestVarDep2 <- function(cmemLearner, effects, node, causes, causesNotPerm, causesPerm, msr, numPerms){
  print("enters permTestVarDep2")
  
  cmemLearnerAux <- cmemLearner$learn(cmemLearner)
  mesr <- do.call(msr, list(learner=cmemLearnerAux))
  
  n <- nrow(causes)
  
  set.seed(12)
  rperms <- sapply(1:numPerms, function(i) sample(n))
  
  print("calculating under null: the tested parent not true parent")
  msrPerm <- apply(rperms, 2, function(col){
    # i<-1; col <- rperms[,i]
    causesPerms <- causes[col,causesPerm,drop=FALSE]
    causesNotPerms <- causes[,causesNotPerm, drop=FALSE]
    trainDataO <- constructData(cbind(causesNotPerms, causesPerms), effects)
    cmemLearnerAux <- cmemLearner
    cmemLearnerAux$hyperParams$trainData <- trainDataO
    cmemLearnerAux <- cmemLearnerAux$learn(cmemLearnerAux)
    mesr <- do.call(msr, list(learner=cmemLearnerAux))
    return(mesr)
  })
  
  #hist(msrPerm, 30)
  #abline(v=mesr, col="red")
  
  #hist(msrPerm, xlim=range(c(msrPerm, mesr)), main=paste("node: ", node, " prnt test take-off: ", causesPerm))
  #abline(v=mesr, lwd=2, col="red")
  
  pval <- sum(mesr<msrPerm)/numPerms
  return(pval)
  print("exits permTestVarDep")
}


# CO: causal ordering; PD: parent deletion
minimalDAG_cmem <- function(x, cmemLearner, msrCO, msrPD, numPerms, co, alpha=0.05){
  print("enters minimalDAG_cmem")
  d <- ncol(x)
  print("calculating causal ordering")
  
  ord <- causalOrdering_cmem(x, cmemLearner, msr=msrCO)
  
  if(is.null(ord)){
    print("no consitent DAGs")
    return(NULL)
  }  
  
  if(co){
    parents <- lapply(1:d, function(j){
      node <- ord[j]
      if(j == 1){
        prnts_node <- as.numeric()
      } else{
        prnts_node <- ord[1:(j-1)]
      }
    })
  } else{
  
  totRegs <- d*(d-1)/2
  #count <- 0
  
  print("calculating parents")
  # for each node starts with all the parents and removes certain parents if we keep
  # independence of residuals with reduced inputs
  #pm0 <- proc.time()
  parents <- lapply(1:d, function(j){
    # j <- 2
    # print(paste("j: ", j, sep=""))
    node <- ord[j]
    print(paste("node: ", node))
    if(j  == 1){
      prnts_node <- as.numeric()
      #print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    } else{
      
      prnts_node <- ord[1:(j-1)]
      causes <- as.matrix(x[,prnts_node,drop=F])
      effects <- x[,node, drop=F]
      # we will try and take off as many inputs as we can maintaining independence
      for(k in (j-1):1){
        # k <- 1
        #print(paste("k: ", k, sep=""))
        print(paste("for node ", node, "testing taking off input: ", ord[k]))
        #count <<- count + 1
        notParent_node_try <- ord[k]
        prnts_node_try <- setdiff(prnts_node, notParent_node_try)
        
        # to test if we shd eliminate a parent we apply measure to sample and compare,
        # measure when applying it to "unordered" (in the parent var) samples

        #causesPerm <- notParent_node_try
        
        
        pval <- permTestVarDep(cmemLearner, effects, causes, node=node, causesNotPerm=prnts_node_try, causesPerm=notParent_node_try, msr=msrPD, numPerms)
        
        # if we maintain independence of residuals and inputs (p-value big) while not
        # using input k we take it off the parents list
        print(paste("pval: ", pval))
        if(pval >= alpha){
          print(paste("for node ", node,"taking off input: ", notParent_node_try))
          prnts_node <- setdiff(prnts_node, notParent_node_try)
        }
      }
      
    }
    print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    return(prnts_node)
  })
  
  }
  
  # make into a matrix and plot
  V <- ord
  parents2 <- sapply(parents, function(el) (V %in% el)*1)
  
  rownames(parents2) <- V
  colnames(parents2) <- V
  
  parents2 <- parents2[,colnames(x)]
  parents2 <- parents2[colnames(x),]
  
  parents3 <- getGraph(parents2)
  #plot(parents3)
  
  print("exits minimalDAG_cmem")
  return(parents2)
}


# generalte all dags
exhaustive <- function(data, trueDAG=NULL, pars=NULL){
  p <- ncol(data)
  DAGset <- genAllDAGS(p)
  return(DAGset)
}

# generalte 2 dags from a dag with 6 nodes: x_t, y_t, x_{t-1}, y_{t-1}
# z_xy_{t-1} z_yx_{t-1} where z_xy and z_yx are estimates of z
timeSeriesXYZ_ts <- function(data, trueDAG, pars=NULL){
  
  dagMat <- matrix(0,8,8)
  colnames(dagMat) <- rownames(dagMat) <- c("x","y","x0","y0","z_xy","z_yx","z_xx","z_yy")
  dagMat_xy <- dagMat_yx <- dagMat
  indxMat_xy <- matrix(c(c(3,2), c(4,2), c(5,2), c(3,1),c(7,1)),5,2,byrow=T)
  indxMat_yx <- matrix(c(c(3,1), c(4,1), c(6,1), c(4,2), c(8,2)),5,2,byrow=T)
  dagMat_xy[indxMat_xy] <- 1
  dagMat_yx[indxMat_yx] <- 1
  dagArr <- abind(dagMat_xy, dagMat_yx, along=3)
  
  dgNms <- getHypID(dagArr)$id
  dimnames(dagArr) <- list(from=dimnames(dagMat)[[1]], to=dimnames(dagMat)[[2]], dag=dgNms)
  
  return(dagArr)
}
timeSeriesXYZ_cmpr <- function(data, trueDAG, pars=NULL){
  
  dagMat <- matrix(0,8,8)
  colnames(dagMat) <- rownames(dagMat) <- c("x","y","x0","y0","z_xy","z_yx","z_xx","z_yy")
  dagMat_xy_simp <- dagMat_yx_simp <- dagMat
  indxMat_xy_simp <- matrix(c(1,2),1,2,byrow=T)
  indxMat_yx_simp <- matrix(c(2,1),1,2,byrow=T)
  dagMat_xy_simp[indxMat_xy_simp] <- 1
  dagMat_yx_simp[indxMat_yx_simp] <- 1
  dagArr <- abind(dagMat_xy_simp, dagMat_yx_simp, along=3)
  
  dgNms <- getHypID(dagArr)$id
  dimnames(dagArr) <- list(from=dimnames(dagMat)[[1]], to=dimnames(dagMat)[[2]], dag=dgNms)
  
  return(dagArr)
}

# narrows it down to the true markov equivalence class
oracleME <- function(data=NULL, trueDAG, pars=NULL){
  # trueDAG <- matrix(c(0,0,1,0),2,2)
  MEdagTrue <- getMarkovEquivClass(trueDAG)
  #p <- ncol(trueDAG)
  #trueCPDAG <- dag2cpdag(trueDAG)
  #MEdagTrue <- pdag2allDags(trueCPDAG)$dags
  #MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
  return(MEdagTrue)
}

pcSkeleton <- function(data, trueDAG=NULL, pars){
  V <- colnames(data)
  uDAG <- pcalg:::skeleton(suffStat=list(data=data, ic.method=pars$ic.method), indepTest=eval(parse(text=pars$indepTest)), alpha=pars$alpha, labels=V)
  uDAG <- amat(as.bn(uDAG))
  DAGset <- udag2dags(uDAG)
  return(DAGset)
}

pcMarkovEquiv <- function(data, trueDAG=NULL, pars){
  V <- colnames(data)
  pDAG <- pc(suffStat=list(data=data, ic.method=pars$ic.method), indepTest=eval(parse(text=pars$indepTest)), alpha=pars$alpha, labels=V)
  pDAG <- amat(as.bn(pDAG))
  DAGset <- pdag2dags(pDAG)
  return(DAGset)
}


# I COULD FIX THIS SO THAT THE DATA REGIMENT (HOLDOUT OR RECYCLE) USES TRAIN PREDICTIONS ON WHOLE SAMPLE OR
# CROSS VALIDATED OUT OF SAMPLE PREDICTIONS INSTEAD OF DIVIDING INTO TWO
#mooij et al 2009 minimal dag algo with ANMs
minDAG <- function(data, trueDAG=NULL, pars){
  dataReg <- dataRegime(data, type=pars$dataRegime)
  xTrain <- dataReg$train
  xTest <- dataReg$test
  minDAG <- minimalDAG(Xtr=xTrain, Xte=xTest, learner=eval(parse(text=pars$learner)), method=pars$method, co=pars$resolution=="co")
  # obtain generating set of "super-DAGs" (those obtained by adding edges but not taking away)
  # 1. identify non-present edges
  # 2. obtain all super-DAGs
  
  #DAGset <- amat(as.bn(minDAG))
  DAGset <- minDAG
  if(pars$resolution=="addEdges"){
    DAGset <- mindag2dags(DAGset)
  } else if(pars$resolution %in% c("dag","co")){
    dim(DAGset) <- c(dim(DAGset),1)
  }
  return(DAGset)
}

#adaptation of above for cmems (including hsic)
minDAG_cmem <- function(data, trueDAG=NULL, pars){
  
  minDAG <- minimalDAG_cmem(x=data, cmemLearner=pars$cmemLearner, msrCO=pars$msrCO, msrPD=pars$msrPD, numPerms=pars$numPerms, co=pars$resolution=="co")
  # obtain generating set of "super-DAGs" (those obtained by adding edges but not taking away)
  # 1. identify non-present edges
  # 2. obtain all super-DAGs
  
  DAGset <- minDAG #amat(as.bn(minDAG))
  if(pars$resolution=="addEdges"){
    DAGset <- mindag2dags(DAGset)
  } else if(pars$resolution %in% c("dag","co")){
    dim(DAGset) <- c(dim(DAGset),1)
  }
  return(DAGset)
}


getKnownEdgesFromMEC <- function(MEC){
  
  
  knownEdge    <- apply(MEC, c(1,2), function(vec) prod(vec))
  knownNonEdge <- apply(MEC, c(1,2), function(vec) prod(!vec)*1)
  unknownEdge <- (!knownEdge+knownNonEdge)*1
  
  # edges we know to exist in all dags in MEC
  
  
  knownEdge <- melt(knownEdge)
  #knownEdge[,1] <- as.character(knownEdge[,1])
  #knownEdge[,2] <- as.character(knownEdge[,2])
  knownEdge <- knownEdge[which(knownEdge$value==1),]
  knownEdge <- knownEdge[,c(1,2)]
  colnames(knownEdge) <- c("from","to")
  knownEdge$bin <- rep(1, nrow(knownEdge))
  
  
  # edges we know not to exist in all dags in MEC
  knownNonEdge <- melt(knownNonEdge)
  #knownNonEdge[,1] <- as.character(knownNonEdge[,1])
  #knownNonEdge[,2] <- as.character(knownNonEdge[,2])
  knownNonEdge <- knownNonEdge[which(knownNonEdge$value==1),]
  knownNonEdge <- knownNonEdge[,c(1,2)]
  colnames(knownNonEdge) <- c("from","to")
  knownNonEdge$bin <- 0
  
  knownEdge <- rbind(knownEdge, knownNonEdge)
  

  all(is.na(cast(knownEdge, from~to, value="bin")[,2:(dim(MEC)[1]+1)])*1==unknownEdge)
  
  # edges we don't know if they are ->, <- or non-edge
  direction <- upper.tri(unknownEdge)*1+lower.tri(unknownEdge)*-1
  unknownEdge <- cbind(melt((unknownEdge | t(unknownEdge))*1), direction=melt(direction)[,3])
  unknownEdge <- unknownEdge[which(unknownEdge$value==1 & unknownEdge$direction==1), c(1,2)]
  colnames(unknownEdge) <- c("i","j")
  
  knownEdges <- list(known=knownEdge, unknown=unknownEdge)
  return(knownEdges)
}


pairwise_cmem <- function(data, trueDAG=NULL, pars){
  
  # pars <- list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=FALSE)
  if(!is.null(trueDAG) & pars$withinMEC){
    DAGset <- pairwise_cmem_withinMEC(x=data, trueDAG, cmemLearner=pars$cmemLearner, msrCO=pars$msrCO)
    
    
  } else{
    DAGset <- pairwise_cmem_notWithinMEC(x=data, cmemLearner=pars$cmemLearner, msrCO=pars$msrCO, msrPD=pars$msrPD, 
                                             numPerms=pars$numPerms, resolution=pars$resolution)
  }
  dim(DAGset) <- c(dim(DAGset),1)
  return(DAGset)
}


pairwise_cmem_withinMEC <- function(x, trueDAG, cmemLearner, msrCO){
  print("enters pairwise_cmem_withinMEC")
  MEdag <- getMarkovEquivClass(trueDAG)
  knownEdges <- getKnownEdgesFromMEC(MEC=MEdag)
  
  
  unknownEdges <- knownEdges$unknown
  knownEdges <- knownEdges$known
  
  dag <- cast(knownEdges, from~to, value="bin", add.missing=T)
  dag <- dag[1:ncol(x),]
  nms <- dag[,1]
  dag <- dag[,2:(ncol(x)+1)]  
  dag <- as.matrix(as.matrix(dag))
  rownames(dag) <- colnames(dag) <- as.character(nms)
  
  
  S <- colnames(x)
  d <- ncol(x)
  diag(dag) <- 0
  
  cmemLearnerAux <- eval(parse(text=cmemLearner))
  
  
  # A. Pairwise setting of unknown edges 
  
  if(nrow(unknownEdges)>0 ){
    
    unknownEdges$S12 <- NA
    unknownEdges$S21 <- NA
    
    # calculate measures of unknown edges
    for(i in 1:nrow(unknownEdges)){
      # i <- 1; unknownEdges[i,]
      causes12 <- setdiff(S, unknownEdges[i,2])
      effect12 <- unknownEdges[i,2]
      causes21 <- setdiff(S, unknownEdges[i,1])
      effect21 <- unknownEdges[i,1]
      
      trainData12 <- constructData(x[,causes12, drop=F], x[,effect12, drop=F])
      cmemLearner12 <- setParams(learner=cmemLearnerAux, trainData=trainData12)
      cmemLearner12 <- cmemLearner12$learn(cmemLearner12)
      mesr12 <- do.call(msrCO, list(learner=cmemLearner12))
      
      trainData21 <- constructData(x[, causes21, drop=F], x[,effect21, drop=F])
      cmemLearner21 <- setParams(learner=cmemLearnerAux, trainData=trainData21)
      cmemLearner21 <- cmemLearner21$learn(cmemLearner21)
      mesr21 <- do.call(msrCO, list(learner=cmemLearner21))
      
      unknownEdges[i,"S12"] <- mesr12
      unknownEdges[i,"S21"] <- mesr21
     
    }
    
    
    unknownEdges$dif <- unknownEdges$S12 - unknownEdges$S21 
    unknownEdges <- unknownEdges[order(unknownEdges$dif, decreasing=F),]
    
    
    markovEquiv <- function(MEdag, dag){
      #mat <- MEdag[,,5]
      res <- apply(MEdag, "dag", function(mat) all(dag==mat, na.rm=T))
      #res <- any(res)
      return(res)
    }
    
    
    # we make a binary tree of ME class where each level corresponds to each unknown edge
    # level 1 corresponds to node with smallest s12-s21 (smallest "dif"), level 2 to second smallest etc
    # we represent node tree with lists list(node=, dags=, left=list(), right())
    
    MEdagAux <- MEdag
    
    j <- 0
    while(dim(MEdagAux)[3]>1){
      j <- j + 1
      dagAux12 <- dagAux21 <- dag
      dagAux12[unknownEdges[j,1],unknownEdges[j,2]] <- 1
      dagAux12[unknownEdges[j,2],unknownEdges[j,1]] <- 0
      dagAux21[unknownEdges[j,1],unknownEdges[j,2]] <- 0
      dagAux21[unknownEdges[j,2],unknownEdges[j,1]] <- 1
      
      indx12 <- markovEquiv(MEdag, dagAux12)
      indx21 <- markovEquiv(MEdag, dagAux21)
      MEdagAux12  <- MEdag[,,which(indx12), drop=F] 
      MEdagAux21 <-  MEdag[,,which(indx21), drop=F]
      
      if(dim(MEdagAux12)[3]==0 & dim(MEdagAux21)[3]==0){
        stop("Error!")
      } else if(dim(MEdagAux12)[3]==0 & dim(MEdagAux21)[3]>0){
        MEdagAux <- MEdagAux21
        dag <- dagAux21
      } else if(dim(MEdagAux12)[3]>0 & dim(MEdagAux21)[3]==0){
        MEdagAux <- MEdagAux12
        dag <- dagAux12
      } else{
        # this is the dim(MEdagAux12)[3]>0 & dim(MEdagAux21)[3]>0 cse
        # here we check whether S12 or S21 is smaller
        if(unknownEdges$dif[j]<=0){
          MEdagAux <- MEdagAux12
          dag <- dagAux12
        } else{
          MEdagAux <- MEdagAux21
          dag <- dagAux21
        }
      }
    }
    
    if(dim(MEdagAux)[3]!=1) stop("Error2 !")
  
    dag <- MEdagAux[,,1]
  }
  
  print("exits pairwise_cmem_withinMEC")
  return(dag)
}


pairwise_cmem_notWithinMEC <- function(x, cmemLearner, msrCO, msrPD, numPerms, resolution, alpha=0.05){
  unknownEdges <- t(apply(combinations(ncol(x), 2), 1, function(row) colnames(x)[row]))
  unknownEdges <- as.data.frame(unknownEdges)
  colnames(unknownEdges) <- c("i","j")
  
  dag <- diag(ncol(x))*NA
  diag(dag) <- 0  
  colnames(dag) <- rownames(dag) <- colnames(x)
  
  cmemLearnerAux <- eval(parse(text=cmemLearner))
  
  # A. Pairwise setting of unknown edges 
  S <- colnames(x) 
  d <- ncol(x)
  
  for(i in 1:nrow(unknownEdges)){
      # i <- 1; unknownEdges[i,]
      # print(paste("i:", i))
      causes12 <- setdiff(S, unknownEdges[i,2])
      effect12 <- unknownEdges[i,2]
      causes21 <- setdiff(S, unknownEdges[i,1])
      effect21 <- unknownEdges[i,1]
    
      trainData12 <- constructData(x[,causes12, drop=F], x[,effect12, drop=F])
      cmemLearner12 <- setParams(learner=cmemLearnerAux, trainData=trainData12)
      cmemLearner12 <- cmemLearner12$learn(cmemLearner12)
      mesr12 <- do.call(msrCO, list(learner=cmemLearner12))
    
      trainData21 <- constructData(x[, causes21, drop=F], x[,effect21, drop=F])
      cmemLearner21 <- setParams(learner=cmemLearnerAux, trainData=trainData21)
      cmemLearner21 <- cmemLearner21$learn(cmemLearner21)
      mesr21 <- do.call(msrCO, list(learner=cmemLearner21))
    
      dagAux <- dag
      if(mesr12 < mesr21){
        dagAux[unknownEdges[i,1],unknownEdges[i,2]] <- 1
        dagAux[unknownEdges[i,2],unknownEdges[i,1]] <- 0
      } else{
        dagAux[unknownEdges[i,1],unknownEdges[i,2]] <- 0
        dagAux[unknownEdges[i,2],unknownEdges[i,1]] <- 1
      }
      dag <- dagAux
    }
  
    
  # B. Check that it is a valid graph and elimnate cycles if not
  # dag <- matrix(c(0,1,0,1,0, 0,0,1,0,1, 1,0,0,0,0, 0,0,0,0,1, 1,0,0,1,0), p, p); dimnames(dag) <- list(from=letters[1:p], to=letters[1:p]) ; plot(as(dag, "graphNEL"))
  if(!isValidGraph(dag, "dag")){
    stop("cycles generated!!")
    print("cycles generated!!")
    # identify cycles
    
    edges <- melt(dag)
    edges <- edges[which(edges$value==1),1:2]
    
    edges_g <- graph(t(as.matrix(edges)))
    #plot(edges_g)
    
    # cycles expressed in terms of V(edges_g) !!!, not in terms of (dag)
    cycles <- lapply(1:p, function(n) get_n_cycles_directed_E( as(edges_g, "igraph"), n, list_cycles=T))
    V_edges_g <- names(as.factor(V(edges_g)))
    
    # break each cycle by 
    breakLinks <- lapply(cycles, function(cycle){
      # cycle <- cycles[[3]]
      if(cycle$count>=1){
        cycs <- cycle$cycles
        cycs <- cbind(cycs, cycs[,1])
        
        breakLinks <- lapply(1:cycle$count, function(i){
          # i <- 1
          msrs <- sapply(1:(ncol(cycs)-1), function(j){
            # j <- 1
            effect <- V_edges_g[cycs[i,j+1]]
            causes <- rownames(dag)[which(dag[,effect]==1)]
            
            trainData0 <- constructData(x[, causes, drop=F], x[, effect, drop=F])
            cmemLearnerij <- setParams(learner=cmemLearnerAux, trainData=trainData0)
            cmemLearnerij <- cmemLearnerij$learn(cmemLearnerij)
            mesr <- do.call(msrCO, list(learner=cmemLearnerij))  
            
          })
          indx.break <- which.max(msrs)
          breakLink <- V_edges_g[cycs[i,indx.break:(indx.break+1)]]
          return(breakLink)
        })
      } else{
        res <- NULL
      }
    })
    
    breakLinks <- unlist(breakLinks)
    breakLinks <- match(breakLinks, colnames(dag))
    breakLinks <- matrix(breakLinks, nrow=length(breakLinks)/2, ncol=2, byrow=T)
    
    dag[breakLinks] <- 0
    #plot(as(dag, "graphNEL"))
    
  }
  
  
  
  # C. Take off unnecessary parents using msrPD (if resolution=="max")
  if(resolution=="max"){
    
    ord <- node.ordering(as.bn(as(as.matrix(dag), "graphNEL")))
    # see igraph:::topo_sort  also
    
    parents <- lapply(1:d, function(j){
      # j <- 2
      # print(paste("j: ", j, sep=""))
      node <- ord[j]
      print(paste("node: ", node))
      if(j  == 1){
        prnts_node <- as.numeric()
        #print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
      } else{
        
        prnts_node <- rownames(dag)[which(dag[,node]==1)]
        causes <- as.matrix(x[,prnts_node,drop=F])
        effects <- x[,node, drop=F]
        # we will try and take off as many inputs as we can maintaining independence
        for(notParent_node_try in prnts_node){
          # notParent_node_try <- prnts_node[1]
          print(paste("for node ", node, "testing taking off input: ", notParent_node_try))
          
          prnts_node_try <- setdiff(prnts_node, notParent_node_try)
          
          # to test if we shd eliminate a parent we apply measure to sample and compare,
          # measure when applying it to "unordered" (in the parent var) samples
          
          #causesPerm <- notParent_node_try
          
          
          pval <- permTestVarDep(cmemLearner, effects, causes, node=node, causesNotPerm=prnts_node_try, causesPerm=notParent_node_try, msr=msrPD, numPerms)
          
          # if we maintain independence of residuals and inputs (p-value big) while not
          # using input k we take it off the parents list
          print(paste("pval: ", pval))
          if(pval >= alpha){
            print(paste("for node ", node,"taking off input: ", notParent_node_try))
            prnts_node <- setdiff(prnts_node, notParent_node_try)
          }
        }
        
      }
      print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
      return(prnts_node)
    })
   
    V <- ord
    dag <- sapply(parents, function(el) (V %in% el)*1)
    
    rownames(dag) <- V
    colnames(dag) <- V
    
    dag <- dag[,colnames(x)]
    dag <- dag[colnames(x),] 
    
  } 
  
  return(dag)
  
}

pairwise_cmemComp  <- function(data, trueDAG=NULL, pars){
  
  # pars <- list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=FALSE)

    DAGset <- pairwise_cmemComp_notWithinMEC(x=data, hypScorer=pars$hypScorer, ppTab=pars$ppTab, msrCO=pars$msrCO,
                                             msrPD=pars$msrPD, numPerms=pars$numPerms, resolution=pars$resolution,
                                             dataNm=pars$dataNm, folderSave=pars$folderSave)
  
  dim(DAGset) <- c(dim(DAGset),1)
  return(DAGset)
}

pairwise_cmemComp_notWithinMEC <- function(x, hypScorer, ppTab, msrCO, msrPD, numPerms, resolution, dataNm, folderSave, alpha=0.1){
  
  hypScorer <- eval(parse(text=hypScorer))
  cmemLearner1 <- hypScorer$cmemLearner1
  cmemLearner2 <- hypScorer$cmemLearner2
  
  nodes <- colnames(x)
  p <- length(nodes)
  hypArray <- sapply(nodes, function(nodeTo){
    # nodeTo <- nodes[2]
    #print(nodeTo)
    res <- matrix(0, p, p)
    colnames(res) <- rownames(res) <- 1:p
    res[,nodeTo] <- 1
    res[nodeTo,nodeTo] <- 0
    dimnames(res) <- list(nodeFrom=1:p, nodeTo=1:p)
    return(res)
  }, simplify="array")
  
  names(dimnames(hypArray))[3] <- "dag"
            
  scores <- cmem_hypScorer_comp(x, hypArray, cmemLearner1, cmemLearner2, ppTab=ppTab, plot=FALSE, dataNm=dataNm, folderSave=folderSave)
  
  # the lower the score the better fit that regression is meaning the node is
  # more of an "effect" so we should order form highest to lowest
  
  ord <- as.character(order(scores[,msrCO], decreasing=T))
  
  parents <- lapply(1:p, function(j){
    node <- ord[j]
    if(j == 1){
      prnts_node <- as.numeric()
    } else{
      prnts_node <- ord[1:(j-1)]
    }
    return(prnts_node)
  })
  
  V <- ord
  parents2 <- sapply(parents, function(el) (V %in% el)*1)
  
  rownames(parents2) <- V
  colnames(parents2) <- V
  
  parents2 <- parents2[,nodes]
  parents2 <- parents2[nodes,]
  
  dag <- parents2
  
  
  if(resolution=="max"){
    # Take off unnecessary parents using msrPD (if resolution=="max")
    
    #ord2 <- node.ordering(as.bn(as(as.matrix(dag), "graphNEL")))
    # see igraph:::topo_sort  also
    
    parents <- lapply(1:p, function(j){
    # j <- 2
    # print(paste("j: ", j, sep=""))
    node <- ord[j]
    print("***************************")
    print(paste("iter: ", j))
    print(paste("node: ", node))
    if(j  == 1){
        prnts_node <- as.numeric()
        #print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
    } else{
        
      prnts_node <- rownames(dag)[which(dag[,node]==1)]
      #causes <- as.matrix(x[,prnts_node,drop=F])
      causes <- as.matrix(x[,setdiff(nodes, node),drop=F])
      effects <- x[,node, drop=F]
      # we will try and take off as many inputs as we can maintaining independence
      for(notParent_node_try in prnts_node){
        # notParent_node_try <- prnts_node[1]
        print(paste("for node ", node, "testing taking off input: ", notParent_node_try))
          
        #prnts_node_try <- setdiff(prnts_node, notParent_node_try)
        prnts_node_try <- setdiff(nodes, c(node, notParent_node_try))
          
        # to test if we shd eliminate a parent we apply measure to sample and compare,
        # measure when applying it to "unordered" (in the parent var) samples
          
        #causesPerm <- notParent_node_try
          
          # Use params cmemLearner  for effects on everything else
        regressorsChar <- paste(setdiff(nodes, node), collapse="-")
        regressionChar <- paste(node, "on", regressorsChar, sep="")
        nm <- paste(dataNm, cmemLearner2, "2",sep="_")
        fileSave <- paste(nm, "_", regressionChar, ".RData", sep="")
        load(file=paste(folderSave, fileSave, sep=""))
        
        pval <- permTestVarDep2(cmemLearner=cmemLearnerAux, effects, node, causes, causesNotPerm=prnts_node_try, causesPerm=notParent_node_try, msr=msrPD, numPerms)
          
          # if we maintain independence of residuals and inputs (p-value big) while not
          # using input k we take it off the parents list
          print(paste("pval: ", pval))
          if(pval >= alpha){
            print(paste("for node ", node,"taking off input: ", notParent_node_try))
            prnts_node <- setdiff(prnts_node, notParent_node_try)
          }
        }
        
      }
      print(paste("parents of node ", node, " are ", paste(prnts_node, collapse=", ")))
      return(prnts_node)
    })
    
    V <- ord
    dag <- sapply(parents, function(el) (V %in% el)*1)
    
    rownames(dag) <- V
    colnames(dag) <- V
    
    dag <- dag[,nodes]
    dag <- dag[nodes,] 
    
    #check, there shouldn't be any edges on this dag that are not on non-edge deleted 
    # one (parents2)
    any(dag & !parents2) 
    
  }
  
  return(dag)
  
}
