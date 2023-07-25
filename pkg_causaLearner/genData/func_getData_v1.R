# obtain data on which to test causal learner package

library(pcalg) #
library(unifDAG) #unifDAG
library(gRbase) # topoSort
library(purrr) # rbernoulli
library(bnlearn) # amat in (simRandAddSEM)
library(graph) # nodes in (simRandAddSEM)
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_kernels/func_kernel_pkg.R")

# get Tubingen Cause Effect Pairs

createTCEPList <- function(pairs, meta, folder){
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- pairs
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  
  
  
  
  
  xyCols <- lapply(pairs, function(pair){
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- pairs
  
  dags <- lapply(fileFolders, function(el) dag) 
  xs <- fileFolders
  
  ns <-xyCols
  res <- list(dags=dags, xs=xs, noiss=ns, names=pairs)
  return(res)
}

numDataPairs <- function(pairs, folder){
  num <- sapply(pairs, function(pair){ 
    
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    dat <- read.csv(fileFolder, sep="", header=F)
    n <- nrow(dat)
    return(n)
  })
  return(num)
}

summaryDataPairs <- function(meta, pairs, folder){
  sapply(pairs, function(pair){ 
    indxPair <- which(meta$pairNumber==pair)
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    dat <- read.csv(fileFolder, sep="", header=F)
    print("********************")
    print(paste("pair:", pair))
    print("summary cause:")
    print(summary(dat[,meta[indxPair, "firstCauseCol"]]))
    print("summary effect:")
    print(summary(dat[,meta[indxPair,"firstEffectCol"]]))
    print("dim(x)")
    print(dim(dat))
    return(NULL)
  })
  return(NULL)
}


# Simulate Additive Noise SEMs

# simulate from a given sem defined as a dag, set of functions and noise distributions
# defined as n = a*dist(n_iid, pars)^b, dag could be obtained from functions but we
# allow for this redundancy

simSEM <- function(sem){
  if(all(names(sem$simPars$nodes)!=nodes(sem$dag))) stop("names (sem,ns) dont match")
  if(all(names(sem$funcs)!=paste("f",nodes(sem$dag),sep=""))) stop("names (nodes, funcs) dont match")
  
  n <- sem$simPars$n
  p <- length(sem$simPars$nodes)
  
  
  
  ns <- mapply(FUN=function(dist, pars, a, b){
    
    a*(do.call(dist, c(list(n=n), pars)))^b
    
    
  }, dist=lapply(sem$simPars$nodes, function(el) el$dist), 
  pars=lapply(sem$simPars$nodes, function(el) el$pars),
  a=lapply(sem$simPars$nodes, function(el) el$a),
  b=lapply(sem$simPars$nodes, function(el) el$a))
  
  colnames(ns) <- nodes(sem$dag)
  
  dagAmat <- amat(as.bn(sem$dag))
  Nodes <- nodes(sem$dag)
  topoNodes <- topoSort(sem$dag)
  
  x <- matrix(NA, n, p)
  colnames(x) <- Nodes
  indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
  x[,Nodes[indxNoParents]] <- ns[,Nodes[indxNoParents]]
  
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    
    indxParents <- which(dagAmat[,nd]==1)
    Parents <- Nodes[indxParents]
    argsNd <- cbind(x[,Parents], ns[,nd])
    argsNd <- as.list(as.data.frame(argsNd))
    names(argsNd) <- c(Parents,"n")
    
    x[,nd] <- do.call(sem$funcs[[paste("f",nd,sep="")]], argsNd)
  }
  
  return(list(dag=dagAmat,x=x,n=ns))
}

simSEMs <- function(q, sem){
  
  dagMat <- amat(as.bn(sem$dag))
  dags <- lapply(1:q, function(i) dagMat)
  
  sims <- lapply(1:q, function(i) simSEM(sem))
  
  nms <- 1:q
  
  xs <- lapply(sims, function(el) el$x)
  ns <- lapply(sims, function(el) el$n)
  
  names(xs) <- names(ns)  <- names(dags) <- nms
  
  
  
  return(list(dags=dags, xs=xs, noiss=ns, names=nms))  
}

# simulates a random function from hilbert space, indexed by (z,alpha), where z are the 
# inducing points
simfxH <- function(numVars, pctHidVars=1, sigma, sigmaErr=0, dist, distPars, geU=function(y, n, scale, const) y){
  
  m <- 100 #number of inducing points
  numHiddenVars <- floor(pctHidVars*numVars)
  numTotVars <- numVars + numHiddenVars
  
  # inducing points - non- hidden
  z <- matrix(rnorm(m*(numVars)), m, numVars)  
  # dummy points to get an idea of bias and variance of simulated points
  # we will evaluate fz() on w to obtain y
  w <- matrix(runif(m*(numVars)), m, numVars)
  
  
  distPars2 <- c(list(n=m*numHiddenVars), distPars)
  
  # inducing points - hidden/latent
  zh <- matrix(do.call(dist, distPars2), m, numHiddenVars)
  wh <- matrix(apply(zh, 2, mean), m, numHiddenVars, byrow=T)
  
  
  
  alpha <- rnorm(m, 0 , 1) 
  
  Kx <-  kern_rbf2(z, zh, w, wh, sigma, sigmaErr)
  y <- Kx%*%alpha
  return(list(z=z, zh=zh, alpha=alpha, const=min(y)+1, scale=(max(y)-min(y))/2, 
              numHiddenVars=numHiddenVars, sigma=sigma, sigmaErr=sigmaErr, dist=dist, distPars=distPars, geU=geU))  
}

# applies random hilbert function to obtain f_{z,alpha}(x)
applySimfxH <- function(x, simfx){
  z <- simfx$z
  zh <- simfx$zh
  alpha <- simfx$alpha
  sigma  <- simfx$sigma
  sigmaErr <- simfx$sigmaErr
  const <- simfx$const
  scale <- simfx$scale
  geU <- simfx$geU
  
  n <- nrow(x)
  dist <- simfx$dist
  distPars <- list()
  for(i in 1:length(simfx$distPars)) distPars <- c(distPars, list(simfx$distPars[[i]]))
  names(distPars) <- names(simfx$distPars)
  numHiddenVars <- simfx$numHiddenVars
  distPars2 <- c(distPars, n=numHiddenVars*n)
  xh <- do.call(dist,distPars2)
  xh <- matrix(xh, n, numHiddenVars)
  Kx <- kern_rbf2(x, xh, z, zh, sigma, sigmaErr)
  res <- as.numeric(Kx %*% alpha)
  
  print(paste("min train:", simfx$const, ", min test:", min(res)))
  
  distPars2 <- c(distPars, n=n)
  no <- do.call(dist,distPars2)
  res <- do.call(geU, list(y=res, const=const, scale=scale, n=no)) 
  return(res)
}

# writes f(z,alpha,x,n) so that z, alpha are written internally as parameters and we may have a 
# function of the type f(x,n) for use in a SEM
writefxH <- function(simfx){
  
  nmsDistPars <- names(distPars)
  valsDistPars <- unlist(distPars)
  strDistPars <- paste(paste(nmsDistPars, valsDistPars, sep="="), collapse=", ")
  
  
  
  fx <- eval(parse(text=
                     paste("fx <- function(x,n){z <- matrix(c(", 
                           paste(as.numeric(simfx$z), collapse=", "), "),", nrow(simfx$z), ", ", ncol(simfx$z), "); ",  
                           "zh <- matrix(c(", paste(as.numeric(simfx$zh), collapse=", "),  "),", nrow(simfx$zh), ", ", ncol(simfx$zh), "); ",  
                           "alpha <- c(", paste(simfx$alpha, collapse=", "), "); ",
                           "sigma <- ", simfx$sigma, "; ",
                           "sigmaErr <- ", simfx$sigmaErr, "; ",
                           "const <- ", simfx$const, "; ",
                           "scale <- ", simfx$scale, "; ", 
                           "geU <- ", paste(deparse(simfx$geU), collapse=""), "; ",
                           "N <- nrow(x); ",
                           "dist <- ", simfx$dist, "; ",
                           "distPars <- list(", strDistPars, "); ",
                           "numHiddenVars <- ", simfx$numHiddenVars, "; ",
                           "distPars2 <- c(distPars, list(n=numHiddenVars*N)); ",
                           "xh <- do.call(dist,distPars2); ",
                           "xh <- matrix(xh, N, numHiddenVars); ",
                           "Kx <- kern_rbf2(x, xh, z, zh, sigma, sigmaErr); ",
                           "res <- as.numeric(Kx %*% alpha); ", 
                           "distPars2 <- c(distPars, n=N);",
                           "no <- do.call(dist,distPars2);",
                           "res <- do.call(geU, list(y=res, const=const, scale=scale, n=no)); ", 
                           "return(res)}", sep="")))
  
  return(fx)
}


simRandAddSEM <- function(p, n, nodes, sigma, dagMat=NULL, markovEquiv=FALSE){
  # p <- 4; n <- 100 
  # nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
  # nodes <- rep(list(nodes),4)
  
  m <- 100 #number of inducing points
  
  if(!is.null(dagMat)){
    dag <- as(dagMat, "graphNEL")
  } else{
    dag <- unifDAG(p)
  }
  
  if(markovEquiv){
    MEdagTrue <- getMarkovEquivClass(dagMat)
    #cpdag <- dag2cpdag(dag)
    #cpdagMat <- as(cpdag, "matrix")
    #MEdagTrue <- pdag2allDags(cpdagMat)$dags
    #MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
    dagMat <- (MEdagTrue[,,sample(1:dim(MEdagTrue)[3], size=1)])*1
    dag <- as(dagMat, "graphNEL")
  }
  
  
  nois <- mapply(FUN=function(dist, pars, a, b){
    
    a*(do.call(dist, c(list(n=n), pars)))^b
    
    
  }, dist=lapply(nodes, function(el) el$dist), 
  pars=lapply(nodes, function(el) el$pars),
  a=lapply(nodes, function(el) el$a),
  b=lapply(nodes, function(el) el$a))
  
  
  
  colnames(nois) <- nodes(dag)
  
  
  dagAmat <- amat(as.bn(dag))
  Nodes <- nodes(dag)
  topoNodes <- topoSort(dagAmat)
  
  
  x <- matrix(NA, n, p)
  colnames(x) <- Nodes
  indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
  x[,Nodes[indxNoParents]] <- nois[,Nodes[indxNoParents]]
 
  
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # (nd <- setdiff(topoNodes, Nodes[indxNoParents])[1])
    indxParents <- which(dagAmat[,nd]==1)
    Parents <- Nodes[indxParents]
    
    K <- kernelMatrix("kern_rbf", x=x[,Parents, drop=F], 
                               y=matrix(rnorm(m*length(Parents), 0, 1), m, length(Parents)), 
                      pars=list(sigma=sigma))
    alpha <- rnorm(m, 0 , 1)    
    fx <- as.numeric(K %*% alpha)
    
    # plot(x[,Parents[2]],fx)
    
    x[,nd] <- fx + nois[,nd]
    
    
  }
  
  
  
  return(list(dag=dagAmat, x=x, nois=nois))
}

simRandAddSEMs <- function(q, ps, ns, nodess, nms=1:q, dagMat=NULL, markovEquiv=FALSE){
  
  
  res <- mapply(FUN=function(p, n, nodes){
    # i <- 1; p <- ps[i]; n <- ns[i]; nodes <- nodess[[i]]
    simRandAddSEM(p, n, nodes, dagMat=dagMat, markovEquiv=markovEquiv)
  }, p=ps, n=ns, nodes=nodess, SIMPLIFY=F)
  
  dags <- lapply(res, function(el) el$dag)
  xs <- lapply(res, function(el) el$x)
  noiss <- lapply(res, function(el) el$nois)
  
  names(dags) <- nms
  names(xs) <- nms
  names(noiss) <- nms
  
  return(list(dags=dags, xs=xs, noiss=noiss, names=nms))
  
}


# Simulate Time Series from given SEM



# simulate a directed graph
simDG <- function(p, pDiag, pOffDiag){
  ps <- c(pOffDiag, pDiag)[diag(p)+1]
  dg <- matrix(rbernoulli(p^2, ps)*1, p, p)
  return(dg)
}

# obtain time indices wrt a given window (relative indices) from cycle and lag 
getT <- function(cycleLag, sizeWindow){
  aux <- strsplit(cycleLag, split="\\.")
  cycle <- as.numeric(sapply(aux, function(el) el[[1]]))
  lag <- as.numeric(sapply(aux, function(el) el[[2]]))
  res <- sizeWindow - cycle*lag + 1
  return(res)
}

# sim time dag
simTimeDag <- function(p, C, L, pDiag, pOffDiag, instantaneous=TRUE){
  #simulate time dag
  
  # lagged relationships
  dgs <- lapply(L, function(l) lapply(1:l, function(i) simDG(p, pDiag, pOffDiag)))
  dgs <- do.call(abind, c(unlist(dgs, recursive=F), along=3))
  cycleLags <- unlist(mapply(function(cycle, lags) sapply(1:lags, function(lag) paste(cycle, lag, sep=".")), 
                             cycle=C, lags=L, SIMPLIFY=FALSE))
  dimnames(dgs) <- list(from=1:p, to=1:p, cycleLag=cycleLags)
  
  if(instantaneous){
    # instantanious relationships
    dag <- as(unifDAG(p),"matrix")
    
  } else{
    dag <- matrix(0, p, p)
  }
  
  dim(dag) <- c(dim(dag),1 )
  dimnames(dag) <- list(from=1:p, to=1:p, cycleLags="0.0")
  timeDag <- abind(dag, dgs, along=3)
  
  names(dimnames(timeDag)) <- c("from","to","cycleLag")
  
  
  return(timeDag)
}

# create according to user input
createTimeDag <- function(C, L, matList){
  
  if(length(matList)-1 != length(C) | length(C) != length(L)) stop("Error! length of cycle (C) and lag (L) vector should 
                                          be one less than number of matrices!")
  
  if(!isValidGraph(matList[[1]], type ="dag")) stop("Error! first matrix corresponds to instantaneous
                                                    relationships and must be a dag.")
  
  timeDag <- do.call(abind, c(matList, list(along=3)))
  p <- dim(timeDag)[2]
  
  cycleLags <- unlist(mapply(function(cycle, lags) sapply(1:lags, function(lag) paste(cycle, lag, sep=".")), 
                             cycle=C, lags=L))
  
  cycleLags <- c("0.0", cycleLags)
  dimnames(timeDag) <- list(from=1:p, to=1:p, cycleLag=cycleLags)
  
  return(timeDag)
}

# should simulate a function per process (using simfxH and wirtefxH)
# functions should be of the form  f <- function(x, n){ y<-stuff(x,n)} return(y)
# where x is a matrix of n x p (p the number of incoming edges in time dag) 
#       n is a matrix of n x 1 is the noise matrix
simFunsSem <- function(timeDag, pctHidVars=1, sigma, sigmaErr=0, distNs, distNsPars, geU=function(y,n) y){
  # simulate functions for timeSEM
  # obtain number of parameters for each node
  numVarss <- apply(timeDag, "to", function(mat) sum(mat)) 
  # simulate functions
  
  funcs <- lapply(numVarss, function(numVars) simfxH(numVars, pctHidVars, sigma, sigmaErr, distNs, distNsPars))
  
  funcs <- lapply(funcs, function(func) fx <- writefxH(func))
  
  names(funcs) <- dimnames(timeDag)$to
  
  return(funcs)
}

# should simulate a function per process (using simfxH and wirtefxH)
# functions should be of the form  g <- function(y, n){ y<-f(x); y <- g(y, n)} return(y)
# where y is a matrix of n x 1 (the function will be appied post-hoc g(f(x), n)) 
#       n is a matrix of n x 1 is the noise matrix
simFunsErr <- function(timeDag, pctHidVars, sigma, sigmaErr, distErr, distErrPars){
  funcs <- lapply(1:dim(timeDag)[2], function(numVars) simfxH(numVars=1, pctHidVars, sigma, sigmaErr, distErr, distErrPars))
  
  funcs <- lapply(funcs, function(func) fx <- writefxH(func))
  
  names(funcs) <- dimnames(timeDag)$to
  return(funcs)
}

# simulates n time steps of time series specified by a DAG (in the form of a timeDag)
simTimeSEM <- function(n, burnin, timeDag, funcsSem, funcsErr, distNs, distNsPars, distErr, distErrPars){
  
  p <- length(dimnames(timeDag)$from)
  aux <- strsplit(dimnames(timeDag)$cycleLag, split="\\.")
  C <- as.numeric(sapply(aux, function(el) el[[1]]))
  L <- as.numeric(sapply(aux, function(el) el[[2]]))
  
  indx <- which(C==max(C))
  sizeWindow <- max(C)*L[indx][which.max(L[indx])]
  
  distNsPars2 <- c(list(n=n), distNsPars)
  ns <- do.call(distNs, distNsPars2)
  
  # initialize sim
  var_ts <- array(NA, dim=c(n+burnin, dim(timeDag)[1]))
  dimnames(var_ts) <- c(list(t=1:(n+burnin)), dimnames(timeDag)[1])
  names(dimnames(var_ts))[2] <- "node"
  var_ts[1:sizeWindow,] <- 0
  
  # simulate sequentially
  
  t <- sizeWindow+1
  nodes <- dimnames(timeDag)$from
  topoNodes <- topoSort(timeDag[,,"0.0"])
  
  numVarss <- apply(timeDag, "to", function(mat) sum(mat)) 
  
  while(t <= (n+burnin)){
    
    if(t%%100==0) print(paste("t: ", t))  
    
    for(node in topoNodes){
      # i <- 1; node <- topoNodes[i]
      # print(paste("node: ", node))
      indxMat <- timeDag[,node,]==1
      indxMat2 <- matrix(FALSE, sizeWindow+1, p) 
      cycleLagToT <- getT(cycleLag=colnames(indxMat), sizeWindow)
      indxMat2[cycleLagToT,] <- t(indxMat)
      # sum(indxMat); sum(indxMat2); numVarss[node]
      # a <- matrix(seq(12),3, 4)
      # set.seed(2); b <- matrix(sample(c(TRUE,FALSE)),3,4))
      # a;b; a[b]
      # args are passed 1st node all times, 2nd node all times, etc times go from oldest to newest
      args <- var_ts[(t-sizeWindow):(t),][indxMat2]
      if(numVarss[node]>0){
        val <- do.call(funcsSem[[node]], list(x=matrix(args,1,length(args)),n=ns[t]))
      } else{
        val <- ns[t]
      }
      var_ts[t,node] <- val  
    }
    
    t <- t+1
    
  }
  
  # add measurement error
  
  distErrPars2 <- c(list(n=n*p), distErrPars)
  ns <- do.call(distErr, distErrPars2)
  ns <- matrix(ns, n, p) 
  
  varErr_ts <- mapply(function(x,n,f){
    ns <- do.call(distErr, distErrPars2)
    xErr <- do.call(f, list(x=matrix(x, length(x), 1),n=n))
  }, x=as.list(as.data.frame(var_ts)), n=as.list(as.data.frame(ns)), f=funcsErr, SIMPLIFY="array")
  
  #var_ts + matrix(do.call(distErr, c(list(n=(n+burnin)*p), distErrPars)), n+burnin, p)
  
  # return timeDag, summaryDAG, timeSeries, timeSeries with error
  
  return(list(timeSeries=var_ts[(burnin+1):(burnin+n),], timeSeriesErr=varErr_ts[(burnin+1):(burnin+n),]))
}

# obtains the lags from a multi-dimensional time series specified by timeDag or by (C,L(c))=(Cycles, LagsPerCycle)
getTimeDF <- function(timeSeries, timeDag, funcsSem, C=NULL, L=NULL){
  
  
  if(is.null(C) | is.null(L)){  
    
    cycleLag <- dimnames(timeDag)$cycleLag
  } else{
    cycleLag <- unlist(mapply(function(cycle, lags) sapply(1:lags, function(lag) paste(cycle, lag, sep=".")), 
                              cycle=C, lags=L))
  }
  
  aux <- strsplit(cycleLag, split="\\.")
  cycle <- as.numeric(sapply(aux, function(el) el[[1]]))
  lag <- as.numeric(sapply(aux, function(el) el[[2]]))
  
  
  indx <- which(cycle==max(cycle))
  sizeWindow <- max(cycle)*lag[indx][which.max(lag[indx])]
  
  n <- nrow(timeSeries)
  
  df <- lapply(colnames(timeSeries), function(nodeDep){
    #nodeDep <- "2"
    #print(paste("nodeDep: ", nodeDep))
    
    if(is.null(C) | is.null(L)){  
      nodeIndeps <- timeDag[,nodeDep,]
      nodeIndeps <- reshape2:::melt.array(nodeIndeps, as.is=TRUE)
      nodeIndeps <- nodeIndeps[which(nodeIndeps$value>0),]
      cycleLag <- as.character(nodeIndeps$cycleLag)
      nodeIndeps <- nodeIndeps$from
    } else{
      cycleLag <- unlist(mapply(function(cycle, lags) sapply(1:lags, function(lag) paste(cycle, lag, sep=".")), 
                                cycle=C, lags=L, SIMPLIFY=FALSE))
      nodeIndeps <- colnames(timeSeries)
      aux <- expand.grid(nodeIndep=nodeIndeps, cycleLag=cycleLag, stringsAsFactors=FALSE)
      cycleLag <- aux$cycleLag
      nodeIndeps <- aux$nodeIndep
    }
    
    aux <- strsplit(cycleLag, split="\\.")
    cycle <- as.numeric(sapply(aux, function(el) el[[1]]))
    lag <- as.numeric(sapply(aux, function(el) el[[2]]))
    
    if(length(lag)>0){
    
      df <- data.frame(nodeDep=nodeDep, depVal=timeSeries[(sizeWindow+1):n, nodeDep])
    
      dfs <- mapply(function(cyc, lg, nd){
        # i <- 1; cyc <- cycle[i]; lg <- lag[i]; nd <- nodeIndeps$from[i]
        #print(paste("cycle: ", cyc, " lag: ", lg, " node: ", nd, sep=""))
        LAG <- cyc*lg
        df <- data.frame(smTmSem$timeSeries[(sizeWindow+1-LAG):(n-LAG),nd])
        return(df)
      }, cyc=cycle, lg=lag, nd=nodeIndeps, SIMPLIFY=FALSE)
      dfs <- do.call(cbind, dfs)
      colnames(dfs) <- paste(nodeIndeps, cycle, lag, sep=".")
    
      LAG <- cycle*lag
      aux <- data.frame(as.numeric(nodeIndeps), LAG)
      indx <- order(aux[,1],aux[,2], decreasing=c(FALSE, TRUE))
      dfs <- dfs[,indx,drop=FALSE]
      #dfs$y <- do.call(funcsSem[[nodeDep]], list(x=as.matrix(dfs), n=matrix(0, nrow(dfs), 1)))
      df <- cbind(df, dfs)
      #plot(df$y, df$depVal)  
      #abline(a=0, b=1, col="red")
    
      df <- reshape2:::melt(df, id.vars=c("nodeDep","depVal"), factorsAsStrings=FALSE)
      df$variable <- as.character(df$variable)
    } else{
      df <- data.frame(nodeDep=character(), depVal=numeric(), variable=character(), value=numeric())
    }
    return(df)
  })
  df <- do.call(rbind, df)
  aux <- strsplit(df$variable, split="\\.")
  df$nodeIndep <- sapply(aux, function(el) el[[1]])
  df$cycle <- sapply(aux, function(el) el[[2]])
  df$lag <- sapply(aux, function(el) el[[3]])
  df$cycleLag <- paste(df$cycle, df$lag, sep=".")
  return(df)
}
