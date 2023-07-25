# Conditional Mean Embeddding Measure (CMEM) Hypothesis Scorer functions

library(psych) #geometric.mean
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_kernels/func_kernel_pkg.R")

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


cmem_hypScorer <- function(x, hypArray, kernelPack, measurePack, ppTab=NULL, plot=FALSE){
  
  
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
  cmemSet <- getMeasuresGivenDAGSet(uniqueRegsList, x, kernelPack, measurePack)
  
  # calculate HSIC measures for unique residuals (variables with no parents)
  
  
  # calculate CME measures for each node against all other nodes - measure of independence in same scale as above
  # for each DAG well take measures corresponding to that DAGs no parents nodes. Some will already be in cmemSet.
  # We can see this by comparing unique Regs List and eachAgainstRestList. We can fill in those that we already
  # have form cmemSet and calculate the others -> actually better yet we can make sure these regs are in uniqueRegsList
  # before and tan-tan.
  
  
  #hsicSet <- getHSICGivenDAGSet(uniqueNoParentsList, x, kernelPack)
  
  # for each hypothesis obtain conditional mean embedding measures  encoded in hypothesis
  print(paste("scoring all hypotheses"))
  
  
  
  scores <- assembleDagScores(dags=hypArray, uniqueRegsList, cmemSet, measurePack, prnt=FALSE)
  
  scores <- melt(scores)
  scores <- cast(scores, dag~measures+kernels)
  dagNms <- scores$dag
  scores <- scores[,-1]
  scoreNms <- colnames(scores)
  scores <- as.matrix(scores)
  dimnames(scores) <- list(dag=dagNms, score=scoreNms)
  
  return(scores)
}


getMeasuresGivenDAGSet <- function(uniqueRegsList, x, kernelPack, measurePack){
  
  kernelPack <- eval(parse(text=kernelPack))
  measurePack <- eval(parse(text=measurePack))
  
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
  
  kernelPerms <- expand.grid(x=1:length(kernelPack$kernelsX), y=1:length(kernelPack$kernelsY))
  
  #count <- 0
  #pm0 <- proc.time()
  measureList <- lapply(nodes,  function(nodeTo){
    # nodeTo <- "1"
    print("*********************")
    print(paste(nodeTo, " regressions: (", length(numRegs[[nodeTo]]), ", total)", sep=""))
    res <- lapply(numRegs[[nodeTo]], function(numReg){
      # numReg <- "1"
      #print("******")
      #print(paste(nodeTo,", reg # ", numReg, sep=""))
      indxPreds <- which(uniqueRegsList[[nodeTo]][,numReg]==1)
      
      #count <- count + 1
      
      dataX <- as.matrix(x[,nodes[indxPreds]])
      dataY <- as.matrix(x[,nodeTo])
      
      cmems <- sapply(measurePack, function(msr) apply(kernelPerms, 1, function(perm){
        # i <- 2; j <- 1; msr <- measurePack[[i]]; perm <- as.matrix(kernelPerms)[j,]
        #print(paste("msr: ", msr$func, ", kernelX: ", kernelPack$kernelsX[[perm["x"]]]$func, ", kernelY: ", kernelPack$kernelsY[[perm["y"]]]$func))
        parsMsr <- list(pars=msr$pars)
        parsMsr$x <- dataX
        parsMsr$y <- dataY
        parsMsr$kernelX <- kernelPack$kernelsX[[perm["x"]]]
        parsMsr$kernelY <- kernelPack$kernelsY[[perm["y"]]]
        # x <- parsMsr$x; y <- parsMsr$y; kernelX <- parsMsr$kernelX; kernelY <- parsMsr$kernelY; pars <- parsMsr$pars
        mesr <- do.call(msr$func, parsMsr)
        
        
        return(mesr)
      }), simplify="array") 
      
      rownames(cmems) <- paste(names(kernelPack$kernelsX)[kernelPerms$x], names(kernelPack$kernelsY)[kernelPerms$y], sep="_") 
      
      
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

getHSICGivenDAGSet <- function(uniqueNoParentsList, x, kernelPack){
  
  kernelPack <- eval(parse(text=kernelPack))
  
  n <- nrow(x)
  
  hsics <- sapply(uniqueNoParentsList, function(noParents){
    # i <- 1; noParents <- uniqueNoParentsList[[i]]
    data <- x[,noParents]
    hsic <- dhsic.test(data, matrix.input=TRUE, pairwise=FALSE, method="gamma")$statistic
    
    set.seed(12345)            
    rperms <- sapply(1:kernelPack$numPerms, function(i) sapply(1:length(noParents), function(j) sample(n), simplify="array"), simplify="array")
    rhsic <- mean(apply(rperms, 3, function(rperm){
      # k <- 1; rperm <- rperms[,,k]
      data <- mapply(function(col, m) x[col,noParents[m]], col=as.data.frame(rperm), m=1:length(noParents), SIMPLIFY="array")
      rhsic <- dhsic.test(data, matrix.input=TRUE, pairwise=FALSE, method="gamma")$statistic
      return(rhsic)
    }))
    q_hsic <- hsic/rhsic
    
    return(q_hsic)
  })
  
  return(hsics)
}

assembleDagScores <- function(dags, uniqueRegsList, cmemSet, measurePack, prnt=FALSE){
  
  measurePack <- eval(parse(text=measurePack))
  
  nodes <- dimnames(dags)[[1]]
  
  numDags <- dim(dags)[3]
  count <- 0
  pm0 <- proc.time()
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  names(dagsList) <- dimnames(dags)$dag
  scores <- sapply(dagsList, function(dag){
    #get residuals for the four regressions implied in each column
    # dag <- dags[,,1]
    count <<- count + 1
    if(prnt){
      print("*******************")
      print(paste("dag # ",count, sep=""))
    }
    
    # obtain residuals/vars to evaluate for each dag
    facs <- sapply(nodes, function(nodeTo){
      # (nodeTo <- nodes[1])
      #print(paste("nodeTo: ", nodeTo))
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
    
    names(dimnames(facs)) <- c("kernels","measures","nodeFactors")
    
    aggFuncs <- sapply(measurePack, function(el) el$aggFunc)
    
    score <- sapply(1:length(aggFuncs), function(i){
      # i <- 1
      apply(facs[,i,], 1, aggFuncs[i], na.rm=T)
    })
    
    colnames(score) <- dimnames(facs)$measures
    names(dimnames(score)) <- names(dimnames(facs))[1:2]  
    
    
    if(prnt){
      print("estimated time until completion")
      print((proc.time()-pm0)*(numDags-count)/count)
    }
    return(score)
  }, simplify="array")
  if(prnt) print(proc.time()-pm0) #
  
  
  names(dimnames(scores))[3] <- "dag"
  
  
  scores <- aperm(scores, c(2,1,3))
  
  return(scores)
}


# Get unique (index) list of nodes with no parents in a set of dags
getUniqueNoParentsList <- function(dags){
  dagsList <- lapply(1:dim(dags)[3], function(i) dags[,,i])
  noParents <- lapply(dagsList, function(mat) which(apply(mat, 2, function(col) all(col==0))))
  uniqueNoParents <- unique(noParents)
  return(uniqueNoParents)
}


# set parameter functions
medSig <- function(x){
  sigma0 <- 1/median(as.numeric(dist(x)^2))
  return(list(sigma=sigma0))
}
deg1 <- function(x) return(list(degree=1))
off1 <- function(x) return(list(offset=1))
noPars <- function(x) return(list())

# measure functions: shd take in x, y, kernelX, kernelY (the lists from kernel Pack)
# and pars (the list from measurePack)
# x- alleged causes
# y- alleged effects
KCDC <- function(x, y, kernelX, kernelY, pars){
  
  parsX <- do.call(kernelX$setParams, list(x=x)) 
  parsY <- do.call(kernelY$setParams, list(x=y))
  L  <- kernelMatrix(kernelX$func, x, pars=parsX) 
  K  <- kernelMatrix(kernelY$func, y, pars=parsY)
  
  lambda <- pars$lambda
  
  n <- nrow(x)
  I <- diag(n)
  
  Blambda <- solve(L+n*lambda*I)
  Alambda <- Blambda%*%K%*%Blambda
  
  LAL <- L%*%Alambda%*%L
  
  b <- sum(diag(LAL))
  c <- sum(diag(LAL^2))
  res <- (c/n) - (b/n)^2
  
  return(res)
}

KCDCrel <- function(x, y, kernelX, kernelY, pars){
  
  n <- nrow(x)
  set.seed(12345)            
  rperms <- sapply(1:pars$numPerms, function(i) sample(n))
  
  pars <- pars[-which(names(pars)=="numPerms")]
  
  mesr <- KCDC(x, y, kernelX, kernelY, pars)
  
  rmesr <- mean(apply(rperms, 2, function(rperm){
    # k <- 1; rperm <- rperms[,k]
    
    res <- KCDC(x, as.matrix(y[rperm,]), kernelX, kernelY, pars) 
    return(res)
  }))
  qmesr <- mesr/rmesr
  
  
  return(qmesr)
}

HSIC_cmem <- function(x, y, kernelX, kernelY, pars){
  
  parsX <- do.call(kernelX$setParams, list(x=x)) 
  parsY <- do.call(kernelY$setParams, list(x=y))
  L  <- kernelMatrix(kernelX$func, x, pars=parsX) 
  K  <- kernelMatrix(kernelY$func, y, pars=parsY)
  
  Ks <- vector("list", 2)
  Ks[[1]] <- K
  Ks[[2]] <- L

  res <- dhsic(K=Ks)$dHSIC
  
  return(res)
}

