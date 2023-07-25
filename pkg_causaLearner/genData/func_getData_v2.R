# obtain data on which to test causal learner package

library(pcalg) # this gives problem installing because of V8 dependency  
# install.packages("V8", configure.vars='INCLUDE_DIR=/usr/include LIB_DIR=/usr/lib64') (including path to v8 in server did the job)

library(unifDAG) #unifDAG
library(gRbase) # topoSort
library(purrr) # rbernoulli
library(bnlearn) # amat in (simRandAddSEM)
library(graph) # nodes in (simRandAddSEM)
library(R.matlab) # readMat
library(Compositional) # mkde
library(pdfCluster) #kepdf
library(vegan)
source("./pkg_kernels/func_kernel_pkg.R")
source("./pkg_causaLearner/utilities/func_dagStuff.R")

# get Tubingen Cause Effect Pairs

# dont read data just keep file name for later
createTCEPList1 <- function(pairs, meta, folder){
  
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
  
  ns <- xyCols
  res <- list(dags=dags, xs=xs, noiss=ns, names=as.character(pairs))
  return(res)
}

# read in data, sample n fixed - one direction x->y
createTCEPList2 <- function(pairs, meta, folder, n){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[1]]
    #print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    set.seed(1234)
    x <- x[sample(1:nrow(x), min(nrow(x), n)),1:2]
    return(x)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_cluster <- function(pairs, meta, folder){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[6]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    #plot(x[,1], x[,2])
    
   
    clust <- kmeans(x, centers=10)
    probs <- table(clust$cluster)
    probs <- 1/probs  
    
    probs <- probs[match(clust$cluster, names(probs))]
    probs <- probs/sum(probs)
    
    set.seed(1234)
    xnew <- x[sample(1:nrow(x), size=min(nrow(x), 100), prob=probs, replace=T),]
    
    #plot(xnew[,1], xnew[,2])
    
    return(xnew)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_density <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    #plot(x[,1], x[,2])
    mod <- mkde(x)
    
    
    if(plt){
      x2 <- cbind(x, mod)
      p <- ggplot(as.data.frame(x2))
      p <- p + geom_point(aes(x=x, y=y, colour=mod))
      p
    }
    probs <- 1/mod  
    probs <- probs/sum(probs)
    
    if(plt){
      x3 <- cbind(x2, probs)
      p <- ggplot(as.data.frame(x3))
      p <- p + geom_point(aes(x=x, y=y, colour=probs))
      p
    }
    
    
    set.seed(1234)
    xnew <- x[sample(1:nrow(x), size=n, prob=probs, replace=T),]
    #plot(x[,1],x[,2])
    #plot(xnew[,1], xnew[,2])
    #hist(xnew[,1])
    #hist(xnew[,2])
    
    if(plt){ 
      p <- ggplot(as.data.frame(xnew))
      p <- p + geom_point(aes(x=x, y=y))
      p
    }
    
    #plot(xnew[,1], xnew[,2])
    
    return(xnew)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_denseSquare <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    #plot(x[,1], x[,2])
    mod <- mkde(x)
    indxMax <- which.max(mod)
    #lines(x[indxMax,1],x[indxMax,2], cex=2, col="red", type="p")
    
    
    indices <- knn.index(x, k=100)
    indices <- indices[indxMax,]
    #lines(x[indices,1], x[indices,2],col="green")
    #plot(x[indices,1], x[indices,2])
    
    set.seed(1234)
    xnew <- x[indices,]
    #plot(x[,1],x[,2])
    #plot(xnew[,1], xnew[,2])
    #hist(xnew[,1])
    #hist(xnew[,2])
    

    return(xnew)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

# read in data, sample n var - one direction x->y
createTCEPList_all <- function(pairs, meta, folder){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    set.seed(1234)
    x <- x[sample(1:nrow(x), nrow(x)),1:2]
    return(x)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_copula <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    #plot(x[,1], x[,2])
    xnew <- apply(x, 2, function(col){
      o <- order(col)
      xnew <- (1:length(col))/length(col)
      xnew <- xnew[order(o)]
      return(xnew)
    })
    # plot(xnew[,1], xnew[,2])
    # hist(xnew[,1])
    # hist(xnew[,2])
    
    set.seed(1234)
    xnew <- xnew[sample(1:nrow(xnew), min(nrow(xnew), n)),1:2]
    
    return(xnew)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_copulaCause <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    #plot(x[,1], x[,2])
    xnew <- x
    o <- order(xnew[,1])
    aux <- (1:length(o))/length(o)
    aux <- aux[order(o)]
    xnew[,1] <- aux
    # plot(xnew[,1], xnew[,2])
    # hist(xnew[,1])
    # hist(xnew[,2])
    
    set.seed(1234)
    xnew <- xnew[sample(1:nrow(xnew), min(nrow(xnew), n)),1:2]
    
    return(xnew)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  names(xs) <- nms
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_densityCauseIter <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    
    return(x)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  
  count <- 0
  xs <- lapply(xs, function(x){
    # x <- xs[[65]]
    count <<- count + 1
    print(paste("count: ", count))
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    xnew <- x
    #plot(x[,1], x[,2])
    #hist(xnew[,1])
    for(i in 1:10){
      print(paste("i: ", i))
      mod <- kepdf(xnew[,1], eval.points = xnew[,1], kernel = "gaussian", bwtype = "adaptive")@estimate
      
      probs <- 1/mod  
      probs <- probs/sum(probs)
      
      xnew <- xnew[sample(1:nrow(xnew), size=10000, prob=probs, replace=T),]
      #plot(x[,1],x[,2])
      #plot(xnew[,1], xnew[,2])
      #hist(x[,1])
      #hist(xnew[,1])
      #pval <- ks.test(xnew[,1], y="punif", min=min(xnew[,1]), max=max(xnew[,1]))$p.value
      #print(paste("pval: ", pval))
    }
    
    
    set.seed(1234)
    xnew <- xnew[sample(1:nrow(xnew), size=n),]
    
    #plot(x[,1],x[,2])
    #plot(xnew[,1], xnew[,2])
    #hist(x[,1])
    #hist(xnew[,1])
    
    return(xnew)
  })
  
  names(xs) <- nms
  
  
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createTCEPList_densityEffectIter <- function(pairs, meta, folder, n=100){
  
  nms <- as.character(pairs)
  
  fileFolders <- lapply(pairs, function(pair){
    fileName <- paste("pair", formatC(pair, width=4, flag=0), ".txt", sep="")
    fileFolder <- paste(folder, fileName, sep="")
    return(fileFolder)
  })
  names(fileFolders) <- nms
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  plt <- FALSE
  xs <- lapply(fileFolders, function(file){
    # file <- fileFolders[[65]]
    print(file)
    x <- read.csv(file, sep="", header=F)
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    
    return(x)
  })
  
  xyCols <- lapply(pairs, function(pair){
    # pair <- 1
    indx <- which(meta$pairNumber==pair)  
    res <- c(meta[indx,"firstCauseCol"], meta[indx,"firstEffectCol"])
    return(res)
  })
  names(xyCols) <- nms
  
  xs <- mapply(FUN=function(x, indx){
    # x <- xs[[1]]; indx <- xyCols[[1]]
    res <- x[,indx]
    colnames(res) <- colnames(dag)  
    return(res)
  }, x=xs, indx=xyCols, SIMPLIFY=FALSE)
  
  count <- 0
  xs <- lapply(xs, function(x){
    # x <- xs[[65]]
    count <<- count + 1
    print(paste("count: ", count))
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
    xnew <- x
    #plot(x[,1], x[,2])
    #hist(xnew[,1])
    for(i in 1:10){
      print(paste("i: ", i))
      mod <- kepdf(xnew[,2], eval.points = xnew[,2], kernel = "gaussian", bwtype = "adaptive")@estimate
      
      probs <- 1/mod  
      probs <- probs/sum(probs)
      
      xnew <- xnew[sample(1:nrow(xnew), size=10000, prob=probs, replace=T),]
      #plot(x[,1],x[,2])
      #plot(xnew[,1], xnew[,2])
      #hist(x[,1])
      #hist(xnew[,1])
      #pval <- ks.test(xnew[,1], y="punif", min=min(xnew[,1]), max=max(xnew[,1]))$p.value
      #print(paste("pval: ", pval))
    }
    
    
    set.seed(1234)
    xnew <- xnew[sample(1:nrow(xnew), size=n),]
    
    #plot(x[,1],x[,2])
    #plot(xnew[,1], xnew[,2])
    #hist(x[,1])
    #hist(xnew[,1])
    
    return(xnew)
  })
  
  names(xs) <- nms
  
  
  
  dags <- lapply(fileFolders, function(el) dag) 
  names(dags) <- nms
  
  ns <- rep(NA, length(pairs))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}


createEmulatedList <- function(dataRepos, n, oneDir=FALSE){
  
  
  files <- dir(dataRepos)
  
  aux <- strsplit(files, "_")
  #table(sapply(aux, function(el) el[1])) # p1 or p2
  #table(sapply(aux, function(el) el[3])) # 1-13 band
  #table(sapply(aux, function(el) el[4])) # 1-7 biopar
  
  cause <- sapply(aux, function(el) as.numeric(substr(el[1],2,2)))
  band <- sapply(aux, function(el) as.numeric(substr(el[3],1,nchar(el[3])-6)))
  biopar <- sapply(aux, function(el) as.numeric(substr(el[4],1,1)))
  
  pm <- proc.time()
  dat <- lapply(files, function(file) readMat(paste(dataRepos, file, sep="")))
  proc.time() - pm # 3 mins
 
  # bnd <- 2
  # bpr <- 4 
  # (indx <- which(band==bnd & biopar==bpr))
  # smpl <- sample(50000, 100)
  # plot(dat[[indx[1]]]$x[smpl,], dat[[indx[2]]]$y[smpl,])
  # plot(dat[[indx[1]]]$y[smpl,], dat[[indx[2]]]$x[smpl,])
  
  # length(dat) # 182 data sets
  # names(dat)
  # names(dat[[1]])
  # table(sapply(dat, function(el) el$Y))
  # table(sapply(dat, function(el) nrow(el$x))) #500,000 obs per dataset
  # table(sapply(dat, function(el) ncol(el$x))) #1d dataset
  # table(sapply(dat, function(el) nrow(el$y))) #500,000 obs per dataset
  # table(sapply(dat, function(el) ncol(el$y))) #1d dataset
  # i <- 92
  # plot(dat[[i]]$x[indx,], dat[[i]]$y[indx,])
  # plot(dat[[i]]$y[indx,], dat[[i]]$x[indx,])
  
  
  nms <- paste(biopar, band, cause, sep="_")
  dag1 <- matrix(c(0, 0, 1, 0), 2, 2)
  dag2 <- matrix(c(0, 1, 0, 0), 2, 2)
  colnames(dag1) <- rownames(dag1) <- colnames(dag2) <- rownames(dag2) <- c("x", "y")
  dags <- abind(dag1, dag2, along=3)
  dags <- lapply(cause, function(c) dags[,,c]) 
  if(oneDir) dags <- lapply(1:length(cause), function(el) dag1)
  names(dags) <- nms
  
  
  xs <- mcmapply(function(dt,caus){
    # i <- 1; dt <- dat[[i]]; caus <- cause[i]
    indx <- sample(1:nrow(dt$x), n)
    res <- cbind(dt$x[indx,], dt$y[indx,])
    if(oneDir) res <- res[,c(caus, setdiff(1:2, caus))]
    res <- as.matrix(res)
    colnames(res) <- c("x","y")
    return(res)
  }, dt=dat, caus=cause, SIMPLIFY=FALSE)
  names(xs) <- nms
  ns <- rep(NA, length(cause))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(res)
}

createProsailList <- function(dataRepos, n){
  
  
  files <- dir(dataRepos)
  
  aux <- strsplit(files, "_")
  #table(sapply(aux, function(el) el[1])) # p1 or p2
  #table(sapply(aux, function(el) el[3])) # 1-13 band
  
  cause <- sapply(aux, function(el) as.numeric(substr(el[1],2,2)))
  band <- sapply(aux, function(el) as.numeric(substr(el[3],1,1)))
  
  pm <- proc.time()
  dat <- lapply(files, function(file) readMat(paste(dataRepos, file, sep="")))
  proc.time() - pm # 0.2 secs
  
  # bnd <- 1
  # (indx <- which(band==bnd))
  # plot(dat[[indx[1]]]$x, dat[[indx[2]]]$y)
  # plot(dat[[indx[1]]]$y, dat[[indx[2]]]$x)
  
  # length(dat) # 12 data sets
  # names(dat)
  # names(dat[[1]])
  # table(sapply(dat, function(el) el$Y))
  # table(sapply(dat, function(el) nrow(el$x))) #1000 obs per dataset
  # table(sapply(dat, function(el) ncol(el$x))) #1d dataset
  # table(sapply(dat, function(el) nrow(el$y))) #1000 obs per dataset
  # table(sapply(dat, function(el) ncol(el$y))) #1d dataset
  # i <- 10
  # plot(dat[[i]]$x, dat[[i]]$y)
  # plot(dat[[i]]$y, dat[[i]]$x)
  
  
  nms <- paste(band, cause, sep="_")
  dag1 <- matrix(c(0, 0, 1, 0), 2, 2)
  dag2 <- matrix(c(0, 1, 0, 0), 2, 2)
  colnames(dag1) <- rownames(dag1) <- colnames(dag2) <- rownames(dag2) <- c("x", "y")
  dags <- abind(dag1, dag2, along=3)
  dags <- lapply(cause, function(c) dags[,,c]) 
  names(dags) <- nms
  
  
  xs <- lapply(dat, function(el){
    indx <- sample(1:nrow(el$x), n)
    res <- cbind(el$x[indx,], el$y[indx,])
    res <- as.matrix(res)
    colnames(res) <- c("x","y")
    return(res)
  })
  names(xs) <- nms
  ns <- rep(NA, length(cause))
  names(ns) <- nms
  
  res <- list(dags=dags, xs=xs, noiss=ns, names=nms)
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

# ANLSMN data from bQCD paper

readQNLSM <- function(folder, nm){
  # folder <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/data/ANLSMN/LS-s/"
  file_gt <- "pairs_gt.txt"
  resp <- read.csv(paste(folder, file_gt, sep=""), header=F)
  resp <- resp[,1]
  
  filePairs <- dir(folder)
  filePairs <- filePairs[-grep("gt", filePairs)]
  num <- as.numeric(sapply(strsplit(sapply(strsplit(filePairs,"_"), function(el) el[2]), "\\."), function(el) el[1]))
  filePairs <- filePairs[order(num)]
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- c("x", "y")
  rownames(dag) <- c("x", "y")
  
  xs <- mcmapply(function(file, rsp){
    #i <- 50; file <- filePairs[i]; rsp <- resp[i]
    #print(file)
    
    x <- as.matrix(read.csv(paste(folder, file, sep=""), row.names=1))
    dir <- (!rsp)*1
    dir <- dir + 1
    dir <- c(dir, setdiff(c(1,2), dir))
    x <- x[,dir]
    colnames(x) <- colnames(dag)
    return(x)
  }, file=filePairs, rsp=resp, mc.cores=1, SIMPLIFY=FALSE)  
  
  dags <- lapply(filePairs, function(el) dag) 
  ns <- lapply(1:length(xs), function(i) NA)
  nms <- paste(nm, 1:length(filePairs), sep=".")
  names(xs) <- names(ns) <- names(dags) <- nms
  dataList <- list(dags=dags, xs=xs, noiss=ns, names=nms)
  return(dataList)
}
# dataList <- readQNLSM(folder, nm="LS-s")

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
simfxH <- function(numVars, pctHidVars=1, sigma, sigmaErr=0, dist, distPars, geU=function(y, nois, scale, const) y, seed=NULL){
  
  m <- 100 #number of inducing points
  n <- 100 #number of points on which to evaluate to get approx. constant and scale
  numHiddenVars <- floor(pctHidVars*numVars)
  numTotVars <- numVars + numHiddenVars
  set.seed(seed)
  # inducing points - non- hidden
  z <- matrix(rnorm(m*(numVars)), m, numVars)  
  # dummy points to get an idea of bias and variance of simulated points
  # we will evaluate fz() on w to obtain y
  
  #w <- matrix(runif(n*(numVars)), n, numVars)
  w <- matrix(rnorm(n*(numVars)), n, numVars)
  
  distPars2 <- c(list(n=m*numHiddenVars), distPars)
  
  # inducing points - hidden/latent
  zh <- matrix(do.call(dist, distPars2), m, numHiddenVars)
  wh <- matrix(apply(zh, 2, mean), n, numHiddenVars, byrow=T)
  
  
  
  alpha <- rnorm(m, 0 , 1) 
  
  Kx <-  kern_rbf2(w, wh,z, zh, sigma, sigmaErr)
  y <- Kx%*%alpha
  
  distPars2 <- c(list(n=n), distPars)
  nois <- do.call(dist, distPars2)
  y <- do.call(geU, list(y=y, nois=nois, scale=1, const=0))
  return(list(z=z, zh=zh, alpha=alpha, const=min(y), scale=(max(y)-min(y)), 
              numHiddenVars=numHiddenVars, sigma=sigma, 
              sigmaErr=sigmaErr, dist=dist, distPars=distPars, geU=geU))  
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
  #print("summary(Kx)")
  #print(summary(as.numeric(Kx)))
  #print("summary(alpha)")
  #print(summary(as.numeric(alpha)))
  #print("summary(res)")
  #print(summary(res))
  #print(paste("min train:", simfx$const, ", min test:", min(res)))
  
  distPars2 <- c(distPars, n=n)
  no <- do.call(dist,distPars2)
  res <- do.call(geU, list(y=res, const=const, scale=scale, nois=no)) 
  return(list(res=res,nois=no))
}

# writes f(z,alpha,x,n) so that z, alpha are written internally as parameters and we may have a 
# function of the type f(x,n) for use in a SEM.
# Basically this function allows us to simulate a dag and then simulate randomly functions (possibly non-linear and non-additive) for this
# dag with simfxH in the form of (z, zh, alpha, const, scale, numHiddenVars, sigma, sigmaErr, dist, distPars, geU) and then writes all that 
# internally in an R function. 
writefxH <- function(simfx){
  
  nmsDistPars <- names(distPars)
  valsDistPars <- unlist(distPars)
  strDistPars <- paste(paste(nmsDistPars, valsDistPars, sep="="), collapse=", ")
  
  
  
  fx <- eval(parse(text=
                     paste("fx <- function(x,nois){z <- matrix(c(", 
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
                           "res <- do.call(geU, list(y=res, const=const, scale=scale, nois=no)); ", 
                           "return(res)}", sep="")))
  
  return(fx)
}


simRandSEM <- function(p, n, nodes, nois=NULL, sigma, sigmaErr=0, dagMat=NULL, markovEquiv=FALSE, geU=function(y, nois, scale, constant) y + nois, seed=NULL){
  # p <- 4; n <- 100 
  # nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
  # nodes <- rep(list(nodes),4)
  
  m <- 100 #number of inducing points
  
  if(!is.null(dagMat)){
    dag <- as(dagMat, "graphNEL")
  } else{
    dag <- unifDAG(p)
    dagMat <- amat(as.bn(dag))
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
  
  
  if(is.null(nois)){
    nois <- mapply(FUN=function(dist, pars, a, b){
    a*(do.call(dist, c(list(n=n), pars)))^b
    }, dist=lapply(nodes, function(el) el$dist), 
      pars=lapply(nodes, function(el) el$pars),
      a=lapply(nodes, function(el) el$a),
      b=lapply(nodes, function(el) el$b))
      colnames(nois) <- nodes(dag)
  }
  
  dagAmat <- amat(as.bn(dag))
  Nodes <- nodes(dag)
  topoNodes <- topoSort(dagAmat)
  
  
  x <- matrix(NA, n, p)
  colnames(x) <- Nodes
  indxNoParents <- which(apply(dagAmat,2, function(col) all(col==0)))
  x[,Nodes[indxNoParents]] <- nois[,Nodes[indxNoParents]]
  
  
  
  for(nd in setdiff(topoNodes, Nodes[indxNoParents])){
    # (nd <- setdiff(topoNodes, Nodes[indxNoParents])[1])
    #print(paste("var: ", nd))
    indxParents <- which(dagAmat[,nd]==1)
    Parents <- Nodes[indxParents]
    
    simf <- simfxH(numVars=length(Parents), pctHidVars=1, 
                   sigma=sigma, sigmaErr=sigmaErr, dist=nodes[[nd]]$dist, 
                   distPars=nodes[[nd]]$pars, geU=geU, seed=seed)
    fx <- applySimfxH(x=x[,Parents, drop=F], simf)
    
    x[,nd] <- fx$res
    nois[,nd] <- fx$nois  
    #K <- kernelMatrix("kern_rbf", x=x[,Parents, drop=F], 
    #                           y=matrix(rnorm(m*length(Parents), 0, 1), m, length(Parents)), 
    #                  pars=list(sigma=sigma))
    #alpha <- rnorm(m, 0 , 1)    
    #fx <- as.numeric(K %*% alpha)
    
    # plot(x[,Parents[2]],fx)
    
    #x[,nd] <- fx + nois[,nd]
    
    
  }
  
  
  
  return(list(dag=dagAmat, x=x, nois=nois))
}

simRandSEMs <- function(q, ps, ns, nodess, noiss=NULL, sigma, sigmaErr=0, nms=as.character(1:q), dagMat=NULL, markovEquiv=FALSE, geU=function(y, nois, scale, constant) y + nois, seeds=NULL){
  
  if(is.null(noiss)) noiss <- lapply(1:q, function(el) NULL)
  if(is.null(seeds)) seeds <- lapply(1:q, function(el) NULL)
  
  res <- mapply(FUN=function(p, n, nodes, nois, seed){
    # i <- 1; p <- ps[i]; n <- ns[i]; nodes <- nodess[[i]]; nois <- noiss[[i]]
    res <- simRandSEM(p, n, nodes, nois=nois, sigma=sigma, sigmaErr=sigmaErr, dagMat=dagMat, markovEquiv=markovEquiv, geU=geU, seed=seed)
    #res$dag
    #plot(res$x[,9], res$nois[,9]); abline(a=0, b=1, col="red")
    #summary(res$x[,9]-res$nois[,9])
    
    return(res)
  }, p=ps, n=ns, nodes=nodess, nois=noiss, seed=seeds, SIMPLIFY=F)
  
  dags <- lapply(res, function(el) el$dag)
  xs <- lapply(res, function(el) el$x)
  noiss <- lapply(res, function(el) el$nois)
  
  names(dags) <- nms
  names(xs) <- nms
  names(noiss) <- nms
  
  return(list(dags=dags, xs=xs, noiss=noiss, names=nms))
  
}

dataJoin <- function(dataListList){
  dags <- lapply(dataListList, function(el) el$dags)
  dags <- do.call(c, dags)
  xs <- lapply(dataListList, function(el) el$xs)
  xs <- do.call(c, xs)
  noiss <- lapply(dataListList, function(el) el$nois)
  noiss <- do.call(c, noiss)
  names <- lapply(dataListList, function(el) el$names)
  names <- lapply(1:length(names), function(i) paste(names(names)[i],names[[i]],sep="."))
  names <- do.call(c, names)
  names(names) <- NULL
  names(xs) <- names(noiss) <- names(dags) <- names
  res <- list(dags=dags, xs=xs, noiss=noiss, names=names)
  return(res)
}

# Simulate "methods and benchmarks(Mooij et al 2016)"-style  GP data

rGP_rbf <- function(x, tau, xInd=NULL, inducing=FALSE, sigma=NULL, minDR=NULL, maxDR=NULL, seed=NULL){
  #if((ncol(x) !=length(sigma) & length(sigma)!=1) & (ncol(x) !=length(cutoff) & length(cutoff)!=1)) error("one sigma per column or one sigma for all columns!")
  #if(is.null(sigma) & is.null(cutoff)) error("one of sigma and error shd be provided")
  n <- nrow(x)
  if(inducing){
    x2 <- xInd
    removeDiag <- FALSE
  } else{
    x2 <- x
    removeDiag <- TRUE
  }
  
  if(any(is.na(sigma))){
    
    #print("minDR"); print(minDR)
    #print("maxDR"); print(maxDR)
    sigmaAux <- sapply(which(is.na(sigma)), function(i){
      # i <- 1
      # fit sigma
      # print(paste("i: ", i))
      ord <- 10
      sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F], sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
        res <- var(distsX)
        return(res)
      })
      indxMaxVar <- which.max(varsHsics.sig)
      sigmaVar <- sigmas1[indxMaxVar]
      
      # obtain dynamic range of sigma
      spl <- spline(log(sigmas1,10), varsHsics.sig)
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
      
      sigmas1 <- (10^seq(DR.sig[1],DR.sig[2], length.out=20)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F],sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
        res <- var(distsX)
        return(res)
      })
      indxMaxVar <- which.max(varsHsics.sig)
      sigmaVar <- sigmas1[indxMaxVar]
      
      # obtain dynamic range of sigma
      spl <- spline(log(sigmas1,10), varsHsics.sig, n=100*length(sigmas1))
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
      
      #minDR <- 0.2
      #maxDR <- 0.4
      
      indx <- which(spl$x>=DR.sig[1] & spl$x<DR.sig[2])
      xx <- spl$x[indx]
      yy <- spl$y [indx]
      indxMax <- which.max(yy)
      xx2 <- rep(NA, length(xx))
      xx2[1:indxMax] <- yy[1:indxMax]/(2*yy[indxMax])
      xx2[(indxMax+1):length(xx)] <- 0.5+(0.5-yy[(indxMax+1):length(xx)]/(2*max(yy[(indxMax+1):length(xx)])))
      #plot(xx2, xx)
      
      splf <- splinefun(c(0,xx2,1), c(xx[1],xx,xx[length(xx)]), method="monoH.FC")
      print(c(minDR[i], maxDR[i]))
      rngSigs <- splf(c(minDR[i], maxDR[i]))
      if(abs(rngSigs[2]-rngSigs[1])<1e-3){
        rngSigs[1] <- mean(rngSigs)-1e-3
        rngSigs[2] <- mean(rngSigs)+1e-3
      }
      
      res <- runif(1, min=rngSigs[1], max=rngSigs[2])
      res <- 10^res
      
      if(FALSE){
        plot(log(sigmas1,10), varsHsics.sig, col="red", type="p")
        lines(spl)
        abline(v=DR.sig, col="purple")
        abline(v=log(sigmaVar,10), col="blue")
        abline(h=max(varsHsics.sig), col="green")
        abline(v=rngSigs, col="orange", type="p")
      }
      
      return(res)
    })
    sigma[which(is.na(sigma))] <- sigmaAux 
  }
  
  
  sigma2 <- matrix(0, ncol(x), ncol(x))
  diag(sigma2) <- 1/sigma
  x2 <- x %*% sigma2
  x2 <- x2/sqrt(2)
  
  
  if(inducing){
    
    xInd2 <- xInd %*% sigma2
    xInd2 <- xInd2/sqrt(2)
    K <- kern_rbf(x2, xInd2, sigma=1)
    #alternative K calculation
    if(FALSE){
      
      K_alt <- sapply(1:length(sigma), function(i) kern_rbf(x=x[,i,drop=F], y=xInd[,i,drop=F], sigma=1/(2*sigma[i]^2)), simplify="array")
      K_alt <- apply(K_alt, c(1,2), prod)
      smpl1 <- sample(length(K), 100)
      plot(c(K)[smpl1], c(K_alt)[smpl1]); abline(a=0, b=1, col="red")
    }
    svdK <- svd(K) 
    set.seed(seed)
    y <- rnorm(nrow(xInd), sd=sqrt(svdK$d)+tau)
    set.seed(NULL)
    y <- matrix(y, length(y), 1)
    if(ncol(svdK$u) == nrow(y)){
      y <- svdK$u %*%  y
    } else{
      y <- t(svdK$v) %*% y 
    }
    
  } else{
    K <- kern_rbf(x2, sigma=1) 
    n <- nrow(K)
    I <- diag(n)
    C <- K + tau^2*I
    if(!is.positive.definite(C)) C <- make.positive.definite((C))
    C  <- chol(C)
    set.seed(seed)
    y <- rnorm(nrow(x))
    set.seed(NULL)
    y <- t(C) %*% matrix(y,n,1)
    
  }
  
  # we want cov(y, y') = K(y, y') 
  # say y=t(C) %*% x and y'= t(C') %*% x', x,x'~ N(0,I) 
  # then cov(y, y') = cov(t(C)%*%x, t(C')%*%x')
  #                 = 
  
  
  return(y)  
}

RDdist_GP <- function(n, tau, sigma=NA, minDR=NULL, maxDR=NULL, seedDist=NULL, seedFun=NULL){
  # n <- 100
  # sigma  <- rgamma(1, shape=5, scale=0.1)
  # tau <- 1e-4
  set.seed(seedDist)
  x <- matrix(rnorm(n))
  smpl <-sample(nrow(x))
  set.seed(NULL)
  F <- rGP_rbf(x, tau, sigma=sigma, minDR=minDR, maxDR=maxDR, seed=seedFun)
  eF <- exp(F)
  o <- order(x)
  G <- cumtrapz(x[o], eF[o])
  G <- matrix(G[smpl])
  return(G)
}

SIMdist_GP_sigPar <- function(n,  tau, sig_RDx=NA, sig_RDy=NA, sig_Ex=NA, sig_Ey=NA, sig_x=NA,  sig_RDz=NA, sig_Ez=NA, 
                              min_RDx=NULL, min_RDy=NULL, min_Ex=NULL, min_Ey=NULL, min_x=NULL, min_RDz=NULL, min_Ez=NULL,
                              max_RDx=NULL, max_RDy=NULL, max_Ex=NULL, max_Ey=NULL, max_x=NULL, max_RDz=NULL, max_Ez=NULL,
                              inducing=FALSE, numInducing=n,
                              sig_RDxInd=sig_RDx, sig_RDyInd=sig_RDy, sig_ExInd=sig_Ex, sig_EyInd=sig_Ey, sig_xInd=sig_x, sig_RDzInd=sig_RDz, sig_EzInd=sig_Ez, 
                              min_RDxInd=min_RDx, min_RDyInd=min_RDy, min_ExInd=min_Ex, min_EyInd=min_Ey, min_xInd=min_x, min_RDzInd=min_RDz, min_EzInd=min_Ez,
                              max_RDxInd=max_RDx, max_RDyInd=max_RDy, max_ExInd=max_Ex, max_EyInd=max_Ey, max_xInd=max_x, max_RDzInd=max_RDz, max_EzInd=max_Ez,
                              sig_nois_x=NULL, sig_nois_y=NULL,
                              addNoise=TRUE, 
                              seedDist=NULL, seedFun=NULL, seedNois=NULL){
  # n <- 1000
  # sig_RDz=NULL; sig_Ez=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  # min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL
  # max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL
  if(is.null(sig_RDz)!=is.null(sig_Ez) & is.null(min_RDz)!=is.null(min_Ez) ) stop("if there is a confounder both sig_RDz and sig_Ez shoud be given otherwise both should be NULL")
  set.seed(seedDist)
  seedDist2 <- sample(1:10000, 7)
  set.seed(NULL)
  set.seed(seedNois)
  seedNois2 <- sample(1:10000, 8)
  set.seed(NULL)
  set.seed(seedFun)
  seedFun2 <- sample(10000,2)
  set.seed(NULL)
  
  Ex <- RDdist_GP(n, tau, sigma=sig_RDx, minDR=min_RDx, maxDR=max_RDx, seedDist=seedDist2[1], seedFun=seedDist2[2])
  Ex <- apply(Ex, 2, stdrize)
  Ey <- RDdist_GP(n, tau, sigma=sig_RDy, minDR=min_RDy, maxDR=max_RDy, seedDist=seedNois2[1], seedFun=seedNois2[2])
  Ey <- apply(Ey, 2, stdrize)
  
  Ex2 <- Ex
  noisNms <- c("Ex")
  sig_Ex2 <- sig_Ex
  min_Ex2 <- min_Ex
  max_Ex2 <- max_Ex
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    #print("conf")
    Ez <- RDdist_GP(n, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[3], seedFun=seedNois2[4])
    Ez <- apply(Ez, 2, stdrize)
    Ex2 <- cbind(Ex2, Ez)
    noisNms <- c(noisNms,"Ez")
    sig_Ex2 <- c(sig_Ex2, sig_Ez)
    min_Ex2 <- c(min_Ex2, min_Ez)
    max_Ex2 <- c(max_Ex2, max_Ez)
  }
  nois <- cbind(Ex2, Ey)
  noisNms <- c(noisNms, "Ey")
  
  x <- rGP_rbf(Ex2, tau, sigma=sig_Ex2, minDR=min_Ex2, maxDR=max_Ex2, seed=seedDist2[3])
  x <- apply(x, 2, stdrize)
  x2 <- cbind(x, Ey)
  
  sig_x2 <- c(sig_x, sig_Ey)
  min_x2 <- c(min_x, min_Ey)
  max_x2 <- c(max_x, max_Ey)
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    x2 <- cbind(x2, Ez)
    sig_x2 <- c(sig_x2, sig_Ez)
    min_x2 <- c(min_x2, min_Ez)
    max_x2 <- c(max_x2, max_Ez)
  }
  
  if(inducing){
    #numInducing <- 10
    ExInd <- RDdist_GP(numInducing, tau, sigma=sig_RDxInd, minDR=min_RDxInd, maxDR=max_RDxInd, seedDist=seedDist2[4], seedFun=seedDist2[5])
    ExInd <- apply(ExInd, 2, stdrize)
    
    EyInd <- RDdist_GP(numInducing, tau, sigma=sig_RDyInd, minDR=min_RDyInd, maxDR=max_RDyInd, seedDist=seedNois2[5], seedFun=seedNois2[6])
    EyInd <- apply(EyInd, 2, stdrize)
    
    Ex2Ind <- ExInd
    noisNmsInd <- c("ExInd")
    sig_Ex2Ind <- sig_ExInd
    min_Ex2Ind <- min_ExInd
    max_Ex2Ind <- max_ExInd
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      #print("conf")
      EzInd <- RDdist_GP(numInducing, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[7], seedFun=seedNois2[8])
      EzInd <- apply(Ez, 2, stdrize)
      Ex2Ind <- cbind(Ex2Ind, EzInd)
      noisNmsInd <- c(noisNmsInd,"EzInd")
      sig_Ex2Ind <- c(sig_Ex2Ind, sig_EzInd)
      min_Ex2Ind <- c(min_Ex2Ind, min_EzInd)
      max_Ex2Ind <- c(max_Ex2Ind, max_EzInd)
    }
    noisInd <- cbind(Ex2Ind, EyInd)
    noisNmsInd <- c(noisNmsInd, "EyInd")
    
    xInd <- rGP_rbf(Ex2Ind, tau, sigma=sig_Ex2Ind, minDR=min_Ex2Ind, maxDR=max_Ex2Ind, seed=seedDist2[6])
    xInd <- apply(xInd, 2, stdrize)
    x2Ind <- cbind(xInd, EyInd)
    
    sig_x2Ind <- c(sig_xInd, sig_EyInd)
    min_x2Ind <- c(min_xInd, min_EyInd)
    max_x2Ind <- c(max_xInd, max_EyInd)
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      x2Ind <- cbind(x2Ind, EzInd)
      sig_x2Ind <- c(sig_x2Ind, sig_EzInd)
      min_x2Ind <- c(min_x2Ind, min_EzInd)
      max_x2Ind <- c(max_x2Ind, max_EzInd)
    }
  } else{
    xInd2 <- NULL
  }
  
  
  #x=x2; tau; sigma=sig_x2; minDR=min_x2; maxDR=max_x2; seed=seedFun2[1]; seedInd=seedFun2[2]
  y <- rGP_rbf(x=x2, tau, xInd=x2Ind, inducing=inducing, sigma=sig_x2, minDR=min_x2, maxDR=max_x2, seed=seedFun2[1])
  set.seed(NULL)
  
  y <- apply(y, 2, stdrize)
  
  # add noise
  if(addNoise){
    set.seed(seedDist2[7])
    addNoiseX <- rnorm(n, sd=sig_nois_x)
    x <- x + addNoiseX
    set.seed(NULL)
    addNoiseY <- rnorm(n, sd=sig_nois_y) 
    y <- y + addNoiseY
    nois <- cbind(nois, addNoiseX, addNoiseY)
    noisNms <- c(noisNms, "addNoiseX","addNoiseY")
  }
  #plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y, nois)
  colnames(res) <- c("x","y", noisNms)
  res <- list(res=res)
  if(inducing){
    resInd <- cbind(xInd, noisInd)
    colnames(resInd) <- c("xInd", noisNmsInd)
    res <- c(res, list(resInd=resInd))
  }
  
  return(res)
}



SIMdist_GP_sigPar_wrapper <- function(q=100, n=1000, seed=NULL, seed_rep=1, tau=NULL, type=c("SIM", "SIMc","SIMG","SIMln"), 
                                      calcStats=FALSE, nms=NULL, ...){
  #q=100; n=1000; seed=NULL; tau=NULL; sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL; sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  type <- match.arg(type, choices=c("SIM", "SIMc","SIMG","SIMln"))
  
  # pars <- list()
  pars <- list(...)
  pars$n <- n
  #n <- 100
  #type <- "SIM"
  
  pars <- switch(type, 
                 SIM={
                   set.seed(seed) 
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                   
                 },
                 SIMc={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDz"]]) & is.null(pars[["min_RDz"]])) pars$sig_RDz <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_Ez"]]) & is.null(pars[["min_Ez"]])) pars$sig_Ez <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                 },
                 SIMG={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                 },
                 SIMln={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=1.5*200), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                   
                   
                 })
  
  print("names(pars)")
  print(names(pars))
  
  pars <- lapply(1:q, function(i) lapply(pars, function(el) if(length(el)==1) el else el[i]))
  
  
  # par <- pars[[601]]
  if(FALSE){
    sig_RDx=NA; sig_RDy=NA; sig_Ex=NA; sig_Ey=NA; sig_x=NA;  sig_RDz=NA; sig_Ez=NA; 
    min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL;
    max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL;
    inducing=FALSE; numInducing=n;
    sig_nois_x=NULL; sig_nois_y=NULL;
    addNoise=TRUE; 
    seedDist=NULL; seedFun=NULL; seedNois=NULL
    indx <- setdiff(1:length(par),grep("ind",names(pars)))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
    sig_RDxInd=sig_RDx; sig_RDyInd=sig_RDy; sig_ExInd=sig_Ex; sig_EyInd=sig_Ey; sig_xInd=sig_x; sig_RDzInd=sig_RDz; sig_EzInd=sig_Ez; 
    min_RDxInd=min_RDx; min_RDyInd=min_RDy; min_ExInd=min_Ex; min_EyInd=min_Ey; min_xInd=min_x; min_RDzInd=min_RDz; min_EzInd=min_Ez;
    max_RDxInd=max_RDx; max_RDyInd=max_RDy; max_ExInd=max_Ex; max_EyInd=max_Ey; max_xInd=max_x; max_RDzInd=max_RDz; max_EzInd=max_Ez;
    indx <- grep("ind",names(pars))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
  }
  # dat <- do.call("SIMdist_GP_sigPar", par)
  
  
  # plot(dat[,"x"], dat[,"y"])
  
  if(!is.null(nms)) nms <- paste(type, 1:q, nms, sep=".") else nms <- paste(type, 1:q, sep=".")
  
  mc_cores <- detectCores()/2
  #mc_cores <- 1
  pm <- proc.time()
  xs <- mcmapply(function(par, nm){
    # par <- pars[[64]]
    print(paste("nm: ", nm))
    res <- do.call("SIMdist_GP_sigPar", par)
    return(res)
  }, par=pars, nm=nms, mc.cores=mc_cores, SIMPLIFY=FALSE)
  proc.time() - pm 
  # 1.241 secs for q=100, n=100 and 4 cores with confounder (3 gps)
  # 7.15 mins for q=100, n=1000 and 4 cores with confounder (3 gs)
  xs2 <- lapply(xs, function(el) el$res[,c("x","y")])
  ind_xs <- lapply(xs, function(el) el$resInd)
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs, function(el) el$res[,-which(colnames(el$res) %in% c("x","y"))])
  names(xs2) <- names(ns) <- names(dags) <- names(ind_xs) <- nms
  dataList <- list(xs=xs2, ns=ns, dags=dags, names=nms, ind_xs=ind_xs)
  
  if(calcStats){
    print("calc stats")
    pm <- proc.time()
    count <- 0
    descr_indices <- sapply(xs2, function(xs){
      count <<- count + 1
      # count <- 1;  xs <- xs2[[count]]
      print("*****************")
      print(count)
      pm <- proc.time()
      num_obs <- 100
      num_reps <- floor(nrow(xs))/num_obs
      pm <- proc.time()
      res <- mcmapply(function(j){
        #j <- 1
        #print(paste("j: ",j))
        
        # plot(xs[,"x"], xs[,"y"])
        # hist(xs[,"x"])
        # hist(xs[,"y"])
        smpl <- sample(1:nrow(xs), num_obs)
        df <- as.data.frame(xs[smpl,])
        mod <- lm(y~x,df)
        #plot(df$x, residuals(mod))
        #hist(residuals(mod))
        xs_2 <- apply(xs[smpl,], 2, norml)
        #plot(xs_2[,"x"], xs_2[,"y"])
        trDat <- constructData(x=as.matrix(xs_2[,"x"]), y=as.matrix(xs_2[,"y"]))
        krrAux <- setParams(learner=krr1, trainData=trDat)
        krrAux <- krrAux$learn(krrAux)
        predKrr <- pred.CV(krrAux, trDat)$test
        #o <- order(trDat$x); plot(trDat$x, trDat$y); lines(trDat$x[o], predKrr$gyh_class[o], col="red")
        residsLin <- residuals(mod)
        # plot(xs[,"x"], residsLin)
        residsNonLin <- krrAux$resids(krrAux, predKrr)[,"resid"]
        #plot(predKrr$x_class, residsNonLin)
        #plot(o, residsNonLin[o])
        res1 <- dhsic.test(X=df[,"x",drop=F],Y=matrix(residsLin,ncol=1))$p.value
        res2 <- dhsic.test(X=predKrr$x_class,Y=residsNonLin)$p.value
        #o <- order(predKrr$x_class)
        #res2b <- dhsic.test(X=o,Y=residsNonLin[o])$p.value
        res3 <- ks.test(x=residsLin, y="pnorm", mean=mean(residsLin), sd=sd(residsLin))$p.value
        res4 <- ks.test(x=residsNonLin, y="pnorm", mean=mean(residsNonLin), sd=sd(residsNonLin))$p.value
        res5 <- ks.test(x=xs_2[,"x"], y="pnorm", mean=mean(xs_2[,"x"]), sd=sd(xs_2[,"x"]))$p.value
        res6 <- ks.test(x=xs_2[,"x"], y="punif", min=min(xs_2[,"x"]), max=max(xs_2[,"x"]))$p.value
        res7 <- Shannon_KDP(xs_2[,"x"])
        
        modDens <- kepdf(xs_2[,"x"], eval.points = xs_2[,"x"], kernel = "gaussian", bwtype = "adaptive")
        #hist(xs_2[,"x"], prob=T); o <- order(modDens@eval.points); lines(modDens@eval.points[o], modDens@estimate[o], col="red")
        res8 <- max(modDens@estimate)-min(modDens@estimate)
        res <- c(res1, res2, res3, res4, res5, res6, res7, res8)
        names(res) <- c("lm_indep","add_indep","lm_gauss","add_gauss",
                        "cause_gauss","cause_unif", "cause_ent","cause_rngPdf")
        #print(res)
        proc.time() - pm 
        return(res)
      }, SIMPLIFY="array", j=1:num_reps, mc.cores=detectCores()/2)
      proc.time()- pm #37.5  with 1 core, 22.5 with 4 cores (100pts)
      res <- apply(res, 1, mean)
      return(res)
    })
    proc.time() - pm #34 mins
    descr_indices <- t(descr_indices)
    dataList <- c(dataList, list(descr_indices=descr_indices))
  }
  
  
  #plotPairsList(dataList)
  return(dataList)  
}


rKalpha_rbf <- function(x, tau, xInd=NULL, inducing=FALSE, sigma=NULL, minDR=NULL, maxDR=NULL, seed=NULL, mean_alpha=0){
  #if((ncol(x) !=length(sigma) & length(sigma)!=1) & (ncol(x) !=length(cutoff) & length(cutoff)!=1)) error("one sigma per column or one sigma for all columns!")
  #if(is.null(sigma) & is.null(cutoff)) error("one of sigma and error shd be provided")
  n <- nrow(x)
  if(inducing){
    x2 <- xInd
    removeDiag <- FALSE
  } else{
    x2 <- x
    removeDiag <- TRUE
  }
  
  if(any(is.na(sigma))){
    
    #print("minDR"); print(minDR)
    #print("maxDR"); print(maxDR)
    sigmaAux <- sapply(which(is.na(sigma)), function(i){
      # i <- 2
      # fit sigma
      # print(paste("i: ", i))
      ord <- 10
      sigmas1 <- (10^seq(-ord,ord,1)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        #print(paste("sd: ", sd1))
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F], sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
        res <- var(distsX)
        return(res)
      })
      indxMaxVar <- which.max(varsHsics.sig)
      sigmaVar <- sigmas1[indxMaxVar]
      
      # obtain dynamic range of sigma
      spl <- spline(log(sigmas1,10), varsHsics.sig)
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
      
      sigmas1 <- (10^seq(DR.sig[1],DR.sig[2], length.out=20)) # *sigma0
      varsHsics.sig <- sapply(sigmas1, function(sd1){
        # i <- 1; sd1 <- sigmas1[i]
        Kxs <- kern_rbf(x=x[,i,drop=F], x2[,i,drop=F],sigma=1/(2*sd1^2))
        if(removeDiag) distsX <- (Kxs)[lower.tri(Kxs)] else distsX <- c(Kxs)
        res <- var(distsX)
        return(res)
      })
      indxMaxVar <- which.max(varsHsics.sig)
      sigmaVar <- sigmas1[indxMaxVar]
      
      # obtain dynamic range of sigma
      spl <- spline(log(sigmas1,10), varsHsics.sig, n=100*length(sigmas1))
      splf <- splinefun(log(sigmas1,10), varsHsics.sig)
      dVar.sig_dlog.sig <- splf(log(sigmas1,10), deriv=1)
      tol <- 1e-3
      DR.sig <- sigmas1[which(abs(dVar.sig_dlog.sig)>tol)]
      DR.sig <- range(log(DR.sig,10))
      
      #minDR <- 0.2
      #maxDR <- 0.4
      
      indx <- which(spl$x>=DR.sig[1] & spl$x<DR.sig[2])
      xx <- spl$x[indx]
      yy <- spl$y [indx]
      indxMax <- which.max(yy)
      xx2 <- rep(NA, length(xx))
      xx2[1:indxMax] <- yy[1:indxMax]/(2*yy[indxMax])
      xx2[(indxMax+1):length(xx)] <- 0.5+(0.5-yy[(indxMax+1):length(xx)]/(2*max(yy[(indxMax+1):length(xx)])))
      #plot(xx2, xx)
      
      splf <- splinefun(c(0,xx2,1), c(xx[1],xx,xx[length(xx)]), method="monoH.FC")
      print(c(minDR[i], maxDR[i]))
      rngSigs <- splf(c(minDR[i], maxDR[i]))
      if(abs(rngSigs[2]-rngSigs[1])<1e-3){
        rngSigs[1] <- mean(rngSigs)-1e-3
        rngSigs[2] <- mean(rngSigs)+1e-3
      }
      
      res <- runif(1, min=rngSigs[1], max=rngSigs[2])
      res <- 10^res
      
      if(FALSE){
        plot(log(sigmas1,10), varsHsics.sig, col="red", type="p")
        lines(spl)
        abline(v=DR.sig, col="purple")
        abline(v=log(sigmaVar,10), col="blue")
        abline(h=max(varsHsics.sig), col="green")
        abline(v=rngSigs, col="orange", type="p")
      }
      
      return(res)
    })
    sigma[which(is.na(sigma))] <- sigmaAux 
  }
  
  
  sigma2 <- matrix(0, ncol(x), ncol(x))
  diag(sigma2) <- 1/sigma
  x2 <- x %*% sigma2
  x2 <- x2/sqrt(2)
  
  
  if(inducing){
    
    xInd2 <- xInd %*% sigma2
    xInd2 <- xInd2/sqrt(2)
    K <- kern_rbf(x2, xInd2, sigma=1)
    #alternative K calculation
    if(FALSE){
      
      K_alt <- sapply(1:length(sigma), function(i) kern_rbf(x=x[,i,drop=F], y=xInd[,i,drop=F], sigma=1/(2*sigma[i]^2)), simplify="array")
      K_alt <- apply(K_alt, c(1,2), prod)
      smpl1 <- sample(length(K), 100)
      plot(c(K)[smpl1], c(K_alt)[smpl1]); abline(a=0, b=1, col="red")
    }
    
    
  } else{
    K <- kern_rbf(x2, sigma=1) 
    n <- nrow(K)
    I <- diag(n)
    #C <- K + tau^2*I
    
    
  }
  
  #alpha <- matrix(abs(rnorm(ncol(K))), ncol(K), 1)
  alpha <- matrix(rnorm(ncol(K), mean=mean_alpha), ncol(K), 1)
  y <- K %*% alpha 
  
  # we want cov(y, y') = K(y, y') 
  # say y=t(C) %*% x and y'= t(C') %*% x', x,x'~ N(0,I) 
  # then cov(y, y') = cov(t(C)%*%x, t(C')%*%x')
  #                 = 
  
  
  return(y)  
}



SIMdist_Kalpha_sigPar <- function(n,  tau, sig_RDx=NA, sig_RDy=NA, sig_Ex=NA, sig_Ey=NA, sig_x=NA,  sig_RDz=NA, sig_Ez=NA, 
                                  min_RDx=NULL, min_RDy=NULL, min_Ex=NULL, min_Ey=NULL, min_x=NULL, min_RDz=NULL, min_Ez=NULL,
                                  max_RDx=NULL, max_RDy=NULL, max_Ex=NULL, max_Ey=NULL, max_x=NULL, max_RDz=NULL, max_Ez=NULL,
                                  inducing=FALSE, numInducing=n, inducing_mean_x=0,
                                  sig_RDxInd=sig_RDx, sig_RDyInd=sig_RDy, sig_ExInd=sig_Ex, sig_EyInd=sig_Ey, sig_xInd=sig_x, sig_RDzInd=sig_RDz, sig_EzInd=sig_Ez, 
                                  min_RDxInd=min_RDx, min_RDyInd=min_RDy, min_ExInd=min_Ex, min_EyInd=min_Ey, min_xInd=min_x, min_RDzInd=min_RDz, min_EzInd=min_Ez,
                                  max_RDxInd=max_RDx, max_RDyInd=max_RDy, max_ExInd=max_Ex, max_EyInd=max_Ey, max_xInd=max_x, max_RDzInd=max_RDz, max_EzInd=max_Ez,
                                  sig_nois_x=NULL, sig_nois_y=NULL,
                                  addNoise=TRUE, 
                                  seedDist=NULL, seedFun=NULL, seedNois=NULL){
  # n <- 1000
  # sig_RDz=NULL; sig_Ez=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  # min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL
  # max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL
  if(is.null(sig_RDz)!=is.null(sig_Ez) & is.null(min_RDz)!=is.null(min_Ez) ) stop("if there is a confounder both sig_RDz and sig_Ez shoud be given otherwise both should be NULL")
  set.seed(seedDist)
  seedDist2 <- sample(1:10000, 7)
  set.seed(NULL)
  set.seed(seedNois)
  seedNois2 <- sample(1:10000, 8)
  set.seed(NULL)
  set.seed(seedFun)
  seedFun2 <- sample(10000,2)
  set.seed(NULL)
  
  Ex <- RDdist_GP(n, tau, sigma=sig_RDx, minDR=min_RDx, maxDR=max_RDx, seedDist=seedDist2[1], seedFun=seedDist2[2])
  Ex <- apply(Ex, 2, stdrize)
  Ey <- RDdist_GP(n, tau, sigma=sig_RDy, minDR=min_RDy, maxDR=max_RDy, seedDist=seedNois2[1], seedFun=seedNois2[2])
  Ey <- apply(Ey, 2, stdrize)
  
  Ex2 <- Ex
  noisNms <- c("Ex")
  sig_Ex2 <- sig_Ex
  min_Ex2 <- min_Ex
  max_Ex2 <- max_Ex
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    #print("conf")
    Ez <- RDdist_GP(n, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[3], seedFun=seedNois2[4])
    Ez <- apply(Ez, 2, stdrize)
    Ex2 <- cbind(Ex2, Ez)
    noisNms <- c(noisNms,"Ez")
    sig_Ex2 <- c(sig_Ex2, sig_Ez)
    min_Ex2 <- c(min_Ex2, min_Ez)
    max_Ex2 <- c(max_Ex2, max_Ez)
  }
  nois <- cbind(Ex2, Ey)
  noisNms <- c(noisNms, "Ey")
  
  x <- rGP_rbf(Ex2, tau, sigma=sig_Ex2, minDR=min_Ex2, maxDR=max_Ex2, seed=seedDist2[3])
  x <- apply(x, 2, stdrize)
  x2 <- cbind(x, Ey)
  
  sig_x2 <- c(sig_x, sig_Ey)
  min_x2 <- c(min_x, min_Ey)
  max_x2 <- c(max_x, max_Ey)
  
  if(!is.na(sig_RDz) | !is.null(min_RDz)){
    x2 <- cbind(x2, Ez)
    sig_x2 <- c(sig_x2, sig_Ez)
    min_x2 <- c(min_x2, min_Ez)
    max_x2 <- c(max_x2, max_Ez)
  }
  
  if(inducing){
    #numInducing <- 10
    ExInd <- RDdist_GP(numInducing, tau, sigma=sig_RDxInd, minDR=min_RDxInd, maxDR=max_RDxInd, seedDist=seedDist2[4], seedFun=seedDist2[5])
    ExInd <- apply(ExInd, 2, stdrize)
    
    EyInd <- RDdist_GP(numInducing, tau, sigma=sig_RDyInd, minDR=min_RDyInd, maxDR=max_RDyInd, seedDist=seedNois2[5], seedFun=seedNois2[6])
    EyInd <- apply(EyInd, 2, stdrize)
    
    Ex2Ind <- ExInd
    noisNmsInd <- c("ExInd")
    sig_Ex2Ind <- sig_ExInd
    min_Ex2Ind <- min_ExInd
    max_Ex2Ind <- max_ExInd
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      #print("conf")
      EzInd <- RDdist_GP(numInducing, tau, sigma=sig_RDz, minDR=min_RDz, maxDR=max_RDz, seedDist=seedNois2[7], seedFun=seedNois2[8])
      EzInd <- apply(Ez, 2, stdrize)
      Ex2Ind <- cbind(Ex2Ind, EzInd)
      noisNmsInd <- c(noisNmsInd,"EzInd")
      sig_Ex2Ind <- c(sig_Ex2Ind, sig_EzInd)
      min_Ex2Ind <- c(min_Ex2Ind, min_EzInd)
      max_Ex2Ind <- c(max_Ex2Ind, max_EzInd)
    }
    noisInd <- cbind(Ex2Ind, EyInd)
    noisNmsInd <- c(noisNmsInd, "EyInd")
    
    xInd <- rGP_rbf(Ex2Ind, tau, sigma=sig_Ex2Ind, minDR=min_Ex2Ind, maxDR=max_Ex2Ind, seed=seedDist2[6])
    xInd <- apply(xInd, 2, stdrize)
    xInd <- xInd + inducing_mean_x
    x2Ind <- cbind(xInd, EyInd)
    
    sig_x2Ind <- c(sig_xInd, sig_EyInd)
    min_x2Ind <- c(min_xInd, min_EyInd)
    max_x2Ind <- c(max_xInd, max_EyInd)
    
    if(!is.na(sig_RDzInd) | !is.null(min_RDzInd)){
      x2Ind <- cbind(x2Ind, EzInd)
      sig_x2Ind <- c(sig_x2Ind, sig_EzInd)
      min_x2Ind <- c(min_x2Ind, min_EzInd)
      max_x2Ind <- c(max_x2Ind, max_EzInd)
    }
  } else{
    xInd2 <- NULL
  }
  
  
  #x=x2; tau; xInd=x2Ind; inducing=inducing; sigma=sig_x2; minDR=min_x2; maxDR=max_x2; seed=seedFun2[1]
  y <- rKalpha_rbf(x=x2, tau, mean_alpha=-inducing_mean_x,xInd=x2Ind, inducing=inducing, sigma=sig_x2, minDR=min_x2, maxDR=max_x2, seed=seedFun2[1])
  set.seed(NULL)
  
  y <- apply(y, 2, stdrize)
  
  # add noise
  if(addNoise){
    set.seed(seedDist2[7])
    addNoiseX <- rnorm(n, sd=sig_nois_x)
    x <- x + addNoiseX
    set.seed(NULL)
    addNoiseY <- rnorm(n, sd=sig_nois_y) 
    y <- y + addNoiseY
    nois <- cbind(nois, addNoiseX, addNoiseY)
    noisNms <- c(noisNms, "addNoiseX","addNoiseY")
  }
  #plot(x,y)
  #plot(Ey, y)
  res <- cbind(x,y, nois)
  colnames(res) <- c("x","y", noisNms)
  res <- list(res=res)
  if(inducing){
    resInd <- cbind(xInd, noisInd)
    colnames(resInd) <- c("xInd", noisNmsInd)
    res <- c(res, list(resInd=resInd))
  }
  
  return(res)
}


SIMdist_Kalpha_sigPar_wrapper <- function(q=100, n=1000, seed=NULL, seed_rep=1, tau=NULL, type=c("SIM", "SIMc","SIMG","SIMln"), 
                                          calcStats=FALSE, nms=NULL, ...){
  #q=100; n=1000; seed=NULL; tau=NULL; sig_RDx=NULL; sig_RDy=NULL; sig_RDz=NULL; sig_Ex=NULL; sig_Ey=NULL; sig_Ez=NULL; sig_x=NULL; sig_nois_x=NULL; sig_nois_y=NULL; addNoise=TRUE; seedDist=NULL; seedFun=NULL; seedNois=NULL
  type <- match.arg(type, choices=c("SIM", "SIMc","SIMG","SIMln"))
  
  # pars <- list()
  pars <- list(...)
  print("extra pars: ")
  print(names(pars))
  pars$n <- n
  #n <- 100
  #type <- "SIM"
  
  pars <- switch(type, 
                 SIM={
                   set.seed(seed) 
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                   
                 },
                 SIMc={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDz"]]) & is.null(pars[["min_RDz"]])) pars$sig_RDz <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_Ez"]]) & is.null(pars[["min_Ez"]])) pars$sig_Ez <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   set.seed(NULL)
                   pars
                 },
                 SIMG={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=1e6, scale=1e-3), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=15), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.1), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                 },
                 SIMln={
                   set.seed(seed)
                   if(is.null(pars[["sig_RDx"]]) & is.null(pars[["min_RDx"]])) pars$sig_RDx <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   if(is.null(pars[["sig_RDy"]]) & is.null(pars[["min_RDy"]])) pars$sig_RDy <- rep(rgamma(q/seed_rep, shape=5, scale=0.1), seed_rep)
                   pars$sig_RDz <- NULL
                   if(is.null(pars[["sig_Ex"]]) & is.null(pars[["min_Ex"]])) pars$sig_Ex <- rep(rgamma(q/seed_rep, shape=2, scale=1.5), seed_rep)
                   if(is.null(pars[["sig_x"]]) & is.null(pars[["min_x"]])) pars$sig_x <- rep(rgamma(q/seed_rep, shape=2, scale=3), seed_rep)
                   if(is.null(pars[["sig_Ey"]]) & is.null(pars[["min_Ey"]])) pars$sig_Ey <- rep(rgamma(q/seed_rep, shape=2, scale=1.5*200), seed_rep)
                   pars$sig_Ez <- NULL
                   if(is.null(pars[["sig_nois_x"]])) pars[["sig_nois_x"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["sig_nois_y"]])) pars[["sig_nois_y"]] <- rep(rgamma(q/seed_rep, shape=2, scale=0.01), seed_rep)
                   if(is.null(pars[["tau"]])) pars$tau <- 1e-4
                   if(is.null(pars[["addNoise"]])) pars$addNoise=TRUE
                   
                   set.seed(NULL)
                   pars
                   
                   
                 })
  
  print("names(pars)")
  print(names(pars))
  
  pars <- lapply(1:q, function(i) lapply(pars, function(el) if(length(el)==1) el else el[i]))
  
  
  # par <- pars[[1]]
  if(FALSE){
    sig_RDx=NA; sig_RDy=NA; sig_Ex=NA; sig_Ey=NA; sig_x=NA;  sig_RDz=NA; sig_Ez=NA; 
    min_RDx=NULL; min_RDy=NULL; min_Ex=NULL; min_Ey=NULL; min_x=NULL; min_RDz=NULL; min_Ez=NULL;
    max_RDx=NULL; max_RDy=NULL; max_Ex=NULL; max_Ey=NULL; max_x=NULL; max_RDz=NULL; max_Ez=NULL;
    inducing=FALSE; numInducing=n; inducing_mean_x=0
    sig_nois_x=NULL; sig_nois_y=NULL;
    addNoise=TRUE; 
    seedDist=NULL; seedFun=NULL; seedNois=NULL
    indx <- setdiff(1:length(par),grep("ind",names(pars)))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
    sig_RDxInd=sig_RDx; sig_RDyInd=sig_RDy; sig_ExInd=sig_Ex; sig_EyInd=sig_Ey; sig_xInd=sig_x; sig_RDzInd=sig_RDz; sig_EzInd=sig_Ez; 
    min_RDxInd=min_RDx; min_RDyInd=min_RDy; min_ExInd=min_Ex; min_EyInd=min_Ey; min_xInd=min_x; min_RDzInd=min_RDz; min_EzInd=min_Ez;
    max_RDxInd=max_RDx; max_RDyInd=max_RDy; max_ExInd=max_Ex; max_EyInd=max_Ey; max_xInd=max_x; max_RDzInd=max_RDz; max_EzInd=max_Ez;
    indx <- grep("ind",names(pars))
    if(length(indx)>0) for(i in indx) assign(names(par)[i], par[[i]])
  }
  # dat <- do.call("SIMdist_Kalpha_sigPar", par)
  
  
  # plot(dat$res[,"x"], dat$res[,"y"])
  
  if(!is.null(nms)) nms <- paste(type, 1:q, nms, sep=".") else nms <- paste(type, 1:q, sep=".")
  
  mc_cores <- detectCores()/2
  #mc_cores <- 1
  print("xs")
  pm <- proc.time()
  xs <- mcmapply(function(par, nm){
    # par <- pars[[64]]
    print(paste("nm: ", nm))
    res <- do.call("SIMdist_Kalpha_sigPar", par)
    return(res)
  }, par=pars, nm=nms, mc.cores=mc_cores, SIMPLIFY=FALSE)
  proc.time() - pm 
  # 1.241 secs for q=100, n=100 and 4 cores with confounder (3 gps)
  # 7.15 mins for q=100, n=1000 and 4 cores with confounder (3 gs)
  print("xs2")
  xs2 <- lapply(xs, function(el) el$res[,c("x","y")])
  print("ind_xs")
  ind_xs <- lapply(xs, function(el) el$resInd)
  print("dags,ns")
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs, function(el) el$res[,-which(colnames(el$res) %in% c("x","y"))])
  names(xs2) <- names(ns) <- names(dags) <- names(ind_xs) <- nms
  print("dagList")
  dataList <- list(xs=xs2, ns=ns, dags=dags, names=nms, ind_xs=ind_xs)
  
  
  
  #plotPairsList(dataList)
  return(dataList)  
}


simProsailLike <- function(n, c, aMin, aMax, maxNoise, seed=NULL){
  set.seed(seed)
  x <- runif(n)
  fx <- exp(-c*(x))
  #plot(x, fx)
  
  g1 <- aMin+(aMax-aMin)*(1-exp(-c*x))
  shape2b <- 5#runif(1,0,5)
  shape1b <- 2#runif(1,0,shape2b)
  notDone <- TRUE
  bMin <- aMin
  bMax <- aMax
  while(notDone){
    bMin <- bMin*rbeta(1, shape1=shape1b, shape2=shape2b) #bMin*0.99#
    bMax <- bMax*rbeta(1, shape1=shape1b, shape2=shape2b) #bMax*0.99#
    #print(paste("bmin: ", bMin, ", bmax: ", bMax))
    #b2 <- 5 
    h <- bMin+(bMax-bMin)*(1-exp(-c*x))
    #plot(x, g1, ylim=range(c(g1,h)))
    #lines(x, h, type="p", col="red")
    #plot(x, h/g1, ylim=range(c(h/g1, 1-g1)))
    #lines(x, 1-g1, type="p", col="red")
    g2 <- (1-g1)/g1
    shape1 <- (1-g1-(1+g2)*h)/(h*(1+g2)^2)
    shape2 <- shape1*g2
    if(any(shape1<=0 | shape2<=0)){ #) | any(shape1 > shape2)
      #print("not Done")
    } else{
      notDone <- FALSE
      #print("done")
    }
  }
  #plot(x, shape1); abline(h=0, col="red")
  #plot(x, shape2); abline(h=0, col="red")
  #plot(shape1, shape2); abline(a=0, b=1, col="red")
  ny <- maxNoise*rbeta(n, shape1=shape1, shape2=shape2)
  #summary(ny)
  #hist(ny)
  #plot(x, ny)
  y <-fx +ny
  #plot(x, y)
  res <- cbind(x,y)
  colnames(res) <- c("x","y")
  return(res)
}

simProsailLike_apply <- function(q, n, c, aMin, aMax, maxNoise, seeds=NULL, nms=1:q){
  
  if(is.null(seeds)) seeds <- sample(1:1000000, q)
  xs <- mcmapply(FUN=function(c, aMin, aMax, maxNoise, seed){
    x <- simProsailLike(n, c, aMin, aMax, maxNoise, seed)
  }, c=c, aMin=aMin, aMax=aMax, maxNoise=maxNoise, seed=seeds, SIMPLIFY=FALSE, mc.cores=1)
  
  dag <- matrix(c(0,0,1,0), 2, 2)
  colnames(dag) <- rownames(dag) <- c("x","y")
  dags <- lapply(1:q, function(el) dag)
  ns <- lapply(xs, function(el) NA)
  names(xs) <- names(ns) <- names(dags) <- nms
  
  return(list(dags=dags, xs=xs, noiss=ns, names=nms))
}

calcStatss <- function(dataList){
  xs2 <- dataList$xs
  print("calc stats")
  pm <- proc.time()
  count <- 0
  descr_indices <- mcmapply(function(x){
    count <<- count + 1
    # count <- 1;  x <- xs2[[count]]
    # x <- runif(100,-1,1); y = x^2; x <- cbind(x,y); plot(x[,1],x[,2])
    print("*****************")
    print(count)
    pm <- proc.time()
    num_obs <- 100
    num_reps <- floor(nrow(x))/num_obs
    x <- apply(x, 2, norml)
    mod <- mkde(x)
    probs <- 1/mod  
    probs <- probs/sum(probs)
    
    pm <- proc.time()
    res <- mcmapply(function(j){
      #j <- 1
      #print(paste("j: ",j))
      
      # plot(xs[,"x"], xs[,"y"])
      # hist(xs[,"x"])
      # hist(xs[,"y"])
      smpl <- sample(1:nrow(x), num_obs, prob=probs)
      xs_2 <- apply(x[smpl,], 2, norml)
      dirs <- list(c("x","y"), c("y","x"))
      if(TRUE){
      res <- sapply(dirs, function(nms){
        # nms <- dirs[[1]]
        xs_3 <- xs_2
        colnames(xs_3) <- nms
        df <- as.data.frame(xs_3)
        mod <- lm(y~x,df)
        #plot(df$x, residuals(mod))
        #hist(residuals(mod))
        #plot(xs_3[,"x"], xs_3[,"y"])
        trDat <- constructData(x=as.matrix(xs_3[,"x"]), y=as.matrix(xs_3[,"y"]))
        krrAux <- setParams(learner=krr1, trainData=trDat)
        krrAux <- krrAux$learn(krrAux)
        predKrr <- pred.CV(krrAux, trDat)$test
        #o <- order(trDat$x); plot(trDat$x, trDat$y); lines(trDat$x[o], predKrr$gyh_class[o], col="red")
        residsLin <- residuals(mod)
        # plot(xs_3[,"x"], residsLin)
        residsNonLin <- krrAux$resids(krrAux, predKrr)[,"resid"]
        #plot(predKrr$x_class, residsNonLin)
        #plot(o, residsNonLin[o])
        res1 <- dhsic.test(X=df[,"x",drop=F],Y=matrix(residsLin,ncol=1))$p.value
        res2 <- dhsic.test(X=predKrr$x_class,Y=residsNonLin)$p.value
        #o <- order(predKrr$x_class)
        #res2b <- dhsic.test(X=o,Y=residsNonLin[o])$p.value
        res3 <- ks.test(x=residsLin, y="pnorm", mean=mean(residsLin), sd=sd(residsLin))$p.value
        res4 <- ks.test(x=residsNonLin, y="pnorm", mean=mean(residsNonLin), sd=sd(residsNonLin))$p.value
        res5 <- ks.test(x=xs_3[,"x"], y="pnorm", mean=mean(xs_3[,"x"]), sd=sd(xs_3[,"x"]))$p.value
        res6 <- ks.test(x=xs_3[,"x"], y="punif", min=min(xs_3[,"x"]), max=max(xs_3[,"x"]))$p.value
        res7 <- Shannon_KDP(xs_3[,"x"])
        
        modDens <- kepdf(xs_3[,"x"], eval.points = xs_3[,"x"], kernel = "gaussian", bwtype = "adaptive")
        #hist(xs_3[,"x"], prob=T); o <- order(modDens@eval.points); lines(modDens@eval.points[o], modDens@estimate[o], col="red")
        res8 <- max(modDens@estimate)-min(modDens@estimate)
        res <- c(res1, res2, res3, res4, res5, res6, res7, res8)
        names(res) <- c("lm_indep","add_indep","lm_gauss","add_gauss",
                        "cause_gauss","cause_unif", "cause_ent","cause_rngPdf")
        return(res)
      }, simplify="array")  
      absDiff <- function(x) max(x) - min(x) 
      funcs <- c(rep("max",6), rep("absDiff",2))
      #res <- sapply(1:nrow(res), function(m) do.call(funcs[m], list(res[m,])))  
      }
      #bijectivity index
      
      #sigmaX <- 1/median(as.numeric(dist(unique(xs_2))^2))
      #K <- kern_rbf(xs_2, sigma=sigmaX)
      #K <- 1-K
      K <- sqrt(matNorm2(xs_2))
      
      #library(cppRouting)
      ktry <- 20
      ksuccess <- FALSE
      #K2 <- K
      #K2 <- apply(K, 2, function(col){
      #  res <- col
      #  res[res>sort(col)[ktry]] <- NA
      #  return(res)
      #})
      #K3 <- melt(K2)
      #K3 <- na.omit(K3)
      #colnames(K3) <- c("from","to","cost")
      #G <- makegraph(K3, directed = FALSE, coords = NULL)
      #K4 <- get_distance_matrix(G, from=1:nrow(K), to=1:nrow(K))
      #plot(c(K), c(K4)); abline(a=0, b=1, col="red")
      
      #diag(K2) <- Inf
      #twoLong <- max(apply(K2,2,min))*1.01
      while(!ksuccess){
        #print(paste("ktry: ", ktry))
        #L <- as.matrix(stepacross(K, toolong=twoLong))
        L <- try(as.matrix(isomapdist(K, k=ktry)), silent=TRUE)
        ktry <- ktry + 1
        if(class(L) != "try-error") ksuccess <- TRUE
      }
      L <- as.matrix(L)
      #L <- apply(L, 2, norml)+0.5
      #K2 <- apply(K, 2, norml)+0.5
      #i <- 1
      #plot(c(L[,i]), c(K2[,i])); abline(a=0,b=1, col="red"); i <- i + 1
      
      #apply(L, 2, which.max)
      #L2 <- as.matrix(isomapdist(K, k=nrow(K)))
      #L2 <- as.matrix(stepacross(K, toolong=max(K)+0.0001))
      res9 <- log(L/K) # quantile(, na.rm=T, probs=0.5)
      #print(summary(c(res9)))
      #plot(c(K), c(L), main=count); abline(a=0, b=1, col="red")
      #plot(c(L2), c(L)); abline(a=0, b=1, col="red")
      res9 <- mean(res9, na.rm=T)
      names(res9) <- "bij"
      #res9 <- apply(res9,2,min, na.rm=T)
      res <- list(bi=res, uni=res9)
      #names(res) <- c("lm_indep","add_indep","lm_gauss","add_gauss",
      #                "cause_gauss","cause_unif", "cause_ent","cause_rngPdf", "bij")
      #res <- res9
      #print(res)
      proc.time() - pm 
      return(res)
    }, SIMPLIFY=FALSE, j=1:num_reps, mc.cores=1)#
    proc.time()- pm #37.5  with 1 core, 22.5 with 4 cores (100pts)
    #res <- apply(res,1,mean)
    res <- list(bi=sapply(res, function(el) el$bi, simplify="array"), 
                uni=sapply(res, function(el) el$uni, simplify="array"))
    return(res)
  }, x=xs2, SIMPLIFY=FALSE, mc.cores=detectCores()/2)
  proc.time() - pm #34 mins
  
  descr_indices <- list(uni=sapply(descr_indices, function(el) el$uni, simplify="array"),
                        bi=sapply(descr_indices, function(el) el$bi, simplify="array"))
  
  #alpha <- as.numeric(sapply(strsplit(sapply(strsplit(names(descr_indices), "_"), function(el) el[[2]]),"\\."), function(el) el[[1]]))
  #alpha <- as.numeric(sapply(strsplit(sapply(strsplit(colnames(descr_indices), "_"), function(el) el[[2]]),"\\."), function(el) el[[1]]))
  #alpha <- rep(alpha, rep(nrow(descr_indices), length(alpha)))
  #boxplot(c(descr_indices)~alpha)
  #agg <- aggregate(c(descr_indices[which(!is.infinite(descr_indices))]), list(alpha[which(!is.infinite(descr_indices))]), FUN=mean, na.rm=T)
  #plot(agg[,1], agg[,2])
  
  #descr_indices <- t(descr_indices)
  dataList <- c(dataList, list(descr_indices=descr_indices))
  return(dataList)
}

getStatArrs <- function(stats){
  # 1. bidrectional - 
  statsArrBi <- sapply(stats$bi, function(matList){
    # matList <- stats$bi[[1]]
    
    # a) aggregate by direction: max, min, absDif
    agResMat <- sapply(matList, function(mat){
      # mat <- matList[[1]]
      agRes <- apply(mat, 1, function(row) c(min(row), max(row), max(row)-min(row)))
      agRes <- t(agRes)
      colnames(agRes) <- c("min","max","absDif")
      return(agRes)
    }, simplify="array")
    
    # b) aggregate by boostrap: max, min, mean, median, sd, absDif
    dim(agResMat)
    agAgResMat <- apply(agResMat, c(1,2), function(col) c(min=min(col), max=max(col), absDif=max(col)-min(col), mean=mean(col), median=median(col), sd=sd(col) ))  
    names(dimnames(agAgResMat)) <- c("aggBoot", "msr","aggDir")
    return(agAgResMat)
  }, simplify="array")
  names(dimnames(statsArrBi))[4] <- "names"  
  dim(statsArrBi)
  # 2. unidirectional
  statsArrUni <- sapply(stats$uni, function(matList){
    # matList <- stats$uni[[1]]
    
    # a) aggregate by boostrap: max, min, mean, median, sd, absDif
    agResVec <- unlist(matList)
    agAgResVec <- c(min=min(agResVec), 
                    max=max(agResVec), 
                    absDif=max(agResVec)-min(agResVec), 
                    mean=mean(agResVec), median=median(agResVec), 
                    sd=sd(agResVec))
    return(agAgResVec)
  }, simplify="array")
  names(dimnames(statsArrUni)) <- c("aggBoot", "names")  
  dim(statsArrUni)
  return(list(uni=statsArrUni, bi=statsArrBi))
}

getStatsLong <- function(statArrs){  
  statsBiDf <- melt(statArrs$bi)
  statsUniDf <- melt(statArrs$uni)
  
  statsUniDf$msr <- "bij"
  statsUniDf$aggDir <- "only"
  statsUniDf <- statsUniDf[,colnames(statsBiDf)]
  statsDf <- rbind(statsBiDf, statsUniDf)
  return(statsDf)
}

plotStats_alphaProg <- function(statLongDF, msrs, aggBoots, aggDirs){
  aux <- strsplit( sapply(strsplit(as.character(statLongDF$name), "\\."), function(el) el[length(el)]), "_")
  statLongDF$alpha <- as.numeric(sapply(aux, function(el) el[2]))
  statLongDF$rep <- as.numeric(sapply(aux, function(el) el[1]))
  indx <- with(statLongDF, which(msr %in% msrs &  aggBoot %in% aggBoots &
                                   aggDir %in% aggDirs))
  
  p <- ggplot(statLongDF[indx,])
  p <- p + geom_boxplot(aes(x=factor(alpha), y=value, colour=aggBoot))
  p <- p + facet_wrap(aggDir~msr, scales="free")
  print(p)
}

getStatsWide <- function(statLongDF){
  statLongDF$var <- with(statLongDF, paste(msr, aggDir, aggBoot, sep="_"))
  res <- cast(statLongDF, names~var, value="value")
  return(res)
}


# visualize outputof simRandSEMs when all(ps==2)
plotPairsList <- function(dataList, hit=NULL, sizePts=0.5, shapePts=4, alpha=0.2, legend=TRUE, bg=TRUE){
  # hit=NULL; sizePts=0.5; shapePts=4; alpha=0.2; legend=TRUE; bg=TRUE
  trueDAGid <-  getHypID(sapply(dataList$dags, function(dag) dag, simplify="array"))$id
  data <- melt(dataList$xs)
  colnames(data) <- c("numData","xY","val","dataset")
  data <- cast(data, dataset+numData~xY, value="val")
  colnames(data) <- c("dataset","numData","x","y")
  data$id <- as.factor(trueDAGid[match(data$dataset, dataList$names)])
  data$dataset2 <- as.numeric(factor(data$dataset))
  
  if(!is.null(hit)){
    data$hit <- hit[match(data$dataset, names(hit))]
    data$hit <- as.factor(c("miss","hit")[(data$hit)*1+1])
  }
  
  if(length(unique(trueDAGid))==2 | !is.null(hit)){
    clrs <- c("darkmagenta", "steelblue4")
  } else{
      clrs <- "steelblue4"
    }
  
  
  p <- ggplot(data)
  if(!is.null(hit)){
    p <- p + geom_point(aes(x, y, colour=hit), alpha=alpha,  size=sizePts, shape=shapePts)
  } else{
    p <- p + geom_point(aes(x, y, colour=id), alpha=alpha,  size=sizePts, shape=shapePts) 
  }
  # make numeric so that it orders it int he right way and
  # you can identify them on the plot
  p <- p + facet_wrap(as.numeric(dataset2)~., scales="free")
  p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
                 axis.text=element_blank(), axis.ticks=element_blank())
  p <- p + theme(axis.title.x=element_blank(),  axis.title.y=element_blank())
  if(!legend) p <- p + theme(legend.position= "none")
  if(!bg) p <- p + theme(panel.background = element_blank())
   p <- p + scale_color_manual(values=clrs)
  print(p)
  return(p)
}

plotDataList <- function(dataList){
  data <- melt(dataList$xs)
  colnames(data) <- c("numData","var","val","dataset")
  vars <- unique(data$var)
  data <- cast(data, dataset+numData~var, value="val")
  
  combs <- combn(vars, 2)
  
  dataL <- lapply(1:ncol(combs), function(i){
    # i <- 1
    res <- data[,c("dataset","numData")]
    res$V1 <- combs[1,i]
    res$V2 <- combs[2,i]
    res <- cbind(res, data[,combs[,i]+2])
    colnames(res) <- c("dataset","numData","V1","V2","val1","val2")
    return(res)
  })
  data <- do.call(rbind, dataL)
  data$combi <- paste(data$V1, data$V2, sep="-")
  
  p <- ggplot(data)
  p <- p + geom_point(aes(val1, val2, colour=combi), alpha=0.2, size=0.5, shape=4)
  p <- p + facet_wrap(dataset+combi~., scales="free")
  p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
                 axis.text=element_blank(), axis.ticks=element_blank())
  p <- p + scale_color_brewer(type = 'qual', palette = "Dark2", direction = 1)
  print(p)
  return(p)
}



plotPairsList_alphaProg <- function(dataList, reps=1:10){
  data <- melt(dataList$xs)
  colnames(data) <- c("numData","xY","val","dataset")
  data <- cast(data, dataset+numData~xY, value="val")
  colnames(data) <- c("dataset","numData","x","y")
  
  aux <- strsplit( sapply(strsplit(data$dataset, "\\."), function(el) el[length(el)]), "_")
  data$alpha <- as.numeric(sapply(aux, function(el) el[2]))
  data$rep <- as.numeric(sapply(aux, function(el) el[1]))
  # table(data$alpha); table(data$rep)
  
  # scatter
  sizePts=0.5; shapePts=4; alpha=0.2
  p <- ggplot(data[which(data$rep %in% reps),])
  p <- p + geom_point(aes(x, y), colour="steelblue4", alpha=alpha,  size=sizePts, shape=shapePts)
  p <- p + facet_grid(alpha~rep, scales="free")
  #p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
  #p <- p + theme(axis.text=element_blank(), axis.ticks=element_blank())
  p <- p + theme(axis.title.x=element_blank(),  axis.title.y=element_blank())
  p <- p + theme(legend.position= "none")
  p <- p + theme(panel.background = element_blank())
  #p <- p + scale_color_manual(values=clrs)
  print(p)
}

plotPairsList_alphaProg2 <- function(dataList, varRep, varProg, reps, alphas){
  data <- melt(dataList$xs)
  colnames(data) <- c("numData","xY","val","dataset")
  data <- cast(data, dataset+numData~xY, value="val")
  colnames(data) <- c("dataset","numData","x","y")
  
  aux <- strsplit(data$dataset, "\\.")
  indxRep <- match(varRep, aux[[1]])+1
  indxProg <- match(varProg, aux[[1]])+1
  data$alpha <- as.numeric(sapply(aux, function(el) el[indxProg]))
  data$rep <- as.numeric(sapply(aux, function(el) el[indxRep]))
  # table(data$alpha); table(data$rep)
  
  # scatter
  sizePts=0.5; shapePts=4; alpha=0.2
  p <- ggplot(data[which(data$rep %in% reps & data$alpha %in% alphas),])
  p <- p + geom_point(aes(x, y), colour="steelblue4", alpha=alpha,  size=sizePts, shape=shapePts)
  if(length(reps)==1){
    p <- p + facet_wrap(alpha~., scales="free")
  } else{
    p <- p + facet_grid(alpha~rep, scales="free")
  }
  p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
  p <- p + theme(axis.text=element_blank(), axis.ticks=element_blank())
  p <- p + theme(axis.title.x=element_blank(),  axis.title.y=element_blank())
  p <- p + theme(legend.position= "none")
  p <- p + theme(panel.background = element_blank())
  #p <- p + scale_color_manual(values=clrs)
  print(p)
}




addLinDesc <- function(dataList){
  
  desc <- sapply(dataList$xs, function(x){
    #x <- dataList$xs[[1]]
    print("******************")
    
    lin <- dhsic.test(x[,"x"],residuals(lm(y~x, data=as.data.frame(x))), method="gamma")$p.value
    krr1_1 <- setParams(krr1, trainData=constructData(x=x[,"x",drop=F], y=x[,"y"]), mc_cores=5)  
    krr1_1 <- krr1_1$learn(krr1_1)
    add <- dhsic.test(x[,"x"],krr1_1$resids(krr1_1, krr1_1$predict(krr1_1, constructData(x=x[,"x",drop=F], y=x[,"y"])))[,"resid"], method="gamma")$p.value
    res <- c(lin=log(lin,10), add=log(add,10))
    plot(x[,"x"],x[,"y"], 
         main=paste("lin: ", round(res["lin"],2),"add: ", round(res["add"],2) ))
    print(res)
    return(res)
  }, simplify="array")
  desc <- t(desc)
  names(dimnames(desc)) <- c("nm","msr")
  dataList$desc <- desc
  
  return(dataList)
}


# Simulate Time Series from given SEM



# simulate a directed graph
simDG <- function(px, py=px, pDiag, pOffDiag){
  ps <- c(pOffDiag, pDiag)[diag(nrow=px, ncol=py)+1]
  dg <- matrix(rbernoulli(px*py, ps)*1, px, py)
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

####################################*
# Two functions to get time dags
# A- simulatate a time dag
# B- create from specific matrices
####################################*
# A- sim time dag
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

# B- create according to user input
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
