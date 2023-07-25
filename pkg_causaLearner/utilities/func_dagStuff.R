# dag utility functions

library(igraph) #V(), graph in get_n_cycles_directed_E

# Generate all dags implied in an undirected graph
udag2dags <- function(uDAG, checkValid=TRUE){
  numEdges <- sum(uDAG)/2
  numDags <- 2^numEdges
  numNodes <- nrow(uDAG)
  V <- colnames(uDAG)
  
  aux <- uDAG  
  aux[upper.tri(aux)] <- 0
  aux <- melt(aux, stringsAsFactors=F)
  aux <- aux[which(aux$value>0),1:2]
  aux$X1 <- as.character(aux$X1)
  aux$X2 <- as.character(aux$X2)
  
  
  lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
  upperOnes <- (lowerOnes==0)*1
  
  lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
  upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
  
  matOnes <- mapply(FUN=function(a,b){
    aux1 <- a
    colnames(aux1) <- c("x1","x2")
    aux2 <- b
    colnames(aux2) <- c("x2","x1")
    base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
  
  res <- matrix(0, numNodes, numNodes)
  rownames(res) <- V
  colnames(res) <- V
  
  dags <- sapply(matOnes, function(mat){
    # mat <- matOnes[[1]]
    mat2 <- mat[which(mat[,3]==1),]
    indxRow <- match(mat2[,1], V)
    indxCol <- match(mat2[,2], V)
    indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
    
    res2 <- res
    res2[indxMat] <-   1
    return(res2)
    
  }, simplify="array")
  
  # check for and delete graphs with cycles
  
  if(checkValid){
    indxValid <- which(apply(dags, c(3), function(mat){ 
      isValidGraph(mat, "dag")
    }))
    dags <- dags[,,indxValid]
  }
  
  return(dags)
  
}

# Generate all dags implied in a partially directed graph: ¿is this not thes same as finding the
# markov equivalence graph and getting all the graphs in that? (there are functions for that in library(pcalg) -> pdag2allDags & dag2cpdag)
# basically permute undirected edges
pdag2dags <- function(pDAG){
  
  # undirected edges adjacency matrix
  uDAG <- (pDAG==1 & t(pDAG)==1)*1
  
  # number of undirected edges
  numUEdges <- sum(uDAG)/2
  
  # number of dags we can get from cpdag
  numDags <- 2^numUEdges
  numNodes <- nrow(pDAG)
  V <- colnames(pDAG)
  
  # if there is at least one undirected edge 
  # we permute those edges
  
  if(numDags>1){
    
    aux <- uDAG  
    aux[upper.tri(aux)] <- 0
    aux <- melt(aux, stringsAsFactors=F)
    aux <- aux[which(aux$value>0),1:2]
    aux$X1 <- as.character(aux$X1)
    aux$X2 <- as.character(aux$X2)
    
    
    lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
    upperOnes <- (lowerOnes==0)*1
    
    lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
    upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
    
    matOnes <- mapply(FUN=function(a,b){
      aux1 <- a
      colnames(aux1) <- c("x1","x2")
      aux2 <- b
      colnames(aux2) <- c("x2","x1")
      base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
    
    res <- pDAG - uDAG
    
    
    dags <- sapply(matOnes, function(mat){
      # mat <- matOnes[[1]]
      mat2 <- mat[which(mat[,3]==1),]
      indxRow <- match(mat2[,1], V)
      indxCol <- match(mat2[,2], V)
      indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
      
      res2 <- res
      res2[indxMat] <-   1
      return(res2)
      
    }, simplify="array")
  } else{
    dags <- pDAG
    dim(dags) <- c(dim(dags),1)
  }
  
  
  
  return(dags)
  
}

# generate all dags implied in a "minimal dag": basically all possible ways of adding edges without creating cycles. 
mindag2dags <- function(minDAG, checkValid=TRUE){
  
  noLinks <- (minDAG==0 & t(minDAG)==0)*1
  diag(noLinks) <- 0
  
  numEdges <- sum(noLinks)/2
  numDags <- 3^numEdges
  numNodes <- nrow(minDAG)
  V <- colnames(minDAG)
  
  # if there is at least one pair of nodes with no edges
  # we permute those
  if(numDags>1){
    
    aux <- noLinks  
    aux[upper.tri(aux)] <- 0
    aux <- melt(aux, stringsAsFactors=F)
    aux <- aux[which(aux$value>0),1:2]
    aux$X1 <- as.character(aux$X1)
    aux$X2 <- as.character(aux$X2)
    
    
    lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
    lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
    
    
    res <-  minDAG*0
    
    noLinks2 <- sapply(lowerOnes, function(mat){
      # mat <- matOnes[[1]]
      mat2 <- mat[which(mat[,3]==1),]
      indxRow <- match(mat2[,1], V)
      indxCol <- match(mat2[,2], V)
      indxMat1 <- matrix(c(indxRow, indxCol), length(indxRow), 2)
      indxMat2 <- matrix(c(indxCol, indxRow), length(indxRow), 2)
      
      res2 <- res
      res2[indxMat1] <-   1
      res2[indxMat2] <-   1
      return(res2)
      
    }, simplify="array")
    
    dagss <- sapply(1:(dim(noLinks2)[3]), function(i){
      noLinks <- noLinks2[,,i]
      if(sum(noLinks)==0){
        dags <- minDAG
        dim(dags) <- c(dim(dags),1)
      } else{
        aux <- noLinks  
        aux[upper.tri(aux)] <- 0
        aux <- melt(aux, stringsAsFactors=F)
        aux <- aux[which(aux$value>0),1:2]
        aux$X1 <- as.character(aux$X1)
        aux$X2 <- as.character(aux$X2)
        
        
        lowerOnes <- as.matrix(do.call(expand.grid, lapply(1:nrow(aux), function(i) c(0,1))))
        upperOnes <- (lowerOnes==0)*1
        
        lowerOnes <-  apply(lowerOnes, 1, function(row) cbind(aux, row))
        upperOnes <-  apply(upperOnes, 1, function(row) cbind(aux, row))
        
        matOnes <- mapply(FUN=function(a,b){
          aux1 <- a
          colnames(aux1) <- c("x1","x2")
          aux2 <- b
          colnames(aux2) <- c("x2","x1")
          base:::rbind(aux1, aux2[,c(2,1,3)], stringsAsFactors=F)}, a=lowerOnes, b=upperOnes, SIMPLIFY=F)
        
        res <- minDAG
        
        dags <- sapply(matOnes, function(mat){
          # mat <- matOnes[[1]]
          mat2 <- mat[which(mat[,3]==1),]
          indxRow <- match(mat2[,1], V)
          indxCol <- match(mat2[,2], V)
          indxMat <- matrix(c(indxRow, indxCol), length(indxRow), 2)
          
          res2 <- res
          res2[indxMat] <-   1
          return(res2)
          
        }, simplify="array")
      }
      return(dags)
      
    }, simplify="array")
    
    dagss <- do.call(what=abind, args=c(dagss, along=3) )
    
    if(checkValid){
      indxValid <- which(apply(dagss, c(3), function(mat){ 
        isValidGraph(mat, "dag")
      }))
      dagss <- dagss[,,indxValid]
    }
  } else{
    dagss <- minDAG
    dim(dagss) <- c(dim(dagss),1)
  }
  
  return(dagss)
  
}

#test.x <- 10^seq(-20,20); ctch <- sapply(test.x, function(x) tryCatch(intToBits(x-1), warning=function(w) w)); rslt <- sapply(ctch, function(el) class(el)[1]); data.frame(test.x, rslt)
# intToBitss only works with x < 10^10 so, for dags in p >= 6 we have problems
intToBitss <- function(x, maxInt) paste(rev(strsplit(substr(paste(as.integer(intToBits(x-1)), collapse=""),1, ceiling(log(maxInt,2))), split="")[[1]]), collapse="")

#x <- 6; maxInt <- 6000
# intToBitss(x,maxInt)
# nchar(intToBitss(x,maxInt))

BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+1

dagIDtoMatrix <- function(dagID, p){
  if(p<=5){
    dagBin <- intToBitss(x=dagID, maxInt=2^(p^2))
    dag <- matrix(as.numeric(sapply(1:nchar(dagBin), function(i) substr(dagBin, i, i))), p , p)
  } else{
    dag <- dagIDtoMatrix_large(dagID, p)
  }
  return(dag)
}

dagIDtoMatrix_large <- function(dagID, p){
  dag <- rep(0, p^2)
  indx <- as.numeric(strsplit(as.character(dagID), "0")[[1]])
  dag[indx] <- 1
  dag <- matrix(dag, p, p)
  return(dag)
}


# obtain binary id string for set of dags. allHyps is a 3d array, diff hyps along 3rd dim
getHypID <- function(allHyps){
  
  if(dim(allHyps)[1]>5){
    hypIDs <- getHypID_large(allHyps)
  } else{
    hypIDs <- data.frame(bin=apply(allHyps, 3, function(mat) paste(as.character(as.numeric(mat)), collapse="")))
    hypIDs$id <- sapply(hypIDs$bin, BinToDec)
    hypIDs <- hypIDs[,c("id","bin")]
  }
  
  return(hypIDs)
}  

# for larger dags lets j
getHypID_large <- function(allHyps){
  
  hypIDs <- data.frame(bin=apply(allHyps, 3, function(mat) paste(which(c(mat)==1), collapse=".")), stringsAsFactors = F)
  hypIDs$id <- sapply(hypIDs$bin, function(el) as.numeric(paste(strsplit(el, "\\.")[[1]], collapse="0")))
  hypIDs <- hypIDs[,c("id","bin")]
  
  return(hypIDs)
}

# get markov equivalent class corresponding to a certain dag
getMarkovEquivClass <- function(dagMat){
  p <- nrow(dagMat)
  dag <- as(dagMat, "graphNEL")
  cpdag <- dag2cpdag(dag)
  cpdagMat <- as(cpdag, "matrix")
  MEdagTrue <- pdag2allDags(t(cpdagMat))$dags
  MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=F), simplify="array")
  dgNms <- getHypID(MEdagTrue)$id
  dimnames(MEdagTrue) <- list(from=dimnames(dagMat)[[1]], to=dimnames(dagMat)[[2]], dag=dgNms)
  return(MEdagTrue)
}


# generate all the DAGS with m-nodes
genAllDAGS <- function(m){
  # possible topological orderings
  topOrd <- permutations(m,m)
  numTopOrds <- nrow(topOrd)
  
  # obtain all permutations of DAG matrix for one topological ordering
  maxEdges <- m*(m-1)/2
  numEdgePerms = 2^maxEdges
  
  # create binary strings indicating whether a certain directed edge is switched on or not for all possible configs
  edgePerms <- sapply(0:(numEdgePerms-1), function(i) as.numeric(strsplit(substr(paste(as.integer(intToBits(i)), collapse=""),1,maxEdges),"")[[1]]))
  
  # now we put in adjacency matrix form by placing in 
  dags1topOrd <- array(0, dim=c(m,m,numEdgePerms))
  indMat <- matrix(1:m^2,m,m)
  dags1topOrd[array(upper.tri(indMat), dim=c(m,m,numEdgePerms))] <- as.numeric(edgePerms)
  
  # we now create for all topological orders by permuting
  dagsAlltopOrd <- sapply(1:numTopOrds, function(i){
    perm <-  topOrd[i,]
    aux <- dags1topOrd[,perm,]
    return(aux[perm,,])
  }, simplify="array")
  # stack adjacency matrices in one direction
  
  dim(dagsAlltopOrd) <- c(m,m,numEdgePerms*numTopOrds) 
  
  # remove duplicates
  dagsAlltopOrd <- unique(dagsAlltopOrd, MARGIN=3)
  return(dagsAlltopOrd)
  
  
}



# DAG distances 
edgeDist <- function(dagChosen, dagTrue){
  difs <- dagChosen - dagTrue
  difs[which(difs<0)] <- 0
  res <- sum(difs)
  return(res)
}    
nonEdgeDist <- function(dagChosen, dagTrue){
  difs <- dagTrue - dagChosen
  difs[which(difs<0)] <- 0
  res <- sum(difs)
  return(res)
}
totalEdgeDist <- function(dagChosen, dagTrue) sum(abs(dagTrue -  dagChosen))



# normalized structural hamming distance (SHD)

SHD <- function(dagChosen, dagTrue){
  # p <- 4; dagChosen <- unifDAG(p); dagTrue <- unifDAG(p)
  pcalg:::shd(as(dagChosen, "graphNEL"),as(dagTrue, "graphNEL"))
}

# Functions taken from github jakob-fischer/jrnf_R_tools to identify cycles in a directed graph

 # tools folder. Converts an igraph graph to an adjacency matrix...
graph_to_amatrix <- function(g) {
  a <- get.adjacency(g)
  return(matrix(as.numeric(a), ncol=ncol(a)))
}
# tools.R file 
# This function removes all the duplicates / permutations if one checks for 
# subisomorphisms with a directed ring with the igraph method 
# "graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))".
# This is done by simply requiring that the first edge in the sumisomorphism
# is the one with the lowest id.
subisomorphism_rm_permutation <- function(si) {
  # method that identifies the right subisomorphisms to keep
  is_first_min <- function(x)  {  return(x[1] == min(x))  }
  sel <- lapply(si, is_first_min)    # select by list apply
  return(si[unlist(sel)])            # return selected elements from list si
} 
# cyles.R file

# The function calculates the number of cycles of length <n> in the directed 
# graph <g> and the number of occurences of each edge in this cycles by using
# the graph subisomorphism method of igraph. Because the subisomorphism is 
# defined in terms of nodes, multiple edges between two nodes don't imply the
# cycle being counted multiple times. The suffix "_V" indicates thus that the
# cycles are identified in respect to the vertices. If <list_cycles> is set
# then the third part of the returned list is a matrix for which each row 
# indicates the vertices taking part in one cycle.
#
# This function returns the results in a list: 
# The first element is the <count> of cycles followed by a incidence vector <v_incidence>
# that describes how many cycles the corresponding nodes are part of. If <list_cycles> is
# used a matrix <cycles> (each row corresponding to one cycle) and a vector <cycles_m>
# describing how often the respective cycles occur are also returned.
get_n_cycles_directed_E <- function(g, n, list_cycles=F) {
  gc()                              # Try using garbace gollector because function tends to segfault (igraph)
  v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)
  
  if(n == 1) {
    loops <- which(is.loop(g))        # which edges are loops?
    
    for(i in loops) {                 
      v <- ends(g, i, names=F)[1,1]           # For each, the first node is the only one
      v_count[v] <- v_count[v] + 1   # 
    }
    
    if(list_cycles) {
      cycles <- matrix(0, ncol=length(V(g)), nrow=sum(v_count != 0))  
      cycles_m <- rep(0, nrow(cycles))
      
      u_loops_s <- which(v_count != 0)
      cycles <- matrix(u_loops_s, ncol=1)  
      cycles_m <- rep(1, length(u_loops_s))
    }
  } else {
    m <- graph_to_amatrix(g)                               # calculate graph adjacency matrix
    g_s <- simplify(g, remove.multiple=T, remove.loops=T)  # remove multiple edges and loops
    
    # find subisomorphisms and count occurence of nodes
    si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
    si <- subisomorphism_rm_permutation(si)  # remove permutations
    
    cycles_m <- c()    
    
    for(i in si) { 
      k <- 1  # number of variants to traverse the multiple edges for this subisomorphism / cycle
      
      for(j in 1:length(i))    # just multiply all the values of the adjacency matrix in the cycle 
        if(j != length(i))            
          k <- k * m[i[j], i[j+1]]
        else
          k <- k * m[i[length(i)], i[1]]    # edge between first and last node
        
        v_count[i] <- v_count[i] + k    # accumulate number of cycles individual nodes are part of
        cycles_m <- c(cycles_m, k)
    }
    
    if(list_cycles) {
      if(length(si) != 0)
        cycles <- do.call(rbind, si)
      else
        cycles <- matrix(0, ncol=n, nrow=0)   
    }  
  }
  
  if(list_cycles)
    return(list(count=sum(v_count)/n, v_incidence=v_count, cycles=cycles, c_mul=cycles_m))
  else
    return(list(count=sum(v_count)/n, v_incidence=v_count))
}
# get_n_cycles_directed_E fixes the problem of not counting cycles with multiple
# edges multiple times by calculating their effect on the number of cycles after 
# the subisomorphism algorithm has been invoked.
#
# The return format is equivalent to get_n_cycles_directed_V
get_n_cycles_directed_V <- function(g, n, list_cycles=F) { 
  gc()                              # Try using garbace gollector because function tends to segfault (igraph)
  v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)
  
  # if n==1 just count loops
  if(n == 1) {
    loops <- which(is.loop(g))        # which edges are loops?
    
    for(i in loops) {                 
      v <- ends(g, i, names=F)[1,1]           # For each, the first node is the only one
      v_count[v] <- 1                # Count only once for each node (to be consistent with n>1) 
    }
    
    # Listing cycles is simple here, just find the vertices that have nonzero
    # loops (v_count != 0) and add one row for each of them in cycles matrix.
    if(list_cycles) {
      u_loops_s <- which(v_count != 0)
      cycles <- matrix(u_loops_s, ncol=1)  
      cycles_m <- rep(1, length(u_loops_s))
    }
    # handles the general case of n > 1
  } else {
    # find subisomorphisms of directed graph with ring of length n
    si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
    si <- subisomorphism_rm_permutation(si)  # remove permutations
    
    for(i in si)
      v_count[i] <- v_count[i] + 1
    
    if(list_cycles) {
      if(length(si) != 0)
        cycles <- do.call(rbind, si)
      else
        cycles <- matrix(0, ncol=n, nrow=0)    
      
      cycles_m <- rep(1, nrow(cycles))
    }
  }
  
  if(list_cycles)
    return(list(count=sum(v_count)/n, v_incidence=v_count, cycles=cycles, c_mul=cycles_m))
  else
    return(list(count=sum(v_count)/n, v_incidence=v_count))
}






