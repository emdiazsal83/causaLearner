# obtaining the Markov equivalence class of a DAG

remove(list=ls())

library(pcalg)
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner/utilities/func_dagStuff.R")

# pcalg examples - first example

## generate a random DAG
set.seed(42)
g <- randomDAG(3, 0.4)
plot(g)

g01 <- 1*(as(g,"matrix") > 0) # 0/1-version of adjacency matrix
print.table(g01, zero.=".")


p <- 3
# v-structure
dagVmat <- matrix(c(0,1,1,0,0,0,0,0,0), p, p)
plot(getGraph(dagVmat))
# non-vstructure
dagNVmat <- matrix(c(0,0,0,1,0,0,0,1,0),p,p)
plot(getGraph(dagNVmat))


# Get Markov equivalence graphs in same class as dagV
dagV <- as(dagVmat, "graphNEL")
cpdagV <- dag2cpdag(dagV)
plot(cpdagV)
cpdagVmat <- as(cpdagV, "matrix")
MEdagsV <- pdag2allDags(t(cpdagVmat))$dags
MEdagsV <- lapply(1:nrow(MEdagsV),  function(i) matrix(MEdagsV[i,], p, p, byrow=F))
plot(getGraph(MEdagsV[[1]]))
# only graph in markov equivalence class is v-structure as it should be

gm <- cpdagVmat

pdag2allDags <- function (gm) {
  nodeNms <- colnames(gm)
  p <- ncol(gm)
  rownames(gm) <- colnames(gm) <- as.character(1:p)
  res <- allDags.internal(gm = gm, a = gm, tmp = NULL)
  list(dags = res, nodeNms = nodeNms)
}

# A sink is a node with no arrows going out.
# If rows has from and cols to then rowSums must equal zero
# If rows has to and cols from then colSums must equal zero
# however since we input a cpdag we first delete undirected edges

find.sink <- function (gm){
  gm1 <- gm
  arrow <- gm1 == 1
  bothArrowsOrNone <- gm1 == t(gm1)
  undirected <- bothArrowsOrNone & arrow
  gm1[undirected] <- 0
  which(colSums(gm1) == 0)
}

allDags.internal <- function (gm, a, tmp){
  if (sum(a) == 0) {
    tmp2 <- rbind(tmp, c(t(gm)))
    if (all(!duplicated(tmp2))) 
      tmp <- tmp2
  }
  else {
    sinks <- find.sink(a)
    for (x in sinks) {
      gm2 <- gm
      a2 <- a
      if (adj.check(a, x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(rownames(a)[inc.to.x])
          real.x <- as.numeric(rownames(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x, -x]
        tmp <- allDags.internal(gm2, a2, tmp, verbose)
      }
    }
  }
  tmp
}

# Get Markov equivalence graphs in same class as dagNV
dagNV <- as(dagNVmat, "graphNEL")
cpdagNV <- dag2cpdag(dagNV)
plot(cpdagNV)
cpdagNVmat <- as(cpdagNV, "matrix")
MEdagsNV <- pdag2allDags(t(cpdagNVmat))$dags
MEdagsNV <- lapply(1:nrow(MEdagsNV),  function(i) matrix(MEdagsNV[i,], p, p, byrow=F))
plot(getGraph(MEdagsNV[[3]]))
# theres a v-structure (where a 1<-2->3 structure shd be)
# in markov equivalence class! what am I doing wrong???




#For all the dags
# generate all the DAGs with m-nodes and split into their respective Markov Equivalence Classes
p <- 3
dagsMat <- genAllDAGS(p)
dim(dagsMat)

dagsMatList <- lapply(1:dim(dagsMat)[3], function(i) dagsMat[,,i])
MEList <- lapply(1:length(dagsMatList), function(i){
  # i <- 7
  print("************************")
  print(paste("i:", i))
  dagMat <- dagsMatList[[i]]
  dag <- as(dagMat, "graphNEL")
  cpdag <- dag2cpdag(dag)
  
  cpdagMat <- as(cpdag, "matrix")
  MEdagTrue <- pdag2allDags(t(cpdagMat))$dags
  
  MEdagTrue <- lapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=F))
  
  plot(dag, main=paste("dag: ",i, sep=""))
  plot(cpdag, main=paste("cpdag:", i, sep=""))
  for(j in 1:length(MEdagTrue)){
    plot(getGraph(MEdagTrue[[j]]), main=paste("genDag: ",i, ", dag:", j))
  }
  
  return(MEdagTrue)
})

uniMEList <- unique(MEList)

length(uniMEList) # 11 markov equivalence classes
sapply(uniMEList, function(arr) length(arr))
sum(sapply(uniMEList, function(arr) length(arr)))
uniMEList

getMarkovEquivClass <- function(dagMat){
  dag <- as(dagMat, "graphNEL")
  cpdag <- dag2cpdag(dag)
  cpdagMat <- as(cpdag, "matrix")
  MEdagTrue <- pdag2allDags(t(cpdagMat))$dags
  MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=F), simplify="array")
  dgNms <- getHypID(MEdagTrue)$id
  dimnames(MEdagTrue) <- list(from=dimnames(dagMat)[[1]], to=dimnames(dagMat)[[2]], dag=dgNms)
  return(MEdagTrue)
}

rownames(dagVmat) <- colnames(dagVmat) <- rownames(dagNVmat) <- colnames(dagNVmat) <- c("x", "y", "z")
getMarkovEquivClass(dagVmat); dim(getMarkovEquivClass(dagVmat)); dimnames(getMarkovEquivClass(dagVmat)) 
getMarkovEquivClass(dagNVmat); dim(getMarkovEquivClass(dagNVmat)); dimnames(getMarkovEquivClass(dagNVmat))
