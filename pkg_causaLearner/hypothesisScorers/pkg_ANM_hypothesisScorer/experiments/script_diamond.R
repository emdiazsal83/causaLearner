


remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R")
source("./functions_v1.R")


# number of dags as a function of number of nodes

plot(log(count.graphs(type="all-dags", nodes=1:9)))
count.graphs(type="all-dags", nodes=1:8)


dag <- unifDAG(4)
plot(dag)
class(dag)
nodes(dag)
degree(dag)
adj(dag,c("1","4"))
acc(dag,c("2","3"))


plot(as.bn(dag))
amat(as.bn(dag))
getGraph(amat(as.bn(dag)))
class(getGraph(amat(as.bn(dag))))




allDags4 <- genAllDAGS(4)
#check they are all really DAGS
all(apply(allDags4, 3, isValidGraph, type="dag"))

# plot a random one
allDags4[,,127]
plot(getGraph(allDags4[,,127]))

# find and plot diamond dag
edgL <- vector("list", length=4)
names(edgL) <- as.character(1:4)
edgL[[1]] <- list(edges=c(2,3))
edgL[[2]] <- list(edges=c(4))
edgL[[3]] <- list(edges=c(4))
diamond <- graphNEL(nodes=as.character(1:4), edgeL=edgL, edgemode="directed")
plot(diamond)
amat(as.bn(diamond))
indx <- which(apply(allDags4, 3, function(g) all(g==amat(as.bn(diamond)))))
allDags4[,,indx]
plot(getGraph(allDags4[,,indx]))
dimnames(allDags4) <- list(from=nodes(diamond), to=nodes(diamond), dag=seq(dim(allDags4)[3]))


# Diamond problem from Hoyer


# SEM: a list with a DAG as an adjacency matrix and a list of functions where last argument is noise

# create diamond SEM
nodes(diamond) <- c("w","x","y","z")
funcs <- list(fw=function(n) n, fx=function(w,n) w^2 + n, fy=function(w, n) 4*sqrt(abs(w))+n, fz=function(x,y,n) 2*sin(x)+2*sin(y) + n)
sem1 <- list(dag=diamond, funcs=funcs)

# simulate from it
n <- 500
set.seed(2)
ns <- matrix(runif(4*n)-0.5,4,n)
ns <- t(ns*c(6,2,2,2))
colnames(ns) <- nodes(diamond)
xs <- simSEM(sem1, ns)




# obtain residuals for a given graph
midpoint <- ceiling(n/2)
xTrain <- xs[1:midpoint,]
xTest <- xs[(midpoint+1):n,]

rs.krr   <- getResiduals(G=diamond, xTrain, xTest, learner=krr, numFolds=5)
rs.qhsic <- getResiduals(diamond, xTrain, xTest, qhsic, 3)
#rs.hsic  <- getResiduals(diamond, xTrain, xTest, hsic, 3)

rs.krr2 <- getResidualsMat(amat(as.bn(diamond)), xTrain, xTest, krr, 5)

par(mfrow=c(2,2))
plot(rs.krr$test[,1],rs.krr2$test[,1])
plot(rs.krr$test[,2],rs.krr2$test[,2])
plot(rs.krr$test[,3],rs.krr2$test[,3])
plot(rs.krr$test[,4],rs.krr2$test[,4])


pairs(rs.krr$test, lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(rs.qhsic$test, lower.panel = panel.smooth, upper.panel = panel.cor)

#pairs(rs.hsic$test)

# fit a dag to the data

# perform independence test on residuals


dhsic.test(rs.krr$test, matrix.input=T, alpha=0.02, pairwise=F) 
dhsic.test(rs.qhsic$test, matrix.input=T, alpha=0.02, pairwise=F) 
#dhsic.test(rs.hsic$test, matrix.input=T, alpha=0.02, pairwise=F) 


# For a collection of graphs given as  3D array of stacked adjacency matrices
# obtain p-value for joint independence of residuals

# dont really need to do so much! just look at unique columns for each node

allDags4 <- genAllDAGS(4)
dimnames(allDags4) <- list(from=nodes(diamond), to=nodes(diamond), dag=seq(dim(allDags4)[3]))

dim(allDags4)


uniqueRegs <- getUniqueRegs(allDags4)

# krr
uniqueResids.krr <- getUniqueResids(uniqueRegs, xTrain, xTest, learner=krr, numFolds=5)
#takes roughly 30mns mins for krr and 500 train + test samples from hoyers diamond
pvals.krr <- getPvals(allDags4, uniqueRegs, uniqueResids.krr)
#takes roughly 66 mins for krr and 500 train + test samples from hoyers diamond

# qhsic
pm <- proc.time()
uniqueResids.qhsic <- getUniqueResids(uniqueRegs, xTrain, xTest, learner=qhsic, numFolds=5)
proc.time()-pm #
pm <- proc.time()
pvals.qhsic <- getPvals(allDags4, uniqueRegs, uniqueResids.qhsic)
proc.time()-pm #


# save(list=ls(), file="pvalsDiamond_krrAndQhisic.RData")
# load(file="pvalsDiamond_krrAndQhisic.RData")

plot(getGraph(allDags4[,,52]))



# selected graphs krr
indx.krr <- which(pvals.krr[,"test"]>0.001)
length(indx.krr)
par(mfrow=c(3,3))
for(i in indx.krr) plot(getGraph(allDags4[,,i]))
pvals.krr[indx.krr,"test"]	

# selected graphs qhsic
indx.qhsic <- which(pvals.qhsic[,"test"]>0.001)
length(indx.qhsic)
par(mfrow=c(3,3))
for(i in indx.qhsic) plot(getGraph(allDags4[,,i]))
pvals.qhsic[indx.qhsic,"test"]	

indx.krrqhsic <- intersect(indx.krr, indx.qhsic)
length(indx.krrqhsic)
par(mfrow=c(2,2))
for(i in indx.krrqhsic) plot(getGraph(allDags4[,,i]))
pvals.qhsic[indx.qhsic,"test"]	

# run hsic for 4 graphs selected by both krr and qhsic

dags4 <- allDags4[,,indx.krrqhsic]
uniqueRegsList4 <- getUniqueRegsList(dags=dags4)

# hsic

pm <- proc.time()
uniqueResidsList.hsic <- getUniqueResidsList(uniqueRegsList4, xTrain, xTest, learner=hsic, numFolds=5)
proc.time() - pm
# 

pvals4.hsic <- getPvalsList(dags=dags4, uniqueRegsList=uniqueRegsList4, uniqueResidsList=uniqueResidsList.hsic)
#




