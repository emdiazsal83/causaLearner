# Try KCDC on Sparc data

remove(list=ls())
server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2MV-ME2-sparc_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]

repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)
print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

library(FNN) #  (in learn.qhsic, myGpCreate)
# FNN sometimes has a problem: cant find C functions perhaps because of some
# unloading of dll libraries: see https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16446

#########################################################################################################################################3


# load sparc data 

dir("./data/sparc")
sparc <- readMat("./data/sparc/SPARC0304.mat")
class(sparc)
names(sparc)
dim(sparc$X) # reflectance data for 62 wavelengths, 135 sample points
dim(sparc$Y) # Chlorophyl content, LAI and fCover for 135 sample points
dim(sparc$WaveLength) # 62 wavelengths for which there is reflectance info

# true causal hypothesis
# Y -> X

plot(sparc$X[,1], sparc$Y[,2])

# true effects - spectra
X <- apply(sparc$X[,c(1), drop=FALSE], 2, stdrize) # c(1,10,21)
apply(X, 2, mean); apply(X, 2, sd)
pairs(X)
summary(as.numeric(X))

# true causes - vegetation parameters
Y <- apply(sparc$Y[,1,drop=FALSE], 2, stdrize)
apply(Y, 2, mean); apply(Y, 2, sd)
pairs(Y)
summary(as.numeric(Y))

# NEED TO change construct data so that it takes matrices for y data!!!

constructData <- function (x, y){
  stopifnot(is.list(x) || is.vector(x) || is.matrix(x))
  stopifnot(is.list(y) || is.vector(y) || is.matrix(y))
  data = list(x = x, y = y)
  class(data) = "CVST.data"
  return(data)
}

getSubset <- function (data, subset){
  stopifnot(class(data) == "CVST.data")
  x = getX(data, subset)
  if(!is.null(dim(data$x))){
    y = getY(data, subset) 
  } else{
    y = data$y[subset]  
  }
  
  ret = constructData(x = x, y = y)
  return(ret)
}

getY <- function (data, subset = NULL){
  stopifnot(class(data) == "CVST.data")
  if (is.null(subset)) {
    ret = data$y
  }
  else {
    if (is.list(data$y) || is.vector(data$y)) {
      ret = data$y[subset]
    }
    else {
      ret = data$y[subset, , drop = FALSE]
    }
  }
  return(ret)
}

trainDataXY <- constructData(X,Y)
trainDataYX <- constructData(Y,X)


# cross validate lambda, sigma 
# train hyperparameters 
pm <- proc.time()
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_log_quad_L2_none, trainData=trainDataXY, plot=FALSE)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_log_quad_L2_none, trainDataYX, plot=FALSE)
proc.time() - pm # 25 secs
#getHyperPar(cmem_rbf_rbf_L2_lambda_kernParsX_xy, "sigmaX")

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$lambda <- 1e-4
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$lambda <- 1e-4

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$degreeX <- 1
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$degreeX <- 1

# train learn parameters
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)


# calculate measures
cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)

# WRONG causal direction inferred with rbf_rbf and quad quad when we use all 62 spectra!!!
# causal direction with log_quad infered correctly for KCDC for X =c(1,10,21) vs Y = Chl, LAI, fCover


# COULD THERE BE A NORMALIZATION PROBLEM WHEN X AND Y HAVE DIFFERENT DIMENSIONALITY???

px <- 3
py_even <- 3
py_skewed <- 12
p_even <- px + py_even
p_skewed <- px + py_skewed
q <- 100
ps_even <- rep(p_even, q)
ps_skewed <- rep(p_skewed, q)
n <- 100
ns <- rep(n, q)
nodes <- list(dist="runif", pars=list(min=-1.5, max=1.5), a=1, b=1)
nodes_even <- rep(list(nodes),p_even)
nodes_skewed <- rep(list(nodes),p_skewed)
names(nodes_even) <- as.character(seq(p_even))
names(nodes_skewed) <- as.character(seq(p_skewed))
nodess_even <- lapply(ps_even, function(p) nodes_even)
nodess_skewed <- lapply(ps_skewed, function(p) nodes_skewed)

dgMatXY_even <- simDG(px, py_even, pDiag=0.5, pOffDiag=0.5)
dagMat_even <- rbind(cbind(matrix(0,px,px), dgMatXY_even), matrix(0, py_even, p_even))
colnames(dagMat_even) <- rownames(dagMat_even) <- seq(p_even)

dgMatXY_skewed <- cbind(dgMatXY_even, dgMatXY_even, dgMatXY_even, dgMatXY_even)
dagMat_skewed <- rbind(cbind(matrix(0,px,px), dgMatXY_skewed), matrix(0, py_skewed, p_skewed))
colnames(dagMat_skewed) <- rownames(dagMat_skewed) <- seq(p_skewed)


dataList_even <- simRandSEMs(q, ps_even, ns, nodess_even, sigma=1.5, dagMat=dagMat_even) #, geU=function(y, nois, scale, constant) (y + nois - constant)/scale
dataList_skewed <- simRandSEMs(q, ps_skewed, ns, nodess_skewed, sigma=1.5, dagMat=dagMat_skewed)

i <- 10

dagMat_even
pairs(dataList_even$xs[[i]])


pm <- proc.time()
msrs <- mcmapply(function(el1, el2, nm){
  # i <- 9; el <- dataTestList$xs[[i]]; nm <- dataTestList$names[i]
  X_even <- apply(el1, 2, stdrize)
  X_skewed <- apply(el2, 2, stdrize)
  
  print(paste("name: ", nm))
  #print("head(X)")
  #print(head(X))
  #print(apply(X, 2, mean))
  #print(apply(X, 2, sd))
  
  x <- X_even[,1:px, drop=F]
  y <- X_even[,(px+1):p_even]
  trainDataXY_even <- constructData(x, y)
  trainDataYX_even <- constructData(y, x)
  
  x <- X_skewed[,1:px, drop=F]
  y <- X_skewed[,(px+1):p_skewed]
  trainDataXY_skewed <- constructData(x, y)
  trainDataYX_skewed <- constructData(y, x)
  
  
  # train hyperparameters
  cmem_rbf_rbf_L2_none_xy <- setParams(learner=cmem_rbf_rbf_L2_none, trainData=trainDataXY_even, plot=FALSE)
  cmem_rbf_rbf_L2_none_yx <- setParams(learner=cmem_rbf_rbf_L2_none, trainDataYX_even, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataXY_even, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataYX_even, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataXY_even, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataYX_even, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataXY_even, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataYX_even, plot=FALSE)
  
  # train learn parameters
  cmem_rbf_rbf_L2_none_xy <- cmem_rbf_rbf_L2_none_xy$learn(cmem_rbf_rbf_L2_none_xy)
  cmem_rbf_rbf_L2_none_yx <- cmem_rbf_rbf_L2_none_yx$learn(cmem_rbf_rbf_L2_none_yx)
  cmem_rbf_rbf_L2_lambda_xy <- cmem_rbf_rbf_L2_lambda_xy$learn(cmem_rbf_rbf_L2_lambda_xy)
  cmem_rbf_rbf_L2_lambda_yx <- cmem_rbf_rbf_L2_lambda_yx$learn(cmem_rbf_rbf_L2_lambda_yx)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  # calculate measures
  msrsFixXY <- cmem_rbf_rbf_L2_none_xy$calcMsrs(cmem_rbf_rbf_L2_none_xy)
  msrsFixYX <- cmem_rbf_rbf_L2_none_yx$calcMsrs(cmem_rbf_rbf_L2_none_yx)
  msrsOpt1XY <- cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
  msrsOpt1YX <- cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
  msrsOpt2XY <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  msrsOpt2YX <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #msrsOpt3XY <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #msrsOpt3YX <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  
  #KCDC lambda fix
  KCDCq <- msrsFixYX["KCDC"]/msrsFixXY["KCDC"] 
  #KCDCrel lambda fix
  KCDCrelq <- msrsFixYX["KCDCrel"]/msrsFixXY["KCDCrel"]
  #KCDC lambda opt1
  KCDCq_opt1 <- msrsOpt1YX["KCDC"]/msrsOpt1XY["KCDC"]
  #KCDCrel lambda opt1
  KCDCrelq_opt1 <- msrsOpt1YX["KCDCrel"]/msrsOpt1XY["KCDCrel"]
  #KCDC lambda opt2
  KCDCq_opt2 <- msrsOpt2YX["KCDC"]/msrsOpt2XY["KCDC"]
  #KCDCrel lambda opt2
  KCDCrelq_opt2 <- msrsOpt2YX["KCDCrel"]/msrsOpt2XY["KCDCrel"]
  #KCDC lambda opt3
  #KCDCq_opt3 <- msrsOpt3YX["KCDC"]/msrsOpt3XY["KCDC"]
  #KCDCrel lambda opt3
  #KCDCrelq_opt3 <- msrsOpt3YX["KCDCrel"]/msrsOpt3XY["KCDCrel"]
  
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2, KCDCq_opt3, KCDCrelq_opt3)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2","KCDCopt3","KCDCrelopt3")
  res_even <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2)
  names(res_even) <- paste(c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2"), "even", sep="")
  
  # train hyperparameters
  cmem_rbf_rbf_L2_none_xy <- setParams(learner=cmem_rbf_rbf_L2_none, trainData=trainDataXY_skewed, plot=FALSE)
  cmem_rbf_rbf_L2_none_yx <- setParams(learner=cmem_rbf_rbf_L2_none, trainDataYX_skewed, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataXY_skewed, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda, trainDataYX_skewed, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataXY_skewed, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsX, trainDataYX_skewed, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataXY_skewed, plot=FALSE)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- setParams(learner=cmem_rbf_rbf_L2_lambda_kernParsXY, trainDataYX_skewed, plot=FALSE)
  
  # train learn parameters
  cmem_rbf_rbf_L2_none_xy <- cmem_rbf_rbf_L2_none_xy$learn(cmem_rbf_rbf_L2_none_xy)
  cmem_rbf_rbf_L2_none_yx <- cmem_rbf_rbf_L2_none_yx$learn(cmem_rbf_rbf_L2_none_yx)
  cmem_rbf_rbf_L2_lambda_xy <- cmem_rbf_rbf_L2_lambda_xy$learn(cmem_rbf_rbf_L2_lambda_xy)
  cmem_rbf_rbf_L2_lambda_yx <- cmem_rbf_rbf_L2_lambda_yx$learn(cmem_rbf_rbf_L2_lambda_yx)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_xy <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #cmem_rbf_rbf_L2_lambda_kernParsXY_yx <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  # calculate measures
  msrsFixXY <- cmem_rbf_rbf_L2_none_xy$calcMsrs(cmem_rbf_rbf_L2_none_xy)
  msrsFixYX <- cmem_rbf_rbf_L2_none_yx$calcMsrs(cmem_rbf_rbf_L2_none_yx)
  msrsOpt1XY <- cmem_rbf_rbf_L2_lambda_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_xy)
  msrsOpt1YX <- cmem_rbf_rbf_L2_lambda_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_yx)
  msrsOpt2XY <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  msrsOpt2YX <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  #msrsOpt3XY <- cmem_rbf_rbf_L2_lambda_kernParsXY_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_xy)
  #msrsOpt3YX <- cmem_rbf_rbf_L2_lambda_kernParsXY_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsXY_yx)
  
  #KCDC lambda fix
  KCDCq <- msrsFixYX["KCDC"]/msrsFixXY["KCDC"] 
  #KCDCrel lambda fix
  KCDCrelq <- msrsFixYX["KCDCrel"]/msrsFixXY["KCDCrel"]
  #KCDC lambda opt1
  KCDCq_opt1 <- msrsOpt1YX["KCDC"]/msrsOpt1XY["KCDC"]
  #KCDCrel lambda opt1
  KCDCrelq_opt1 <- msrsOpt1YX["KCDCrel"]/msrsOpt1XY["KCDCrel"]
  #KCDC lambda opt2
  KCDCq_opt2 <- msrsOpt2YX["KCDC"]/msrsOpt2XY["KCDC"]
  #KCDCrel lambda opt2
  KCDCrelq_opt2 <- msrsOpt2YX["KCDCrel"]/msrsOpt2XY["KCDCrel"]
  #KCDC lambda opt3
  #KCDCq_opt3 <- msrsOpt3YX["KCDC"]/msrsOpt3XY["KCDC"]
  #KCDCrel lambda opt3
  #KCDCrelq_opt3 <- msrsOpt3YX["KCDCrel"]/msrsOpt3XY["KCDCrel"]
  
  #res <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2, KCDCq_opt3, KCDCrelq_opt3)
  #names(res) <- c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2","KCDCopt3","KCDCrelopt3")
  res_skewed <- c(KCDCq, KCDCrelq, KCDCq_opt1, KCDCrelq_opt1,KCDCq_opt2, KCDCrelq_opt2)
  names(res_skewed) <- paste(c("KCDC","KCDCrel","KCDCopt1","KCDCrelopt1","KCDCopt2","KCDCrelopt2"), "skewed", sep="_")
  
  res <- c(res_even, res_skewed)
  
  return(res)
}, el1=dataList_even$xs, el2=dataList_skewed$xs, nm=dataList$names, 
SIMPLIFY="array", mc.cores=1)
proc.time() - pm # 100 mins


pctRight <- function(x){
  
  res <- sum(x>1)/length(x)*100
  return(res)
}

names(dimnames(msrs)) <- c("measure", "database")
msrsDB <- melt(msrs)

save("msrs", file="./experiments/sparcTests.RData")

cast(msrsDB, measure~., value="value", fun.aggregate="pctRight")


