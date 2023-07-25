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

# WAVELENGTHS

# Infra-reds
# FIR               ~               25/40-200/350 microns
# MIR               ~                     5-25/40 microns
# NIR               ~ 700/1000-5000    nm .7/1-5  microns
# http://www.icc.dur.ac.uk/~tt/Lectures/Galaxies/Images/Infrared/Regions/irregions.html

# Red 	            ~ 700–635 nm 	~ 430–480 THz
# Orange 	          ~ 635–590 nm 	~ 480–510 THz
# Yellow 	          ~ 590–560 nm 	~ 510–540 THz
# Green 	          ~ 560–520 nm 	~ 540–580 THz
# Cyan            	~ 520–490 nm 	~ 580–610 THz
# Blue 	            ~ 490–450 nm 	~ 610–670 THz
# Violet or Purple 	~ 450–400 nm
# https://turing.manhattan.edu/~mrifkind01/Homework_2/Colors.html

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
X <- apply(sparc$X[,c(25), drop=FALSE], 2, norml) # c(1,10,21)
apply(X, 2, mean); apply(X, 2, sd)
pairs(X)
summary(as.numeric(X))

# true causes - vegetation parameters
Y <- apply(sparc$Y[,c(3),drop=FALSE], 2, norml) #stdrize, norml
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

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$lambda <- 1e-8
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$lambda <- 1e-8

cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$degreeX <- 0.95
cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$degreeX <- 0.95

# train learn parameters
cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)



# calculate measures
(KCDCxy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy))
(KCDCyx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx))

KCDCyx/KCDCxy

# WRONG causal direction inferred with rbf_rbf and quad quad when we use all 62 spectra!!!
# causal direction with log_quad infered correctly for KCDC for X =c(1,10,21) vs Y = Chl, LAI, fCover

pm <- proc.time()
KCDCs <- sapply(1:ncol(sparc$X), function(i) sapply(1:ncol(sparc$Y), function(j){
  # i <- 1; j <- 1
  print(paste("i-j=", paste(i,j, sep="-")))
  # true effects - spectra
  X <- apply(sparc$X[,i, drop=FALSE], 2, norml) 
  # true causes - vegetation parameters
  Y <- apply(sparc$Y[,j,drop=FALSE], 2, norml)
  trainDataXY <- constructData(X,Y)
  trainDataYX <- constructData(Y,X)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- setParams(learner=cmem_log_quad_L2_none, trainData=trainDataXY, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- setParams(learner=cmem_log_quad_L2_none, trainDataYX, plot=FALSE)
  cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$lambda <- c(1e-4, 1e-9, 1e-8)[j]
  cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$lambda <- c(1e-4, 1e-9, 1e-8)[j]
  cmem_rbf_rbf_L2_lambda_kernParsX_xy$hyperParams$data$non_optimizable$degreeX <- c(1, 0.98, 0.95)[j]
  cmem_rbf_rbf_L2_lambda_kernParsX_yx$hyperParams$data$non_optimizable$degreeX <- c(1, 0.98, 0.95)[j]
  # train learn parameters
  cmem_rbf_rbf_L2_lambda_kernParsX_xy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$learn(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  cmem_rbf_rbf_L2_lambda_kernParsX_yx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$learn(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  # calculate measures
  KCDCxy <- cmem_rbf_rbf_L2_lambda_kernParsX_xy$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_xy)
  KCDCyx <- cmem_rbf_rbf_L2_lambda_kernParsX_yx$calcMsrs(cmem_rbf_rbf_L2_lambda_kernParsX_yx)
  res <- cbind(KCDCxy, KCDCyx)
  return(res)
}, simplify="array"), simplify="array")
proc.time() - pm # 4.5 mins

dim(KCDCs)
dimnames(KCDCs) <- list(msr=c("KCDC","KCDCrel"), dag=c("xy","yx"), par=c("Chl","LAI","fCover"), spectra=1:62)


nameInt <- c("Violet", "Blue", "Cyan", "Green", "Yellow", "Orange", "Red","NearInfraRed","InfraRed")
clrs <- nameInt
posts <- c(400,450,490,520,560, 590,635,700,850,1000)
nameInterval <- nameInt[findInterval(sparc$WaveLength, posts)]

plot(sparc$WaveLength[,1], KCDCs["KCDC","yx","Chl",]/KCDCs["KCDC","xy","Chl",], type="l")
abline(h=1, v=posts, col="red")
plot(sparc$WaveLength[,1], KCDCs["KCDCrel","yx","Chl",]/KCDCs["KCDCrel","xy","Chl",], type="l")
abline(h=1, v=posts, col="red")


plot(sparc$WaveLength[,1], KCDCs["KCDC","yx","LAI",]/KCDCs["KCDC","xy","LAI",], type="l")
abline(h=1, v=posts, col="red")
plot(sparc$WaveLength[,1], KCDCs["KCDCrel","yx","LAI",]/KCDCs["KCDCrel","xy","LAI",], type="l")
abline(h=1, v=posts, col="red")

plot(sparc$WaveLength[,1], KCDCs["KCDC","yx","fCover",]/KCDCs["KCDC","xy","fCover",], type="l")
abline(h=1, v=posts, col="red")
plot(sparc$WaveLength[,1], KCDCs["KCDCrel","yx","fCover",]/KCDCs["KCDCrel","xy","fCover",], type="l")
abline(h=1, v=posts, col="red")

