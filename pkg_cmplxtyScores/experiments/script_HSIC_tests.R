# Compare different HSIC implementations

# 1. dHSIC package
# 2. https://rdrr.io/github/rikenbit/FUCHIKOMA/src/R/HSIC.R
# 3. kpcalg package
# 4. mooij via matlab and HSICMooij function

remove(list=ls())
library(R.matlab)
library(dHSIC)
library(kpcalg)
#devtools::install_github("rikenbit/FUCHIKOMA") # HSIC
library(fuchikoma)

setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner")
source("./func_causaLearners_v1.R", echo=FALSE)


#Simulate 100 pairs

q <- 100
set.seed(4)
(p <- rep(2, q))
N <- 100
(n <- rep(N,q))
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- lapply(p, function(p) rep(list(nodes),p))


# make list of data with corresponding ground truth DAG
# lets sim 50 data matrices with, 50-200 pts each and 2-4 vars in each case
set.seed(5)
dataTestList <- simRandAddSEMs(q, p, n, nodes, forceEdge=TRUE)



#####################################################
# Compare mooij RWrapper2C vs Matlab implementation 
#####################################################

# help(Matlab)

options(matlab="/usr/local/MATLAB/R2014a/R2014a/bin/matlab")

Matlab$startServer(port=9998)
matlabSession <- Matlab(port=9998)
setVerbose(matlabSession, threshold=200000)
open(matlabSession)
# close(matlabSession)

#Matlab mooij
#set path to mooijs fasthsic function
evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/code_mooijs/fasthsic'))")

# pass sample to matlab
i <- 1
data <- dataTestList$xs[[i]]
setVariable(matlabSession, x=data[,1])
setVariable(matlabSession, y=data[,2])

# calculate a variable in matlab
evaluate(matlabSession, "[pval, res]=fasthsic(x,y);")

# Bring back to R
res <- as.numeric(getVariable(matlabSession, c("res"))$res)
HSICMooij(X=data[,1], Y=data[,2], method="gamma")$statistic


res <- sapply(dataTestList$xs, function(data){
  setVariable(matlabSession, x=data[,1])
  setVariable(matlabSession, y=data[,2])
  evaluate(matlabSession, "[pval, res]=fasthsic(x,y);")
  
  res_mooij_c <- HSICMooij(data[,1], data[,2], method="gamma")$statistic
  res_mooij_m <- as.numeric(getVariable(matlabSession, c("res"))$res)
  return(c(mooij_c=res_mooij_c, mooij_m=res_mooij_m))
})

plot(res[1,], res[2,], xlab="mooij C++", ylab="mooij matlab")
abline(a=0, b=1, col="red")

close(matlabSession)

########################################
# compare statistic   PAIRS
########################################
i <- 1
data <- dataTestList$xs[[i]]
dhsic.test(data[,1], data[,2], method="gamma")$statistic

hsic.test(data[,1], data[,2], p=1, hsic.method="gamma", sig=sqrt(median(as.numeric(dist(data[,2])^2))))$statistic #uses same sigma for both!!!
HSICMooij(X=data[,1], Y=data[,2], method="gamma")$statistic
Kx <- kernelMatrix(rbfdot(sigma= 1/median(as.numeric(dist(data[,1])^2))), data[,1])
Ky <- kernelMatrix(rbfdot(sigma= 1/median(as.numeric(dist(data[,2])^2))), data[,2])
HSIC(Kx, Ky, shrink=FALSE, type="gamma")$HSIC


res <- sapply(dataTestList$xs, function(data){
  res_dhsic <- dhsic.test(data[,1], data[,2], method="gamma")$statistic/N
  res_kpcalg <- hsic.test(data[,1], data[,2], p=1, hsic.method="gamma", sig=sqrt(median(as.numeric(dist(data[,2])^2))))$statistic
  names(res_kpcalg) <- NULL
  res_mooij <- HSICMooij(data[,1], data[,2], method="gamma")$statistic
  Kx <- kernelMatrix(rbfdot(sigma= 1/median(as.numeric(dist(data[,1])^2))), data[,1])
  Ky <- kernelMatrix(rbfdot(sigma= 1/median(as.numeric(dist(data[,2])^2))), data[,2])
  res_fuchi <- HSIC(Kx, Ky, shrink=FALSE, type="gamma")$HSIC
  return(c(dhsic=res_dhsic, kpcalg=res_kpcalg, mooij=res_mooij, fuchi=res_fuchi))
})


pairs(t(res), lower.panel = panel.xy, upper.panel = panel.cor)


#################################################
# compare gamma pvalue implementations 
#################################################

#########
# PAIRS
#########

# values

# coverage, power, level

# p-val dist under NULL (independence) - shd be uniform!

###############
# Multivariate 
# - use max 
# pairwise 
# pval for all 
# except dhsic
###############

# values

# coverage, power, level

# p-val dist under NULL (independence) - shd be uniform!

###################################################
# compare permutation pvalue implementations PAIRS
###################################################

#########
# PAIRS
#########

# values

# coverage, power, level

# p-val dist under NULL (independence) - shd be uniform!

###############
# Multivariate 
# - use max 
# pairwise 
# pval for all 
# except dhsic
###############

# values

# coverage, power, level

# p-val dist under NULL (independence) - shd be uniform!
