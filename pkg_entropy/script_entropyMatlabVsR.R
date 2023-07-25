# Test entropy functions comparing to matlab version
remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_entropy")
source("./func_entropy_v1.R")

options(matlab="/usr/local/MATLAB/R2014a/R2014a/bin/matlab")

Matlab$startServer(port=9999)
matlabSession <- Matlab(port=9999)
setVerbose(matlabSession, threshold=200000)
open(matlabSession)
# close(matlabSession)

#Load ITE library into Matlab session
evaluate(matlabSession, "addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/matlab'))")


n <- 100
x <- rnorm(n) # continuous data A
#x <- sample(1:10, n, replace=T) # discrete data A


# Entropy functions

# univariate entropy functions

# 1sp
Shannon_1sp(x)

# gaussian
Shannon_Gauss(x)

# Spacing_V - cant deal with discrete data A (matlab or R)
Shannon_spacing_V(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_V")

# Spacing_Vb - cant deal with discrete data A (matlab or R)
Shannon_spacing_Vb(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vb")

# Spacing_Vp_const - cant deal with discrete data A (matlab or R)
Shannon_spacing_Vpconst(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vpconst")

# Spacing_Vplin - cant deal with discrete data A (matlab or R)
Shannon_spacing_Vplin(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vplin")

# Spacing_Vplin2- cant deal with discrete data A (matlab or R)
Shannon_spacing_Vplin2(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_Vplin2")

# VKDE 
Shannon_spacing_VKDE(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_VKDE")

# spacing_LL
Shannon_spacing_LL(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_spacing_LL")

# PSD_SzegoT 
Shannon_PSD_SzegoT(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_PSD_SzegoT")

# MaxEnt1
Shannon_MaxEnt1(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_MaxEnt1")

# MaxEnt2
Shannon_MaxEnt2(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_MaxEnt2")

###################################################################################333

# Multivariate entropy measures

n <- 100
p <- 10
y <- matrix(rnorm(p*n), n, p)


# expF 
Shannon_expF(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_expF")
Shannon_expF(y) 
genericMatlabEntropy(y, matlabSession, type="Shannon_expF", mult=p)


# kNN_k 
k <- 3
Shannon_kNN_k(x, k)
genericMatlabEntropy(x, matlabSession, type="Shannon_kNN_k" , k=k)
Shannon_kNN_k(y, k)
genericMatlabEntropy(y, matlabSession, type="Shannon_kNN_k" , k=k)


# EdgeWorth 
Shannon_Edgeworth(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_Edgeworth")

Shannon_Edgeworth(y)
genericMatlabEntropy(y, matlabSession, type="Shannon_Edgeworth")


# KDP
Shannon_KDP(x)
genericMatlabEntropy(x, matlabSession, type="Shannon_KDP")

Shannon_KDP(y)
genericMatlabEntropy(y, matlabSession, type="Shannon_KDP")

close(matlabSession)

genericMatlabEntropy(x, matlabSession, type=c("Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", 
                                                          "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", 
                                                          "Shannon_spacing_LL", "Shannon_KDP", "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1", 
                                                          "Shannon_MaxEnt2", "Shannon_expF", "Shannon_vME") , ...)


