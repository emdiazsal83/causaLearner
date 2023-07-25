# correlation in RKHS function

remove(list=ls())

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
print("loading causal learners functions")
hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

n <- 100
x <- runif(2*n,-pi,pi)
y <- sin(x)+rnorm(2*n,sd=0.3)
plot(x,y)

set.seed(1234)
smpl_tr <- sample(2*n,n, replace=T)
smpl_te <- setdiff(1:(2*n), unique(smpl_tr))
nTr <- length(smpl_tr)
nTe <- length(smpl_te)

xTr <- matrix(x[smpl_tr],nTr,1)
xTe <- matrix(x[smpl_te],nTe,1)
yTr <- matrix(y[smpl_tr],nTr,1)
yTe <- matrix(y[smpl_te],nTe,1)

sigmax <- 1
L <- kern_rbf(xTr, sigma=sigmax)
Lte_tr <- kern_rbf(xTe,xTr, sigma=sigmax)
Lte <- kern_rbf(xTe, sigma=sigmax)
lambda <- 0.1
I <- diag(n)
Blambda <- solve(L+n*lambda*I)
sigmay <- 3

corrRKHS <- function(sigmay, B, Lte_tr, ytr, Lte=NULL, yte=NULL){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- Lte_tr%*%B
  
  if(is.null(yte)){
    Lte <- Lte_tr
    Ktr_te <- K 
    }else{
      Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
    }
  
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  res <- cor(c(Lte),c(LB%*%Ktr_te), method="spearman")
  return(res)
}
corrRKHS_col <- function(sigmay, B, Lte_tr, ytr, Lte=NULL, yte=NULL){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- Lte_tr%*%B
  
  if(is.null(yte)){
    Lte <- Lte_tr
    Ktr_te <- K 
  }else{
    Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
  }
  
  LBK <- LB%*%Ktr_te
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    res <- cor(Lte_col,LBK_col, method="spearman")
    #print(res)
    return(res)
    }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- which.max(corrs)
  #plot(Lte_tr[,indx], LBK[,indx])
  res <- mean(corrs)
  return(res)
}
wcorrRKHS_col <- function(sigmay, B, Lte_tr, ytr, Lte=NULL, yte=NULL, method="pearson"){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- Lte_tr%*%B
  
  if(is.null(yte)){
    Lte <- Lte_tr
    Ktr_te <- K 
  }else{
    Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
  }
  
  LBK <- LB%*%Ktr_te
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    #res <- cor(Lte_col,LBK_col, method="spearman")
    res <- weightedCorr(Lte_col,LBK_col, method=method, weights=Lte_col/sum(Lte_col))
    #print(res)
    return(res)
  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- match(sort(corrs,decreasing=T),corrs)[3]
  # indx <- 2
  #plot(Lte[indx,], LBK[indx,], cex=Lte[indx,])
  #res <- mean(corrs)
  res <- corrs
  return(res)
}

pm <- proc.time()
corrRKHS(sigmay, B=Blambda, Lte_tr, ytr=yTr)
corrRKHS(sigmay, B=Blambda, Lte_tr, ytr=yTr, Lte=Lte, yte=yTe)
proc.time() - pm # 0.073

pm <- proc.time()
corrRKHS_col(sigmay, B=Blambda, L, ytr=yTr)
corrRKHS_col(sigmay, B=Blambda, Lte_tr=Lte_tr, ytr=yTr, Lte=Lte, yte=yTe)
proc.time() - pm # 0.062

library(wCorr)
pm <- proc.time()
wcorrRKHS_col(sigmay, B=Blambda, L, ytr=yTr)
wcorrRKHS_col(sigmay, B=Blambda, Lte_tr=Lte_tr, ytr=yTr, Lte=Lte, yte=yTe)
proc.time() - pm # 0.062


sigmasy <- 10^seq(-16,16,1)
corrsRKHS_tr <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda, L, yTr))/(nTr^2)
corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda, Lte_tr, ytr=yTr, Lte=Lte, yte=yTe ))/(nTe^2)

plot(log(sigmasy, 10), corrsRKHS_tr, ylim=range(corrsRKHS_tr, corrsRKHS_te, na.rm=T),type="l")
lines(log(sigmasy, 10), corrsRKHS_te, col="purple", type="l")
indxMax <- which(abs(corrsRKHS_te-max(corrsRKHS_te))/max(corrsRKHS_te)<0.001)
indxMax <- indxMax[length(indxMax)]
indxMax
abline(v=log(sigmasy,10)[indxMax], col="red")

# For a certain dataset and learner load grid of NCE loss for 
# sigmax, sigmay and lambda. Maybe for a certain subspace of 
# lambda and sigmax with optimal NegCE we can fine tune parameter
# selection

# use debug structure to load the grid

experimentName <- "dag2-ME2-SIMdens_cmem-comp-4-l2"
#experimentName <-  "dag2-ME2-Cmplx-SinPlus-SinTimes_cmem-comp-4-l2"
dataName <- strsplit(experimentName, "_")[[1]][1]

block <- 1
pm <- proc.time()
#load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
load(file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name
dag <- dat$dag

Xtr <- dat$x
Xtr <- apply(Xtr, 2, norml)
Xte <- dat$x
Xte <- apply(Xte, 2, norml)
#plot(Xtr[,"x"],Xtr[,"y"])
plot(Xtr[,1],Xtr[,2])

hyp_scorers_cmem <- c(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1,
                      #cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)
#cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)

cmem_learners <- lapply(hyp_scorers_cmem, function(el){
  hypScorer <- eval(parse(text=el))
  cmem_learner1 <- hypScorer$cmemLearner1
  cmem_learner2 <- hypScorer$cmemLearner2
  #print("***************")
  #print(cmfm_learner1)
  #print(cmfm_learner2)
  #return(NULL)
  return(list(cmem_learner1=cmem_learner1, cmem_learner2=cmem_learner2))
})

cmem_learners1 <- sapply(cmem_learners, function(el) el$cmem_learner1)
cmem_learners2 <- sapply(cmem_learners, function(el) el$cmem_learner2)

indx_lrn <- 9
cmem_learner1 <- cmem_learners1[indx_lrn]
cmem_learner2 <- cmem_learners2[indx_lrn]

expType <- "CMEM"
folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/learnersNCE_SIMdens_f_getRbf4mod2_200/",sep="")
#folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/learnersNCE_CMPLXdens_k_getRbf4/",sep="")
dir(folderLearners)[1]
folderSave <- folderLearners

cmemLearner1 <- cmem_learner1
aux <- strsplit(cmemLearner1, "_")[[1]]
indxDEL1 <- which(aux=="DEL1")
indxDEL2 <- which(aux=="DEL2")
cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
nm1 <- paste(dataNm, cmemLearner1_save, "1",sep="_")
nm <- nm1
regressors <- "x"; nodeTo <- "y"
#regressors <- "1"; nodeTo <- "2"
regressorsChar <- paste(regressors, collapse="-")
regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
fileSaveCmemLearner %in% dir(folderSave)
load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
grid <- cmemLearnerAux$hyperParams$data$grid
dmns <- dimnames(grid)[[2]]
summary(grid["test",,1])
ini <- regexpr("lambda",dmns)
fin <- regexpr("sigma.rbf.X",dmns)
lambda <- as.numeric(substr(dmns, ini+7, fin-2))
table(lambda)
ini <- regexpr("sigma.rbf.X",dmns)
fin <- regexpr("sigma.rbf.Y",dmns)
sigma.rbf.X <- as.numeric(substr(dmns, ini+12, fin-2))
table(sigma.rbf.X)
ini <- regexpr("sigma.rbf.Y",dmns)
fin <- nchar(dmns)
sigma.rbf.Y <- as.numeric(substr(dmns, ini+12, fin))
table(sigma.rbf.Y)
negCE <- grid["test",,1] 
names(negCE) <- NULL

randNegCE <- LogLoss(runif(1e7),rbernoulli(1e7)*1)


# 1. for a given lambda, sigmax, compare NCE loss to corrRKHS for diff sigmay



df <- data.frame(lambda=lambda, sigma.rbf.X=sigma.rbf.X, sigma.rbf.Y=sigma.rbf.Y, negCE=negCE)
df$logLambda <- round(log(df$lambda, 10),2)
df$logSigma.rbf.X <- round(log(df$sigma.rbf.X,10),2)

p <- ggplot(df)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=negCE))
p <- p + geom_hline(yintercept=randNegCE, colour="red")
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=negCE))
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
#p <- p + ylim(c(0.25,1))
p

sigmasy <- 10^seq(-16,16,1)

x <- cmemLearnerAux$hyperParams$trainData$x
xTe <- cmemLearnerAux$hyperParams$trainData$x[51:100,,drop=F]
xTr <- cmemLearnerAux$hyperParams$trainData$x[1:50,,drop=F]
y <- as.matrix(cmemLearnerAux$hyperParams$trainData$y)
yTe <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[51:100])
yTr <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[1:50])
nTr <- nrow(xTr)
nTe <- nrow(xTe)
n <- nrow(x) 
Itr <- diag(nTr)  
I <- diag(n)  

corrsRKHS <- sapply(unique(lambda), function(lam){ 
  sapply(unique(sigma.rbf.X), function(sigmax){
    # i <- 2; j <- 3; lam <- unique(lambda)[i]; sigmax <- unique(sigma.rbf.X)[j]
    L <- kern_rbf(xTr, sigma=sigmax)
    Lte_tr <- kern_rbf(xTe,xTr, sigma=sigmax)
    Lte <- kern_rbf(xTe, sigma=sigmax)
    Blambda <- solve(L+nTr*lam*Itr)
    corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda, Lte_tr, ytr=yTr, Lte=Lte,yte=yTe ))
    # plot(sigmasy, corrsRKHS_te)
  return(corrsRKHS_te)
  
}, simplify="array")}, simplify="array")
dim(corrsRKHS)
dimnames(corrsRKHS) <- list(sigma.rbf.Y=sigmasy, sigma.rbf.X=unique(sigma.rbf.X), lambda=unique(lambda)) 

dfCorrRKHS <- melt(corrsRKHS)
colnames(dfCorrRKHS)[4] <- "corrRKHS"

dfCorrRKHS$logLambda <- round(log(dfCorrRKHS$lambda, 10),2)
dfCorrRKHS$logSigma.rbf.X <- round(log(dfCorrRKHS$sigma.rbf.X,10),2)

p <- ggplot(dfCorrRKHS)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS))
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS))
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p

dfNorm <- df[which(df$negCE<=2),]
dfNorm$negCE <- norml(dfNorm$negCE)+0.5
dfCorrRKHSNorm <- dfCorrRKHS[which(!is.nan(dfCorrRKHS$corrRKHS)),]
dfCorrRKHSNorm$corrRKHS <- norml(dfCorrRKHSNorm$corrRKHS)+0.5

p <- ggplot()
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm, size=0.5)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS),colour="red",data=dfCorrRKHSNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS),colour="red",data=dfCorrRKHSNorm, size=0.5)
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p

# correlation is saturating wich means that the y features corresponding
# to low sigma_y values are too easy too project onto (constant?).
# The best negCE seems to be in change of phase zone. 


# Lets obtain p-values for the corrRKHS against a distribution of permuted values
x <- cmemLearnerAux$hyperParams$trainData$x
xTe <- cmemLearnerAux$hyperParams$trainData$x[51:100,,drop=F]
xTr <- cmemLearnerAux$hyperParams$trainData$x[1:50,,drop=F]
y <- as.matrix(cmemLearnerAux$hyperParams$trainData$y)
yTe <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[51:100])
yTr <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[1:50])
nTr <- nrow(xTr)
nTe <- nrow(xTe)
n <- nrow(x) 
Itr <- diag(nTr)  
I <- diag(n)  


numPerms <- 100
set.seed(1234)
rperms <- sapply(1:numPerms, function(i) sample(nTe), simplify="matrix")
yperms <- sapply(1:numPerms, function(i) yTe[rperms[,i]], simplify="matrix")

corrRKHS_dep <- function(sigmay, B, L, ytr, yte=NULL){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- L%*%B
  LBK <- LB%*%K
  if(is.null(yte)){
    Ktr_te <- K 
    LBKtr_te <- LBK
  }else{
    Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
    LBKtr_te <- LB%*%Ktr_te
  }
  
  res <- (L*(LBKtr_te))/diag(LBK%*%B%*%t(L))
  return(sum(res))
}

norml2 <- function(x){
  res <- (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) - 0.5 # centered at 0 
  return(res)
}

sigmasy <- 10^seq(-16,16,1)
i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_difRef <- sapply(unique(lambda), function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(unique(sigma.rbf.X), function(sigmax){
    # i <- 1; j <- 1; lam <- unique(lambda)[i]; sigmax <- unique(sigma.rbf.X)[j]
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    Lte <- kern_rbf(xTe, sigma=sigmax)
    Ltr <- kern_rbf(xTr, sigma=sigmax)
    Lte_tr <- kern_rbf(xTe,xTr, sigma=sigmax)
    Blambda_tr <- solve(Ltr+nTr*lam*Itr)    
    
    corrRKHS(1, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(yTe))/(nTe^2)
    corrRKHS(1, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(yperms[,3]))/(nTe^2)

    corrRKHS_dep(1, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2)
    corrRKHS_dep(1, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yperms[,3]))/(nTe^2)
    
    sigmasy <- 10^seq(-16,16,1)  
    pm <- proc.time()
    corrsRKHS_ref <- apply(yperms, 2, function(col) sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(col))))
    proc.time() - pm # 133secs, 33 sigmas, 100 data points, 1000 permutations
    if(FALSE){
      LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr, as.matrix(yperms[,1]), sigma=sigmasy[which.max(difRef)])
      plot(c(Lte), c(LBK))
      plot(Lte[,21], LBK[,21])
    }
    
    nboots <- 100
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(i){
      #i <- 1
      smpl_tr <- sample(nTr, replace=T)
      smpl_te <- sample(setdiff(1:n, unique(smpl_tr)), nTe, replace=T)
      xTe_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_te,,drop=F]
      xTr_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_tr,,drop=F]
      yTe_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_te])
      yTr_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_tr])
      nTr_b <- nrow(xTr_b)
      nTe_b <- nrow(xTe_b)
      Itr_b <- diag(nTr_b)  
      
      Ltr <- kern_rbf(xTr_b, sigma=sigmax)
      Lte <- kern_rbf(xTe_b, sigma=sigmax)
      Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr*lam*Itr)
      if(FALSE){
        LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[which.max(difRef)])
        plot(c(Lte), c(LBK)) 
        plot(c(Lte[,21]), c(LBK[,21])) 
      }
      
      corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, B=Blambda_tr, L=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b ))
    }, simplify="matrix")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    corrsRKHS_ref_mu <- apply(corrsRKHS_ref,1,mean)
    corrsRKHS_te_boot_mu <- apply(corrsRKHS_te_boot,1,mean)
    if(any(!is.na(corrsRKHS_te_boot_mu))){
      plot(log(sigmasy,10), corrsRKHS_ref_mu, col="red",ylim=range(corrsRKHS_ref_mu,corrsRKHS_te_boot_mu,na.rm=T), type="l")
      lines(log(sigmasy,10), corrsRKHS_te_boot_mu, col="blue")
    }
    difRef <- corrsRKHS_te_boot_mu - corrsRKHS_ref_mu
    corrsRKHS_dep <- sapply(sigmasy, function(sigmay) corrRKHS_dep(sigmay, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2))
    if(any(!is.na(difRef))){
      plot(log(sigmasy,10), difRef, type="l")
      plot(log(sigmasy,10), norml2(difRef), type="l")
      lines(log(sigmasy,10), norml2(corrsRKHS_dep), col="red")
    }
    
    if(FALSE){
    # identify peak
    xx <- log(sigmasy, 10)
    yy <- difRef
    difRef_sp <- splinefun(xx, yy)
    xxx <- seq(log(min(sigmasy),10), log(max(sigmasy),10), length.out=1000)
    yyy <- difRef_sp(xxx)
    plot(xx, yy,ylim=range(yy))
    lines(xxx, yyy, col="red")
    dyyy_dxxx <- difRef_sp(xxx,deriv=1)
    #hist(dyyy_dxxx)
    which(abs(dyyy_dxxx)>0.1)
    indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), 0.95))*1+1
    lines(xxx, c(mean(yy),mean(yy)*2)[indxAux], col=c("blue","green")[indxAux], type="p")
    ps <- seq(0.8,1,length.out=1000)
    propTog <- sapply(ps, function(p){
      indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), p))*1+1
      indxAux2 <- which(indxAux==2)
      res <- length(indxAux2)/(max(indxAux2)-min(indxAux2)+1)
      return(res)
    })
    plot(ps, propTog)
    indxAux3 <- which(propTog<1)
    pMax <- ps[indxAux3][which.max(propTog[indxAux3])]
    plot(xx, yy,ylim=range(yy))
    lines(xxx, yyy, col="red")
    indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), pMax))*1+1
    lines(xxx, c(mean(yy),mean(yy)*2)[indxAux], col=c("blue","green")[indxAux], type="p")
    rngY <- range(xxx[which(indxAux==2 & yyy>0)])
    
    sigmasy <- 10^seq(rngY[1],rngY[2], length.out=20)  
    pm <- proc.time()
    corrsRKHS_ref <- apply(yperms, 2, function(col) sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(col)))/(nTe^2))
    proc.time() - pm # 133secs, 33 sigmas, 100 data points, 1000 permutations
    nboots <- 100
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(i){
      #i <- 1
      smpl_tr <- sample(nTr, replace=T)
      smpl_te <- sample(setdiff(1:n, unique(smpl_tr)), nTe, replace=T)
      xTe_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_te,,drop=F]
      xTr_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_tr,,drop=F]
      yTe_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_te])
      yTr_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_tr])
      nTr_b <- nrow(xTr_b)
      nTe_b <- nrow(xTe_b)
      Itr_b <- diag(nTr_b)  
      
      Ltr <- kern_rbf(xTr_b, sigma=sigmax)
      Lte <- kern_rbf(xTe_b, sigma=sigmax)
      Ltr_te <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr*lam*Itr)
      corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, B=Blambda_tr, L=Ltr_te, ytr=yTr_b, Lte=Lte,yte=yTe_b ))/(nTe_b^2)
    }, simplify="matrix")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    corrsRKHS_ref_mu <- apply(corrsRKHS_ref,1,mean)
    corrsRKHS_te_boot_mu <- apply(corrsRKHS_te_boot,1,mean)
    plot(log(sigmasy,10), corrsRKHS_ref_mu, col="red",ylim=range(corrsRKHS_ref_mu,corrsRKHS_te_boot_mu), type="l")
    lines(log(sigmasy,10), corrsRKHS_te_boot_mu, col="blue")
    difRef <- corrsRKHS_te_boot_mu - corrsRKHS_ref_mu
    plot(log(sigmasy,10), difRef, type="l")
    corrsRKHS_dep <- sapply(sigmasy, function(sigmay) corrRKHS_dep(sigmay, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2))
    plot(log(sigmasy, 10), corrsRKHS_dep)
    plot(log(sigmasy,10), norml(difRef), type="l")
    lines(log(sigmasy,10), norml(corrsRKHS_dep), col="red")
    # identify when diff in correlation is a certain % more than the perm mean
    
    which.max(difRef)
    pm <- proc.time()
    pvals_ks <- sapply(1:length(sigmasy), function(i){
      # i <- 8
      pvals_corr <- sapply(1:nboots, function(j) sum(corrsRKHS_ref[i,]<corrsRKHS_te_boot[i,j])/numPerms)
      # hist(pvals_corr)
      pval_ks <- ks.test(pvals_corr, "runif")$statistic
      return(pval_ks)
    })
    proc.time() -pm# 0.067 secs1000 perms, 100 boots
    
    plot(log(sigmasy,10), pvals_ks,type="l")
    plot(log(sigmasy,10), norml(difRef), type="l")
    lines(log(sigmasy,10), norml(pvals_ks), col="red")
    }
    
    # plot(sigmasy, corrsRKHS_te)
    return(difRef)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins

dim(corrsRKHS_difRef)
dimnames(corrsRKHS_difRef) <- list(sigma.rbf.Y=sigmasy, sigma.rbf.X=unique(sigma.rbf.X), lambda=unique(lambda)) 

dfCorrRKHS_difRef <- melt(corrsRKHS_difRef)
colnames(dfCorrRKHS_difRef)[4] <- "corrRKHS_difRef"

dfCorrRKHS_difRef$logLambda <- round(log(dfCorrRKHS_difRef$lambda, 10),2)
dfCorrRKHS_difRef$logSigma.rbf.X <- round(log(dfCorrRKHS_difRef$sigma.rbf.X,10),2)

dfCorrRKHS_difRefNorm <- dfCorrRKHS_difRef[which(!is.na(dfCorrRKHS_difRef$corrRKHS_difRef)),]
dfCorrRKHS_difRefNorm$corrRKHS_difRef <- norml2(dfCorrRKHS_difRefNorm$corrRKHS_difRef)+0.5

p <- ggplot(dfCorrRKHS_difRef)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef))
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef))
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p

p <- ggplot()
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef),colour="red",data=dfCorrRKHS_difRefNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef),colour="red",data=dfCorrRKHS_difRefNorm)
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p


# Lets repeat but with column to column correlation

sigmasy <- 10^seq(-16,16,1)
i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_difRef <- sapply(unique(lambda), function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(unique(sigma.rbf.X), function(sigmax){
    # i <- 1; j <- 1; lam <- unique(lambda)[i]; sigmax <- unique(sigma.rbf.X)[j]
    # lam <- 0.1; sigmax <- 1
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    Lte <- kern_rbf(xTe, sigma=sigmax)
    Ltr <- kern_rbf(xTr, sigma=sigmax)
    Lte_tr <- kern_rbf(xTe,xTr, sigma=sigmax)
    Blambda_tr <- solve(Ltr+nTr*lam*Itr)    
    
    corrRKHS_col(1, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(yTe))/(nTe^2)
    corrRKHS_col(1, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(yperms[,3]))/(nTe^2)
    
    corrRKHS_dep(1, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2)
    corrRKHS_dep(1, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yperms[,3]))/(nTe^2)
    
    difYs <- sapply(sigmasy, function(sigmay){
      K <- kern_rbf(yTr, sigma=sigmay)
      min(apply(K,2, function(col) length(unique(col))))
    })/nTr
    
    pm <- proc.time()
    corrsRKHS_ref <- apply(yperms, 2, function(col) sapply(sigmasy, function(sigmay) mean(wcorrRKHS_col(sigmay, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(col)))))
    proc.time() - pm # 133secs, 33 sigmas, 100 data points, 1000 permutations
    if(FALSE){
      LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr, as.matrix(yperms[,1]), sigma=sigmasy[which.max(difRef)])
      plot(c(Lte), c(LBK))
      plot(Lte[21,], LBK[21,])
    }
    
    nboots <- 100
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(i){
      #i <- 1
      smpl_tr <- sample(nTr, replace=T)
      smpl_te <- sample(setdiff(1:n, unique(smpl_tr)), nTe+1, replace=T)
      xTe_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_te,,drop=F]
      xTr_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_tr,,drop=F]
      yTe_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_te])
      yTr_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_tr])
      nTr_b <- nrow(xTr_b)
      nTe_b <- nrow(xTe_b)
      Itr_b <- diag(nTr_b)  
      
      Ltr <- kern_rbf(xTr_b, sigma=sigmax)
      Lte <- kern_rbf(xTe_b, sigma=sigmax)
      Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr*lam*Itr)
      if(FALSE){
        LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
        plot(c(Lte), c(LBK)) 
        plot(c(Lte[21,]), c(LBK[21,])) 
        cor(c(Lte[21,]), c(LBK[21,]))
      }
      
      corrsRKHS_te <- sapply(sigmasy, function(sigmay) mean(wcorrRKHS_col(sigmay, B=Blambda_tr, Lte_tr=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b )), simplify="array")
      
      return(corrsRKHS_te)
    }, simplify="matrix")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    corrsRKHS_ref_mu <- apply(corrsRKHS_ref,1,mean)
    corrsRKHS_te_boot_mu <- apply(corrsRKHS_te_boot,1,mean)
    if(any(!is.na(corrsRKHS_te_boot_mu))){
      plot(log(sigmasy,10), corrsRKHS_ref_mu, col="red",ylim=range(corrsRKHS_ref_mu,corrsRKHS_te_boot_mu,1,na.rm=T), type="l")
      lines(log(sigmasy,10), corrsRKHS_te_boot_mu, col="blue")
      lines(log(sigmasy,10), difYs, col="purple")
      abline(v=log(sigmasy,10)[which.max(difYs)], col="green")
      abline(v=log(sigmasy,10)[rev(which(difYs==max(difYs)))[1]], col="green")
    }
    difRef <- corrsRKHS_te_boot_mu - corrsRKHS_ref_mu
    corrsRKHS_dep <- sapply(sigmasy, function(sigmay) corrRKHS_dep(sigmay, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2))
    if(any(!is.na(difRef))){
      plot(log(sigmasy,10), difRef, type="l")
      plot(log(sigmasy,10), norml2(difRef), type="l")
      lines(log(sigmasy,10), norml2(corrsRKHS_dep), col="red")
    }
    
    if(FALSE){
      # identify peak
      xx <- log(sigmasy, 10)
      yy <- difRef
      difRef_sp <- splinefun(xx, yy)
      xxx <- seq(log(min(sigmasy),10), log(max(sigmasy),10), length.out=1000)
      yyy <- difRef_sp(xxx)
      plot(xx, yy,ylim=range(yy))
      lines(xxx, yyy, col="red")
      dyyy_dxxx <- difRef_sp(xxx,deriv=1)
      #hist(dyyy_dxxx)
      which(abs(dyyy_dxxx)>0.1)
      indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), 0.95))*1+1
      lines(xxx, c(mean(yy),mean(yy)*2)[indxAux], col=c("blue","green")[indxAux], type="p")
      ps <- seq(0.8,1,length.out=1000)
      propTog <- sapply(ps, function(p){
        indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), p))*1+1
        indxAux2 <- which(indxAux==2)
        res <- length(indxAux2)/(max(indxAux2)-min(indxAux2)+1)
        return(res)
      })
      plot(ps, propTog)
      indxAux3 <- which(propTog<1)
      pMax <- ps[indxAux3][which.max(propTog[indxAux3])]
      plot(xx, yy,ylim=range(yy))
      lines(xxx, yyy, col="red")
      indxAux <- (abs(difRef_sp(xxx,deriv=1))>quantile(difRef_sp(xxx,deriv=1), pMax))*1+1
      lines(xxx, c(mean(yy),mean(yy)*2)[indxAux], col=c("blue","green")[indxAux], type="p")
      rngY <- range(xxx[which(indxAux==2 & yyy>0)])
      
      sigmasy <- 10^seq(rngY[1],rngY[2], length.out=20)  
      pm <- proc.time()
      corrsRKHS_ref <- apply(yperms, 2, function(col) sapply(sigmasy, function(sigmay) corrRKHS(sigmay, Blambda_tr, Lte_tr, ytr=yTr, Lte=Lte, yte=as.matrix(col)))/(nTe^2))
      proc.time() - pm # 133secs, 33 sigmas, 100 data points, 1000 permutations
      nboots <- 100
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:nboots, function(i){
        #i <- 1
        smpl_tr <- sample(nTr, replace=T)
        smpl_te <- sample(setdiff(1:n, unique(smpl_tr)), nTe, replace=T)
        xTe_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_te,,drop=F]
        xTr_b <- cmemLearnerAux$hyperParams$trainData$x[smpl_tr,,drop=F]
        yTe_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_te])
        yTr_b <- as.matrix(cmemLearnerAux$hyperParams$trainData$y[smpl_tr])
        nTr_b <- nrow(xTr_b)
        nTe_b <- nrow(xTe_b)
        Itr_b <- diag(nTr_b)  
        
        Ltr <- kern_rbf(xTr_b, sigma=sigmax)
        Lte <- kern_rbf(xTe_b, sigma=sigmax)
        Ltr_te <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
        Blambda_tr <- solve(Ltr+nTr*lam*Itr)
        corrsRKHS_te <- sapply(sigmasy, function(sigmay) corrRKHS(sigmay, B=Blambda_tr, L=Ltr_te, ytr=yTr_b, Lte=Lte,yte=yTe_b ))/(nTe_b^2)
      }, simplify="matrix")
      proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
      corrsRKHS_ref_mu <- apply(corrsRKHS_ref,1,mean)
      corrsRKHS_te_boot_mu <- apply(corrsRKHS_te_boot,1,mean)
      plot(log(sigmasy,10), corrsRKHS_ref_mu, col="red",ylim=range(corrsRKHS_ref_mu,corrsRKHS_te_boot_mu), type="l")
      lines(log(sigmasy,10), corrsRKHS_te_boot_mu, col="blue")
      difRef <- corrsRKHS_te_boot_mu - corrsRKHS_ref_mu
      plot(log(sigmasy,10), difRef, type="l")
      corrsRKHS_dep <- sapply(sigmasy, function(sigmay) corrRKHS_dep(sigmay, Blambda_tr, Lte_tr, ytr=yTr,  yte=as.matrix(yTe))/(nTe^2))
      plot(log(sigmasy, 10), corrsRKHS_dep)
      plot(log(sigmasy,10), norml(difRef), type="l")
      lines(log(sigmasy,10), norml(corrsRKHS_dep), col="red")
      # identify when diff in correlation is a certain % more than the perm mean
      
      which.max(difRef)
      pm <- proc.time()
      pvals_ks <- sapply(1:length(sigmasy), function(i){
        # i <- 8
        pvals_corr <- sapply(1:nboots, function(j) sum(corrsRKHS_ref[i,]<corrsRKHS_te_boot[i,j])/numPerms)
        # hist(pvals_corr)
        pval_ks <- ks.test(pvals_corr, "runif")$statistic
        return(pval_ks)
      })
      proc.time() -pm# 0.067 secs1000 perms, 100 boots
      
      plot(log(sigmasy,10), pvals_ks,type="l")
      plot(log(sigmasy,10), norml(difRef), type="l")
      lines(log(sigmasy,10), norml(pvals_ks), col="red")
    }
    
    # plot(sigmasy, corrsRKHS_te)
    return(difRef)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins

dim(corrsRKHS_difRef)
dimnames(corrsRKHS_difRef) <- list(sigma.rbf.Y=sigmasy, sigma.rbf.X=unique(sigma.rbf.X), lambda=unique(lambda)) 

dfCorrRKHS_difRef <- melt(corrsRKHS_difRef)
colnames(dfCorrRKHS_difRef)[4] <- "corrRKHS_difRef"

dfCorrRKHS_difRef$logLambda <- round(log(dfCorrRKHS_difRef$lambda, 10),2)
dfCorrRKHS_difRef$logSigma.rbf.X <- round(log(dfCorrRKHS_difRef$sigma.rbf.X,10),2)

dfCorrRKHS_difRefNorm <- dfCorrRKHS_difRef[which(!is.na(dfCorrRKHS_difRef$corrRKHS_difRef)),]
dfCorrRKHS_difRefNorm$corrRKHS_difRef <- norml2(dfCorrRKHS_difRefNorm$corrRKHS_difRef)+0.5

p <- ggplot(dfCorrRKHS_difRef)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef))
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef))
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p

p <- ggplot()
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=negCE),data=dfNorm)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef),colour="red",data=dfCorrRKHS_difRefNorm)
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=corrRKHS_difRef),colour="red",data=dfCorrRKHS_difRefNorm)
p <- p + facet_grid(logLambda~logSigma.rbf.X, scales="free")
p

# Lets keep all the bootstrap samples and relate optimal sigmay to y
# to see if we could vary sigmay... this perhaps cd also be done
# for sigmax but since x will in general be multidim it might get hairy
# we'll need to take bootstrap samples out

# non hits: 1,2,10,16,18, 24,26, 34, 39, 40,41,42, 45-48,50
# 53, 62, 64, 66, 68-70, 72, 75, 81-83, 85-86, 88-90, 92
#97, 99-100

#experimentName <- "dag2-ME2-SIMdens-200_cmem-comp-4-l2"
#experimentName <-  "dag2-ME2-Cmplx-SinPlus-SinTimes_cmem-comp-4-l2"
#dataName <- strsplit(experimentName, "_")[[1]][1]
dataName <- "dag2-ME2-SIM-1000"

block <- 33
pm <- proc.time()
#load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
load(file=paste("./data/TCEPs/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name
Xtr <- dat$x
Xtr <- apply(Xtr, 2, norml)
Xte <- dat$x
Xte <- apply(Xte, 2, norml)
plot(Xtr[,1],Xtr[,2])

#mod <- mkde(Xtr)
#probs <- 1/mod  
#probs <- probs/sum(probs)
#smplDens <- sample(1:1000, 10000, prob=probs, replace=T)
#Xtr2 <- Xtr[smplDens,]
#hist(Xtr[,"x"]); hist(Xtr2[,"x"])
#hist(Xtr[,"y"]); hist(Xtr2[,"y"])
#plot(Xtr2[,"x"],Xtr2[,"y"])

X2 <- Xtr[,"x", drop=F]
Y2 <- Xtr[,"y",drop=F]
RNGversion("3.5.0")
#RNGversion("3.6.1")
set.seed(1234)
smpl <- sample(1:nrow(X2))
X2 <- X2[smpl,,drop=F]
Y2 <- Y2[smpl,,drop=F]
plot(X2, Y2)
n <- nrow(X2)
hist(residuals(lm(Y2~X2)))


numPerBoot <- 100
nboots <- 100
RNGversion("3.5.0")
#RNGversion("3.6.1")
set.seed(1234)
smpl_tr2 <- sapply(1:nboots, function(i) sample(1:n,numPerBoot, replace=F), simplify="matrix")
smpl_te2 <- sapply(1:nboots, function(i) sample(setdiff(1:n, unique(smpl_tr2[,i])), numPerBoot, replace=F), simplify="matrix")

dists2_xy <- as.numeric(dist(X2)^2)
dists2_xy <- dists2_xy[which(dists2_xy>0)]
dists2_yx <- as.numeric(dist(Y2)^2)
dists2_yx <- dists2_yx[which(dists2_yx>0)]

sigma1_xy <- 1/quantile(dists2_xy, 0.99)  
sigma2_xy <- 1/quantile(dists2_xy, 0.2)
sigma1_yx <- 1/quantile(dists2_yx, 0.99)  
sigma2_yx <- 1/quantile(dists2_yx, 0.2)

sigmaxSeq_xy <- 10^seq(log(sigma1_xy,10), log(sigma2_xy,10), length.out=10)
sigmaxSeq_yx <- 10^seq(log(sigma1_yx,10), log(sigma2_yx,10), length.out=10)
sigmaxSeq <- unique(sort(c(sigmaxSeq_xy, sigmaxSeq_yx)))
lambdaSeq <- 10^seq(-5,-1,1)
#sigmaxSeq <- rep(sigmaxSeq[1],2)

sigmasy <- 10^seq(-1,3,length.out=10)

wcorrRKHS_col <- function(sigmay, B, Lte_tr, ytr, Lte=NULL, yte=NULL, method="pearson"){
  K <- kern_rbf(ytr, sigma=sigmay)
  LB <- Lte_tr%*%B
  ntr <- nrow(ytr)
  Itr <- diag(ntr)
  Htr <- Itr-matrix(1/ntr,ntr,ntr)
  K <- Htr%*%K%*%Htr
  
  if(is.null(yte)){
    Lte <- Lte_tr
    Ktr_te <- K 
    nte <- ntr
  }else{
    nte <- nrow(yte)
    Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
    Ite <- diag(nte)
    Hte <- Ite-matrix(1/nte,nte,nte)
    Ktr_te <- Htr%*%Ktr_te%*%Hte
  }
  
  LBK <- LB%*%Ktr_te
  #res <- (Lte*(LB%*%Ktr_te))/diag(LB%*%K%*%t(LB))
  corrs <- mcmapply(function(Lte_col,LBK_col){
    #res <- cor(Lte_col,LBK_col, method="spearman")
    res <- weightedCorr(Lte_col,LBK_col, method=method, weights=Lte_col/sum(Lte_col))
    #print(res)
    return(res)
  }, Lte_col=as.list(as.data.frame(t(Lte))), LBK_col=as.list(as.data.frame(t(LBK))), mc.cores=1)
  #indx <- match(sort(corrs,decreasing=T),corrs)[3]
  # indx <- 2
  #plot(Lte[indx,], LBK[indx,], cex=Lte[indx,])
  #res <- mean(corrs)
  res <- corrs
  return(res)
}

i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_xy <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_xy , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    # L <- kern_rbf(X, sigma=sigmax)
    # weights <- apply(L, 2, function(col) col/sum(col))
    # wmeansX <- apply(weights, 2, function(col) weighted.mean(c(X), w=col))
    # wmeansY <- apply(weights, 2, function(col) weighted.mean(c(Y), w=col))
    # wsdsX <- apply(weights, 2, function(col) sqrt(weighted.var(c(X), w=col)))
    # wsdsY <- apply(weights, 2, function(col) sqrt(weighted.var(c(Y), w=col)))
    # Xnorm <- (X-wmeansX)/wsdsX
    # Ynorm <- (Y-wmeansY)/wsdsY
    
    difYs <- sapply(sigmasy, function(sigmay){
      K <- kern_rbf(Y2, sigma=sigmay)
      min(apply(K,2, function(col) length(unique(col))))
    })/n
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(k){
      #k <- 1
      
      
      
      xTe_b2 <- X2[smpl_te2[,k],,drop=F]
      xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)  
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      if(FALSE){
        LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[rev(which(difYs==max(difYs)))[1]])
        plot(c(Lte), c(LBK)) 
        plot(c(Lte[21,]), c(LBK[21,])) 
        cor(c(Lte[21,]), c(LBK[21,]))
      }
      
      corrsRKHS_te2 <- sapply(sigmasy, function(sigmay) wcorrRKHS_col(sigmay, B=Blambda_tr2, Lte_tr=Lte_tr2, ytr=yTr_b2, Lte=Lte2,yte=yTe_b2 ), simplify="array")
      
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), sigmay=sigmasy, boot=1:nboots)
    
    #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins
i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_yx <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_yx , function(sigmax){
    # i <- 1; j <- 1; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_yx)[j]
    # lam <- 0.1; sigmax <- 1
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    # L <- kern_rbf(Y, sigma=sigmax)
    # weights <- apply(L, 2, function(col) col/sum(col))
    # wmeansX <- apply(weights, 2, function(col) weighted.mean(c(X), w=col))
    # wmeansY <- apply(weights, 2, function(col) weighted.mean(c(Y), w=col))
    # wsdsX <- apply(weights, 2, function(col) sqrt(weighted.var(c(X), w=col)))
    # wsdsY <- apply(weights, 2, function(col) sqrt(weighted.var(c(Y), w=col)))
    # Xnorm <- (X-wmeansX)/wsdsX
    # Ynorm <- (Y-wmeansY)/wsdsY
    
    difYs <- sapply(sigmasy, function(sigmay){
      K <- kern_rbf(X2, sigma=sigmay)
      min(apply(K,2, function(col) length(unique(col))))
    })/n
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(k){
      #k <- 1
      
      xTe_b <- Y2[smpl_te2[,k],,drop=F]
      xTr_b <- Y2[smpl_tr2[,k],,drop=F]
      yTe_b <- as.matrix(X2[smpl_te2[,k],])
      yTr_b <- as.matrix(X2[smpl_tr2[,k],])
      nTr_b <- nrow(xTr_b)
      nTe_b <- nrow(xTe_b)
      Itr_b <- diag(nTr_b)  
      
      Ltr <- kern_rbf(xTr_b, sigma=sigmax)
      Lte <- kern_rbf(xTe_b, sigma=sigmax)
      Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
      if(FALSE){
        LBK <- Lte_tr %*% Blambda_tr %*% kern_rbf(yTr_b, yTe_b, sigma=sigmasy[11])
        plot(xTr_b,yTr_b)
        plot(xTe_b,yTe_b)
        l <- 2; lines(xTe_b[l], yTe_b[l], cex=2, col="red",type="p")
        plot(c(Lte), c(LBK)) 
        plot(c(Lte[2,]), c(LBK[2,])) 
        cor(c(Lte[2,]), c(LBK[2,]), method="kendall")
      }
      
      corrsRKHS_te <- sapply(sigmasy, function(sigmay) wcorrRKHS_col(sigmay, B=Blambda_tr, Lte_tr=Lte_tr, ytr=yTr_b, Lte=Lte,yte=yTe_b ), simplify="array")
      
      return(corrsRKHS_te)
    }, simplify="array")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), sigmay=sigmasy, boot=1:nboots)
    
    #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins


i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_xy_norml <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_xy , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(k){
      #k <- 1
      
      
      
      xTe_b2 <- X2[smpl_te2[,k],,drop=F]
      xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)  
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      LB <- Lte_tr2%*%Blambda_tr2
      
      ytrNormlL <- mcmapply(function(Lte_col){
        # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
        ws <- Lte_col/sum(Lte_col)
        res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
        #weighted.mean(res, w=ws); weighted.var(res, w=ws)
        return(res)
      }, Lte_col=as.list(as.data.frame(t(Lte_tr2))), mc.cores=1, SIMPLIFY=FALSE)
      yteNormlL <- mcmapply(function(Lte_col){
        # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
        ws <- Lte_col/sum(Lte_col)
        res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
        #weighted.mean(res, w=ws); weighted.var(res, w=ws)
        return(res)
      }, Lte_col=as.list(as.data.frame(t(Lte2))), mc.cores=1, SIMPLIFY=FALSE)
      
      
      corrsRKHS_te2 <- sapply(sigmasy, function(sigmay){ 
        #sigmay <- sigmasy[1]
        
        
        corrs <- mcmapply(function(Lte_col,LB_col, ytrNorml, yteNorml){
          # p <- 1; Lte_col <- as.list(as.data.frame(t(Lte2)))[[p]]; LB_col<-as.list(as.data.frame(t(LB)))[[p]];ytrNorml<-ytrNormlL[[p]];yteNorml<-yteNormlL[[p]] 
          Ktr_te <- kern_rbf(ytrNorml, yteNorml, sigma=sigmay)
          Ktr_te <- Htr %*% Ktr_te %*% Hte
          LBK_col <- t(t(LB_col)%*%Ktr_te)
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          #print(res)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte2))), LB_col=as.list(as.data.frame(t(LB))), 
           ytrNorml=ytrNormlL, yteNorml=yteNormlL, mc.cores=1)
        return(corrs)
      },simplify="array")
      
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 81 secs, 10 sigmas, 100 data points, 100 bootstrap samples
    
    dim(corrsRKHS_te_boot)
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), sigmay=sigmasy, boot=1:nboots)
    
    #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins

i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_yx_norml <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_yx , function(sigmax){
    # i <- 1; j <- 1; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_yx)[j]
    # lam <- 0.1; sigmax <- 1
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:nboots, function(k){
      #k <- 1
      
      xTe_b <- Y2[smpl_te2[,k],,drop=F]
      xTr_b <- Y2[smpl_tr2[,k],,drop=F]
      yTe_b <- as.matrix(X2[smpl_te2[,k],])
      yTr_b <- as.matrix(X2[smpl_tr2[,k],])
      nTr_b <- nrow(xTr_b)
      nTe_b <- nrow(xTe_b)
      Itr_b <- diag(nTr_b)
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr <- kern_rbf(xTr_b, sigma=sigmax)
      Lte <- kern_rbf(xTe_b, sigma=sigmax)
      Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
      Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
       
      LB <- Lte_tr%*%Blambda_tr
      
      ytrNormlL <- mcmapply(function(Lte_col){
        # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
        ws <- Lte_col/sum(Lte_col)
        res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
        #weighted.mean(res, w=ws); weighted.var(res, w=ws)
        return(res)
      }, Lte_col=as.list(as.data.frame(t(Lte_tr))), mc.cores=1, SIMPLIFY=FALSE)
      yteNormlL <- mcmapply(function(Lte_col){
        # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
        ws <- Lte_col/sum(Lte_col)
        res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
        #weighted.mean(res, w=ws); weighted.var(res, w=ws)
        return(res)
      }, Lte_col=as.list(as.data.frame(t(Lte))), mc.cores=1, SIMPLIFY=FALSE)
      
      
      corrsRKHS_te <- sapply(sigmasy, function(sigmay){ 
        
        corrs <- mcmapply(function(Lte_col,LB_col, ytrNorml, yteNorml){
          Ktr_te <- kern_rbf(ytrNorml, yteNorml, sigma=sigmay)
          Ktr_te <- Htr %*% Ktr_te %*% Hte
          LBK_col <- t(t(LB_col)%*%Ktr_te)
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          #print(res)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte))), LB_col=as.list(as.data.frame(t(LB))), 
        ytrNorml=ytrNormlL, yteNorml=yteNormlL, mc.cores=1)
        return(corrs)
      },simplify="array")
      
      return(corrsRKHS_te)
    }, simplify="array")
    proc.time() - pm # 4.9 secs, 33 sigmas, 100 data points, 100 bootstrap samples
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), sigmay=sigmasy, boot=1:nboots)
    
    #corrsRKHS_te_boot[,which(difYs<max(difYs)),] <- NA
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm #21 mins


dim(corrsRKHS_boot_xy)
dim(corrsRKHS_boot_yx)
dimnames(corrsRKHS_boot_xy)[4:5] <- list(sigmax=sigmaxSeq_xy,lambda=lambdaSeq)
dimnames(corrsRKHS_boot_yx)[4:5] <- list(sigmax=sigmaxSeq_yx,lambda=lambdaSeq)
names(dimnames(corrsRKHS_boot_xy))[4:5] <- c("sigmax","lambda")
names(dimnames(corrsRKHS_boot_yx))[4:5] <- c("sigmax","lambda")

dim(corrsRKHS_boot_xy_norml)
dim(corrsRKHS_boot_yx_norml)
dimnames(corrsRKHS_boot_xy_norml)[4:5] <- list(sigmax=sigmaxSeq_xy,lambda=lambdaSeq)
dimnames(corrsRKHS_boot_yx_norml)[4:5] <- list(sigmax=sigmaxSeq_yx,lambda=lambdaSeq)
names(dimnames(corrsRKHS_boot_xy_norml))[4:5] <- c("sigmax","lambda")
names(dimnames(corrsRKHS_boot_yx_norml))[4:5] <- c("sigmax","lambda")

#save(list=ls(),file="block33.RData")

# NON-NORMALIZED
{
indxOpt_xy <- apply(corrsRKHS_boot_xy, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_xy)[[1]] <- c("sigmay","sigmax","lambda")
names(dimnames(indxOpt_xy))[1] <- "indx"
apply(indxOpt_xy, "indx", table)

indxOpt_yx <- apply(corrsRKHS_boot_yx, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_yx[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_yx)[[1]] <- c("sigmay","sigmax","lambda")
names(dimnames(indxOpt_yx))[1] <- "indx"
apply(indxOpt_yx, "indx", table)

indxOther_xy <- sapply(dimnames(indxOpt_yx)[["testPt"]], function(testPt) sapply(dimnames(indxOpt_yx)[["boot"]], function(boot){
  # i<-1; j<-1; testPt <- dimnames(indxOpt_xy)[["testPt"]][i]; boot <- dimnames(indxOpt_xy)[["boot"]][j]
  indxSigmay <- indxOpt_yx["sigmay",testPt,boot]
  mat <- corrsRKHS_boot_xy[testPt,indxSigmay,boot,,]
  maxMat <- mat==max(mat,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
  indxs <- c(indxSigmay, indxs)
  return(indxs)
}, simplify="array"),simplify="array")
indxOther_xy <- aperm(indxOther_xy, c(1,3,2))
dim(indxOpt_xy); dim(indxOther_xy)
dimnames(indxOther_xy)[[1]] <- c("sigmay","sigmax","lambda")
names(dimnames(indxOther_xy)) <- c("indx", "testPt","boot")
dimnames(indxOpt_xy); dimnames(indxOther_xy)
apply(indxOther_xy, "indx", table)

indxOther_yx <- sapply(dimnames(indxOpt_xy)[["testPt"]], function(testPt) sapply(dimnames(indxOpt_xy)[["boot"]], function(boot){
  # i<-5; j<-1; testPt <- dimnames(indxOpt_yx)[["testPt"]][i]; boot <- dimnames(indxOpt_yx)[["boot"]][j]
  indxSigmay <- indxOpt_xy["sigmay",testPt,boot]
  mat <- corrsRKHS_boot_yx[testPt,indxSigmay,boot,,]
  maxMat <- mat==max(mat,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
  indxs <- c(indxSigmay, indxs)
  return(indxs)
}, simplify="array"),simplify="array")
indxOther_yx <- aperm(indxOther_yx, c(1,3,2))
dim(indxOpt_yx); dim(indxOther_yx)
dimnames(indxOther_yx)[[1]] <- c("sigmay","sigmax","lambda")
names(dimnames(indxOther_yx)) <- c("indx", "testPt","boot")
dimnames(indxOpt_yx); dimnames(indxOther_yx)
apply(indxOther_yx, "indx", table)

df_xy <- melt(indxOpt_xy)
df_xy <- cast(df_xy, testPt+boot~indx, value="value")
df_xy$corr <- corrsRKHS_boot_xy[as.matrix(df_xy[,names(dimnames(corrsRKHS_boot_xy))])]
colnames(df_xy) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
df_xy$node <- "y"
df_xy$numReg <- "1"
df_xy$lambda <- lambdaSeq[df_xy$indxLambda]
df_xy$sigmax <- sigmaxSeq_xy[df_xy$indxSigmax]
df_xy$sigmay <- sigmasy[df_xy$indxSigmay]

indxMat <- matrix(c(df_xy$testPt,df_xy$boot),nrow(df_xy),2)
df_xy$x <- X2[smpl_te2[indxMat]]
df_xy$y <- Y2[smpl_te2[indxMat]]
plot(df_xy$x, log(df_xy$sigmay,10))
o <- order(X2)
mod_xy <- sapply(X2[o,], function(x) mean(log(df_xy$sigmay[which(df_xy$x==x)],10)))
lines(X2[o,], mod_xy, col="red")
mod_xy2 <- loess(mod_xy~X2[o,], span=0.25)
optSigmay_xy <- predict(mod_xy2, X2[o,])
lines(X2[o,], optSigmay_xy, col="green")
plot(X2,Y2)
lines(X2[o,], norml(mod_xy), col="red")
lines(X2[o,], norml(optSigmay_xy), col="green")
optSigmay_xy <- 10^optSigmay_xy
optSigmay_xy <- optSigmay_xy[order(o)]
plot(X2, log(optSigmay_xy,10))


i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_xy2 <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_xy , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
      #k <- 1
      
      xTe_b2 <- X2[smpl_te2[,k],,drop=F]
      optSigmay_xy_b <- optSigmay_xy[smpl_te2[,k]]
      xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      
      
      LB <- Lte_tr2%*%Blambda_tr2
      
      corrsRKHS_te2 <- mcmapply(function(sigmay, pt){ 
        # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
        
        
        Ktr_te <- kern_rbf(yTr_b2, yTe_b2, sigma=sigmay)
        Ktr_te <- Htr %*% Ktr_te %*% Hte
        LBK <- LB%*%Ktr_te
        Lte_col <- Lte2[pt,]
        LBK_col <- LBK[pt,]
        res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
        return(res)
      }, sigmay=optSigmay_xy_b, pt=1:length(optSigmay_xy_b), SIMPLIFY="array")
        
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 17 secs
    
    dim(corrsRKHS_te_boot)
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
    
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm # 15 mins mins
dim(corrsRKHS_boot_xy2)
dimnames(corrsRKHS_boot_xy2)[3:4] <- list(sigmaxSeq_xy, lambdaSeq)
names(dimnames(corrsRKHS_boot_xy2))[3:4] <- c("sigmax", "lambda")


indxOpt_xy2 <- apply(corrsRKHS_boot_xy2, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_xy2)[[1]] <- c("sigmax","lambda")
names(dimnames(indxOpt_xy2))[1] <- "indx"


df_xy2 <- melt(indxOpt_xy2)
df_xy2 <- cast(df_xy2, testPt+boot~indx, value="value")
df_xy2$corr <- corrsRKHS_boot_xy2[as.matrix(df_xy2[,names(dimnames(corrsRKHS_boot_xy2))])]
colnames(df_xy2) <- c("testPt","boot","indxLambda","indxSigmax","corr")
df_xy2$node <- "y"
df_xy2$numReg <- "1"
df_xy2$lambda <- lambdaSeq[df_xy2$indxLambda]
df_xy2$sigmax <- sigmaxSeq_xy[df_xy2$indxSigmax]
indxMat <- matrix(c(df_xy2$testPt,df_xy2$boot),nrow(df_xy2),2)
df_xy2$sigmay <- optSigmay_xy[df_xy2$testPt]
df_xy2$x <- X2[smpl_te2[indxMat]]
df_xy2$y <- Y2[smpl_te2[indxMat]]
plot(df_xy2$x, log(df_xy2$sigmax,10))
plot(df_xy2$x, log(df_xy2$lambda,10))

df_othr_xy <- melt(indxOther_xy)
df_othr_xy <- cast(df_othr_xy, testPt+boot~indx, value="value")
df_othr_xy$corr <- corrsRKHS_boot_xy[as.matrix(df_othr_xy[,names(dimnames(corrsRKHS_boot_xy))])]
colnames(df_othr_xy) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
df_othr_xy$node <- "x"
df_othr_xy$numReg <- "1"
df_othr_xy$lambda <- lambdaSeq[df_othr_xy$indxLambda]
df_othr_xy$sigmax <- sigmaxSeq_xy[df_othr_xy$indxSigmax]
df_othr_xy$sigmay <- sigmasy[df_othr_xy$indxSigmay]
df_xy <- rbind(df_xy[,1:11], df_othr_xy)
df_xy <- df_xy[order(df_xy$node, df_xy$numReg,df_xy$testPt, df_xy$boot),]

df_yx <- melt(indxOpt_yx)
df_yx <- cast(df_yx, testPt+boot~indx, value="value")
df_yx$corr <- corrsRKHS_boot_yx[as.matrix(df_yx[,names(dimnames(corrsRKHS_boot_yx))])]
colnames(df_yx) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
df_yx$node <- "x"
df_yx$numReg <- "1"
df_yx$lambda <- lambdaSeq[df_yx$indxLambda]
df_yx$sigmax <- sigmaxSeq_yx[df_yx$indxSigmax]
df_yx$sigmay <- sigmasy[df_yx$indxSigmay]


indxMat <- matrix(c(df_yx$testPt,df_yx$boot),nrow(df_yx),2)
df_yx$x <- Y2[smpl_te2[indxMat]]
df_yx$y <- X2[smpl_te2[indxMat]]
plot(df_yx$x, log(df_yx$sigmay,10))
o <- order(Y2)
mod_yx <- sapply(Y2[o,], function(x) mean(log(df_yx$sigmay[which(df_yx$x==x)],10)))
lines(Y2[o,], mod_yx, col="red")
mod_yx2 <- loess(mod_yx~Y2[o,], span=0.25)
optSigmay_yx <- predict(mod_yx2, Y2[o,])
lines(Y2[o,], optSigmay_yx, col="green")
plot(Y2, X2)
lines(Y2[o,], norml(mod_yx), col="red")
lines(Y2[o,], norml(optSigmay_yx), col="green")
optSigmay_yx <- 10^optSigmay_yx
optSigmay_yx <- optSigmay_yx[order(o)]
plot(Y2, log(optSigmay_yx,10))


i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_yx2 <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_yx , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
      #k <- 1
      
      xTe_b2 <- Y2[smpl_te2[,k],,drop=F]
      optSigmay_yx_b <- optSigmay_yx[smpl_te2[,k]]
      xTr_b2 <- Y2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(X2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(X2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)  
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      
      
      LB <- Lte_tr2%*%Blambda_tr2
      
      corrsRKHS_te2 <- mcmapply(function(sigmay, pt){ 
        # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
        
        
        Ktr_te <- kern_rbf(yTr_b2, yTe_b2, sigma=sigmay)
        Ktr_te <- Htr %*% Ktr_te %*% Hte
        LBK <- LB%*%Ktr_te
        Lte_col <- Lte2[pt,]
        LBK_col <- LBK[pt,]
        res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
        return(res)
      }, sigmay=optSigmay_yx_b, pt=1:length(optSigmay_yx_b), SIMPLIFY="array")
      
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 17 secs
    
    dim(corrsRKHS_te_boot)
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
    
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm # 15 mins mins
dim(corrsRKHS_boot_yx2)
dimnames(corrsRKHS_boot_yx2)[3:4] <- list(sigmaxSeq_yx, lambdaSeq)
names(dimnames(corrsRKHS_boot_yx2))[3:4] <- c("sigmax", "lambda")

indxOpt_yx2 <- apply(corrsRKHS_boot_yx2, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_yx2)[[1]] <- c("sigmax","lambda")
names(dimnames(indxOpt_yx2))[1] <- "indx"

df_yx2 <- melt(indxOpt_yx2)
df_yx2 <- cast(df_yx2, testPt+boot~indx, value="value")
df_yx2$corr <- corrsRKHS_boot_yx2[as.matrix(df_yx2[,names(dimnames(corrsRKHS_boot_yx2))])]
colnames(df_yx2) <- c("testPt","boot","indxLambda","indxSigmax","corr")
df_yx2$node <- "x"
df_yx2$numReg <- "1"
df_yx2$lambda <- lambdaSeq[df_yx2$indxLambda]
df_yx2$sigmax <- sigmaxSeq_yx[df_yx2$indxSigmax]
indxMat <- matrix(c(df_yx2$testPt,df_yx2$boot),nrow(df_yx2),2)
df_yx2$sigmay <- optSigmay_yx[df_yx2$testPt]
df_yx2$x <- Y2[smpl_te2[indxMat]]
df_yx2$y <- X2[smpl_te2[indxMat]]
plot(df_yx2$x, log(df_yx2$sigmax,10))
plot(df_yx2$x, log(df_yx2$lambda,10))


df_othr_yx <- melt(indxOther_yx)
df_othr_yx <- cast(df_othr_yx, testPt+boot~indx, value="value")
df_othr_yx$corr <- corrsRKHS_boot_yx[as.matrix(df_othr_yx[,names(dimnames(corrsRKHS_boot_yx))])]
colnames(df_othr_yx) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
summary(df_othr_yx$corr)
df_othr_yx$node <- "y"
df_othr_yx$numReg <- "1"
df_othr_yx$lambda <- lambdaSeq[df_othr_yx$indxLambda]
df_othr_yx$sigmax <- sigmaxSeq_yx[df_othr_yx$indxSigmax]
df_othr_yx$sigmay <- sigmasy[df_othr_yx$indxSigmay]
df_yx <- rbind(df_yx[,1:11], df_othr_yx)

dim(df_xy); dim(df_yx)
table(df_xy$sigmay); table(df_yx$sigmay) 
indx <- order(df_xy$node, df_xy$numReg,df_xy$testPt, df_xy$boot)
df_xy <- df_xy[indx,]
indx <- order(df_yx$node, df_yx$numReg,df_yx$testPt, df_yx$boot)
df_yx <- df_yx[indx,]
head(df_xy); head(df_yx)
all(df_xy$sigmay==df_yx$sigmay)


indxMat <- matrix(c(df_xy$testPt,df_xy$boot),nrow(df_xy),2)
df_xy$y <- Y2[smpl_te2[indxMat]]
df_xy$x <- X2[smpl_te2[indxMat]]
indxMat <- matrix(c(df_yx$testPt,df_yx$boot),nrow(df_yx),2)
df_yx$y <- X2[smpl_te2[indxMat]]
df_yx$x <- Y2[smpl_te2[indxMat]]

df_xy <- as.data.frame(df_xy)
df_yx <- as.data.frame(df_yx)


i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_xy2_othr <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_xy , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
      #k <- 1
      
      xTe_b2 <- X2[smpl_te2[,k],,drop=F]
      optSigmay_xy_b <- optSigmay_yx[smpl_te2[,k]]
      xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)  
      Itr_b2 <- diag(nTr_b2)  
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      
      
      LB <- Lte_tr2%*%Blambda_tr2
      
      corrsRKHS_te2 <- mcmapply(function(sigmay, pt){ 
        # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
        
        
        Ktr_te <- kern_rbf(yTr_b2, yTe_b2, sigma=sigmay)
        Ktr_te <- Htr %*% Ktr_te %*% Hte
        LBK <- LB%*%Ktr_te
        Lte_col <- Lte2[pt,]
        LBK_col <- LBK[pt,]
        res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
        return(res)
      }, sigmay=optSigmay_xy_b, pt=1:length(optSigmay_xy_b), SIMPLIFY="array")
      
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 17 secs
    
    dim(corrsRKHS_te_boot)
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
    
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm # 15 mins mins
dim(corrsRKHS_boot_xy2_othr)
dimnames(corrsRKHS_boot_xy2_othr)[3:4] <- list(sigmaxSeq_xy, lambdaSeq)
names(dimnames(corrsRKHS_boot_xy2_othr))[3:4] <- c("sigmax", "lambda")


indxOpt_xy2_othr <- apply(corrsRKHS_boot_xy2_othr, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_xy2_othr)[[1]] <- c("sigmax","lambda")
names(dimnames(indxOpt_xy2_othr))[1] <- "indx"


df_xy2_othr <- melt(indxOpt_xy2_othr)
df_xy2_othr <- cast(df_xy2_othr, testPt+boot~indx, value="value")
df_xy2_othr$corr <- corrsRKHS_boot_xy2_othr[as.matrix(df_xy2_othr[,names(dimnames(corrsRKHS_boot_xy2_othr))])]
colnames(df_xy2_othr) <- c("testPt","boot","indxLambda","indxSigmax","corr")
df_xy2_othr$node <- "x"
df_xy2_othr$numReg <- "1"
df_xy2_othr$lambda <- lambdaSeq[df_xy2_othr$indxLambda]
df_xy2_othr$sigmax <- sigmaxSeq_xy[df_xy2_othr$indxSigmax]
indxMat <- matrix(c(df_xy2_othr$testPt,df_xy2_othr$boot),nrow(df_xy2_othr),2)
df_xy2_othr$sigmay <- optSigmay_yx[df_xy2_othr$testPt]
df_xy2_othr$x <- X2[smpl_te2[indxMat]]
df_xy2_othr$y <- Y2[smpl_te2[indxMat]]
plot(df_xy2_othr$x, log(df_xy2_othr$sigmax,10))
plot(df_xy2_othr$x, log(df_xy2_othr$lambda,10))

i <- 0
j <- 0
pm <- proc.time()
corrsRKHS_boot_yx2_othr <- sapply(lambdaSeq, function(lam){ 
  i <<- i + 1
  j <<- 0
  sapply(sigmaxSeq_yx , function(sigmax){
    # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
    # lam <- 0.1; sigmax <- 100
    j <<- j + 1
    print(paste("i-j: ", i, "-",j))
    
    
    pm <- proc.time()
    corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
      #k <- 1
      
      xTe_b2 <- Y2[smpl_te2[,k],,drop=F]
      optSigmay_yx_b <- optSigmay_xy[smpl_te2[,k]]
      xTr_b2 <- Y2[smpl_tr2[,k],,drop=F]
      yTe_b2 <- as.matrix(X2[smpl_te2[,k],])
      yTr_b2 <- as.matrix(X2[smpl_tr2[,k],])
      nTr_b2 <- nrow(xTr_b2)
      nTe_b2 <- nrow(xTe_b2)
      Itr_b2 <- diag(nTr_b2)  
      Ite_b2 <- diag(nTe_b2)
      Htr <- Itr_b2-matrix(1/nTr_b2,nTr_b2,nTr_b2)
      Hte <- Ite_b2-matrix(1/nTe_b2,nTe_b2,nTe_b2)
      
      Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
      Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
      Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
      Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
      
      
      LB <- Lte_tr2%*%Blambda_tr2
      
      corrsRKHS_te2 <- mcmapply(function(sigmay, pt){ 
        # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
        
        
        Ktr_te <- kern_rbf(yTr_b2, yTe_b2, sigma=sigmay)
        Ktr_te <- Htr %*% Ktr_te %*% Hte
        LBK <- LB%*%Ktr_te
        Lte_col <- Lte2[pt,]
        LBK_col <- LBK[pt,]
        res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
        return(res)
      }, sigmay=optSigmay_yx_b, pt=1:length(optSigmay_yx_b), SIMPLIFY="array")
      
      return(corrsRKHS_te2)
    }, simplify="array")
    proc.time() - pm # 17 secs
    
    dim(corrsRKHS_te_boot)
    dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
    
    
    return(corrsRKHS_te_boot)
    
  }, simplify="array")}, simplify="array")
proc.time() - pm # 15 mins mins
dim(corrsRKHS_boot_yx2_othr)
dimnames(corrsRKHS_boot_yx2_othr)[3:4] <- list(sigmaxSeq_yx, lambdaSeq)
names(dimnames(corrsRKHS_boot_yx2_othr))[3:4] <- c("sigmax", "lambda")

indxOpt_yx2_othr <- apply(corrsRKHS_boot_yx2_othr, c("testPt","boot"), function(arr){
  # arr <- corrsRKHS_boot_xy[1,1,,]
  maxArr <- arr==max(arr,na.rm=T)
  indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
  #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
  return(indxs)
})
dimnames(indxOpt_yx2_othr)[[1]] <- c("sigmax","lambda")
names(dimnames(indxOpt_yx2_othr))[1] <- "indx"

df_yx2_othr <- melt(indxOpt_yx2_othr)
df_yx2_othr <- cast(df_yx2_othr, testPt+boot~indx, value="value")
df_yx2_othr$corr <- corrsRKHS_boot_yx2_othr[as.matrix(df_yx2_othr[,names(dimnames(corrsRKHS_boot_yx2_othr))])]
colnames(df_yx2_othr) <- c("testPt","boot","indxLambda","indxSigmax","corr")
df_yx2_othr$node <- "y"
df_yx2_othr$numReg <- "1"
df_yx2_othr$lambda <- lambdaSeq[df_yx2_othr$indxLambda]
df_yx2_othr$sigmax <- sigmaxSeq_yx[df_yx2_othr$indxSigmax]
indxMat <- matrix(c(df_yx2_othr$testPt,df_yx2_othr$boot),nrow(df_yx2_othr),2)
df_yx2_othr$sigmay <- optSigmay_xy[df_yx2_othr$testPt]
df_yx2_othr$x <- Y2[smpl_te2[indxMat]]
df_yx2_othr$y <- X2[smpl_te2[indxMat]]
plot(df_yx2_othr$x, log(df_yx2_othr$sigmax,10))
plot(df_yx2_othr$x, log(df_yx2_othr$lambda,10))

plot(log(optSigmay_xy,10), log(optSigmay_yx,10)); abline(a=0, b=1, col="red")
sum(optSigmay_yx>optSigmay_xy)/length(optSigmay_xy)
sum(optSigmay_xy-optSigmay_yx)
sum(log(optSigmay_xy,10)-log(optSigmay_yx,10))
plot(log(df_xy$sigmay[which(df_xy$node=="y")],10), log(df_yx$sigmay[which(df_xy$node=="x")],10)); abline(a=0, b=1, col="red")
sum(df_xy$sigmay[which(df_xy$node=="y")]<df_yx$sigmay[which(df_xy$node=="x")])/sum(df_xy$node=="x")


plot(X2,Y2, col=c("red","blue")[(optSigmay_yx>optSigmay_xy)*1+1])


plot(X2, Y2)
plot(X2, log(optSigmay_yx,10), ylim=range(log(c(optSigmay_yx,optSigmay_xy),10)))
lines(X2, log(optSigmay_xy,10), type="p", col="red")
plot(Y2, X2)

plot(Y2, log(optSigmay_xy,10), ylim=range(log(c(optSigmay_yx,optSigmay_xy),10)))
lines(Y2, log(optSigmay_yx,10), type="p", col="red")
par(mfrow=c(1,1))
hist(X2[which(X2>-0.1 & X2 < -0)])


plot(X2, Y2)
plot(df_xy2$x[which(df_xy2$node=="y")], log(df_xy2$sigmax[which(df_xy2$node=="y")],10))
plot(Y2, X2)
plot(df_yx2$x[which(df_yx2$node=="x")], log(df_yx2$sigmax[which(df_yx2$node=="x")],10))

df_xy2 <- rbind(df_xy2, df_xy2_othr)
df_yx2 <- rbind(df_yx2, df_yx2_othr)

dim(df_xy2); dim(df_yx2)
indx <- order(df_xy2$node, df_xy2$numReg,df_xy2$testPt, df_xy2$boot)
df_xy2 <- df_xy2[indx,]
indx <- order(df_yx2$node, df_yx2$numReg,df_yx2$testPt, df_yx2$boot)
df_yx2 <- df_yx2[indx,]
head(df_xy2); head(df_yx2)
all(df_xy2$sigmay==df_yx2$sigmay)

indxMat <- matrix(c(df_xy2$testPt,df_xy2$boot),nrow(df_xy2),2)
df_xy2$y <- Y2[smpl_te2[indxMat]]
df_xy2$x <- X2[smpl_te2[indxMat]]
indxMat <- matrix(c(df_yx2$testPt,df_yx2$boot),nrow(df_yx2),2)
df_yx2$y <- X2[smpl_te2[indxMat]]
df_yx2$x <- Y2[smpl_te2[indxMat]]

df_xy2 <- as.data.frame(df_xy2)
df_yx2 <- as.data.frame(df_yx2)

}

# NORMALIZED
{
  indxOpt_xy_norml <- apply(corrsRKHS_boot_xy_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_xy_norml)[[1]] <- c("sigmay","sigmax","lambda")
  names(dimnames(indxOpt_xy_norml))[1] <- "indx"
  apply(indxOpt_xy_norml, "indx", table)
  
  indxOpt_yx_norml <- apply(corrsRKHS_boot_yx_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmay","sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_yx[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_yx_norml)[[1]] <- c("sigmay","sigmax","lambda")
  names(dimnames(indxOpt_yx_norml))[1] <- "indx"
  apply(indxOpt_yx_norml, "indx", table)
  
  indxOther_xy_norml <- sapply(dimnames(indxOpt_yx_norml)[["testPt"]], function(testPt) sapply(dimnames(indxOpt_yx_norml)[["boot"]], function(boot){
    # i<-1; j<-1; testPt <- dimnames(indxOpt_xy)[["testPt"]][i]; boot <- dimnames(indxOpt_xy)[["boot"]][j]
    indxSigmay <- indxOpt_yx_norml["sigmay",testPt,boot]
    mat <- corrsRKHS_boot_xy_norml[testPt,indxSigmay,boot,,]
    maxMat <- mat==max(mat,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
    indxs <- c(indxSigmay, indxs)
    return(indxs)
  }, simplify="array"),simplify="array")
  indxOther_xy_norml <- aperm(indxOther_xy_norml, c(1,3,2))
  dim(indxOpt_xy_norml); dim(indxOther_xy_norml)
  dimnames(indxOther_xy_norml)[[1]] <- c("sigmay","sigmax","lambda")
  names(dimnames(indxOther_xy_norml)) <- c("indx", "testPt","boot")
  dimnames(indxOpt_xy_norml); dimnames(indxOther_xy_norml)
  apply(indxOther_xy_norml, "indx", table)
  
  indxOther_yx_norml <- sapply(dimnames(indxOpt_xy_norml)[["testPt"]], function(testPt) sapply(dimnames(indxOpt_xy_norml)[["boot"]], function(boot){
    # i<-5; j<-1; testPt <- dimnames(indxOpt_yx)[["testPt"]][i]; boot <- dimnames(indxOpt_yx)[["boot"]][j]
    indxSigmay <- indxOpt_xy_norml["sigmay",testPt,boot]
    mat <- corrsRKHS_boot_yx_norml[testPt,indxSigmay,boot,,]
    maxMat <- mat==max(mat,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxMat, margin, any))[1])
    indxs <- c(indxSigmay, indxs)
    return(indxs)
  }, simplify="array"),simplify="array")
  indxOther_yx_norml <- aperm(indxOther_yx_norml, c(1,3,2))
  dim(indxOpt_yx_norml); dim(indxOther_yx_norml)
  dimnames(indxOther_yx_norml)[[1]] <- c("sigmay","sigmax","lambda")
  names(dimnames(indxOther_yx_norml)) <- c("indx", "testPt","boot")
  dimnames(indxOpt_yx_norml); dimnames(indxOther_yx_norml)
  apply(indxOther_yx_norml, "indx", table)
  
  df_xy_norml <- melt(indxOpt_xy_norml)
  df_xy_norml <- cast(df_xy_norml, testPt+boot~indx, value="value")
  df_xy_norml$corr <- corrsRKHS_boot_xy_norml[as.matrix(df_xy_norml[,names(dimnames(corrsRKHS_boot_xy_norml))])]
  colnames(df_xy_norml) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
  df_xy_norml$node <- "y"
  df_xy_norml$numReg <- "1"
  df_xy_norml$lambda <- lambdaSeq[df_xy_norml$indxLambda]
  df_xy_norml$sigmax <- sigmaxSeq_xy[df_xy_norml$indxSigmax]
  df_xy_norml$sigmay <- sigmasy[df_xy_norml$indxSigmay]
  
  indxMat <- matrix(c(df_xy_norml$testPt,df_xy_norml$boot),nrow(df_xy_norml),2)
  df_xy_norml$x <- X2[smpl_te2[indxMat]]
  df_xy_norml$y <- apply(as.data.frame(df_xy_norml), 1, function(row){
    # row <- df_xy_norml[1,]; names(row) <- colnames(df_xy_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- X2[smpl_te2[,boot],,drop=F]
    yTe_b <- Y2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(X2, Y2)
  plot(df_xy_norml$x, df_xy_norml$y)
  plot(df_xy_norml$x, log(df_xy_norml$sigmay,10))
  o <- order(X2)
  mod_xy <- sapply(X2[o,], function(x) mean(log(df_xy_norml$sigmay[which(df_xy_norml$x==x)],10)))
  lines(X2[o,], mod_xy, col="red")
  mod_xy2 <- loess(mod_xy~X2[o,], span=0.25)
  optSigmay_xy_norml <- predict(mod_xy2, X2[o,])
  lines(X2[o,], optSigmay_xy_norml, col="green")
  plot(X2,Y2)
  lines(X2[o,], norml(mod_xy), col="red")
  lines(X2[o,], norml(optSigmay_xy_norml), col="green")
  optSigmay_xy_norml <- 10^optSigmay_xy_norml
  optSigmay_xy_norml <- optSigmay_xy_norml[order(o)]
  plot(X2, log(optSigmay_xy_norml,10))
  plot(X2, log(optSigmay_xy,10))
  
  
  i <- 0
  j <- 0
  pm <- proc.time()
  corrsRKHS_boot_xy2_norml <- sapply(lambdaSeq, function(lam){ 
    i <<- i + 1
    j <<- 0
    sapply(sigmaxSeq_xy , function(sigmax){
      # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
      # lam <- 0.1; sigmax <- 100
      j <<- j + 1
      print(paste("i-j: ", i, "-",j))
      
      
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
        #k <- 1
        
        xTe_b2 <- X2[smpl_te2[,k],,drop=F]
        optSigmay_xy_b <- optSigmay_xy_norml[smpl_te2[,k]]
        xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
        yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
        yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
        nTr_b2 <- nrow(xTr_b2)
        nTe_b2 <- nrow(xTe_b2)
        Itr_b2 <- diag(nTr_b2)  
        
        Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
        Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
        Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
        Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
        
        ytrNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte_tr2))), mc.cores=1, SIMPLIFY=FALSE)
        yteNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte2))), mc.cores=1, SIMPLIFY=FALSE) 
        
        LB <- Lte_tr2%*%Blambda_tr2
        
        corrsRKHS_te2 <- mcmapply(function(sigmay, pt, ytr, yte){ 
          # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l; ytr=ytrNormlL[[l]]; yte=ytrNormlL[[l]]
          
          
          Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
          LBK <- LB%*%Ktr_te
          Lte_col <- Lte2[pt,]
          LBK_col <- LBK[pt,]
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          return(res)
        }, sigmay=optSigmay_xy_b, pt=1:length(optSigmay_xy_b), ytr=ytrNormlL, yte=yteNormlL, SIMPLIFY="array")
        
        return(corrsRKHS_te2)
      }, simplify="array")
      proc.time() - pm # 17 secs
      
      dim(corrsRKHS_te_boot)
      dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
      
      
      return(corrsRKHS_te_boot)
      
    }, simplify="array")}, simplify="array")
  proc.time() - pm # 15 mins mins
  dim(corrsRKHS_boot_xy2_norml)
  dimnames(corrsRKHS_boot_xy2_norml)[3:4] <- list(sigmaxSeq_xy, lambdaSeq)
  names(dimnames(corrsRKHS_boot_xy2_norml))[3:4] <- c("sigmax", "lambda")
  
  
  indxOpt_xy2_norml <- apply(corrsRKHS_boot_xy2_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_xy2_norml)[[1]] <- c("sigmax","lambda")
  names(dimnames(indxOpt_xy2_norml))[1] <- "indx"
  
  
  df_xy2_norml <- melt(indxOpt_xy2_norml)
  df_xy2_norml <- cast(df_xy2_norml, testPt+boot~indx, value="value")
  df_xy2_norml$corr <- corrsRKHS_boot_xy2_norml[as.matrix(df_xy2_norml[,names(dimnames(corrsRKHS_boot_xy2_norml))])]
  colnames(df_xy2_norml) <- c("testPt","boot","indxLambda","indxSigmax","corr")
  df_xy2_norml$node <- "y"
  df_xy2_norml$numReg <- "1"
  df_xy2_norml$lambda <- lambdaSeq[df_xy2_norml$indxLambda]
  df_xy2_norml$sigmax <- sigmaxSeq_xy[df_xy2_norml$indxSigmax]
  indxMat <- matrix(c(df_xy2$testPt,df_xy2_norml$boot),nrow(df_xy2_norml),2)
  df_xy2_norml$sigmay <- optSigmay_xy_norml[df_xy2_norml$testPt]
  df_xy2_norml$x <- X2[smpl_te2[indxMat]]
  df_xy2_norml$y <- apply(as.data.frame(df_xy2_norml), 1, function(row){
    # row <- df_xy_norml[1,]; names(row) <- colnames(df_xy_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- X2[smpl_te2[,boot],,drop=F]
    yTe_b <- Y2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(df_xy_norml$sigmax, df_xy2_norml$sigmax)
  plot(X2, Y2)
  plot(df_xy_norml$x, df_xy_norml$y)
  plot(df_xy2_norml$x, df_xy2_norml$y)
  plot(df_xy2_norml$x, log(df_xy2_norml$sigmax,10))
  plot(df_xy2_norml$x, log(df_xy2_norml$lambda,10))
  
  df_othr_xy_norml <- melt(indxOther_xy_norml)
  df_othr_xy_norml <- cast(df_othr_xy_norml, testPt+boot~indx, value="value")
  df_othr_xy_norml$corr <- corrsRKHS_boot_xy_norml[as.matrix(df_othr_xy_norml[,names(dimnames(corrsRKHS_boot_xy_norml))])]
  colnames(df_othr_xy_norml) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
  df_othr_xy_norml$node <- "x"
  df_othr_xy_norml$numReg <- "1"
  df_othr_xy_norml$lambda <- lambdaSeq[df_othr_xy_norml$indxLambda]
  df_othr_xy_norml$sigmax <- sigmaxSeq_xy[df_othr_xy_norml$indxSigmax]
  df_othr_xy_norml$sigmay <- sigmasy[df_othr_xy_norml$indxSigmay]
  df_xy_norml <- rbind(df_xy_norml[,1:11], df_othr_xy_norml)
  df_xy_norml <- df_xy_norml[order(df_xy_norml$node, df_xy_norml$numReg,df_xy_norml$testPt, df_xy_norml$boot),]
  
  df_yx_norml <- melt(indxOpt_yx_norml)
  df_yx_norml <- cast(df_yx_norml, testPt+boot~indx, value="value")
  df_yx_norml$corr <- corrsRKHS_boot_yx_norml[as.matrix(df_yx_norml[,names(dimnames(corrsRKHS_boot_yx_norml))])]
  colnames(df_yx_norml) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
  df_yx_norml$node <- "x"
  df_yx_norml$numReg <- "1"
  df_yx_norml$lambda <- lambdaSeq[df_yx_norml$indxLambda]
  df_yx_norml$sigmax <- sigmaxSeq_yx[df_yx_norml$indxSigmax]
  df_yx_norml$sigmay <- sigmasy[df_yx_norml$indxSigmay]
  
  
  indxMat <- matrix(c(df_yx_norml$testPt,df_yx_norml$boot),nrow(df_yx_norml),2)
  df_yx_norml$x <- Y2[smpl_te2[indxMat]]
  df_yx_norml$y <- apply(as.data.frame(df_yx_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- Y2[smpl_te2[,boot],,drop=F]
    yTe_b <- X2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(Y2, X2)
  plot(df_yx_norml$x, df_yx_norml$y)
  plot(df_yx_norml$x, log(df_yx_norml$sigmay,10))
  o <- order(Y2)
  mod_yx <- sapply(Y2[o,], function(x) mean(log(df_yx_norml$sigmay[which(df_yx_norml$x==x)],10)))
  lines(Y2[o,], mod_yx, col="red")
  mod_yx2 <- loess(mod_yx~Y2[o,], span=0.25)
  optSigmay_yx_norml <- predict(mod_yx2, Y2[o,])
  lines(Y2[o,], optSigmay_yx_norml, col="green")
  plot(Y2, X2)
  lines(Y2[o,], norml(mod_yx), col="red")
  lines(Y2[o,], norml(optSigmay_yx_norml), col="green")
  optSigmay_yx_norml <- 10^optSigmay_yx_norml
  optSigmay_yx_norml <- optSigmay_yx_norml[order(o)]
  plot(Y2, log(optSigmay_yx_norml,10))
  plot(Y2, log(optSigmay_yx,10))
  
  
  i <- 0
  j <- 0
  pm <- proc.time()
  corrsRKHS_boot_yx2_norml <- sapply(lambdaSeq, function(lam){ 
    i <<- i + 1
    j <<- 0
    sapply(sigmaxSeq_yx , function(sigmax){
      # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
      # lam <- 0.1; sigmax <- 100
      j <<- j + 1
      print(paste("i-j: ", i, "-",j))
      
      
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
        #k <- 1
        
        xTe_b2 <- Y2[smpl_te2[,k],,drop=F]
        optSigmay_yx_b <- optSigmay_yx_norml[smpl_te2[,k]]
        xTr_b2 <- Y2[smpl_tr2[,k],,drop=F]
        yTe_b2 <- as.matrix(X2[smpl_te2[,k],])
        yTr_b2 <- as.matrix(X2[smpl_tr2[,k],])
        nTr_b2 <- nrow(xTr_b2)
        nTe_b2 <- nrow(xTe_b2)
        Itr_b2 <- diag(nTr_b2)  
        
        Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
        Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
        Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
        Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
        
        ytrNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte_tr2))), mc.cores=1, SIMPLIFY=FALSE)
        yteNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte2))), mc.cores=1, SIMPLIFY=FALSE) 
        
        LB <- Lte_tr2%*%Blambda_tr2
        
        corrsRKHS_te2 <- mcmapply(function(sigmay, pt, ytr, yte){ 
          # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
          
          
          Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
          LBK <- LB%*%Ktr_te
          Lte_col <- Lte2[pt,]
          LBK_col <- LBK[pt,]
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          return(res)
        }, sigmay=optSigmay_yx_b, pt=1:length(optSigmay_yx_b), ytr=ytrNormlL, yte=yteNormlL, SIMPLIFY="array")
        
        return(corrsRKHS_te2)
      }, simplify="array")
      proc.time() - pm # 17 secs
      
      dim(corrsRKHS_te_boot)
      dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
      
      
      return(corrsRKHS_te_boot)
      
    }, simplify="array")}, simplify="array")
  proc.time() - pm # 15 mins mins
  dim(corrsRKHS_boot_yx2_norml)
  dimnames(corrsRKHS_boot_yx2_norml)[3:4] <- list(sigmaxSeq_yx, lambdaSeq)
  names(dimnames(corrsRKHS_boot_yx2_norml))[3:4] <- c("sigmax", "lambda")
  
  indxOpt_yx2_norml <- apply(corrsRKHS_boot_yx2_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_yx2_norml)[[1]] <- c("sigmax","lambda")
  names(dimnames(indxOpt_yx2_norml))[1] <- "indx"
  
  df_yx2_norml <- melt(indxOpt_yx2_norml)
  df_yx2_norml <- cast(df_yx2_norml, testPt+boot~indx, value="value")
  df_yx2_norml$corr <- corrsRKHS_boot_yx2_norml[as.matrix(df_yx2_norml[,names(dimnames(corrsRKHS_boot_yx2_norml))])]
  colnames(df_yx2_norml) <- c("testPt","boot","indxLambda","indxSigmax","corr")
  df_yx2_norml$node <- "x"
  df_yx2_norml$numReg <- "1"
  df_yx2_norml$lambda <- lambdaSeq[df_yx2_norml$indxLambda]
  df_yx2_norml$sigmax <- sigmaxSeq_yx[df_yx2_norml$indxSigmax]
  indxMat <- matrix(c(df_yx2_norml$testPt,df_yx2_norml$boot),nrow(df_yx2_norml),2)
  df_yx2_norml$sigmay <- optSigmay_yx_norml[df_yx2_norml$testPt]
  df_yx2_norml$x <- Y2[smpl_te2[indxMat]]
  df_yx2_norml$y <- apply(as.data.frame(df_yx2_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- Y2[smpl_te2[,boot],,drop=F]
    yTe_b <- X2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(Y2, X2)  
  plot(df_yx2_norml$x, df_yx2_norml$y)
  plot(df_yx2_norml$x, log(df_yx2_norml$sigmax,10))
  plot(df_yx2_norml$x, log(df_yx2_norml$lambda,10))
  
  
  df_othr_yx_norml <- melt(indxOther_yx_norml)
  df_othr_yx_norml <- cast(df_othr_yx_norml, testPt+boot~indx, value="value")
  df_othr_yx_norml$corr <- corrsRKHS_boot_yx_norml[as.matrix(df_othr_yx_norml[,names(dimnames(corrsRKHS_boot_yx_norml))])]
  colnames(df_othr_yx_norml) <- c("testPt","boot","indxLambda","indxSigmax","indxSigmay","corr")
  summary(df_othr_yx_norml$corr)
  df_othr_yx_norml$node <- "y"
  df_othr_yx_norml$numReg <- "1"
  df_othr_yx_norml$lambda <- lambdaSeq[df_othr_yx_norml$indxLambda]
  df_othr_yx_norml$sigmax <- sigmaxSeq_yx[df_othr_yx_norml$indxSigmax]
  df_othr_yx_norml$sigmay <- sigmasy[df_othr_yx_norml$indxSigmay]
  df_yx_norml <- rbind(df_yx_norml[,1:11], df_othr_yx_norml)
  
  dim(df_xy_norml); dim(df_yx_norml)
  table(df_xy_norml$sigmay); table(df_yx_norml$sigmay) 
  indx <- order(df_xy_norml$node, df_xy_norml$numReg,df_xy_norml$testPt, df_xy_norml$boot)
  df_xy_norml <- df_xy_norml[indx,]
  indx <- order(df_yx_norml$node, df_yx_norml$numReg,df_yx_norml$testPt, df_yx_norml$boot)
  df_yx_norml <- df_yx_norml[indx,]
  head(df_xy_norml); head(df_yx_norml)
  all(df_xy_norml$sigmay==df_yx_norml$sigmay)
  
  
  indxMat <- matrix(c(df_xy_norml$testPt,df_xy_norml$boot),nrow(df_xy_norml),2)
  df_xy_norml$x <- X2[smpl_te2[indxMat]]
  df_xy_norml$y <- apply(as.data.frame(df_xy_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- X2[smpl_te2[,boot],,drop=F]
    yTe_b <- Y2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  indxMat <- matrix(c(df_yx_norml$testPt,df_yx_norml$boot),nrow(df_yx_norml),2)
  df_yx_norml$x <- Y2[smpl_te2[indxMat]]
  df_yx_norml$y <- apply(as.data.frame(df_yx_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- Y2[smpl_te2[,boot],,drop=F]
    yTe_b <- X2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  
  df_xy_norml <- as.data.frame(df_xy_norml)
  df_yx_norml <- as.data.frame(df_yx_norml)
  
  i <- 0
  j <- 0
  pm <- proc.time()
  corrsRKHS_boot_xy2_othr_norml <- sapply(lambdaSeq, function(lam){ 
    i <<- i + 1
    j <<- 0
    sapply(sigmaxSeq_xy , function(sigmax){
      # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
      # lam <- 0.1; sigmax <- 100
      j <<- j + 1
      print(paste("i-j: ", i, "-",j))
      
      
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
        #k <- 1
        
        xTe_b2 <- X2[smpl_te2[,k],,drop=F]
        optSigmay_xy_b <- optSigmay_yx_norml[smpl_te2[,k]]
        xTr_b2 <- X2[smpl_tr2[,k],,drop=F]
        yTe_b2 <- as.matrix(Y2[smpl_te2[,k],])
        yTr_b2 <- as.matrix(Y2[smpl_tr2[,k],])
        nTr_b2 <- nrow(xTr_b2)
        nTe_b2 <- nrow(xTe_b2)
        Itr_b2 <- diag(nTr_b2)  
        
        Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
        Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
        Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
        Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
        
        ytrNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte_tr2))), mc.cores=1, SIMPLIFY=FALSE)
        yteNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte2))), mc.cores=1, SIMPLIFY=FALSE)
        
        LB <- Lte_tr2%*%Blambda_tr2
        
        corrsRKHS_te2 <- mcmapply(function(sigmay, pt, ytr, yte){ 
          # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
          
          
          Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
          LBK <- LB%*%Ktr_te
          Lte_col <- Lte2[pt,]
          LBK_col <- LBK[pt,]
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          return(res)
        }, sigmay=optSigmay_xy_b, pt=1:length(optSigmay_xy_b),ytr=ytrNormlL,yte=yteNormlL, SIMPLIFY="array")
        
        return(corrsRKHS_te2)
      }, simplify="array")
      proc.time() - pm # 17 secs
      
      dim(corrsRKHS_te_boot)
      dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
      
      
      return(corrsRKHS_te_boot)
      
    }, simplify="array")}, simplify="array")
  proc.time() - pm # 15 mins mins
  dim(corrsRKHS_boot_xy2_othr_norml)
  dimnames(corrsRKHS_boot_xy2_othr_norml)[3:4] <- list(sigmaxSeq_xy, lambdaSeq)
  names(dimnames(corrsRKHS_boot_xy2_othr_norml))[3:4] <- c("sigmax", "lambda")
  
  
  indxOpt_xy2_othr_norml <- apply(corrsRKHS_boot_xy2_othr_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_xy2_othr_norml)[[1]] <- c("sigmax","lambda")
  names(dimnames(indxOpt_xy2_othr_norml))[1] <- "indx"
  
  
  df_xy2_othr_norml <- melt(indxOpt_xy2_othr_norml)
  df_xy2_othr_norml <- cast(df_xy2_othr_norml, testPt+boot~indx, value="value")
  df_xy2_othr_norml$corr <- corrsRKHS_boot_xy2_othr_norml[as.matrix(df_xy2_othr_norml[,names(dimnames(corrsRKHS_boot_xy2_othr_norml))])]
  colnames(df_xy2_othr_norml) <- c("testPt","boot","indxLambda","indxSigmax","corr")
  df_xy2_othr_norml$node <- "x"
  df_xy2_othr_norml$numReg <- "1"
  df_xy2_othr_norml$lambda <- lambdaSeq[df_xy2_othr_norml$indxLambda]
  df_xy2_othr_norml$sigmax <- sigmaxSeq_xy[df_xy2_othr_norml$indxSigmax]
  indxMat <- matrix(c(df_xy2_othr_norml$testPt,df_xy2_othr_norml$boot),nrow(df_xy2_othr_norml),2)
  df_xy2_othr_norml$sigmay <- optSigmay_yx_norml[df_xy2_othr_norml$testPt]
  df_xy2_othr_norml$x <- X2[smpl_te2[indxMat]]
  df_xy2_othr_norml$y <- apply(as.data.frame(df_xy2_othr_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- X2[smpl_te2[,boot],,drop=F]
    yTe_b <- Y2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(X2, Y2)
  plot(df_xy2_othr_norml$x, df_xy2_othr_norml$y)
  plot(df_xy2_othr_norml$x, log(df_xy2_othr_norml$sigmax,10))
  plot(df_xy2_othr_norml$x, log(df_xy2_othr_norml$lambda,10))
  
  i <- 0
  j <- 0
  pm <- proc.time()
  corrsRKHS_boot_yx2_othr_norml <- sapply(lambdaSeq, function(lam){ 
    i <<- i + 1
    j <<- 0
    sapply(sigmaxSeq_yx , function(sigmax){
      # i <- 1; j <- 2; lam <- unique(lambdaSeq)[i]; sigmax <- unique(sigmaxSeq_xy)[j]
      # lam <- 0.1; sigmax <- 100
      j <<- j + 1
      print(paste("i-j: ", i, "-",j))
      
      
      pm <- proc.time()
      corrsRKHS_te_boot <- sapply(1:(nboots), function(k){
        #k <- 1
        
        xTe_b2 <- Y2[smpl_te2[,k],,drop=F]
        optSigmay_yx_b <- optSigmay_xy_norml[smpl_te2[,k]]
        xTr_b2 <- Y2[smpl_tr2[,k],,drop=F]
        yTe_b2 <- as.matrix(X2[smpl_te2[,k],])
        yTr_b2 <- as.matrix(X2[smpl_tr2[,k],])
        nTr_b2 <- nrow(xTr_b2)
        nTe_b2 <- nrow(xTe_b2)
        Itr_b2 <- diag(nTr_b2)  
        
        Ltr2 <- kern_rbf(xTr_b2, sigma=sigmax)
        Lte2 <- kern_rbf(xTe_b2, sigma=sigmax)
        Lte_tr2 <- kern_rbf(xTe_b2,xTr_b2, sigma=sigmax)
        Blambda_tr2 <- solve(Ltr2+nTr_b2*lam*Itr_b2)
        
        ytrNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTr_b2-weighted.mean(yTr_b2, w=ws))/sqrt(weighted.var(yTr_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte_tr2))), mc.cores=1, SIMPLIFY=FALSE)
        yteNormlL <- mcmapply(function(Lte_col){
          # Lte_col <- as.list(as.data.frame(t(Lte2)))[[1]]
          ws <- Lte_col/sum(Lte_col)
          res <- (yTe_b2-weighted.mean(yTe_b2, w=ws))/sqrt(weighted.var(yTe_b2, w=ws))
          #weighted.mean(res, w=ws); weighted.var(res, w=ws)
          return(res)
        }, Lte_col=as.list(as.data.frame(t(Lte2))), mc.cores=1, SIMPLIFY=FALSE)
        
        LB <- Lte_tr2%*%Blambda_tr2
        
        corrsRKHS_te2 <- mcmapply(function(sigmay, pt, ytr, yte){ 
          # l <- 1; sigmay <- optSigmay_xy_b[l]; pt <- l
          
          
          Ktr_te <- kern_rbf(ytr, yte, sigma=sigmay)
          LBK <- LB%*%Ktr_te
          Lte_col <- Lte2[pt,]
          LBK_col <- LBK[pt,]
          res <- weightedCorr(Lte_col,LBK_col, method="pearson", weights=Lte_col/sum(Lte_col))
          return(res)
        }, sigmay=optSigmay_yx_b, pt=1:length(optSigmay_yx_b),ytr=ytrNormlL,yte=yteNormlL,SIMPLIFY="array")
        
        return(corrsRKHS_te2)
      }, simplify="array")
      proc.time() - pm # 17 secs
      
      dim(corrsRKHS_te_boot)
      dimnames(corrsRKHS_te_boot) <- list(testPt=1:(numPerBoot), boot=1:nboots)
      
      
      return(corrsRKHS_te_boot)
      
    }, simplify="array")}, simplify="array")
  proc.time() - pm # 15 mins mins
  dim(corrsRKHS_boot_yx2_othr_norml)
  dimnames(corrsRKHS_boot_yx2_othr_norml)[3:4] <- list(sigmaxSeq_yx, lambdaSeq)
  names(dimnames(corrsRKHS_boot_yx2_othr_norml))[3:4] <- c("sigmax", "lambda")
  
  indxOpt_yx2_othr_norml <- apply(corrsRKHS_boot_yx2_othr_norml, c("testPt","boot"), function(arr){
    # arr <- corrsRKHS_boot_xy[1,1,,]
    maxArr <- arr==max(arr,na.rm=T)
    indxs <- sapply(c("sigmax","lambda"), function(margin) which(apply(maxArr, margin, any))[1])
    #corrsRKHS_boot_xy[matrix(c(1,indxs[1],1,indxs[2:3]),1,5)]
    return(indxs)
  })
  dimnames(indxOpt_yx2_othr_norml)[[1]] <- c("sigmax","lambda")
  names(dimnames(indxOpt_yx2_othr_norml))[1] <- "indx"
  
  df_yx2_othr_norml <- melt(indxOpt_yx2_othr_norml)
  df_yx2_othr_norml <- cast(df_yx2_othr_norml, testPt+boot~indx, value="value")
  df_yx2_othr_norml$corr <- corrsRKHS_boot_yx2_othr_norml[as.matrix(df_yx2_othr_norml[,names(dimnames(corrsRKHS_boot_yx2_othr_norml))])]
  colnames(df_yx2_othr_norml) <- c("testPt","boot","indxLambda","indxSigmax","corr")
  df_yx2_othr_norml$node <- "y"
  df_yx2_othr_norml$numReg <- "1"
  df_yx2_othr_norml$lambda <- lambdaSeq[df_yx2_othr_norml$indxLambda]
  df_yx2_othr_norml$sigmax <- sigmaxSeq_yx[df_yx2_othr_norml$indxSigmax]
  indxMat <- matrix(c(df_yx2_othr_norml$testPt,df_yx2_othr_norml$boot),nrow(df_yx2_othr_norml),2)
  df_yx2_othr_norml$sigmay <- optSigmay_xy_norml[df_yx2_othr_norml$testPt]
  df_yx2_othr_norml$x <- Y2[smpl_te2[indxMat]]
  df_yx2_othr_norml$y <- apply(as.data.frame(df_yx2_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- Y2[smpl_te2[,boot],,drop=F]
    yTe_b <- X2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  plot(Y2, X2)
  plot(df_yx2_othr_norml$x, df_yx2_othr_norml$y)
  plot(df_yx2_othr_norml$x, log(df_yx2_othr_norml$sigmax,10))
  plot(df_yx2_othr_norml$x, log(df_yx2_othr_norml$lambda,10))
  
  plot(log(optSigmay_xy,10), log(optSigmay_yx,10)); abline(a=0, b=1, col="red")
  plot(log(optSigmay_xy_norml,10), log(optSigmay_yx_norml,10)); abline(a=0, b=1, col="red")
  sum(optSigmay_yx>optSigmay_xy)/length(optSigmay_xy)
  sum(optSigmay_yx_norml>optSigmay_xy_norml)/length(optSigmay_xy_norml)
  plot(log(df_xy$sigmay[which(df_xy$node=="y")],10), log(df_yx$sigmay[which(df_xy$node=="x")],10)); abline(a=0, b=1, col="red")
  plot(log(df_xy_norml$sigmay[which(df_xy_norml$node=="y")],10), log(df_yx_norml$sigmay[which(df_xy_norml$node=="x")],10)); abline(a=0, b=1, col="red")
  sum(df_xy$sigmay[which(df_xy$node=="y")]<df_yx$sigmay[which(df_xy$node=="x")])/sum(df_xy$node=="x")
  sum(df_xy_norml$sigmay[which(df_xy_norml$node=="y")]<df_yx_norml$sigmay[which(df_xy_norml$node=="x")])/sum(df_xy_norml$node=="x")
  
  plot(X2,Y2, col=c("red","blue")[(optSigmay_yx>optSigmay_xy)*1+1])
  plot(X2,Y2, col=c("red","blue")[(optSigmay_yx_norml>optSigmay_xy_norml)*1+1])
  

  plot(X2, Y2)
  plot(df_xy_norml$x, df_xy_norml$y)
  plot(X2, log(optSigmay_yx,10), ylim=range(log(c(optSigmay_yx,optSigmay_xy),10)))
  lines(X2, log(optSigmay_xy,10), type="p", col="red")
  plot(X2, log(optSigmay_yx_norml,10), ylim=range(log(c(optSigmay_yx_norml,optSigmay_xy_norml),10)))
  lines(X2, log(optSigmay_xy_norml,10), type="p", col="red")
  plot(Y2, X2)
  plot(Y2, log(optSigmay_xy,10), ylim=range(log(c(optSigmay_yx,optSigmay_xy),10)))
  lines(Y2, log(optSigmay_yx,10), type="p", col="red")
  plot(Y2, log(optSigmay_xy_norml,10), ylim=range(log(c(optSigmay_yx_norml,optSigmay_xy_norml),10)))
  lines(Y2, log(optSigmay_yx_norml,10), type="p", col="red")
  
  
  
  df_xy2_norml <- rbind(df_xy2_norml, df_xy2_othr_norml)
  df_yx2_norml <- rbind(df_yx2_norml, df_yx2_othr_norml)
  
  dim(df_xy2_norml); dim(df_yx2_norml)
  indx <- order(df_xy2_norml$node, df_xy2_norml$numReg,df_xy2_norml$testPt, df_xy2_norml$boot)
  df_xy2_norml <- df_xy2_norml[indx,]
  indx <- order(df_yx2_norml$node, df_yx2_norml$numReg,df_yx2_norml$testPt, df_yx2_norml$boot)
  df_yx2_norml <- df_yx2_norml[indx,]
  head(df_xy2_norml); head(df_yx2_norml)
  all(df_xy2_norml$sigmay==df_yx2_norml$sigmay)
  
  indxMat <- matrix(c(df_xy2_norml$testPt,df_xy2_norml$boot),nrow(df_xy2_norml),2)
  
  df_xy2_norml$x <- X2[smpl_te2[indxMat]]
  df_xy2_norml$y <- apply(as.data.frame(df_xy2_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- X2[smpl_te2[,boot],,drop=F]
    yTe_b <- Y2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  
  indxMat <- matrix(c(df_yx2_norml$testPt,df_yx2_norml$boot),nrow(df_yx2_norml),2)
  df_yx2_norml$x <- Y2[smpl_te2[indxMat]]
  df_yx2_norml$y <- apply(as.data.frame(df_yx2_norml), 1, function(row){
    # row <- df_yx_norml[1,]; names(row) <- colnames(df_yx_norml) 
    sigmax <- as.numeric(row["sigmax"])
    boot <- as.numeric(row["boot"])
    testPt <- as.numeric(row["testPt"])
    xTe_b <- Y2[smpl_te2[,boot],,drop=F]
    yTe_b <- X2[smpl_te2[,boot],,drop=F]
    Lte <- kern_rbf(xTe_b, sigma=sigmax)
    weights <- Lte[,testPt]/sum(Lte[,testPt])
    yTe_b <- (yTe_b - weighted.mean(yTe_b, w=weights))/sqrt(weighted.var(yTe_b, w=weights))
    #mean(yTe_b); var(yTe_b); weighted.mean(yTe_b,w=weights); weighted.var(yTe_b, w=weights)
    res <- yTe_b[testPt]
    return(res)
  })
  df_xy2_norml <- as.data.frame(df_xy2_norml)
  df_yx2_norml <- as.data.frame(df_yx2_norml)
  
}



indx0 <- which(df_xy$node=="x")
indx1 <- which(df_xy$node=="y")


# 2 issues: 
#1-more local sigma adds variability (less data) which we're not interested in
#2- shd we actually normalize effect locally

table(df_xy$sigmay[indx0])
table(df_yx$sigmay[indx0])
table(df_xy$sigmay[indx1])
table(df_yx$sigmay[indx1])

table(df_xy$sigmax)
table(df_yx$sigmax)
table(df_xy$sigmax,df_yx$sigmax)
plot(table(df_xy$sigmax,df_yx$sigmax), xlab="x->y",ylab="x<-y")

table(df_xy$lambda)
table(df_yx$lambda)


plot(df_xy$corr,df_yx$corr); abline(a=0, b=1, col="red")
plot(df_xy2$corr,df_yx2$corr); abline(a=0, b=1, col="red")
plot(df_xy_norml$corr,df_yx_norml$corr); abline(a=0, b=1, col="red")
plot(df_xy2_norml$corr,df_yx2_norml$corr); abline(a=0, b=1, col="red")
mean(df_xy$corr); mean(df_yx$corr)
mean(df_xy2$corr); mean(df_yx2$corr)
mean(df_xy_norml$corr); mean(df_yx_norml$corr)
mean(df_xy2_norml$corr); mean(df_yx2_norml$corr)

sum(df_xy$corr > df_yx$corr)/nrow(df_xy)
sum(df_xy2$corr > df_yx2$corr)/nrow(df_xy2)
sum(df_xy_norml$corr > df_yx_norml$corr)/nrow(df_xy_norml)
sum(df_xy2_norml$corr > df_yx2_norml$corr)/nrow(df_xy2_norml)
sum(df_xy$corr[indx0] > df_yx$corr[indx0])/length(indx0)
sum(df_xy$corr[indx1] > df_yx$corr[indx1])/length(indx1)
sum(df_xy_norml$corr[indx0] > df_yx_norml$corr[indx0])/length(indx0)
sum(df_xy_norml$corr[indx1] > df_yx_norml$corr[indx1])/length(indx1)

# block 65 SIM: 0.66, 0.53 for  x->y, x<-y 
# block 15 SIM: 0.833, 0.7968  for  x->y, x<-y 
# block 81 CMPLX;   0.4934  , 0.5486   for x->y, x<-y
# block 7 SIM: 0.9596129, 0.9607941  x->y, x<-y
# block 8 SIM: 0.3929798, 0.448001  x->y, x<-y
# block 9 SIM: 0.6893967, 0.6050134  x->y, x<-y
# block 11 SIM: 0.9466006, 0.9389972    x->y, x<-y
# block 14 SIM: 0.9912979, 0.9901761    x->y, x<-y
# block 15 SIM: 0.8421221, 0.8670485    x->y, x<-y
# block 16 SIM: 0.7970918, 0.7641301      x->y, x<-y
# block 21 SIM: 0.6653446, 0.6126333      x->y, x<-y
# block  1 SIM: 0.8284874, 0.8082067      x->y, x<-y
# block  33 SIM: 0.9383226, 0.9271018      x->y, x<-y
# block 99 SIM: 0.8001908, 0.7675734
# block 50 SIM: 0.7591844, 0.7962408 

KCRDCs_xy <- apply(df_xy, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  res <- res[testPt]/mean(res)
  return(res)
})
KCRDCs_yx <- apply(df_yx, 1, function(row){
  # row <- as.numeric(df_yx[1,]); names(row) <- colnames(df_yx) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b) 
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  res <- res[testPt]/mean(res)
  return(res)
})
KCRDCs_xy2 <- apply(df_xy2, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b) 
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  res <- res[testPt]/mean(res)
  return(res)
})
KCRDCs_yx2 <- apply(df_yx2, 1, function(row){
  # row <- as.numeric(df_yx[1,]); names(row) <- colnames(df_yx) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  res <- res[testPt]/mean(res)
  return(res)
})


#plot(sort(KCRDCs_xy), sort(cmemSet$y$`1`[,"kcrdc"]))

hist(KCRDCs_xy);hist(KCRDCs_yx) 
plot(KCRDCs_xy, KCRDCs_yx); abline(a=0, b=1, col="red")
plot(KCRDCs_xy[indx0], KCRDCs_yx[indx0]); abline(a=0, b=1, col="red")
plot(KCRDCs_xy[indx1], KCRDCs_yx[indx1]); abline(a=0, b=1, col="red")
sum(KCRDCs_xy<KCRDCs_yx)/length(KCRDCs_xy)
sum(KCRDCs_xy[indx0]<KCRDCs_yx[indx0])/length(indx0)
sum(KCRDCs_xy[indx1]<KCRDCs_yx[indx1])/length(indx1)
sum(KCRDCs_xy-KCRDCs_yx)
sum(KCRDCs_xy[indx0]-KCRDCs_yx[indx0])
sum(KCRDCs_xy[indx1]-KCRDCs_yx[indx1])

plot(KCRDCs_xy2, KCRDCs_yx2); abline(a=0, b=1, col="red")
plot(KCRDCs_xy2[indx0], KCRDCs_yx2[indx0]); abline(a=0, b=1, col="red")
plot(KCRDCs_xy2[indx1], KCRDCs_yx2[indx1]); abline(a=0, b=1, col="red")
sum(KCRDCs_xy2<KCRDCs_yx2)/length(KCRDCs_xy2)
sum(KCRDCs_xy2[indx0]<KCRDCs_yx2[indx0])/length(indx0)
sum(KCRDCs_xy2[indx1]<KCRDCs_yx2[indx1])/length(indx1)
sum(KCRDCs_xy2-KCRDCs_yx2)
sum(KCRDCs_xy2[indx0]-KCRDCs_yx2[indx0])
sum(KCRDCs_xy2[indx1]-KCRDCs_yx2[indx1])

var(KCRDCs_xy); var(KCRDCs_yx)
var(KCRDCs_xy[indx0]); var(KCRDCs_yx[indx0])
var(KCRDCs_xy[indx1]); var(KCRDCs_yx[indx1])
var(KCRDCs_xy2); var(KCRDCs_yx2)
var(KCRDCs_xy2[indx0]); var(KCRDCs_yx2[indx0])
var(KCRDCs_xy2[indx1]); var(KCRDCs_yx2[indx1])


plot(df_xy$x,df_xy$y, col=c("red","blue")[(KCRDCs_xy<KCRDCs_yx)*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(df_xy$corr>df_yx$corr)*1+1])

realKCDCs_xy <- apply(df_xy, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  yTe_b <- as.matrix(Y2[smpl_te2[,boot],])
  nTr_b <- nrow(yTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  nTe_b <- nrow(yTe_b)
  Ite_b <- diag(nTe_b)  
  Hte <- Ite_b-matrix(1/nTe_b, nTe_b, nTe_b)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmax)
  Kte_tr <- Hte %*% Kte_tr %*% Htr
  res <- diag(Kte_tr%*%t(Kte_tr))    
  res <- res^0.5
  res <- res[testPt]
  return(res)
})
realKCDCs_yx <- apply(df_yx, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  yTe_b <- as.matrix(X2[smpl_te2[,boot],])
  nTr_b <- nrow(yTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  nTe_b <- nrow(yTe_b)
  Ite_b <- diag(nTe_b)  
  Hte <- Ite_b-matrix(1/nTe_b, nTe_b, nTe_b)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmax)
  Kte_tr <- Hte %*% Kte_tr %*% Htr
  res <- diag(Kte_tr%*%t(Kte_tr))    
  res <- res^0.5
  res <- res[testPt]
  return(res)
})
realKCDCs_xy2 <- apply(df_xy2, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  yTe_b <- as.matrix(Y2[smpl_te2[,boot],])
  nTr_b <- nrow(yTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  nTe_b <- nrow(yTe_b)
  Ite_b <- diag(nTe_b)  
  Hte <- Ite_b-matrix(1/nTe_b, nTe_b, nTe_b)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmax)
  Kte_tr <- Hte %*% Kte_tr %*% Htr
  res <- diag(Kte_tr%*%t(Kte_tr))    
  res <- res^0.5
  res <- res[testPt]
  return(res)
})
realKCDCs_yx2 <- apply(df_yx2, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  yTe_b <- as.matrix(X2[smpl_te2[,boot],])
  nTr_b <- nrow(yTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  nTe_b <- nrow(yTe_b)
  Ite_b <- diag(nTe_b)  
  Hte <- Ite_b-matrix(1/nTe_b, nTe_b, nTe_b)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmax)
  Kte_tr <- Hte %*% Kte_tr %*% Htr
  res <- diag(Kte_tr%*%t(Kte_tr))    
  res <- res^0.5
  res <- res[testPt]
  return(res)
})

plot(KCDCs_xy, realKCDCs_xy); abline(a=0, b=1, col="red")

KCDCs_xy <- apply(df_xy, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay)
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  #Alambda <- Blambda_tr%*%Ktr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  #res <- diag(Ltr%*%Alambda%*%t(Ltr))
  res <- res^0.5
  #res <- res[testPt]
  modDens <- kepdf(xTe_b, eval.points = xTe_b, kernel = "gaussian", bwtype = "adaptive")
  weights1 <- Lte[,testPt]/sum(Lte[,testPt])
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  res <- weighted.var(res, w=weights)/weighted.mean(res, w=weights)
  return(res)
})
KCDCs_yx <- apply(df_yx, 1, function(row){
  # row <- as.numeric(df_yx[1,]); names(row) <- colnames(df_yx) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  #Alambda <- Blambda_tr%*%Ktr%*%Ktr%*%t(Blambda_tr)
  #res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))
  res <- diag(Ltr%*%Alambda%*%t(Ltr))
  res <- res^0.5
  modDens <- kepdf(xTe_b, eval.points = xTe_b, kernel = "gaussian", bwtype = "adaptive")
  weights1 <- Lte[,testPt]/sum(Lte[,testPt])
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #res <- res[testPt]
  res <- weighted.var(res, w=weights)/weighted.mean(res, w=weights)
  return(res)
})
KCDCs_xy2 <- apply(df_xy2, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b) 
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  #Alambda <- Blambda_tr%*%Ktr%*%Ktr%*%t(Blambda_tr)
  #res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))
  res <- diag(Ltr%*%Alambda%*%t(Ltr))
  res <- res^0.5
  modDens <- kepdf(xTe_b, eval.points = xTe_b, kernel = "gaussian", bwtype = "adaptive")
  weights1 <- Lte[,testPt]/sum(Lte[,testPt])
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #res <- res[testPt]
  res <- weighted.var(res, w=weights)/weighted.mean(res, w=weights)
  return(res)
})
KCDCs_yx2 <- apply(df_yx2, 1, function(row){
  # row <- as.numeric(df_yx[1,]); names(row) <- colnames(df_yx) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Htr <- Itr_b-matrix(1/nTr_b, nTr_b, nTr_b)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Ktr <- Htr %*% Ktr %*% Htr
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  #Alambda <- Blambda_tr%*%Ktr%*%Ktr%*%t(Blambda_tr)
  #res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))
  res <- diag(Ltr%*%Alambda%*%t(Ltr))
  res <- res^0.5
  modDens <- kepdf(xTe_b, eval.points = xTe_b, kernel = "gaussian", bwtype = "adaptive")
  weights1 <- Lte[,testPt]/sum(Lte[,testPt])
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #res <- res[testPt]
  res <- weighted.var(res, w=weights)/weighted.mean(res, w=weights)
  return(res)
})

plot(KCDCs_xy, KCDCs_yx); abline(a=0, b=1, col="red")
plot(KCDCs_xy[indx0], KCDCs_yx[indx0]); abline(a=0, b=1, col="red")
plot(KCDCs_xy[indx1], KCDCs_yx[indx1]); abline(a=0, b=1, col="red")
plot(KCDCs_xy[indx1], KCDCs_yx[indx0]); abline(a=0, b=1, col="red")
sum(KCDCs_xy<KCDCs_yx)/length(KCDCs_xy)
sum(KCDCs_xy[indx0]<KCDCs_yx[indx0])/length(indx0)
sum(KCDCs_xy[indx1]<KCDCs_yx[indx1])/length(indx1)
sum(KCDCs_xy[indx1]<KCDCs_yx[indx0])/length(indx1)
sum(KCDCs_xy[indx0]<KCDCs_yx[indx1])/length(indx1)
sum(KCDCs_xy-KCDCs_yx)
sum(KCDCs_xy[indx0]-KCDCs_yx[indx0])
sum(KCDCs_xy[indx1]-KCDCs_yx[indx1])
sum(KCDCs_xy[indx1]-KCDCs_yx[indx0])

plot(KCDCs_xy2, KCDCs_yx2); abline(a=0, b=1, col="red")
plot(KCDCs_xy2[indx0], KCDCs_yx2[indx0]); abline(a=0, b=1, col="red")
plot(KCDCs_xy2[indx1], KCDCs_yx2[indx1]); abline(a=0, b=1, col="red")
sum(KCDCs_xy2<KCDCs_yx2)/length(KCDCs_xy2)
sum(KCDCs_xy2[indx0]<KCDCs_yx2[indx0])/length(indx0)
sum(KCDCs_xy2[indx1]<KCDCs_yx2[indx1])/length(indx1)
sum(KCDCs_xy2-KCDCs_yx2)
sum(KCDCs_xy2[indx0]-KCDCs_yx2[indx0])
sum(KCDCs_xy2[indx1]-KCDCs_yx2[indx1])

var(KCDCs_xy); var(KCDCs_yx)
var(KCDCs_xy[indx0]); var(KCDCs_yx[indx0])
var(KCDCs_xy[indx1]); var(KCDCs_yx[indx1])
var(KCDCs_xy2); var(KCDCs_yx2)
var(KCDCs_xy2[indx0]); var(KCDCs_yx2[indx0])
var(KCDCs_xy2[indx1]); var(KCDCs_yx2[indx1])


plot(df_xy$x,df_xy$y, col=c("red","blue")[(KCDCs_xy2<KCDCs_yx2)*1+1])


# block 65 SIM 0.02, 0.1253  x->y, x<-y
# block 15 SIM 0.0404, 0.0153  x->y, x<-y
# block 81 CMPLX 0.02088, 0.0110926  x->y, x<-y
# block 7 SIM: 0.001700339, 0.00243697
# block 8 SIM: 0.03218845, 0.01891072
# block 9 SIM: 0.05235915, 0.02263131
# block 11 SIM: 0.008177681 , 0.004877557
# block 14 SIM: 0.0003406048, 0.0005657612
# block 15 SIM: 0.009029916, 0.008595049
# block 16 SIM: 0.02027921, 0.005813614
# block 21 SIM: 0.03015909, 0.03207162
# block 1 SIM: 0.01924026, 0.02259842
# block 33 SIM: 0.009103923, 0.01471695
# block 99 SIM: 0.007616867, 0.02878545
# block 50 SIM: 0.04379022, 0.03299293

wKCRDCs_xy <- apply(df_xy, 1, function(row){
  # row <- as.numeric(dfBestLambda[1,]); names(row) <- colnames(dfBestLambda) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  res <- res[testPt]/weighted.mean(res, w=weights)
  return(res)
})
wKCRDCs_yx <- apply(df_yx, 1, function(row){
  # row <- as.numeric(dfBestLambda[1,]); names(row) <- colnames(dfBestLambda) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  res <- res[testPt]/weighted.mean(res, w=weights)
  return(res)
})
wKCRDCs_xy2 <- apply(df_xy2, 1, function(row){
  # row <- as.numeric(dfBestLambda[1,]); names(row) <- colnames(dfBestLambda) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  res <- res[testPt]/weighted.mean(res, w=weights)
  return(res)
})
wKCRDCs_yx2 <- apply(df_yx2, 1, function(row){
  # row <- as.numeric(dfBestLambda[1,]); names(row) <- colnames(dfBestLambda) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+nTr_b*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay) 
  Alambda <- Blambda_tr%*%Ktr%*%t(Blambda_tr)
  res <- diag(Lte_tr%*%Alambda%*%t(Lte_tr))    
  res <- res^0.5
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  res <- res[testPt]/weighted.mean(res, w=weights)
  return(res)
})

hist(wKCRDCs_xy); hist(wKCRDCs_yx)
plot(wKCRDCs_xy, wKCRDCs_yx); abline(a=0, b=1, col="red")
plot(wKCRDCs_xy[indx0], wKCRDCs_yx[indx0]); abline(a=0, b=1, col="red")
plot(wKCRDCs_xy[indx1], wKCRDCs_yx[indx1]); abline(a=0, b=1, col="red")
sum(wKCRDCs_xy<wKCRDCs_yx)/length(wKCRDCs_xy)
sum(wKCRDCs_xy[indx0]<wKCRDCs_yx[indx0])/length(indx0)
sum(wKCRDCs_xy[indx1]<wKCRDCs_yx[indx1])/length(indx0)
sum(wKCRDCs_xy-wKCRDCs_yx)
sum(wKCRDCs_xy[indx0]-wKCRDCs_yx[indx0])
sum(wKCRDCs_xy[indx1]-wKCRDCs_yx[indx1])

plot(wKCRDCs_xy2, wKCRDCs_yx2); abline(a=0, b=1, col="red")
plot(wKCRDCs_xy2[indx0], wKCRDCs_yx2[indx0]); abline(a=0, b=1, col="red")
plot(wKCRDCs_xy2[indx1], wKCRDCs_yx2[indx1]); abline(a=0, b=1, col="red")
sum(wKCRDCs_xy2<wKCRDCs_yx2)/length(wKCRDCs_xy2)
sum(wKCRDCs_xy2[indx0]<wKCRDCs_yx2[indx0])/length(indx0)
sum(wKCRDCs_xy2[indx1]<wKCRDCs_yx2[indx1])/length(indx0)
sum(wKCRDCs_xy2-wKCRDCs_yx2)
sum(wKCRDCs_xy2[indx0]-wKCRDCs_yx2[indx0])
sum(wKCRDCs_xy2[indx1]-wKCRDCs_yx2[indx1])


var(wKCRDCs_xy); var(wKCRDCs_yx)
var(wKCRDCs_xy[indx0]); var(wKCRDCs_yx[indx0])
var(wKCRDCs_xy[indx1]); var(wKCRDCs_yx[indx1])

var(wKCRDCs_xy2); var(wKCRDCs_yx2)
var(wKCRDCs_xy2[indx0]); var(wKCRDCs_yx2[indx0])
var(wKCRDCs_xy2[indx1]); var(wKCRDCs_yx2[indx1])


plot(df_xy$x,df_xy$y, col=c("red","blue")[(wKCRDCs_xy<wKCRDCs_yx)*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wKCRDCs_xy2<wKCRDCs_yx2)*1+1])



# block 65 SIM ,   x->y, x<-y
# block 15 SIM 0.0161, 0.00867   x->y, x<-y
# block 81 CMPLX  0.007136, 0.004443494  x->y, x<-y
# block 7 SIM: 0.001092074, 0.001202124
# block 8 SIM: 0.006907666, 0.007798203
# block 9 SIM: 0.008293172, 0.01078343
# block 11 SIM: 0.001763947, 0.003696441
# block 14 SIM: 0.0002028876, 0.0003575043
# block 15 SIM: 0.005439272, 0.005394823
# block 16 SIM: 0.006696513, 0.002904347
# block 21 SIM: 0.01101336, 0.0111455 
# block 1 SIM: 0.01138606, 0.008067374  
# block 33 SIM: 0.003968222, 0.009271898  
# block 99 SIM: 0.003297856, 0.006440523
# block 50 SIM: 0.01386305, 0.008220684

wwKCRDCs_xyspot <- sapply(df_xy$x, function(xpt){
  # xpt <- df_xy$x
  var(wKCRDCs_xy[which(df_xy$x==xpt)])
})
wwKCRDCs_yxspot <- sapply(df_yx$x, function(xpt){
  # xpt <- df_xy$x
  var(wKCRDCs_yx[which(df_yx$x==xpt)])
})
wwKCRDCs_xyspot2 <- sapply(df_xy2$x, function(xpt){
  # xpt <- df_xy$x
  var(wKCRDCs_xy2[which(df_xy2$x==xpt)])
})
wwKCRDCs_yxspot2 <- sapply(df_yx2$x, function(xpt){
  # xpt <- df_xy$x
  var(wKCRDCs_yx2[which(df_yx2$x==xpt)])
})

plot(wwKCRDCs_xyspot, wwKCRDCs_yxspot); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xyspot2, wwKCRDCs_yxspot2); abline(a=0, b=1, col="red")
sum(wwKCRDCs_xyspot<wwKCRDCs_yxspot)/length(wwKCRDCs_xyspot)
sum(wwKCRDCs_xyspot2<wwKCRDCs_yxspot2)/length(wwKCRDCs_xyspot2)

library(pdfCluster)
modDens_xy0 <- kepdf(df_xy$x[indx0], eval.points = df_xy$x[indx0], kernel = "gaussian", bwtype = "adaptive")
wwKCRDCs_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_xy$x[indx0], xTr)]
  weights1 <- weights1/sum(weights1)
  
  #hist(df_xy$x[indx0], prob=T)
  #o <- order(modDens_xy0@eval.points)
  #lines(modDens_xy0@eval.points[o], modDens_xy0@estimate[o], col="red")
  weights2 <- 1/modDens_xy0@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #plot(df_xy$x[indx0], weights1); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights2); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights); abline(v=xPt, col="red")
  
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy[indx0],w=weights)
    
  return(res)
})
modDens_yx0 <- kepdf(df_yx$x[indx0], eval.points = df_yx$x[indx0], kernel = "gaussian", bwtype = "adaptive")
wwKCRDCs_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_yx$x[indx0], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_yx0@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_yx[indx0],w=weights)
  
  return(res)
})
modDens_xy1 <- kepdf(df_xy$x[indx1], eval.points = df_xy$x[indx1], kernel = "gaussian", bwtype = "adaptive")
wwKCRDCs_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_xy$x[indx1], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_xy1@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy[indx1],w=weights)
  
  return(res)
})
modDens_yx1 <- kepdf(df_yx$x[indx1], eval.points = df_yx$x[indx1], kernel = "gaussian", bwtype = "adaptive")
wwKCRDCs_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_yx$x[indx1], xTr)]
  weights1 <- weights1/sum(weights1)
  # hist(weights)
  weights2 <- 1/modDens_yx1@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  res <- weighted.var(wKCRDCs_yx[indx1],w=weights)
  
  return(res)
})
wwKCRDCs_xy <- c(wwKCRDCs_xy0, wwKCRDCs_xy1)
wwKCRDCs_yx <- c(wwKCRDCs_yx0, wwKCRDCs_yx1)

plot(wwKCRDCs_xy, wwKCRDCs_yx); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy[indx0], wwKCRDCs_yx[indx0]); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy[indx1], wwKCRDCs_yx[indx1]); abline(a=0, b=1, col="red")
sum(wwKCRDCs_xy<wwKCRDCs_yx)/length(wwKCRDCs_xy)
sum(wwKCRDCs_xy[indx0]<wwKCRDCs_yx[indx0])/length(wwKCRDCs_xy[indx0])
sum(wwKCRDCs_xy[indx1]<wwKCRDCs_yx[indx1])/length(wwKCRDCs_xy[indx1])
sum(wwKCRDCs_xy-wwKCRDCs_yx)
sum(wwKCRDCs_xy[indx0]-wwKCRDCs_yx[indx0])
sum(wwKCRDCs_xy[indx1]-wwKCRDCs_yx[indx1])

max(wwKCRDCs_xy[indx0]); max(wwKCRDCs_yx[indx0])
quantile(wwKCRDCs_xy[indx0], 0.9);quantile(wwKCRDCs_yx[indx0], 0.9)
median(wwKCRDCs_xy[indx0]);median(wwKCRDCs_yx[indx0])

max(wwKCRDCs_xy[indx1]); max(wwKCRDCs_yx[indx1])
quantile(wwKCRDCs_xy[indx1], 0.9);quantile(wwKCRDCs_yx[indx1], 0.9)
median(wwKCRDCs_xy[indx1]);median(wwKCRDCs_yx[indx1])

plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy<wwKCRDCs_yx)*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy[indx0]<wwKCRDCs_yx[indx0])*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy[indx1]<wwKCRDCs_yx[indx1])*1+1])

wwKCRDCs_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy2[i,]); names(row) <- colnames(df_xy2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_xy2$x[indx0], xTr)]
  weights1 <- weights1/sum(weights1)
  
  #hist(df_xy$x[indx0], prob=T)
  #o <- order(modDens_xy0@eval.points)
  #lines(modDens_xy0@eval.points[o], modDens_xy0@estimate[o], col="red")
  weights2 <- 1/modDens_xy0@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #plot(df_xy$x[indx0], weights1); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights2); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights); abline(v=xPt, col="red")
  
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy2[indx0],w=weights)
  
  return(res)
})
wwKCRDCs_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx2[i,]); names(row) <- colnames(df_yx2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_yx2$x[indx0], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_yx0@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_yx2[indx0],w=weights)
  
  return(res)
})
wwKCRDCs_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy2[i,]); names(row) <- colnames(df_xy2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_xy2$x[indx1], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_xy1@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy2[indx1],w=weights)
  
  return(res)
})
wwKCRDCs_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx2[i,]); names(row) <- colnames(df_yx2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights1 <- Lte_tr[,match(df_yx2$x[indx1], xTr)]
  weights1 <- weights1/sum(weights1)
  # hist(weights)
  weights2 <- 1/modDens_yx1@estimate
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  res <- weighted.var(wKCRDCs_yx2[indx1],w=weights)
  
  return(res)
})
wwKCRDCs_xy2 <- c(wwKCRDCs_xy0, wwKCRDCs_xy1)
wwKCRDCs_yx2 <- c(wwKCRDCs_yx0, wwKCRDCs_yx1)

plot(wwKCRDCs_xy2, wwKCRDCs_yx2); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy2[indx0], wwKCRDCs_yx2[indx0]); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy2[indx1], wwKCRDCs_yx2[indx1]); abline(a=0, b=1, col="red")
sum(wwKCRDCs_xy2<wwKCRDCs_yx2)/length(wwKCRDCs_xy2)
sum(wwKCRDCs_xy2[indx0]<wwKCRDCs_yx2[indx0])/length(wwKCRDCs_xy2[indx0])
sum(wwKCRDCs_xy2[indx1]<wwKCRDCs_yx2[indx1])/length(wwKCRDCs_xy2[indx1])
sum(wwKCRDCs_xy2-wwKCRDCs_yx2)
sum(wwKCRDCs_xy2[indx0]-wwKCRDCs_yx2[indx0])
sum(wwKCRDCs_xy2[indx1]-wwKCRDCs_yx2[indx1])

max(wwKCRDCs_xy2[indx0]); max(wwKCRDCs_yx2[indx0])
quantile(wwKCRDCs_xy2[indx0], 0.9);quantile(wwKCRDCs_yx2[indx0], 0.9)
median(wwKCRDCs_xy2[indx0]);median(wwKCRDCs_yx2[indx0])

max(wwKCRDCs_xy2[indx1]); max(wwKCRDCs_yx2[indx1])
quantile(wwKCRDCs_xy2[indx1], 0.9);quantile(wwKCRDCs_yx2[indx1], 0.9)
median(wwKCRDCs_xy2[indx1]);median(wwKCRDCs_yx2[indx1])

plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy2<wwKCRDCs_yx2)*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy2[indx0]<wwKCRDCs_yx2[indx0])*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy2[indx1]<wwKCRDCs_yx2[indx1])*1+1])


wwKCRDCs_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx0_2 <- which(df_xy$boot[indx0]==boot)
  
  weights1 <- Lte_tr[,match(df_xy$x[indx0][indx0_2], xTr)]
  weights1 <- weights1/sum(weights1)
  
  #hist(df_xy$x[indx0], prob=T)
  #o <- order(modDens_xy0@eval.points)
  #lines(modDens_xy0@eval.points[o], modDens_xy0@estimate[o], col="red")
  
  weights2 <- 1/modDens_xy0@estimate[indx0_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #plot(df_xy$x[indx0], weights1); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights2); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights); abline(v=xPt, col="red")
  
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy[indx0_2],w=weights)
  
  return(res)
})
wwKCRDCs_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx0_2 <- which(df_yx$boot[indx0]==boot)
  weights1 <- Lte_tr[,match(df_yx$x[indx0][indx0_2], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_yx0@estimate[indx0_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_yx[indx0_2],w=weights)
  
  return(res)
})
wwKCRDCs_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx1_2 <- which(df_xy$boot[indx1]==boot)
  weights1 <- Lte_tr[,match(df_xy$x[indx1][indx1_2], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_xy1@estimate[indx1_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy[indx1_2],w=weights)
  
  return(res)
})
wwKCRDCs_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx1_2 <- which(df_yx$boot[indx1]==boot)
  weights1 <- Lte_tr[,match(df_yx$x[indx1][indx1_2], xTr)]
  weights1 <- weights1/sum(weights1)
  # hist(weights)
  weights2 <- 1/modDens_yx1@estimate[indx1_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  res <- weighted.var(wKCRDCs_yx[indx1_2],w=weights)
  
  return(res)
})
wwKCRDCs_xy_2 <- c(wwKCRDCs_xy0, wwKCRDCs_xy1)
wwKCRDCs_yx_2 <- c(wwKCRDCs_yx0, wwKCRDCs_yx1)

plot(wwKCRDCs_xy_2, wwKCRDCs_yx_2); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy_2[indx0], wwKCRDCs_yx_2[indx0]); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy_2[indx1], wwKCRDCs_yx_2[indx1]); abline(a=0, b=1, col="red")
sum(wwKCRDCs_xy_2<wwKCRDCs_yx_2)/length(wwKCRDCs_xy_2)
sum(wwKCRDCs_xy_2[indx0]<wwKCRDCs_yx_2[indx0])/length(indx0)
sum(wwKCRDCs_xy_2[indx1]<wwKCRDCs_yx_2[indx1])/length(indx1)
sum(wwKCRDCs_xy_2-wwKCRDCs_yx_2)
sum(wwKCRDCs_xy_2[indx0]-wwKCRDCs_yx_2[indx0])
sum(wwKCRDCs_xy_2[indx1]-wwKCRDCs_yx_2[indx1])

wwKCRDCs_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy2[i,]); names(row) <- colnames(df_xy2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx0_2 <- which(df_xy2$boot[indx0]==boot)
  
  weights1 <- Lte_tr[,match(df_xy2$x[indx0][indx0_2], xTr)]
  weights1 <- weights1/sum(weights1)
  
  #hist(df_xy$x[indx0], prob=T)
  #o <- order(modDens_xy0@eval.points)
  #lines(modDens_xy0@eval.points[o], modDens_xy0@estimate[o], col="red")
  
  weights2 <- 1/modDens_xy0@estimate[indx0_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  #plot(df_xy$x[indx0], weights1); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights2); abline(v=xPt, col="red")
  #plot(df_xy$x[indx0], weights); abline(v=xPt, col="red")
  
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy2[indx0_2],w=weights)
  
  return(res)
})
wwKCRDCs_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx2[i,]); names(row) <- colnames(df_yx2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx0_2 <- which(df_yx2$boot[indx0]==boot)
  weights1 <- Lte_tr[,match(df_yx2$x[indx0][indx0_2], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_yx0@estimate[indx0_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_yx2[indx0_2],w=weights)
  
  return(res)
})
wwKCRDCs_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy2[i,]); names(row) <- colnames(df_xy2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx1_2 <- which(df_xy2$boot[indx1]==boot)
  weights1 <- Lte_tr[,match(df_xy2$x[indx1][indx1_2], xTr)]
  weights1 <- weights1/sum(weights1)
  weights2 <- 1/modDens_xy1@estimate[indx1_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  # hist(weights)
  res <- weighted.var(wKCRDCs_xy2[indx1_2],w=weights)
  
  return(res)
})
wwKCRDCs_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx2[i,]); names(row) <- colnames(df_yx2)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  indx1_2 <- which(df_yx2$boot[indx1]==boot)
  weights1 <- Lte_tr[,match(df_yx2$x[indx1][indx1_2], xTr)]
  weights1 <- weights1/sum(weights1)
  # hist(weights)
  weights2 <- 1/modDens_yx1@estimate[indx1_2]
  weights2 <- weights2/sum(weights2)
  weights <- weights1*weights2
  weights <- weights/sum(weights)
  res <- weighted.var(wKCRDCs_yx2[indx1_2],w=weights)
  
  return(res)
})
wwKCRDCs_xy2_2 <- c(wwKCRDCs_xy0, wwKCRDCs_xy1)
wwKCRDCs_yx2_2 <- c(wwKCRDCs_yx0, wwKCRDCs_yx1)

plot(wwKCRDCs_xy2_2, wwKCRDCs_yx2_2); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy2_2[indx0], wwKCRDCs_yx2_2[indx0]); abline(a=0, b=1, col="red")
plot(wwKCRDCs_xy2_2[indx1], wwKCRDCs_yx2_2[indx1]); abline(a=0, b=1, col="red")
sum(wwKCRDCs_xy2_2<wwKCRDCs_yx2_2)/length(wwKCRDCs_xy2_2)
sum(wwKCRDCs_xy2_2[indx0]<wwKCRDCs_yx2_2[indx0])/length(indx0)
sum(wwKCRDCs_xy2_2[indx1]<wwKCRDCs_yx2_2[indx1])/length(indx1)
sum(wwKCRDCs_xy2_2-wwKCRDCs_yx2_2)
sum(wwKCRDCs_xy2_2[indx0]-wwKCRDCs_yx2_2[indx0])
sum(wwKCRDCs_xy2_2[indx1]-wwKCRDCs_yx2_2[indx1])








df_xy$dif_wwKCRDCs <- wwKCRDCs_yx-wwKCRDCs_xy
p <- ggplot()
p <- p + geom_jitter(aes(x=x, y=y, colour=dif_wwKCRDCs), data=df_xy)
p <- p + facet_wrap(node~.)
p
p <- ggplot()
p <- p + geom_jitter(aes(x=y, y=x, colour=sigmax), data=df_yx)
p <- p + facet_wrap(node~.)
p

# block 65 SIM,   x->y, x<-y
# block 15 SIM max 0.0740, 0.0255  , median 0.00541226, 0.009147 x->y, x<-y
# block 81 CMPLX max 0.032244, 0.0154601   , median 0.00576, 0.003942976  x->y, x<-y
# block 7 SIM: max 0.004086159,0.003112318 ;q90 ,0.001652293 ;median  0.001021424, 0.001098591 
# block 8 SIM: max 0.05628, 0.0784189 ;q90 0.006641959 ,0.006904969 ;median  0.004220819, 0.005449294 
# block 9 SIM: max 0.05435728, 0.0171446  ;q90 0.01015099 ,0.01193194 ;median  0.007011652,  0.01002084 
# block 11 SIM: max 0.003769709,0.009168329   ;q90  0.002835336,0.004113942 ;median 0.001850111 ,0.002040604
# block 14 SIM: max 0.001603582,0.001095048   ;q90  0.0002449282,0.0003937126 ;median 7.19475e-05,0.0002612308
# block 15 SIM: max 0.02135701, 0.009884263  ;q90 0.008993436 , 0.006836715 ;median 0.004157616, 0.005775074
# block 16 SIM: max 0.02685977,  0.00649411 ;q90 0.007898834  ,0.003705647  ;median 0.00601549, 0.002696601
# block 21 SIM: max 0.05071291, 0.02031206  ;q90 0.01721862 , 0.0111145  ;median 0.008398663, 0.01001411
# block 1 SIM: max 0.04217174, 0.01781195   ;q90 0.01368515 , 0.008502712  ;median 0.01012215 , 0.007101338
# block 33 SIM: max 0.01024512, 0.02634627    ;q90 0.004205429  ,0.01451006 ;median  0.00331792 , 0.008797866
# block 99 SIM: max 0.008345049, 0.009037073    ;q90  0.004153131 , 0.007045233 ;median  0.002829577 , 0.006392713
# block 50 SIM: max 0.03807169, 0.05308692    ;q90 0.01696671 , 0.0174476  ;median   0.0133127, 0.00725461


wVarCorrs_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_xy$x[indx0], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  mu <- weighted.mean(df_xy$corr[indx0],w=weights)
  res <- weighted.var(df_xy$corr[indx0],w=weights)
  
  return(res/mu)
})
wVarCorrs_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_yx$x[indx0], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  mu <- weighted.mean(df_yx$corr[indx0],w=weights)
  res <- weighted.var(df_yx$corr[indx0],w=weights)
  
  return(res/mu)
})

wVarCorrs_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_xy$x[indx1], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  mu <- weighted.mean(df_xy$corr[indx1],w=weights)
  res <- weighted.var(df_xy$corr[indx1],w=weights)
  
  return(res/mu)
})
wVarCorrs_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_yx)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_yx$x[indx1], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  mu <- weighted.mean(df_yx$corr[indx1],w=weights)
  res <- weighted.var(df_yx$corr[indx1],w=weights)
  
  return(res/mu)
})

wVarCorrs_xy <- c(wVarCorrs_xy0, wVarCorrs_xy1)
wVarCorrs_yx <- c(wVarCorrs_yx0, wVarCorrs_yx1)
plot(wVarCorrs_xy, wVarCorrs_yx, col=c("red","blue")[rep(c(1,2), rep(length(wVarCorrs_xy0),2))]); abline(a=0, b=1, col="red")

plot(wVarCorrs_xy, wVarCorrs_yx, xlim=range(wVarCorrs_xy, wVarCorrs_yx), ylim=range(wVarCorrs_xy, wVarCorrs_yx), col=c("red","blue")[rep(c(1,2), rep(length(wVarCorrs_xy0),2))]); abline(a=0, b=1, col="red"); 
plot(wVarCorrs_xy0, wVarCorrs_yx0); abline(a=0, b=1, col="red"); abline(a=0, b=1, col="red")
plot(wVarCorrs_xy1, wVarCorrs_yx1); abline(a=0, b=1, col="red"); abline(a=0, b=1, col="red")

plot(df_xy$x,df_xy$y, col=c("red","blue")[(wwKCRDCs_xy<wwKCRDCs_yx)*1+1])
plot(df_xy$x,df_xy$y, col=c("red","blue")[(wVarCorrs_xy<wVarCorrs_yx)*1+1])

hist(wVarCorrs_xy-wVarCorrs_yx)

sum(wVarCorrs_xy<wVarCorrs_yx)/length(wVarCorrs_xy)
sum(wVarCorrs_xy0<wVarCorrs_yx0)/length(wVarCorrs_xy0)
sum(wVarCorrs_xy1<wVarCorrs_yx1)/length(wVarCorrs_xy1)

sum(wVarCorrs_xy[indx1]<wVarCorrs_yx[indx0])/length(indx0)
plot(wVarCorrs_xy[indx1],wVarCorrs_yx[indx0])

mean(wVarCorrs_yx-wVarCorrs_xy)
mean(wVarCorrs_yx0-wVarCorrs_xy0)
mean(wVarCorrs_yx1-wVarCorrs_xy1)

plot(wwKCRDCs_yx-wwKCRDCs_xy, wVarCorrs_yx-wVarCorrs_xy); abline(h=0, v=0, col="red")

indxPos <- which(wVarCorrs_yx-wVarCorrs_xy >0)
indxNeg <- which(wVarCorrs_yx-wVarCorrs_xy <0)
sum((wwKCRDCs_yx>wwKCRDCs_xy)[indxPos])/length(indxPos)
sum((wwKCRDCs_yx<wwKCRDCs_xy)[indxNeg])/length(indxNeg)
sum((wwKCRDCs_yx-wwKCRDCs_xy)[indxPos])
sum((wwKCRDCs_yx-wwKCRDCs_xy)[indxNeg])


wL2_f_xy <- apply(df_xy, 1, function(row){
  # row <- df_xy[1,]; names(row) <- colnames(df_xy) 
  
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- X2[smpl_te2[,boot],,drop=F]
  xTr_b <- X2[smpl_tr2[,boot],,drop=F]
  yTe_b <- as.matrix(Y2[smpl_te2[,boot],])
  yTr_b <- as.matrix(Y2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+numPerBoot*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay)
  Kte <- kern_rbf(yTe_b, sigma=sigmay)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmay)
  
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  
  resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
  #l2rel <- sqrt(resids[testPt]/mean(diag(Kte)^2))
  wl2rel <- sqrt(resids[testPt]/weighted.mean(diag(Kte)^2,w=weights))
  
  
  return(wl2rel)
})
wL2_f_yx <- apply(df_yx, 1, function(row){
  # row <- as.numeric(dfBestLambda_yx[1,]); names(row) <- colnames(dfBestLambda_yx) 
  lam <- as.numeric(row["lambda"])
  sigmax <- as.numeric(row["sigmax"])
  sigmay <- as.numeric(row["sigmay"])
  boot <- as.numeric(row["boot"])
  testPt <- as.numeric(row["testPt"])
  xTe_b <- Y2[smpl_te2[,boot],,drop=F]
  xTr_b <- Y2[smpl_tr2[,boot],,drop=F]
  yTe_b <- as.matrix(X2[smpl_te2[,boot],])
  yTr_b <- as.matrix(X2[smpl_tr2[,boot],])
  nTr_b <- nrow(xTr_b)
  Itr_b <- diag(nTr_b)  
  Lte <- kern_rbf(xTe_b, sigma=sigmax)
  Ltr <- kern_rbf(xTr_b, sigma=sigmax)
  Lte_tr <- kern_rbf(xTe_b,xTr_b, sigma=sigmax)
  Blambda_tr <- solve(Ltr+numPerBoot*lam*Itr_b)
  Ktr <- kern_rbf(yTr_b, sigma=sigmay)
  Kte <- kern_rbf(yTe_b, sigma=sigmay)
  Kte_tr <- kern_rbf(yTe_b,yTr_b, sigma=sigmay)
  
  weights <- Lte[testPt,]
  weights <- weights/sum(weights)
  
  resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
  resids <- diag(Kte) + diag( Lte_tr%*%Blambda_tr%*%(-2*t(Kte_tr) + Ktr%*%t(Blambda_tr)%*%t(Lte_tr)))
  #l2rel <- sqrt(resids[testPt]/mean(diag(Kte)^2))
  wl2rel <- sqrt(resids[testPt]/weighted.mean(diag(Kte)^2,w=weights))
  
  
  return(wl2rel)
})

wwL2_f_xy0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_xy$x[indx0], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  wVarL2rels <- weighted.var(wL2_f_xy[indx0], w=weights)
  return(wVarL2rels)
})
wwL2_f_yx0 <- sapply(indx0, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_yx$x[indx0], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  wVarL2rels <- weighted.var(wL2_f_xy[indx0], w=weights)
  return(wVarL2rels)
})

wwL2_f_xy1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_xy[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- X2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_xy$x[indx1], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  wVarL2rels <- weighted.var(wL2_f_xy[indx1], w=weights)
  return(wVarL2rels)
})
wwL2_f_yx1 <- sapply(indx1, function(i){
  # i <- 47
  row <- as.numeric(df_yx[i,]); names(row) <- colnames(df_xy)
  lam <- row["lambda"]
  sigmax <- row["sigmax"]
  sigmay <- row["sigmay"]
  boot <- row["boot"]
  testPt <- row["testPt"]
  xPt <- matrix(row["x"],1,1)
  xTr <- Y2
  Lte_tr <- kern_rbf(xPt,xTr, sigma=sigmax)
  weights <- Lte_tr[,match(df_yx$x[indx1], xTr)]
  weights <- weights/sum(weights)
  # hist(weights)
  wVarL2rels <- weighted.var(wL2_f_xy[indx1], w=weights)
  return(wVarL2rels)
})

wwL2_f_xy <- c(wwL2_f_xy0, wwL2_f_xy1)
wwL2_f_yx <- c(wwL2_f_yx0, wwL2_f_yx1)

plot(wL2_f_xy, wL2_f_yx); abline(a=0, b=1, col="red")
plot(wL2_f_xy[indx0], wL2_f_yx[indx0]); abline(a=0, b=1, col="red")
plot(wL2_f_xy[indx1], wL2_f_yx[indx1]); abline(a=0, b=1, col="red")

sum(wL2_f_xy<wL2_f_yx)/length(wL2_f_xy)
sum(wL2_f_xy[indx0]<wL2_f_yx[indx0])/length(indx0)
sum(wL2_f_xy[indx1]<wL2_f_yx[indx1])/length(indx1)

plot(df_xy$x, df_xy$y, cex=(norml(wL2_f_yx)+0.6)*1, col=c("red","blue")[(wL2_f_xy<wL2_f_yx)*1+1])
plot(df_xy$x, df_xy$y, cex=(norml(wwKCRDCs_yx)+0.6)*1.3, col=c("red","blue")[(wwKCRDCs_xy<wwKCRDCs_yx)*1+1])

plot(wwKCRDCs_yx-wwKCRDCs_xy, wL2_f_yx-wL2_f_xy); abline(h=0, v=0, col="red")
plot(wVarCorrs_yx-wVarCorrs_xy, wL2_f_yx-wL2_f_xy); abline(h=0, v=0, col="red")

plot(wwL2_f_xy, wwL2_f_yx); abline(a=0, b=1, col="red")
plot(wwL2_f_xy0, wwL2_f_yx0); abline(a=0, b=1, col="red")
plot(wwL2_f_xy1, wwL2_f_yx1); abline(a=0, b=1, col="red")

sum(wwL2_f_xy<wwL2_f_yx)/length(wwL2_f_xy)
sum(wwL2_f_xy0<wwL2_f_yx0)/length(indx0)
sum(wwL2_f_xy1<wwL2_f_yx1)/length(indx1)

plot(df_xy$x, df_xy$y, cex=(norml(wwL2_f_yx)+0.6)*1, col=c("red","blue")[(wwL2_f_xy<wwL2_f_yx)*1+1])
plot(df_xy$x, df_xy$y, cex=(norml(wwKCRDCs_yx)+0.6)*1.3, col=c("red","blue")[(wwKCRDCs_xy<wwKCRDCs_yx)*1+1])

plot(wwKCRDCs_yx-wwKCRDCs_xy, wwL2_f_yx-wwL2_f_xy); abline(h=0, v=0, col="red")
plot(wVarCorrs_yx-wVarCorrs_xy, wwL2_f_yx-wwL2_f_xy); abline(h=0, v=0, col="red")




sum(wwKCRDCs_xy<wwKCRDCs_yx)/length(wwKCRDCs_xy)
sum((wwKCRDCs_xy <wwKCRDCs_yx)[indx])/length(indx)
cor((wwKCRDCs_yx-wwKCRDCs_xy)[unique(indx)], (wVarCorrs_yx-wVarCorrs_xy)[unique(indx)], method="spearman")
sum(wwKCRDCs_xy-wwKCRDCs_yx)
indxR <- which(wVarCorrs_yx>wVarCorrs_xy)
indxW <- which(wVarCorrs_yx<wVarCorrs_xy)
sum((wwKCRDCs_xy <wwKCRDCs_yx)[indxR])/length(indxR) - sum((wwKCRDCs_xy >wwKCRDCs_yx)[indxW])/length(indxW)

# prop pos, prop pos indx, corr indx, sum difs , prop_pos_varCorrPos, sum_pos_varCorrPos
# block 7 SIM: 0.5017822,0.6680672, -0.5737559
# block 8 SIM: 0.5932673, 0.5534591, -0.9352747 
# block 9 SIM: 0.8023762, 0.3504098,  0.01756754,  
# block 11 SIM: 0.6557426, 0.8795455, -0.9063106
# block 14 SIM: 0.790099,0.75,-0.5297447
# block 15 SIM: 
# block 16 SIM: 
# block 21 SIM: 
# block 1 SIM: 0.529703, 0.07949791, -0.7212306
# block 33 SIM: 0.7485149, 0.1157895 , 0.3137019
# block 99 SIM: 0.9432673, 0.894958 , -0.06304228, -33.54947, 0.5149985
# block 50 SIM: 0.2047525, 0.1181619 , -0.1693706, 40.80764,  -0.5743416

summary(wwKCRDCs_yx-wwKCRDCs_xy)
summary(wVarCorrs_yx-wVarCorrs_xy)
table(wwKCRDCs_yx-wwKCRDCs_xy>0)
table(wVarCorrs_yx-wVarCorrs_xy>0)
table(wwKCRDCs_yx-wwKCRDCs_xy>0, wVarCorrs_yx-wVarCorrs_xy>0)

plot(dfBestLambda$x, log(dfBestLambda$sigmay,10))
plot(dfBestLambda$x, log(dfBestLambda$sigmax,10))
plot(dfBestLambda$x, log(dfBestLambda$lambda,10))
plot(X, Y)
set.seed(1234)
smpl <- sample(nrow(dfBestLambda),200)

# sigmay
plot(dfBestLambda$x[smpl], log(dfBestLambda$sigmay[smpl],10))
trainDataAux <- constructData(x=as.matrix(dfBestLambda$x[smpl]), y=log(dfBestLambda$sigmay[smpl],10))
krr1_aux <- setParams(krr1, trainData=trainDataAux, mc_cores=5)
krr1_aux <- krr1_aux$learn(krr1_aux)
pred_aux <- krr1_aux$predict(krr1_aux, constructData(x=as.matrix(dfBestLambda$x), y=log(dfBestLambda$sigmay,10)))
plot(dfBestLambda$x, log(dfBestLambda$sigmay,10))
xx <- pred_aux$x_class
yy <- pred_aux$gyh_class
o <- order(xx)
lines(xx[o], yy[o], col="red",lwd=2)

# sigmax
trainDataAux <- constructData(x=as.matrix(dfBestLambda$x[smpl]), y=log(dfBestLambda$sigmax[smpl],10))
krr1_aux <- setParams(krr1, trainData=trainDataAux, mc_cores=5)
krr1_aux <- krr1_aux$learn(krr1_aux)
pred_aux <- krr1_aux$predict(krr1_aux, constructData(x=as.matrix(dfBestLambda$x), y=log(dfBestLambda$sigmax,10)))
plot(dfBestLambda$x, log(dfBestLambda$sigmax,10))
xx <- pred_aux$x_class
yy <- pred_aux$gyh_class
o <- order(xx)
lines(xx[o], yy[o], col="red",lwd=2)

# lambda
trainDataAux <- constructData(x=as.matrix(dfBestLambda$x[smpl]), y=log(dfBestLambda$lambda[smpl],10))
krr1_aux <- setParams(krr1, trainData=trainDataAux, mc_cores=5)
krr1_aux <- krr1_aux$learn(krr1_aux)
pred_aux <- krr1_aux$predict(krr1_aux, constructData(x=as.matrix(dfBestLambda$x), y=log(dfBestLambda$lambda,10)))
plot(dfBestLambda$x, log(dfBestLambda$lambda,10))
xx <- pred_aux$x_class
yy <- pred_aux$gyh_class
o <- order(xx)
lines(xx[o], yy[o], col="red",lwd=2)


