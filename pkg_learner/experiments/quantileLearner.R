
remove(list=ls())


server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

reposOptimus <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
reposErc <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/pkg_causaLearner/", sep="")
if(server=="optimus.uv.es") repos <- reposOptimus
if(server=="erc.uv.es") repos <- reposErc

dir(repos)



source(paste(repos, "pkg_causaLearner/genData/func_getData_v2.R", sep=""))
#library(dHSIC)
source(paste(repos, "pkg_dHSIC/dHSIC.R", sep=""))
source(paste(repos, "pkg_causaLearner/dataTreatments/func_dataTreatments.R", sep=""))
source(paste(repos, "pkg_learner/func_learners_v3.R", sep=""))
library(reshape)
library(ggplot2)
library(opera) # loss(x,y, loss.type=list(name="pinball", tau=0.95))
library(kernlab)

# Quantile regression

# chosen function
edgL <- vector("list", length=2)
names(edgL) <- c("x","y")
edgL[["x"]] <- list(edges=c("y"))
XtoY <- graphNEL(nodes=c("x","y"), edgeL=edgL, edgemode="directed")
plot(XtoY)
# nodesUnif <-  list(x=list(dist="runif", pars=list(min=-pi, max=pi), a=1, b=1),
#                   y=list(dist="runif", pars=list(min=-0.5, max=0.5), a=1, b=1))
# funcs <- list(fx=function(n) n, fy=function(x, n) sin(x)+n)
# sem <- list(dag=XtoY, funcs=funcs, simPars=list(n=100, nodes=nodesUnif))
# nodesUnif <-  list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1),
#                    y=list(dist="runif", pars=list(min=-1, max=1), a=1, b=1))
# funcs <- list(fx=function(n) n, fy=function(x, n) x*n)
# sem <- list(dag=XtoY, funcs=funcs, simPars=list(n=100, nodes=nodesUnif))
# 
# dataTestList <- simSEMs(q=100, sem)
# plotPairsList(dataList=dataTestList)


folder <- "datos/pairs/"
fileName <- "pairmeta.txt"
fileFolder <- paste(repos, folder, fileName, sep="")
meta <- read.csv(fileFolder, sep="", header=F)
# From README file
colnames(meta) <- c("pairNumber", "firstCauseCol", "lastCauseCol", "firstEffectCol", "lastEffectCol", "dataSetWeight")
head(meta)
table(meta$firstCauseCol)
table(meta$firstEffectCol)
table(meta$firstCauseCol, meta$firstEffectCol)
# which one is cause and which one effect?
indxUni <- which(with(meta, (lastCauseCol-firstCauseCol)==0 & (lastEffectCol-firstEffectCol)==0))
length(indxUni)
pairs <- meta$pairNumber[indxUni]
summary(meta[indxUni,])
# how many data points do they have?
(nPerPair <- numDataPairs(pairs, folder))
# get range for x and y for each pair
summaryDataPairs(meta, pairs, folder)
# create data test
dataTestList <- createTCEPList(pairs, meta, folder)
i <- 1; dataTestList$xs[[i]]; dataTestList$dags[[i]]; dataTestList$noiss[[i]]; dataTestList$names[i]
# read in files
dataTestList$xs <- lapply(dataTestList$xs, function(el){
  # el <- dataTestList$xs[[1]]
  #print(el)
  res <- read.csv(el, sep="", header=F)
  res <- as.matrix(res)
  return(res)
})
# delete non-univariate pairs
indx <- which(sapply(dataTestList$xs, ncol)==2)
dataTestList$xs <- dataTestList$xs[indx]
dataTestList$dags <- dataTestList$dags[indx]
dataTestList$noiss <- dataTestList$noiss[indx]
dataTestList$names <- dataTestList$names[indx]
#put in right order
dataTestList$xs <- mapply(FUN=function(xs, indx) xs[,indx], xs=dataTestList$xs, indx=dataTestList$noiss)

# add colnames
dataTestList$xs <- lapply(dataTestList$xs, function(el){
  res <- el
  colnames(res) <- c("x","y")
  return(res)
})
#take only 100 max per pair
dataTestList$xs <- lapply(dataTestList$xs, function(el) el[sample(1:nrow(el), size=min(100, nrow(el))),])
plotPairsList(dataList=dataTestList)


i <- 5

X <- dataTestList$xs[[i]]
X <- apply(X, 2, stdrize)
#print("head(X)")
#print(head(X))
#print(apply(X, 2, mean))
#print(apply(X, 2, sd))

xs <- X[,1]
ys <- X[,2]

par(mfrow=c(2,2))
plot(xs, ys)
plot(ys, xs)
plot(xs, ys)
plot(ys, xs)
par(mfrow=c(1,1))

n <- 100
pm <- proc.time()
qrm <- kqr(xs, ys, tau = 0.5 , kernel = "rbfdot", kpar= list(sigma=0.1), C=100)
proc.time() - pm # 0.33 secs

ytest <- predict(qrm, xs)
plot(xs, ys)
ord <- order(xs)
lines(xs[ord], ytest[ord], col="red")


estimateCDF <- function(qs, ps){
  
  # remove qs outside support
  
  indxOutsideSupport <- which(findInterval(qs[2:(length(qs)-1)], qs[c(1,length(qs))]) %in% c(0,2))
  if(length(indxOutsideSupport)>0){
    qs <- c(qs[1], qs[2:(length(qs)-1)][-indxOutsideSupport], qs[length(qs)])
    ps <- c(ps[1], ps[2:(length(ps)-1)][-indxOutsideSupport], ps[length(ps)])
  }
  # remove NAs from these pairs
  indxNA <- which(!is.na(qs))
  qs <- qs[indxNA]
  ps <- ps[indxNA]
  
  # the cdf is non decreasing so we  delete such points - if we dont take appropriate measures
  # we could end up with only one point!! 
  # possible measures:
  # 1. use very smooth quantile corresponding to p=0 and p=1 such that they will almost always be q0<q1 and you wd at least have 2 pts
  # 2. correct q0 as q0 = min(qs) and q1 as q1 = max(qs)  before going on - this corresponds to a uniform in range of qs
  # 3. for training first and last quantiles use whole sample (we wd still need to do 1 or 2 but it might help)
  # in the spirit of letting the noise and causal direction to punish/reward the parameters i prefer option 2
  
  #qs[1] <- min(qs)
  #qs[length(qs)] <- max(qs)
  
  
  
  difqs <- diff(qs)
  negDif <- which(difqs<0)
  while(any(negDif)){
    qs[negDif+1] <- NA
    indxNA <- which(!is.na(qs))
    qs <- qs[indxNA]
    ps <- ps[indxNA]
    difqs <- diff(qs)
    negDif <- which(difqs<0)
  }
  #qs <- sort(qs)
  # itnerpolate to estimate cdf
  splf <- splinefun(qs, ps, method="monoH.FC") 
  # plot(seq(min(qs), max(qs), length.out=100), splf(seq(min(qs), max(qs), length.out=100)), type="l"); lines(qs, ps, type="p")
  return(list(qs=qs, ps=ps, CDF=splf))
}

estimateCDF2 <- function(qs, ps){
  
  # remove qs outside support
  
  
  # remove NAs from these pairs
  indxNA <- which(!is.na(qs))
  qs <- qs[indxNA]
  ps <- ps[indxNA]
  
  # the cdf is non decreasing so we  delete such points - if we dont take appropriate measures
  # we could end up with only one point!! 
  # possible measures:
  # 1. use very smooth quantile corresponding to p=0 and p=1 such that they will almost always be q0<q1 and you wd at least have 2 pts
  # 2. correct q0 as q0 = min(qs) and q1 as q1 = max(qs)  before going on - this corresponds to a uniform in range of qs
  # 3. for training first and last quantiles use whole sample (we wd still need to do 1 or 2 but it might help)
  # in the spirit of letting the noise and causal direction to punish/reward the parameters i prefer option 2
  qsAux <- qs
  
  qs[1] <- min(qsAux)
  qs[length(qs)] <- max(qsAux)
  
  
  
  difqs <- diff(qs)
  negDif <- which(difqs<0)
  while(any(negDif)){
    qs[negDif+1] <- NA
    indxNA <- which(!is.na(qs))
    qs <- qs[indxNA]
    ps <- ps[indxNA]
    difqs <- diff(qs)
    negDif <- which(difqs<0)
  }
  #qs <- sort(qs)
  # itnerpolate to estimate cdf
  splf <- splinefun(qs, ps, method="monoH.FC") 
  # plot(seq(min(qs), max(qs), length.out=100), splf(seq(min(qs), max(qs), length.out=100)), type="l"); lines(qs, ps, type="p")
  return(list(qs=qs, ps=ps, CDF=splf))
}


predQ <- function(xTr, yTr, xTe, yTe, sig, lam, taus, trainTest){
  
  if(trainTest=="train"){
    xTot <- xTr
    yTot <- yTr
  } else{
    xTot <- rbind(xTe, xTr)
    yTot <- rbind(yTe,yTr)
  }
  
  Cmax <- 1e5
  qrm <- NULL
  
  preds <- mcmapply(FUN=function(t){
    # t <- 0
    #print(paste("t: ", t))
    if(t==0){
      Caux <- Cmax
      while(class(qrm)=="try-error"|is.null(qrm)){
        qrm <- try(kqr(xTot, yTot, tau = 1/(10000*length(xTr)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux))
        #qrm <- try(kqr(xTr, yTr, tau = 3/(length(xTr)), kernel = "rbfdot", kpar= list(sigma=sig), C=lam))
        Caux <- Caux*0.1
      }
    } else if(t==1){
      Caux <- Cmax
      while(class(qrm)=="try-error"|is.null(qrm)){
        qrm <- try(kqr(xTot, yTot, tau = 1-1/(10000*length(xTr)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux))
        #qrm <- try(kqr(xTr, yTr, tau = 1-3/(length(xTr)), kernel = "rbfdot", kpar= list(sigma=sig), C=lam))
        Caux <- Caux*0.1
      }
    } else{ 
      qrm <- try(kqr(xTr, yTr, tau = t, kernel = "rbfdot", kpar= list(sigma=sig), C=lam))
    }
    if(class(qrm)=="try-error"){
      pred <- matrix(NA, length(xTe),1)
    } else{
      pred <- predict(qrm, xTe)
    }
    return(pred)
  }, t=c(0,   taus, 1), mc.cores=1) 
  return(preds)
}

predQ2 <- function(xTr, yTr, xTe, yTe, sig, lam, taus, trainTest){
  
 
 
  preds <- mcmapply(FUN=function(t){
    # t <- 0
    #print(paste("t: ", t))
    
    qrm <- try(kqr(xTr, yTr, tau = t, kernel = "rbfdot", kpar= list(sigma=sig), C=lam))
    
    if(class(qrm)=="try-error"){
      pred <- matrix(NA, length(xTe),1)
    } else{
      pred <- predict(qrm, xTe)
    }
    return(pred)
  }, t=taus, mc.cores=1) 
  return(preds)
}


lossesQ <- function(preds, xs, ys, sig, taus){
  
  #plot(xs, ys, ylim=range(preds))
  #ord <- order(xs)
  #for(i in 1:(ncol(preds))) lines(xs[ord], preds[ord, i], col=i)
  #abline(v=xs[8], col="red")
  
  res <- sapply(1:length(ys), function(i){
    # i <- 8
    #print(paste("i: ", i))
    
    # CDF estimation
    # get tau quantiles and cdf pairs (q,p) from our quantile reg model
    qs_i <- preds[i, ]
    ps_i <- c(0,  taus,  1) #c(0,taus,1)
    #plot(qs_i, ps_i)
    CDFi <- estimateCDF(qs=qs_i, ps=ps_i)
    qs_i <- CDFi$qs
    ps_i <- CDFi$ps
    CDFi <- CDFi$CDF
    
    #plot(qs_i, CDFi(qs_i))
    #arg <- seq(min(qs_i), max(qs_i), length.out=100); plot(arg, CDFi(arg), type="l"); lines(qs_i, CDFi(qs_i), type="p", col="red")
    
    # we are interested in evaluating the cdf and pdf py(.|x_i) evaluating at y_i so we get set pt <- y_i
    y_i <- ys[i]
    # we evaluate py(y_i|x_i) and Py(y_i|x_i) (pdf and cdf)
    lik <- CDFi(y_i, deriv=1)
    cdf <- CDFi(y_i, deriv=0)
    
    # for testing out of training sample y_i could fall out of support of py(y_i|x_i) in this case 
    # if there are some y_j for any j that fall within qi0-qi1 we simply take lik=min_j p(y_j|x_i) j such that qi0<=y_j<=qi1
    # else we look for nearby x_j (using kernel) and set lik= p(y_i|x_j) if qj0<=y_i<=qj1 and lik=min_k p(y_k|x_j) if
    # !(qj0<=y_i<=qj1) but there exists at least 1 yk such that qj0<=y_k<=qj1
    if(lik<=0 | findInterval(y_i, range(qs_i), rightmost.closed=T) %in% c(0,2)){
      print("entering no mans land  !!!")
      lik <- NA
      aux0 <- CDFi(ys, deriv=1)
      indxValid <- which(aux0>0 & findInterval(ys, range(qs_i), rightmost.closed=T)==1)
      if(length(indxValid)>0){
        indx <- which.min(aux0[indxValid])
        #lik <- aux0[indxValid][indx]*0.5
        lik <- NA
      } else{
        print("entering no mans land 2 !!!")
        lik <- NA
        # foundOne <- FALSE
        # if(is.null(dim(xs))) Xs <- matrix(xs, length(xs), 1) else Xs <- xs
        # K <- kernelMatrix("kern_rbf", Xs[i,,drop=FALSE], Xs, pars=list(sigma=sig))
        # ord <- order(as.numeric(K), decreasing=TRUE)
        # ord <- ord[-1]
        # j <- ord[1]
        # while(!foundOne & length(ord)>0){
        #   #print(paste("j: ", j))
        #   qs_j <- preds[j,]
        #   ps_j <- c(0, taus, 1)
        #   CDFj <- estimateCDF(qs_j, ps_j)
        #   qs_j <- CDFj$qs
        #   CDFj <- CDFj$CDF
        #   lik <- CDFj(y_i, deriv=1)
        #   if(lik>0 & findInterval(y_i, range(qs_j))==1){
        #     foundOne <- TRUE
        #   } else {
        #     aux0 <- CDFj(ys, deriv=1)
        #     indxValid <- which(aux0>0 & findInterval(ys, range(qs_j))==1)
        #     if(length(indxValid)>0){
        #       foundOne <- TRUE
        #       indx <- which.min(aux0[indxValid])
        #       lik <- aux0[indxValid][indx]*0.5
        #     }
        #     ord <- ord[-1]
        #     j <- ord[1]
        #   }
        # }
        # 
        # if(!foundOne){
        #   lik <- NA
        # }
        
      }
    }
    if(cdf<0) cdf <- 0
    if(cdf>1) cdf <- 1
    
    qs <- findInterval(y_i, qs_i, rightmost.closed=T)
    return(c(lik, cdf, qs))
  }, simplify="array")
  res <- t(res)
  colnames(res) <- c("lik","cdf","q")
  indxNA <- which(is.na(res[,"lik"]))
  res[indxNA,"lik"] <- min(res[,"lik"], na.rm=T)*0.1 #1e-6
  return(res)
  
}


lossesQ2 <- function(preds, xs, ys, sig, taus){
  
  #plot(xs, ys, ylim=range(preds))
  #ord <- order(xs)
  #for(i in 1:(ncol(preds))) lines(xs[ord], preds[ord, i], col=i)
  #abline(v=xs[8], col="red")
  
  res <- sapply(1:length(ys), function(i){
    # i <- 9
    #print(paste("i: ", i))
    
    # CDF estimation
    # get tau quantiles and cdf pairs (q,p) from our quantile reg model
    qs_i <- preds[i, ]
    ps_i <- taus
    #plot(qs_i, ps_i)
    
    if(sum(!is.na(qs_i))<2){
      res <- c(0,0,1)
      return(res)
    }
    
    CDFi <- estimateCDF2(qs=qs_i, ps=ps_i)
    qs_i <- CDFi$qs
    ps_i <- CDFi$ps
    CDFi <- CDFi$CDF
    
    #plot(qs_i, CDFi(qs_i))
    #arg <- seq(min(qs_i), max(qs_i), length.out=100); plot(arg, CDFi(arg), type="l"); lines(qs_i, CDFi(qs_i), type="p", col="red")
    
    # we are interested in evaluating the cdf and pdf py(.|x_i) evaluating at y_i so we get set pt <- y_i
    y_i <- ys[i]
    # we evaluate py(y_i|x_i) and Py(y_i|x_i) (pdf and cdf)
    lik <- CDFi(y_i, deriv=1)
    cdf <- CDFi(y_i, deriv=0)
    
    # for testing out of training sample y_i could fall out of support of py(y_i|x_i) in this case 
    # if there are some y_j for any j that fall within qi0-qi1 we simply take lik=min_j p(y_j|x_i) j such that qi0<=y_j<=qi1
    # else we look for nearby x_j (using kernel) and set lik= p(y_i|x_j) if qj0<=y_i<=qj1 and lik=min_k p(y_k|x_j) if
    # !(qj0<=y_i<=qj1) but there exists at least 1 yk such that qj0<=y_k<=qj1
    if(lik<=0 | findInterval(y_i, range(qs_i), rightmost.closed=T) %in% c(0,2)){
      #print("entering no mans land  !!!")
      lik <- NA
      aux0 <- CDFi(ys, deriv=1)
      indxValid <- which(aux0>0 & findInterval(ys, range(qs_i), rightmost.closed=T)==1)
      if(length(indxValid)>0){
        indx <- which.min(aux0[indxValid])
        #lik <- aux0[indxValid][indx]*0.5
        lik <- NA
      } else{
        #print("entering no mans land 2 !!!")
        lik <- NA
        # foundOne <- FALSE
        # if(is.null(dim(xs))) Xs <- matrix(xs, length(xs), 1) else Xs <- xs
        # K <- kernelMatrix("kern_rbf", Xs[i,,drop=FALSE], Xs, pars=list(sigma=sig))
        # ord <- order(as.numeric(K), decreasing=TRUE)
        # ord <- ord[-1]
        # j <- ord[1]
        # while(!foundOne & length(ord)>0){
        #   #print(paste("j: ", j))
        #   qs_j <- preds[j,]
        #   ps_j <- c(0, taus, 1)
        #   CDFj <- estimateCDF(qs_j, ps_j)
        #   qs_j <- CDFj$qs
        #   CDFj <- CDFj$CDF
        #   lik <- CDFj(y_i, deriv=1)
        #   if(lik>0 & findInterval(y_i, range(qs_j))==1){
        #     foundOne <- TRUE
        #   } else {
        #     aux0 <- CDFj(ys, deriv=1)
        #     indxValid <- which(aux0>0 & findInterval(ys, range(qs_j))==1)
        #     if(length(indxValid)>0){
        #       foundOne <- TRUE
        #       indx <- which.min(aux0[indxValid])
        #       lik <- aux0[indxValid][indx]*0.5
        #     }
        #     ord <- ord[-1]
        #     j <- ord[1]
        #   }
        # }
        # 
        # if(!foundOne){
        #   lik <- NA
        # }
        
      }
    }
    if(cdf<0) cdf <- 0
    if(cdf>1) cdf <- 1
    
    qs <- findInterval(y_i, qs_i, rightmost.closed=T)
    return(c(lik, cdf, qs))
  }, simplify="array")
  res <- t(res)
  colnames(res) <- c("lik","cdf","q")
  indxNA <- which(is.na(res[,"lik"]))
  res[indxNA,"lik"] <- min(res[,"lik"], na.rm=T)*0.1 #1e-6
  return(res)
  
}



CV.parallel.lik <- function(x, y, numFolds, taus, lambdas, sigmas,  fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  set.seed(12)
  disord <- sample(length(x))
  trainData <- cbind(x=x[disord],y=y[disord])
  
  params <- constructParams(lambda=lambdas, sigma=sigmas)
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), lambdas, sigmas)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    # p <- params[[1]] 
    print(paste("lambda: ", p$lambda, ", sigma: ", p$sigma , sep=""))
    
   
    
    #  quantiles
    losses <- mcmapply(FUN=function(f){
      # f <- 2; p <- params[[1]]
      print(paste("lambda: ", p$lambda, ", sigma: ", p$sigma , " fold: ", f, sep=""))
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      # i <- 1; t <- taus[i]; l <- optLambdas[i]; s <- optSigmas[i]
      # print(paste("tau: ", t, ", lambda: ", l, ", sigma: ", s))
      predsTe <- predQ(xTr=curTrain[,"x", drop=FALSE], yTr=curTrain[,"y", drop=FALSE], xTe=curTest[,"x", drop=FALSE], yTe=curTest[,"y", drop=FALSE], sig=p$sigma, lam=p$lambda, taus, trainTest="test")
      predsTr <- predQ(xTr=curTrain[,"x", drop=FALSE], yTr=curTrain[,"y", drop=FALSE], xTe=curTrain[,"x", drop=FALSE], yTe=curTrain[,"y", drop=FALSE], sig=p$sigma, lam=p$lambda, taus, trainTest="train")
      
      liksTe <- lossesQ(preds=predsTe, xs=curTest[,"x"], ys=curTest[,"y"], sig=p$sigma, taus)
      liksTr <- lossesQ(predsTr, curTrain[,"x"], curTrain[,"y"], p$sigma, taus)  
      
      # plots for debugging
      {
      # plot(curTrain[,"x"], curTrain[,"y"], ylim=range(predsTr))
      # ord <- order(curTrain[,"x"])
      # for(i in 1:(length(taus)+2)) lines(curTrain[ord,"x"], predsTr[ord, i], col=i)
      # 
      
      #plot(curTest[,"x"], curTest[,"y"], ylim=range(predsTe))
      #ord <- order(curTest[,"x"])
      #for(i in 1:(length(taus)+2)) lines(curTest[ord,"x"], predsTe[ord, i], col=i)
      
      
      # 
      # plot(curTest[,"x"], curTest[,"y"], ylim=range(predsTe))
      # ord <- order(curTest[,"x"])
      # for(i in 1:(length(taus)+2)) lines(curTest[ord,"x"], predsTe[ord, i], col=i)
      # indx <- which(lossesTe[ord,"lik"]>100)
      # lines(curTest[ord[indx],"x"], curTest[ord[indx],"y"], cex=2, col="red", type="p")
      # lossesTe[ord,"lik"]
      # log(lossesTe[ord,"lik"])
      # indx
      # 
      # i <- 28
      # (x <- curTest[ord[i],"y"])
      # qs <- predsTe[ord[i], ]
      # ps <- c(0,taus,1)
      # difqs <- diff(qs)
      # negDif <- which(difqs<0)
      # while(any(negDif)){
      #   qs[negDif+1] <- NA
      #   indxNA <- which(!is.na(qs))
      #   qs <- qs[indxNA]
      #   ps <- ps[indxNA]
      #   difqs <- diff(qs)
      #   negDif <- which(difqs<0)
      # }
      # 
      # splf <- splinefun(qs, ps, method="monoH.FC")
      # plot(seq(min(qs), max(qs), length.out=100), splf(seq(min(qs), max(qs), length.out=100)), type="l")
      # lines(qs, ps, type="p")
      # lines(x, splf(x, deriv=0), type="p")
      # lik <- splf(x, deriv=1)
      }
      
      
      
      liksTe <- cbind(liksTe, predsTe)
      liksTr <- cbind(liksTr, predsTr)
      losses <- list(liksTe, liksTr)
      
      
      return(losses) 
    }, f=1:numFolds, mc.cores=5, SIMPLIFY=FALSE)
    
    lossesTe <- lapply(losses, function(el) el[[1]])
    lossesTr <- lapply(losses, function(el) el[[2]])
    lossesTe <- do.call(rbind, lossesTe)
    lossesTr <- do.call(rbind, lossesTr)
    
    predsTe <- lossesTe[,4:ncol(lossesTr)]
    predsTr <- lossesTr[,4:ncol(lossesTr)]
    
    lossesTe <- lossesTe[,1:3]
    lossesTr <- lossesTr[,1:3]
    
    pinballsTe <- sapply(1:(ncol(predsTe)), function(i) {
      # i <- 5
      loss(predsTe[,i], trainData[,"y"], loss.type=list(name="pinball", tau=c(0, taus, 1)[i]))
    }, simplify="array")
    
    
    
    
    lik <- -sum(log(lossesTe[,"lik"]))/nrow(lossesTe)
    cdfDhsic <- dhsic.test(trainData[,"x"], lossesTe[,"cdf"])$statistic
    qDhsic <- dhsic.test(trainData[,"x"], lossesTe[,"q"])$statistic
    pinballs <- apply(pinballsTe, 2, function(col) col/sd(col, na.rm=T))
    pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
    lossesTe <- c(lik, cdfDhsic, qDhsic, pinball)
    
    
    indxTr <- unlist(sapply(1:numFolds, function(f){
      indx <- seq((f-1)*size + 1, min(f*size,n))
      indx <- setdiff(1:n, indx)
      return(indx)
    }, simplify=FALSE))
    
    
    pinballsTr <- sapply(1:(ncol(predsTr)), function(i) {
      # i <- 1
      loss(predsTr[,i], trainData[indxTr,"y"], loss.type=list(name="pinball", tau=c(0,  taus, 1)[i]))
    }, simplify="array")
    
    lik <- -sum(log(lossesTr[,"lik"]))/nrow(lossesTr)
    cdfDhsic <- dhsic.test(trainData[indxTr,"x"], lossesTr[,"cdf"])$statistic
    qDhsic <- dhsic.test(trainData[indxTr,"x"], lossesTr[,"q"])$statistic
    pinballs <- apply(pinballsTr, 2, function(col) col/sd(col, na.rm=T))
    pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
    lossesTr <- c(lik, cdfDhsic, qDhsic, pinball)
    
    losses <- cbind(lossesTe, lossesTr)
    
    return(losses)
  }, p=params, mc.cores=2, SIMPLIFY="array")
  
  
  dimnames(losses) <- list(loss=c("lik", "cdfDhsic", "qDhsic","pinball"), testTrain=c("test","train"),params=names(params))
  
  #print("exits CV parallel")
  return(losses)
}


CV.parallel.lik2 <- function(x, y, numFolds, taus, lambdas, sigmas,  fac=1, verbose=TRUE) {
  
  #print("enters CV parallel")
  
  set.seed(12)
  disord <- sample(length(x))
  trainData <- cbind(x=x[disord],y=y[disord])
  
  params <- constructParams(lambda=lambdas, sigma=sigmas)
  
  nParams <- length(params)
  dimnames <- list(as.character(1:numFolds), lambdas, sigmas)
  
  n <- nrow(trainData)
  size <- ceiling(n / numFolds)
  
  losses <- mcmapply(FUN=function(p){
    # p <- params[[5]] 
    print(paste("lambda: ", p$lambda, ", sigma: ", p$sigma , sep=""))
    
    
    
    #  quantiles
    losses <- mcmapply(FUN=function(f){
      # f <- 3; p <- params[[5]]
      print(paste("lambda: ", p$lambda, ", sigma: ", p$sigma , " fold: ", f, sep=""))
      validationIndex <- seq((f-1)*size + 1, min(f*size,n))
      curTrain <- trainData[setdiff(1:n, validationIndex),]
      curTest <- trainData[validationIndex,]
      # either mean squared error or mean classification error
      
      # i <- 1; t <- taus[i]; l <- optLambdas[i]; s <- optSigmas[i]
      # print(paste("tau: ", t, ", lambda: ", l, ", sigma: ", s))
      predsTe <- predQ2(xTr=curTrain[,"x", drop=FALSE], yTr=curTrain[,"y", drop=FALSE], xTe=curTest[,"x", drop=FALSE], yTe=curTest[,"y", drop=FALSE], sig=p$sigma, lam=p$lambda, taus, trainTest="test")
      predsTr <- predQ2(xTr=curTrain[,"x", drop=FALSE], yTr=curTrain[,"y", drop=FALSE], xTe=curTrain[,"x", drop=FALSE], yTe=curTrain[,"y", drop=FALSE], sig=p$sigma, lam=p$lambda, taus, trainTest="train")
      
      liksTe <- lossesQ2(preds=predsTe, xs=curTest[,"x"], ys=curTest[,"y"], sig=p$sigma, taus)
      liksTr <- lossesQ2(predsTr, curTrain[,"x"], curTrain[,"y"], p$sigma, taus)  
      
      # plots for debugging
      {
        # plot(curTrain[,"x"], curTrain[,"y"], ylim=range(predsTr))
        # ord <- order(curTrain[,"x"])
        # for(i in 1:(length(taus)+2)) lines(curTrain[ord,"x"], predsTr[ord, i], col=i)
        # 
        
        #plot(curTest[,"x"], curTest[,"y"], ylim=range(predsTe))
        #ord <- order(curTest[,"x"])
        #for(i in 1:(length(taus)+2)) lines(curTest[ord,"x"], predsTe[ord, i], col=i)
        
        
        # 
        # plot(curTest[,"x"], curTest[,"y"], ylim=range(predsTe))
        # ord <- order(curTest[,"x"])
        # for(i in 1:(length(taus)+2)) lines(curTest[ord,"x"], predsTe[ord, i], col=i)
        # indx <- which(lossesTe[ord,"lik"]>100)
        # lines(curTest[ord[indx],"x"], curTest[ord[indx],"y"], cex=2, col="red", type="p")
        # lossesTe[ord,"lik"]
        # log(lossesTe[ord,"lik"])
        # indx
        # 
        # i <- 28
        # (x <- curTest[ord[i],"y"])
        # qs <- predsTe[ord[i], ]
        # ps <- c(0,taus,1)
        # difqs <- diff(qs)
        # negDif <- which(difqs<0)
        # while(any(negDif)){
        #   qs[negDif+1] <- NA
        #   indxNA <- which(!is.na(qs))
        #   qs <- qs[indxNA]
        #   ps <- ps[indxNA]
        #   difqs <- diff(qs)
        #   negDif <- which(difqs<0)
        # }
        # 
        # splf <- splinefun(qs, ps, method="monoH.FC")
        # plot(seq(min(qs), max(qs), length.out=100), splf(seq(min(qs), max(qs), length.out=100)), type="l")
        # lines(qs, ps, type="p")
        # lines(x, splf(x, deriv=0), type="p")
        # lik <- splf(x, deriv=1)
      }
      
      
      
      liksTe <- cbind(liksTe, predsTe)
      liksTr <- cbind(liksTr, predsTr)
      losses <- list(liksTe, liksTr)
      
      
      return(losses) 
    }, f=1:numFolds, mc.cores=5, SIMPLIFY=FALSE)
    
    lossesTe <- lapply(losses, function(el) el[[1]])
    lossesTr <- lapply(losses, function(el) el[[2]])
    lossesTe <- do.call(rbind, lossesTe)
    lossesTr <- do.call(rbind, lossesTr)
    
    predsTe <- lossesTe[,4:ncol(lossesTr)]
    predsTr <- lossesTr[,4:ncol(lossesTr)]
    
    lossesTe <- lossesTe[,1:3]
    lossesTr <- lossesTr[,1:3]
    
    pinballsTe <- sapply(1:(ncol(predsTe)), function(i) {
      # i <- 5
      loss(predsTe[,i], trainData[,"y"], loss.type=list(name="pinball", tau=c(0, taus, 1)[i]))
    }, simplify="array")
    
    
    # Test
    
    lik <- -sum(log(lossesTe[,"lik"]))/nrow(lossesTe)
    cdfDhsic <- dhsic.test(trainData[,"x"], lossesTe[,"cdf"])$statistic
    qDhsic <- dhsic.test(trainData[,"x"], lossesTe[,"q"])$statistic
    pinballs <- apply(pinballsTe, 2, function(col) col/sd(col, na.rm=T))
    pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
    # uniformity test
    #hist(lossesTe[,"cdf"])
    unifT <- ks.test(lossesTe[,"cdf"]+rnorm(length(lossesTe[,"cdf"]),mean=0,sd=1e-10), punif)$statistic
    # discretized indep test -> chi2 indep test
    # table(lossesTe[,"q"], findInterval(trainData[,"x"], taus, rightmost.closed=T))
    chi2T <- chisq.test(lossesTe[,"q"], findInterval(trainData[,"x"], taus, rightmost.closed=T) )$statistic
    lossesTe <- c(lik, cdfDhsic, qDhsic, pinball, unifT, chi2T)
    
    # w <- runif(100,-1,1)
    # z <- 5*x+runif(100,-1,1)
    # w <- findInterval(w, c(-0.5,0,0.5))
    # z <- findInterval(z, c(-0.5,0,0.5))
    # table(w, z)
    # chisq.test(w, z)
    # chisq.test(w, z)$statistic
    
    # Train
    
    indxTr <- unlist(sapply(1:numFolds, function(f){
      indx <- seq((f-1)*size + 1, min(f*size,n))
      indx <- setdiff(1:n, indx)
      return(indx)
    }, simplify=FALSE))
    
    
    pinballsTr <- sapply(1:(ncol(predsTr)), function(i) {
      # i <- 1
      loss(predsTr[,i], trainData[indxTr,"y"], loss.type=list(name="pinball", tau=c(0,  taus, 1)[i]))
    }, simplify="array")
    
    lik <- -sum(log(lossesTr[,"lik"]))/nrow(lossesTr)
    cdfDhsic <- dhsic.test(trainData[indxTr,"x"], lossesTr[,"cdf"])$statistic
    qDhsic <- dhsic.test(trainData[indxTr,"x"], lossesTr[,"q"])$statistic
    pinballs <- apply(pinballsTr, 2, function(col) col/sd(col, na.rm=T))
    pinball <- sum(pinballs, na.rm=T)/sum(!is.na(pinballs))
    # uniformity test
    unifT <- ks.test(lossesTr[,"cdf"]+rnorm(length(lossesTr[,"cdf"]),mean=0,sd=1e-10), punif)$statistic
    # discretized indep test -> chi2 indep test
    chi2T <- chisq.test(lossesTr[,"q"], findInterval(trainData[indxTr,"x"], taus, rightmost.closed=T) )$statistic
    lossesTr <- c(lik, cdfDhsic, qDhsic, pinball, unifT, chi2T)
    
    losses <- cbind(lossesTe, lossesTr)
    
    return(losses)
  }, p=params, mc.cores=2, SIMPLIFY="array")
  
  
  dimnames(losses) <- list(loss=c("lik", "cdfDhsic", "qDhsic","pinball", "unif", "chi2"), testTrain=c("test","train"), params=names(params))
  
  #print("exits CV parallel")
  return(losses)
}



taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)
lambdas <- c(0.1, 1,10, 100, 1000)
sigmasXY <-  getRangeRbf2(x=as.matrix(xs), length.out=5)
sigmasXY <- 10^seq(log(sigmasXY[1],10), log(sigmasXY[2],10), length.out=5)
sigmasYX <- getRangeRbf2(as.matrix(ys), length.out=5)
sigmasYX <- 10^seq(log(sigmasYX[1],10), log(sigmasYX[2],10), length.out=5)
sigmas <- sort(c(sigmasXY, sigmasYX))
sigmasXY <- sigmas
sigmasYX <- sigmas
numFolds <- 3



pm <- proc.time()
lossesXY <- CV.parallel.lik2(x=xs, y=ys, numFolds, taus, lambdas, sigmas=sigmasXY)
lossesYX <- CV.parallel.lik2(x=ys, y=xs, numFolds, taus, lambdas, sigmas=sigmasYX)
proc.time() - pm 
# 7.5 secs with 2 sigmas, 2 lambdas, 5 taus
# 27 secs with 4 sigmas, 4 lambdas, 5 taus
# 42 secs with 5 sigmas, 5 lambdas, 5 taus
# 83 secs with 10 sigmas, 5 lambdas, 5 taus

crit <- "cdfDhsic"

# xy
lossesXY_long <- melt(lossesXY[,"test",], stringsAsFactors=FALSE)
lossesXY_loss <- cast(lossesXY_long, params~loss, value="value")
aux <- as.character(lossesXY_loss$params)
lossesXY_loss$lambda <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[1]), function(el) strsplit(el, split="=")), function(el) el[2]))
lossesXY_loss$sigma <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[2]), function(el) strsplit(el, split="=")), function(el) el[2]))

ranks <- apply(lossesXY_loss[,c("cdfDhsic","chi2","lik","pinball","qDhsic","unif")],2, function(col) order(col, decreasing=FALSE))
colnames(ranks) <- c("cdfDhsic","chi2","lik","pinball","qDhsic","unif")
ranks[order(ranks[,"cdfDhsic"]),]

plot(lossesXY_loss$cdfDhsic, lossesXY_loss$chi2)
pairs(lossesXY_loss[,c("cdfDhsic","chi2","lik","pinball","qDhsic","unif")])

pct <- 0.2
indx <- which(lossesXY_loss$lik < quantile(lossesXY_loss$lik, pct) & lossesXY_loss$cdfDhsic < quantile(lossesXY_loss$cdfDhsic, pct)
              & lossesXY_loss$chi2 < quantile(lossesXY_loss$chi2, pct))
lossesXY_loss[indx,]
ranks[indx,]

num <- 1
pct <- num/nrow(lossesXY_loss)
indx <- numeric()
while(length(indx)==0){
  print(num)
  indx <- which(lossesXY_loss$lik < quantile(lossesXY_loss$lik, pct) & lossesXY_loss$cdfDhsic < quantile(lossesXY_loss$cdfDhsic, pct)
                & lossesXY_loss$chi2 < quantile(lossesXY_loss$chi2, pct))
  num <- num +1
  pct <- num/nrow(lossesXY_loss)
}
length(indx)
lossesXY_loss[indx,]
ranks[indx,]

indxOptXY <- indx #which.min(lossesXY_loss[,crit])
parsOptXY <- lossesXY_loss[indxOptXY,]
parsOptXY


# yx
lossesYX_long <- melt(lossesYX[,"test",], stringsAsFactors=FALSE)
lossesYX_loss <- cast(lossesYX_long, params~loss, value="value")
aux <- as.character(lossesYX_loss$params)
lossesYX_loss$lambda <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[1]), function(el) strsplit(el, split="=")), function(el) el[2]))
lossesYX_loss$sigma <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[2]), function(el) strsplit(el, split="=")), function(el) el[2]))
indxOptYX <- which.min(lossesYX_loss[,crit])
parsOptYX <- lossesYX_loss[indxOptYX,]
parsOptYX



# check number of models tried whre x->y was better than y->x
plot(lossesXY_loss$lik, lossesYX_loss$lik, xlim=range(lossesXY_loss$lik, lossesYX_loss$lik), ylim=range(lossesXY_loss$lik, lossesYX_loss$lik))
abline(a=0, b=1, col="red")
(sum(lossesXY_loss$lik < lossesYX_loss$lik)/length(lossesYX_loss$lik) > 0.5)*1

plot(lossesXY_loss$cdfDhsic, lossesYX_loss$cdfDhsic, xlim=range(lossesXY_loss$cdfDhsic, lossesYX_loss$cdfDhsic), ylim=range(lossesXY_loss$cdfDhsic, lossesYX_loss$cdfDhsic))
abline(a=0, b=1, col="red")
(sum(lossesXY_loss$cdfDhsic < lossesYX_loss$cdfDhsic)/length(lossesYX_loss$cdfDhsic) > 0.5)*1

getParsOpt <- function(losses, trainTest, criterion){
  losses_long <- melt(losses[,trainTest,], stringsAsFactors=FALSE)
  losses_loss <- cast(losses_long, params~loss, value="value")
  aux <- as.character(losses_loss$params)
  losses_loss$lambda <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[1]), function(el) strsplit(el, split="=")), function(el) el[2]))
  losses_loss$sigma <- as.numeric(sapply(sapply(sapply(strsplit(aux, " "), function(el) el[2]), function(el) strsplit(el, split="=")), function(el) el[2]))
  
  num <- 1
  pct <- num/nrow(losses_loss)
  indxOpt <- numeric()
  while(length(indxOpt)==0){
    #print(num)
    indxOpt <- which(losses_loss$lik < quantile(losses_loss$lik, pct) & losses_loss$cdfDhsic < quantile(losses_loss$cdfDhsic, pct)
                  & losses_loss$chi2 < quantile(losses_loss$chi2, pct))
    num <- num +1
    pct <- num/nrow(losses_loss)
  }
  indxOpt <- indxOpt[which.min(losses_loss$lik[indxOpt])]
  
  #indxOpt <- which.min(losses_loss[,criterion])
  parsOpt <- losses_loss[indxOpt,]
  return(list(parsOpt=parsOpt, losses_loss=losses_loss))
}


crit <- "lik"
parsOptXYte <- getParsOpt(lossesXY, trainTest="test", criterion=crit)
parsOptYXte <- getParsOpt(lossesYX, trainTest="test", criterion=crit)
parsOptXYte$parsOpt
parsOptYXte$parsOpt

parsOptXYtr <- getParsOpt(lossesXY, trainTest="train", criterion=crit)
parsOptYXtr <- getParsOpt(lossesYX, trainTest="train", criterion=crit)
parsOptXYtr$parsOpt
parsOptYXtr$parsOpt



pm <- proc.time()
predsTe <- sapply(c(0, taus,1), function(tau){
  # tau <- 0.05
  Cmax <- 100000  
  qrmXY <- NULL
  qrmYX <- NULL
  
  print(paste("tau: ", tau))
  if(tau==0){
    Caux <- Cmax
    while(class(qrmXY)=="try-error"|is.null(qrmXY)){
      qrmXY <- kqr(xs, ys, tau = 1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predXY <- predict(qrmXY, xs)
    Caux <- Cmax
    while(class(qrmYX)=="try-error"|is.null(qrmYX)){
      qrmYX <- kqr(ys, xs, tau = 1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  } else if(tau==1){
    Caux <- Cmax
    while(class(qrmXY)=="try-error"|is.null(qrmXY)){
      qrmXY <- kqr(xs, ys, tau = 1-1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predXY <- predict(qrmXY, xs)
    Caux <- Cmax
    while(class(qrmYX)=="try-error"|is.null(qrmYX)){
      qrmYX <- kqr(ys, xs, tau = 1-1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  } else{
    qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptXYte$parsOpt["sigma"])), C=as.numeric(parsOptXYte$parsOpt["lambda"]))
    predXY <- predict(qrmXY, xs)
    qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptYXte$parsOpt["sigma"])), C=as.numeric(parsOptYXte$parsOpt["lambda"]))
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  }
  return(pred)
}, simplify="array")
dimnames(predsTe) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=c(0,  taus, 1))
proc.time() - pm # 2.5 seconds

pm <- proc.time()
predsTr <- sapply(c(0, taus,  1), function(tau){
  # tau <- 0.05
  Cmax <- 100000  
  qrmXY <- NULL
  qrmYX <- NULL
  
  print(paste("tau: ", tau))
  if(tau==0){
    Caux <- Cmax
    while(class(qrmXY)=="try-error"|is.null(qrmXY)){
      qrmXY <- kqr(xs, ys, tau = 1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predXY <- predict(qrmXY, xs)
    Caux <- Cmax
    while(class(qrmYX)=="try-error"|is.null(qrmYX)){
      qrmYX <- kqr(ys, xs, tau = 1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  } else if(tau==1){
    Caux <- Cmax
    while(class(qrmXY)=="try-error"|is.null(qrmXY)){
      qrmXY <- kqr(xs, ys, tau = 1-1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predXY <- predict(qrmXY, xs)
    Caux <- Cmax
    while(class(qrmYX)=="try-error"|is.null(qrmYX)){
      qrmYX <- kqr(ys, xs, tau = 1-1/(10000*length(xs)), kernel = "rbfdot", kpar= list(sigma=1), C=Caux)
      Caux <- Caux*0.1
    }
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  } else{
    qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptXYtr$parsOpt["sigma"])), C=as.numeric(parsOptXYtr$parsOpt["lambda"]))
    predXY <- predict(qrmXY, xs)
    qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptYXtr$parsOpt["sigma"])), C=as.numeric(parsOptYXtr$parsOpt["lambda"]))
    predYX <- predict(qrmYX, ys)
    pred <- cbind(predXY, predYX)
  }
  return(pred)
}, simplify="array")
dimnames(predsTr) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=c(0, taus, 1))
proc.time() - pm # 2.5 seconds



par(mfrow=c(2,2))
plot(xs, ys)
indx <- order(xs)
for(tau in c(0, taus, 1)) lines(xs[indx], predsTe[indx,"xy", as.character(tau)], col=which(tau==c(0,  taus, 1)))
indx <- order(ys)
plot(ys, xs, ylim=range(predsTe[,"yx",]))
for(tau in c(0, taus, 1)) lines(ys[indx], predsTe[indx,"yx", as.character(tau)], col=which(tau==c(0,  taus, 1)))
plot(xs, ys)
indx <- order(xs)
for(tau in c(0, taus, 1)) lines(xs[indx], predsTr[indx,"xy", as.character(tau)], col=which(tau==c(0,  taus, 1)))
indx <- order(ys)
plot(ys, xs)
for(tau in c(0, taus, 1)) lines(ys[indx], predsTr[indx,"yx", as.character(tau)], col=which(tau==c(0, taus,1)))


######################################################################################################################*
# IDEA REGARDING MINIMUM MESSAGE LENGTH OF ZERO LOSS ENCODING -> OPTIMIZING TRAIN 
######################################################################################################################*
{

# HERE I AM LOOKING FOR 2 MODELS, ONE FOR X->Y AND Y->X, WITH SIMILAR LOW LIKELIHOOD
likDif <- sapply(lossesYX_loss$lik, function(lik1) sapply(lossesXY_loss$lik, function(lik2) abs(lik1-lik2)))
likYX <- sapply(lossesYX_loss$lik, function(lik1) sapply(lossesXY_loss$lik, function(lik2) lik1))
likXY <- sapply(lossesYX_loss$lik, function(lik1) sapply(lossesXY_loss$lik, function(lik2) lik2))
indxsYX <- sapply(lossesYX_loss$lik, function(lik1) sapply(lossesXY_loss$lik, function(lik2) which(lik1==lossesYX_loss$lik)))
indxsXY <- sapply(lossesYX_loss$lik, function(lik1) sapply(lossesXY_loss$lik, function(lik2) which(lik2==lossesXY_loss$lik)))
likDif <- as.numeric(likDif)
likXY <- as.numeric(likXY)
likYX <- as.numeric(likYX)
indxsXY <- as.numeric(indxsXY)
indxsYX <- as.numeric(indxsYX)
ord <- order(likDif)
likDif <- likDif[ord]
likXY <- likXY[ord]
likYX <- likYX[ord]
indxsXY <- indxsXY[ord]
indxsYX <- indxsYX[ord]

summary(likDif)


plot(likDif, likXY, col="blue", pch=2)
lines(likDif, likYX, col="red", pch=3, type="p")

indxYX <- which.min(likYX)
indxXY <- which.min(abs(likXY-likYX[indxYX]))
lines(likDif[c(indxYX, indxXY)], c(likXY[indxXY], likYX[indxYX]), col="green", cex=2, type="p")

indxXY_chosen <- indxsXY[indxXY[1]]
indxYX_chosen <- indxsYX[indxYX[1]]

parXY_chosen <-  lossesXY_loss[indxXY_chosen,]
parYX_chosen <-  lossesYX_loss[indxYX_chosen,]
parXY_chosen
parYX_chosen

pm <- proc.time()
predsTe <- sapply(taus, function(tau){
  # tau <- 0.05
  print(paste("tau: ", tau))
  qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parXY_chosen["sigma"])), C=as.numeric(parXY_chosen["lambda"]))
  predXY <- predict(qrmXY, xs)
  qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parYX_chosen["sigma"])), C=as.numeric(parYX_chosen["lambda"]))
  predYX <- predict(qrmYX, ys)
  pred <- cbind(predXY, predYX)
}, simplify="array")
dimnames(predsTe) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=taus)
proc.time() - pm # 2.5 seconds

pm <- proc.time()
predsTr <- sapply(taus, function(tau){
  # tau <- 0.05
  print(paste("tau: ", tau))
  qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptXYtr["sigma"])), C=as.numeric(parsOptXYtr["lambda"]))
  predXY <- predict(qrmXY, xs)
  qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptYXtr["sigma"])), C=as.numeric(parsOptYXtr["lambda"]))
  predYX <- predict(qrmYX, ys)
  pred <- cbind(predXY, predYX)
}, simplify="array")
dimnames(predsTr) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=taus)
proc.time() - pm # 2.5 seconds


par(mfrow=c(2,2))
plot(xs, ys)
indx <- order(xs)
for(tau in taus) lines(xs[indx], predsTe[indx,"xy", as.character(tau)], col=which(tau==taus))
indx <- order(ys)
plot(ys, xs)
for(tau in taus) lines(ys[indx], predsTe[indx,"yx", as.character(tau)], col=which(tau==taus))
plot(xs, ys)
indx <- order(xs)
for(tau in taus) lines(xs[indx], predsTr[indx,"xy", as.character(tau)], col=which(tau==taus))
indx <- order(ys)
plot(ys, xs)
for(tau in taus) lines(ys[indx], predsTr[indx,"yx", as.character(tau)], col=which(tau==taus))

pm <- proc.time()
alphaTe <- sapply(taus, function(tau){
  # tau <- 0.05
  print(paste("tau: ", tau))
  qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptXYte["sigma"])), C=as.numeric(parsOptXYte["lambda"]))
  qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptYXte["sigma"])), C=as.numeric(parsOptYXte["lambda"]))
  res <- cbind(qrmXY@alpha, qrmYX@alpha)
}, simplify="array")
dimnames(alphaTe) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=taus)
proc.time() - pm # 2.5 seconds

pm <- proc.time()
alphaTr <- sapply(taus, function(tau){
  # tau <- 0.05
  print(paste("tau: ", tau))
  qrmXY <- kqr(xs, ys, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptXYtr["sigma"])), C=as.numeric(parsOptXYtr["lambda"]))
  qrmYX <- kqr(ys, xs, tau = tau, kernel = "rbfdot", kpar= list(sigma=as.numeric(parsOptYXtr["sigma"])), C=as.numeric(parsOptYXtr["lambda"]))
  res <- cbind(qrmXY@alpha, qrmYX@alpha)
}, simplify="array")
dimnames(alphaTr) <- list(obs=1:length(xs), direction=c("xy", "yx"), tau=taus)
proc.time() - pm # 2.5 seconds


alphaTe[,"xy",]

apply(alphaTe[,"xy",], 2, function(col) table(col))
apply(alphaTe[,"xy",], 2, function(col) table(round(col,2)))
apply(alphaTe[,"yx",], 2, function(col) table(round(col,2)))

chngsXY <- apply(round(alphaTe[,"xy",],2), 2, function(col) findInterval(col, range(col)))
sum(t(diff(t(chngsXY)))!=0)
chngsYX <- apply(round(alphaTe[,"yx",],2), 2, function(col) findInterval(col, range(col)))
sum(t(diff(t(chngsYX)))!=0)
}
######################################################################################################################*



getLosses <- function(xs, ys, pars, taus){
  
  preds <- predQ2(xTr=xs, yTr=ys, xTe=xs, yTe=ys, sig=pars["sigma"][1,1], lam=pars["lambda"][1,1], taus, trainTest="train")
  
  liks <- lossesQ2(preds, xs, ys, pars["sigma"], taus)
  
  
  return(list(losses=liks, preds=preds))
}

lossesXY <- getLosses(xs, ys, pars=parsOptXY, taus)
lossesYX <- getLosses(ys, xs, pars=parsOptYX, taus)

plotLosses <- function(xs, ys, lossesXY, lossesYX, taus){
  par(mfrow=c(3,2))
  plot(xs, ys)
  indx <- order(xs)
  for(i in 1:(ncol(lossesXY$preds))) lines(xs[indx], lossesXY$preds[indx, i], col=i)
  indx <- order(ys)
  plot(ys, xs)
  for(i in 1:(ncol(lossesYX$preds))) lines(ys[indx], lossesYX$preds[indx, i], col=i)
  plot(xs, lossesXY$losses[,"cdf"])
  plot(ys, lossesYX$losses[,"cdf"])
  plot(xs, lossesXY$losses[,"q"])
  plot(ys, lossesYX$losses[,"q"])
  par(mfrow=c(1,1))
}

plotLosses(xs, ys, lossesXY, lossesYX, taus)


getPerf <- function(losses, x, taus){
  lik <- -sum(log(losses[,"lik"]))/nrow(losses)
  cdfDhsic <- dhsic.test(x, losses[,"cdf"])$statistic
  qDhsic <- dhsic.test(x, losses[,"q"])$statistic
  
  unifT <- ks.test(losses[,"cdf"]+rnorm(length(losses[,"cdf"]),mean=0,sd=1e-10), punif)$statistic
  chi2T <- chisq.test(losses[,"q"], findInterval(x, taus, rightmost.closed=T) )$statistic
  
  
  res <- c(lik, cdfDhsic, qDhsic, unifT, chi2T)
  names(res) <- c("lik", "cdfDhsic", "qDhsic","unif","chi2")
  return(res)
}

getPerf(lossesXY$losses, xs, taus)
getPerf(lossesYX$losses, ys, taus)





lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)

numFolds <- 3
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)
crit <- "cdfDhsic"
sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100)

pm <- proc.time()
msrsQ <- mcmapply(function(el, nm){
  # difficult cases x*n : 41, 44, 57, 62
  # difficult cases sin(x) + n: 84
  # difficult cases tubingen pairs 47
  # i <- 47; el <- dataTestList$xs[[i]]; nm <- dataTestList$names[i]
  X <- apply(el, 2, stdrize)
  print("*********************************************")
  print(paste("name: ", nm))
  #print("head(X)")
  #print(head(X))
  #print(apply(X, 2, mean))
  #print(apply(X, 2, sd))
  
  x <- X[,1]
  y <- X[,2]
  #indx <- which(x>-0.5 & x<0.5 & y>-0.5 & y<0.5)
  #x <- x[indx]
  #y <- y[indx]

  # par(mfrow=c(1,2))
  # plot(x, y)
  # plot(y, x)
  # par(mfrow=c(1,1))

  #trainDataXY <- constructData(as.matrix(x[ord]), y[ord])
  #trainDataYX <- constructData(as.matrix(y[ord]), x[ord])
  
  #sigmasXY <-  getRangeRbf2(x=as.matrix(x))
  #sigmasYX <- getRangeRbf2(x=as.matrix(y))
  
  
  #rng <- range(sigmasXY, sigmasYX)
  
  #sigmasXY <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=5)
  #sigmasYX <- 10^seq(log(rng[1],10), log(rng[2],10), length.out=5)
  
  lossesXY_CV <- CV.parallel.lik2(x=x, y=y, numFolds, taus, lambdas, sigmas=sigmas)
  lossesYX_CV <- CV.parallel.lik2(x=X[,2], y=X[,1], numFolds, taus, lambdas, sigmas=sigmas)
  
  
  
  parsOptXY <- getParsOpt(losses=lossesXY_CV, trainTest="test", criterion=crit)
  parsOptYX <- getParsOpt(lossesYX_CV, "test", criterion=crit)
  
  pars <- rbind(parsOptXY$parsOpt[,c("lambda","sigma")], parsOptYX$parsOpt[,c("lambda","sigma")])
  rownames(pars) <- c("xy","yx")
  
  #plot(parsOptXY$losses_loss$lik, parsOptYX$losses_loss$lik); abline(a=0, b=1, col="red")
  propMods <- sum(parsOptXY$losses_loss[,crit] < parsOptYX$losses_loss[,crit])/length(parsOptXY$losses_loss[,crit])
  
  
  lossesXY <- getLosses(xs=x, ys=y, pars=parsOptXY$parsOpt, taus)
  lossesYX <- getLosses(xs=y, ys=x, pars=parsOptYX$parsOpt, taus)
  
  plotLosses(x, y, lossesXY, lossesYX, taus)
  
  perfXY <- getPerf(losses=lossesXY$losses, x, taus)
  perfYX <- getPerf(lossesYX$losses, y, taus)
  
  msrsCV <- parsOptXY$parsOpt[c("lik","cdfDhsic","qDhsic","unif","chi2")]-parsOptYX$parsOpt[c("lik","cdfDhsic","qDhsic","unif","chi2")]
  msrs <- perfXY-perfYX
  
  print("pars:")
  print(pars)
  print("propMods")
  propMods
  print("msrs:")
  print(msrs)
  print("msrsCV:")
  print(msrsCV)
  
  

  res <- as.matrix(rbind(msrs, msrsCV))
  rownames(res) <- c("all","CV")
  
  return(list(msrs=res, propMods=propMods, pars=pars))
}, el=dataTestList$xs[1:length(dataTestList$xs)], nm=dataTestList$names[1:length(dataTestList$xs)], 
SIMPLIFY=FALSE, mc.cores=1)
proc.time() - pm 
# 12.5 mins with length(taus)=5, length(lambdas)=2, length(sigmas)=2 and numFolds=3, n=100, q=100 for y=sin(x)+n, x~U(-pi,pi), y~U(-0.5,0.5)
# 12.5 mins with length(taus)=5, length(lambdas)=2, length(sigmas)=2 and numFolds=3, n=100, q=100 for y=sin(x)*n, x~N(0,1), y~U(-1,1)
# 81 mins with length(taus)=5, length(lambdas)=5, length(sigmas)=5 and numFolds=3, n=100, q=100 for y=x*n, x~N(0,1), y~U(-1,1)

parsDB <- sapply(msrsQ, function(el) el[[3]], simplify="array")
parsDB <- unlist(parsDB)
parsDB <- data.frame(lambdaXY=parsDB[seq(1, length(parsDB), 4)], lambdaYX=parsDB[seq(2, length(parsDB), 4)], sigmaXY=parsDB[seq(3,length(parsDB),4)], sigmaYX=parsDB[seq(4, length(parsDB),4)])
summary(parsDB)
hist(log(parsDB$sigmaXY, 10))
hist(log(parsDB$sigmaYX, 10))

table(round(log(parsDB$sigmaXY, 10)))
table(round(log(parsDB$sigmaYX, 10)))
table(round(log(parsDB$lambdaXY, 10)))
table(round(log(parsDB$lambdaYX, 10)))

table(round(log(parsDB$sigmaXY, 10)), round(log(parsDB$lambdaXY)))
table(round(log(parsDB$sigmaYX, 10)), round(log(parsDB$lambdaYX)))


propMods <- sapply(msrsQ, function(el) el[[2]])
sum(propMods>0.5) 
#94 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma and trained on all data, setting q[1]=min(qs), q[n]=max(qs)
#99 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma and trained on all data
#100 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma and trained on all data, deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, 
# sigmas set using 0.01-0.99 criterion
#100 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma and trained on all data, deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, 
# sigmas set using 0.1-0.9 criterion
#74 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, depending on lambda and sigma and trained on all data, setting q[1]=min(qs), q[n]=max(qs)
# crit <- "lik", sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
#74 for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, depending on lambda and sigma and trained on all data, setting q[1]=min(qs), q[n]=max(qs)
# crit <- "lik/cdfDhsic/chi2", sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)



#91 for y=sin(x)+n, x~N(-pi,pi), n~U(-0.5,0.5) using q_tau_min, q_tau_max, NOT depending on lambda and sigma and trained on all data, deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, 
# sigmas set using 0.1-0.9 criterion 

# 41 for TCEPs using q_tau_min, q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs
# sigmas <- c(0.001, 0.01, 0.1, 1 ,10), lambdas <- c(0.01, 0.1, 1, 10, 100, 1000),  crit <- "lik"
# 46 for TCEPs using q_tau_min, q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"
# 44 for TCEPs using q_tau_min, q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"
# 48 for TCEPs using q_tau_min, q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik/cdfDhsic/chi2"



msrsQ2 <- sapply(msrsQ, function(el) el[[1]], simplify="array")

msrsDB <- melt(msrsQ2)
head(msrsDB)
colnames(msrsDB) <- c("type","msr","rep","value")
cast(msrsDB, type~msr, fun.aggregate=length)
cast(msrsDB, type~msr, fun.aggregate=function(cell) sum(cell<0)/length(cell)*100)

# 95/96 all/CV for y=sin(x)+n, x~U(-pi,pi), n~U(-0.5,0.5) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# 44/69 for cdfDhsic 74/95 all/CV for y=sin(x)+n, x~U(-pi,pi), n~U(-0.5,0.5) using q_tau_min, q_tau_max, depending on lambda and sigma 
# fixing no increasing qs setting q0=minqs and q1=maxqs, also changed range sigma 
# 44/84 for cdfDhsic 93/96 for lik,  44/86 for qDhsic all/CV for y=sin(x)+n, x~N(-pi,pi), n~U(-0.5,0.5) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, , sigmas set using 0.1-0.9 criterion, used cdfDhsic for CV


# 58/88 for cdfDhsic 85/98 for lik,  all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# 54/82 for cdfDhsic 66/94 for lik,  all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, depending on lambda and sigma 
# fixing no increasing qs setting q0=minqs and q1=maxqs, also changed range sigma 
# 47/90 for cdfDhsic 74/98 for lik,  51/92 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and setting q[1]=min(qs), q[n]=max(qs)
# 53/80 for cdfDhsic 71/97 for lik,  52/78 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and setting q[1]=min(qs), q[n]=max(qs)
# 48/74 for cdfDhsic 58/100 for lik,  49/66 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, , sigmas set using 0.01-0.99 criterion
# 48/74 for cdfDhsic 74/98 for lik,  46/62 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, , sigmas set using 0.1-0.9 criterion
# 53/96 for cdfDhsic 57/67 for lik,  46/81 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, NOT depending on lambda and sigma
# trained on all data, and deleting q[j] if ! q[1]<q[j]<q[m] for j=2,...,m, , sigmas set using 0.1-0.9 criterion, used cdfDhsic for CV
# 70/71 for cdfDhsic 80/90 for lik,  75/69 for qDhsic all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, depending on lambda and sigma
# fixing no increasing qs setting q0=minqs and q1=maxqs, sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"
# 79/76 for cdfDhsic, 66/60 for chi2, 77/87 for lik,  82/68 for qDhsic, 66/52 unif all/CV for y=x*n, x~N(0,1), n~U(-1,1) using q_tau_min, q_tau_max, depending on lambda and sigma
# fixing no increasing qs setting q0=minqs and q1=maxqs, sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"




# 58.6/57.6 for cdfDhsic 43.4/40.4 for lik, 56.6/54.5 for qDhsic, all/CV for TCEPs using q_tau_min, 
#q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs,  
# sigmas <- c(0.001, 0.01, 0.1, 1 ,10), crit <- "lik"
# 57.6/57.6 for cdfDhsic 44.4/39.4 for lik, 56.6/56.6 for qDhsic, all/CV for TCEPs using q_tau_min, 
#q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs,  
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"
# 58.6/57.6 for cdfDhsic 45.5/36.4 for lik, 55.6/53.5 for qDhsic, all/CV for TCEPs using q_tau_min, 
#q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs,  
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik"
# 58.6/54.5 for cdfDhsic , 49.5/52.5 for chi2, 45.5/42.4 for lik, 59.6/49.5 for qDhsic, 46.5,51.5 for unif all/CV for TCEPs using q_tau_min, 
#q_tau_max depending on lambda and sigma fixing no increasing qs setting q0=minqs and q1=maxqs,  
# sigmas <- c(0.0001, 0.001, 0.01, 0.1, 1 ,10, 100), lambdas <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), crit <- "lik/cdfDhsic/chi2"



