# measuring complexity of distributions of data
print("in func_complexity_v1.R")
#library(dHSIC) # dhsic.test
source("./pkg_dHSIC/dHSIC.R")
library(Rcpp) # sourceCpp

# download entropy functions needed for certain complexity functions such as sumMarginalEntropies
source("./pkg_entropy/func_entropy_v1.R")

# Score functions for residuals (model based) or variables (non-parametric)

# Mooij's HSIC function

sourceCpp(file="./pkg_cmplxtyScores/hsic/hsic.cpp")

HSICMooij <- function(X, Y, sigma_x=NA, sigma_y=NA, method=c("gamma","permutation"), unbiased=TRUE, nPerms=1000, dHSICdY=FALSE){
  
  # validate X and Y
  if(class(X) != "matrix" | class(Y) != "matrix"){
    if(class(X) == "numeric" & class(Y) == "numeric"){
      X <- matrix(X, length(X), 1)
      Y <- matrix(Y, length(Y), 1)
    } else{
      stop("X and Y must be numeric vectors or matrices")
    }
  } 
  n <- nrow(X)
  if(n != nrow(Y)) stop("X and Y must contain same number of rows (observations).")
  
  # validate method
  method <- match.arg(method)
  
  # validate sigma_x
  if(!is.na(sigma_x)){
    if(class(sigma_x) != "numeric"){
      stop("sigma_x must be of type numeric")
    } else{
      if(length(sigma_x)>1){
        stop("sigma_x must be of length 1")
      } else{
        if(sigma_x <= 0) stop("sigma_x must be positive")
      }
    }
  } 
  
  # validate sigma_y
  if(!is.na(sigma_y)){
    if(class(sigma_y) != "numeric"){
      stop("sigma_y must be of type numeric")
    } else{
      if(length(sigma_y)>1){
        stop("sigma_y must be of length 1")
      } else{
        if(sigma_y <= 0) stop("sigma_y must be positive")
      }
    }
  }
  
  # validate unbiased
  if(class(unbiased) != "logical"){
    stop("unbiased must be a logical")
  } else{
    if(length(unbiased)>1) stop("unbiased must be of length 1")
  }
  
  # validate nPerms
  if(method=="permutation"){
    if(class(nPerms)!="numeric") stop("nPerms must be numeric")
    if(length(nPerms)>1) stop("nPerms must be of length 1")
    if(nPerms<1 | (nPerms-round(nPerms)) != 0) stop("nPerms must be integer greater or equal to 1")
  }
  
  # validate dHSICdY
  if(class(dHSICdY) != "logical"){
    stop("dHSICdY must be a logical")
  } else{
    if(length(dHSICdY)>1) stop("unbiased must be of length 1")
  }
  
  
  px <- ncol(X)
  py <- ncol(Y)
  x <- as.numeric(X)
  y <- as.numeric(Y)
  
  #sigma_x = bandwidth of RBF kernel for X (0.0 means chosen by a heuristic)
  #sigma_y = bandwidth of RBF kernel for Y (0.0 means chosen by a heuristic)
  if(is.na(sigma_x)) sigma_x <- 0
  if(is.na(sigma_y)) sigma_y <- 0
  
  # nrperm = number of permutations in permutation test 
  #   <  0:  use unbiased HSIC (only calculates fields HSIC, p_value as output), see [2]
  #   == 0:  use approximated gamma distribution instead of permutation test, see [1]
  #   >  0:  use biased HSIC (only calculates fields HSIC, p_value as output), see [1]
  
  if(method=="gamma"){
    nrperm <- 0
  } else{
    if(unbiased){
      nrperm <- -nPerms
    } else{
      nrperm <- nPerms
    }
  }
  
  
  
  res <- calcHSICr(n, px, py, x, y, sigma_x, sigma_y, nrperm, dHSICdY)
  return(res)
}


score_sumLogLik <- function(vars, method){
  
  score <- sum(vars)
  
  return(score)
}


# HSIC - max pairwise HSIC, Using Mooij implementation above
score_pvalHSIC <- function(vars, method){
  # n <- 100; vars <- as.matrix(data.frame(x=rnorm(n), y=rnorm(n), z=rnorm(n))); method="gamma"
  varsCombos <- as.data.frame(combn(colnames(vars), 2))
  
  varsList <-  lapply(varsCombos,  function(col) vars[,col])
  scores <- sapply(varsList, function(el){
    HSICMooij(el[,1], el[,2], method=method)$p.value
  })
  
  # take worse
  score <- 1-min(scores)
  return(score)
}

score_pHSIC <- function(vars, method){
  # n <- 100; vars <- as.matrix(data.frame(x=rnorm(n), y=rnorm(n), z=rnorm(n))); method="gamma"
  varsCombos <- as.data.frame(combn(colnames(vars), 2))
  
  varsList <-  lapply(varsCombos,  function(col) vars[,col])
  scores <- sapply(varsList, function(el){
    HSICMooij(el[,1], el[,2], method=method)$p.value
  })
  
  # take worse
  score <- min(scores)
  score <- -log(score)
  return(score)
}

score_HSIC <- function(vars, method){
  # n <- 100; vars <- as.matrix(data.frame(x=rnorm(n), y=rnorm(n), z=rnorm(n))); method="gamma"
  varsCombos <- as.data.frame(combn(colnames(vars), 2))
  
  varsList <-  lapply(varsCombos,  function(col) vars[,col])
  scores <- sapply(varsList, function(el){
    HSICMooij(el[,1], el[,2], method=method)$statistic
  })
  
  # take worse
  score <- max(scores)
  return(score)
}

score_HSIC_fix <- function(vars, method){
  # n <- 100; vars <- as.matrix(data.frame(x=rnorm(n), y=rnorm(n), z=rnorm(n))); method="gamma"
  varsCombos <- as.data.frame(combn(colnames(vars), 2))
  
  varsList <-  lapply(varsCombos,  function(col) vars[,col])
  scores <- sapply(varsList, function(el){
    # print(el)
    # el <- varsList[[1]]
    # plot(el[,1], el[,2])
    HSICMooij(X=el[,1], Y=el[,2], method=method, sigma_x=1, sigma_y=1)$statistic
  })
  
  # take worse
  score <- max(scores)
  return(score)
}

# dHSIC
score_pvaldHSIC <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- 1-score$p.value
  return(score)
}

score_pdHSIC <- function(vars, method){
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- -log(score$p.value)
  return(score)
}

score_pdHSICgMin <- function(vars, method){
  resids <- vars[,"resid",]
  grps <- vars[,"grp",]
  resMins <- sapply(1:ncol(grps), function(j){
    # j <- 2
    #print(paste("j: ", j))
    grpsUni <- sort(unique(grps[,j]))
    if(length(grpsUni)!=1) grpsUni <- grpsUni[-c(1,length(unique(grps[,j])))]
    resByGrp <- sapply(grpsUni, function(grp){
      # grp <- 6
      #print(paste("grp: ", grp))
      indx <- which(grps[,j]==grp)
      if(length(indx)>1){
        score <- dhsic.test(resids[indx,], matrix.input=TRUE, pairwise=FALSE, method=method)$p.value  
      } else{
        score <- runif(1)
      }
      return(score)
    })
    res <- min(resByGrp)
    return(res)
  })
  res <- min(resMins)
  #score <- -log(res)
  score <- 1-res
  return(score)
}

score_pdHSICgMean <- function(vars, method){
  resids <- vars[,"resid",]
  grps <- vars[,"grp",]
  resMeans <- sapply(1:ncol(grps), function(j){
    # j <- 2
    grpsUni <- sort(unique(grps[,j]))
    if(length(grpsUni)!=1) grpsUni <- grpsUni[-c(1,length(unique(grps[,j])))]
    resByGrp <- sapply(grpsUni, function(grp){
      # grp <- 1
      indx <- which(grps[,j]==grp)
      if(length(indx)>1){
        score <- dhsic.test(resids[indx,], matrix.input=TRUE, pairwise=FALSE, method=method)$p.value  
      } else{
        score <- runif(1)
      }
      return(score)
    })
    res <- mean(resByGrp)
    return(res)
  })
  res <- mean(resMeans)
  #score <- -log(res)
  score <- 1-res
  return(score)
}

score_dHSIC <- function(vars, method){
  n <- nrow(vars)
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
  score <- score$statistic/n
  return(score)
}

score_dHSIC_fix <- function(vars, method){
  n <- nrow(vars)
  score <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method, kernel="gaussian.fixed", bandwidth=1)
  score <- score$statistic/n
  return(score)
}

score_pvalUnifPart <- function(vars, method, numParts){
  
  n <- nrow(vars)

  # n <- 150; numParts <- 8
  numPerPart <- ceil(n/numParts)
  smpl <- rep(1:numParts, numPerPart)[1:n]
  smpl <- smpl[sample(1:n)]
  
  # we sample which residuals to use for each hsic pvalue but we need at least two in each batch
  #smpl <- c(rep(1:numParts, 2) , sample(1:numParts, size=n-numParts*2, replace=TRUE))
  
  
  pvalsPart <- sapply(1:numParts, function(part) {
    indxPart <- which(smpl==part)
    
    
    pval <- dhsic.test(vars[indxPart,], matrix.input=TRUE, pairwise=FALSE, method=method)
    pval <- pval$p.value
    return(pval)
  })
  
  
  
  pvalsPart <- unique(na.omit(pvalsPart))
  if(length(pvalsPart)>2){
    # null is uniformity: high p-values good: dont reject uniform
    ks_test <- ks.test(pvalsPart, punif)  
    score <- ks_test$p.value
  } else{
    # null is independence: high p.values good: dont reject independence
    hsic_test <- dhsic.test(vars, matrix.input=TRUE, pairwise=FALSE, method=method)
    score <- hsic_test$p.value
  } 
  
  
  
  return(-log(score))
  
}

score_pvalUnifBoot <- function(vars, method, numSmpls){
  
  n <- nrow(vars)
  pvalsBoot <- sapply(1:numSmpls, function(i) {
    indxPart <- sample(1:n, replace=T)
    pval <- dhsic.test(vars[indxPart,], matrix.input=TRUE, pairwise=FALSE, method=method)
    pval <- pval$p.value
    return(pval)
  })
  ks_test <- ks.test(pvalsBoot, punif)
  score <- -log(1-ks_test$statistic)
  
  return(score)
  
}




score_sumMarginalEntropies <- function(vars, type=c("Shannon_1sp", "Shannon_Gauss","Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb", 
                                                    "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE", 
                                                    "Shannon_spacing_LL", "Shannon_KDP", "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1", 
                                                    "Shannon_MaxEnt2", "Shannon_expF", "Shannon_vME"), ...){
  

  pars <- list(...)
  
  
  score <- apply(vars, 2, function(col) do.call(type, c(list(x=col), pars)))
  score <- sum(score)
  
  return(score)
  
}

score_entropy <- function(vars, type=c("Shannon_kNN_k", "Shannon_KDP", "Shannon_Edgeworth", "Shannon_expF", "Shannon_vME"), ...){
  
  
  pars <- list(...)
  pars$x <- vars
  
  score <- do.call(type, pars)
  
  
  return(score)
  
}



# Functions to transform scores
# Ranking and probability transformation function  

# function to put each score on gaussian 0,1 scale
gaussianize <- function(x,...){
  
  Fn <- ecdf(x)
  y <- Fn(x)
  nug <- 0.001
  y[which(y==1)] <- 1- nug
  z <- qnorm(y)
  #hist(z)
  
  return(z)
}

# function to express each hypothesis score as a probability
scoreToProb <- function(x,...){
  
  #Fn <- ecdf(x)
  Fn <- ecdf(na.omit(x))
  y <- 1-Fn(x)
  
  #plot(x, y)
  
  nug <- 0.001
  y[which(y==0)] <- nug
  #z <- y / sum(y)
  z <- y / sum(y,na.rm=T)
  # plot(x, z)
  
  return(z)
}

# function to correct a given score for fact many hypothesis are "similar" (share many edges and non edges)
correctScoreToAdd <- function(x, hyps){
  
  p <- dim(hyps)[2]
  offDiag <- !diag(p)
  
  # take mean score for each edge and non edge
  
  unwrappedEdges <- t(apply(hyps, 3, function(col) col[offDiag]))
  unwrappedNonEdges <- t(apply((!hyps)*1, 3, function(col) col[offDiag]))
  
  
  scorePerEdge <- unwrappedEdges*x
  scorePerNonEdge <- unwrappedNonEdges*x
  # Inf*0 produces NaN and well define Inf*0 in this context to be zero
  #scorePerEdge[is.nan(scorePerEdge)] <- 0
  #scorePerNonEdge[is.nan(scorePerNonEdge)] <- 0
  
  #meanScorePerEdge <- apply(scorePerEdge, 2, mean)
  #meanScorePerNonEdge <- apply(scorePerNonEdge, 2, mean)
  meanScorePerEdge <- apply(scorePerEdge, 2, mean, na.rm=T)
  meanScorePerNonEdge <- apply(scorePerNonEdge, 2, mean, na.rm=T)
  
  
  
  # sum score of each edge  and non-edge in hypothesis 
  # to get corrected hypothesis score
  
  corrX <- as.numeric(unwrappedEdges %*% meanScorePerEdge + unwrappedNonEdges %*% meanScorePerNonEdge)
  
  return(corrX)
}


# When there are only two scores- second is true hypothesis, low scores are better
# but these transformations switch it so that we give back good scores are higher

# in rankedDecisions in func_causalLearners.R we will pass the scores in both orders 
# rankFunc(score_h0, score_h1) and rankFunc(score_h1, score_h0)
# in the step before aggregateScores leaves   columns  with c(score_h0, score_h1) and 
# new score wll be rnk_0 <- rankFunc(score_h1, score_h0) so when score_hA < score_hB - > rnk_A > rnk_B
# ie it must flip scores

test_scores_pos <- c(0.5, 1) # first one better
test_scores <- c(-1, 1) # first one better

# BOTH SCORES MUST BE POSITIVE!! then smaller one will become larger one
# suppose score_h0 < score_h1 ->  quot(score_h0, score_h1) < 0 
quot <- function(score_hA, score_hB){
  decision <- c(-1,1)[(score_hA > score_hB)*1+1]
  rnk <- 1/min(score_hA, score_hB)*decision
  return(rnk)
}
new_test_scores_pos <- c(quot(test_scores_pos[2],test_scores_pos[1]), quot(test_scores_pos[1],test_scores_pos[2]))
which.min(test_scores_pos)
which.max(new_test_scores_pos)


# BOTH SCORES MUST BE POSITIVE!! then smaller one will become larger one
addQuot <- function(score_hA, score_hB){
  # my scores are better if they're small so 1 is if we
  # decide B and -1 if we decide A
  decision <- c(-1,1)[(score_hA > score_hB)*1+1] 
  rnk <- 1/(1+min(score_hA, score_hB))*decision
  return(rnk)
}
new_test_scores_pos <- c(addQuot(test_scores_pos[2],test_scores_pos[1]), addQuot(test_scores_pos[1],test_scores_pos[2]))
which.min(test_scores_pos)
which.max(new_test_scores_pos)


differ <- function(score_hA, score_hB){
  res <- NA
  # if both Inf or -Inf then diff is considered 0
  if(!is.na(score_hA) & !is.na(score_hB)){
    res <- 0
    if(score_hA != score_hB) res <- score_hA - score_hB
  }
  return(res)
}
new_test_scores <- c(differ(test_scores[2],test_scores[1]), differ(test_scores[1],test_scores[2]))
which.min(test_scores_pos)
which.max(new_test_scores_pos)

probRD <- function(score_hA, score_hB){
  decision <- c(-1,1)[(score_hA < score_hB)*1+1]
  rnk <- max(score_hA, score_hB)*decision
  return(rnk)
}


