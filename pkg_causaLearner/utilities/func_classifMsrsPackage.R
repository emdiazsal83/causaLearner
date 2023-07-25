# contingency table based measures

contingencyTable <- function(pred, obs, ws=rep(1, length(pred))){
  
  
  ws <- ws/sum(ws)*length(pred)
  tab <- table(factor(pred, levels=c(0,1)), factor(obs, levels=c(0,1)))
  rsums <- apply(tab, 2, sum)
  tab <- rbind(tab, rsums)
  csums <- apply(tab, 1, sum)
  tab <- cbind(tab, csums)
  res <- as.numeric(tab)
  names(res) <- c("tn","fp","f","fn","tp","t","n","p","N")
  
  return(res)
}

correctCR <- function(conTab) (conTab["tn"] + conTab["tp"])/conTab["N"]
misCR <- function(conTab) (conTab["fn"] + conTab["fp"])/conTab["N"]
posPP <- function(conTab){
  if(conTab["p"]==0){
    res <- 1
  } else{
    res <- conTab["tp"]/conTab["p"]
  }
  return(res)
}
negPP <- function(conTab){
  if(conTab["n"]==0){
    res <- 1
  } else{
    res <- conTab["tn"]/conTab["n"]
  }
  return(res)
}
sensitivity <- function(conTab){
  if(conTab["t"]==0){
    # if t==0, tp==0, so you got 100% of trues 
    res <- 1
  } else{
    res <- conTab["tp"]/conTab["t"]
  }
  return(res)
}
specificity <- function(conTab){
  if(conTab["f"]==0){
    # if f==0, tn==0 so you got 100% the falses
    res <- 1
  } else{
    res <- conTab["tn"]/conTab["f"]
  }
  return(res)
}
fpr <- function(conTab){
  if(conTab["f"]==0){
    res <- 0
  } else{
    res <- conTab["fp"]/conTab["f"]
  }
  return(res)
}
fnr <- function(conTab){
  if(conTab["t"]==0){
    res <- 0
  } else{
    res <- conTab["fn"]/conTab["t"]
  }
  return(res)
}
tss <- function(conTab) sensitivity(conTab) + specificity(conTab) - 1
