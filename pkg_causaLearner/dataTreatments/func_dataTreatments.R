# data treatment functions

# function to standarize data

stdrize <- function(x){
  res <- (x - mean(x))/sd(x)
  return(res)
}

norml <- function(x){
  res <- (x - min(x))/(max(x)-min(x)) - 0.5 # centered at 0 
  return(res)
}

gaussianize <- function(x){
  n <- length(x)
  res <- rank(x, ties.method = "random")/(n + 1)
  res <- qnorm(res)
  return(res)
}

iden <- function(x) x

vanillaNorm <- function(x) x

permData <- function(x, seed){
  RNGversion("3.5.0")
  set.seed(seed)
  smpl <- sample(1:nrow(x))
  x <- x[smpl,,drop=F]
  return(x)
}
