remove(list=ls())
# Test KCDC regularized multi-output regression

lossKCDC <- function(alphaVec, Lx, phiy, lambda, n){
  alpha <- matrix(alphaVec, n, ncol(phiy))
  aux <- Lx%*%alpha%*%t(alpha)%*%Lx
  res1 <- -2*sum(diag(Lx%*%alpha%*%t(phiy)))
  res2 <- sum(diag(aux)) 
  
  res3 <- n*lambda*sum(diag(aux*aux)) 
  #res3 <- lambda*n*diag(aux*aux)[3]
  #i<- 2; res3 <- n*diag(aux*aux)[i]
  
  res4 <- -lambda*sum(diag(aux))^2
  
  return(res1+res2+res3+res4)
} 

gradKCDC <- function(alphaVec, Lx, phiy, lambda, n){
  alpha <- matrix(alphaVec, n, ncol(phiy))
  aux <- Lx%*%alpha%*%t(alpha)%*%Lx
  
  #aux2 <- do.call(cbind, lapply(1:n, function(i) Lx[,i, drop=F] %*% Lx[i,, drop=F]))
  
  
  
  In <- diag(n)
  
  res1 <- -2*Lx%*%phiy # right!
  res2 <- 2*Lx%*%Lx%*%alpha #+ right! 
  
  #res3 <- 4*lambda*n*aux2%*%((alpha%*%t(alpha))%x%In)%*%t(aux2) # wrong! 
  #res3 <- 4*lambda*n*(Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha wrong!
  
  # i <- 2
  # Di <- rep(0, n)
  # Di[i] <- 1
  # Di <- diag(Di)
  # res3 <- 4*n*(Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha #right-ish!
  
  pm <- proc.time()
  res3b <- sapply(1:n, function(i){
    Di <- rep(0, n)
    Di[i] <- 1
    Di <- diag(Di)
    res <- (Lx%*%Di%*%Lx)%*%alpha%*%t(alpha)%*%(Lx%*%Di%*%Lx)%*%alpha
  }, simplify="array")
  proc.time() - pm #0.013
  
  pm <- proc.time()
  res3 <- sapply(1:n, function(i){
    
    res <- as.numeric(Lx[i,]%*%alpha%*%t(alpha)%*%Lx[,i])* (Lx[,i,drop=F]%*%Lx[i,,drop=F])%*%alpha
  }, simplify="array")
  proc.time() - pm # 0.011
  
  plot(res3, res3b); abline(a=0, b=1, col="red")
  
  res3 <- 4*n*lambda*apply(res3, c(1,2), sum)
  res3b <- 4*n*lambda*apply(res3b, c(1,2), sum)
  
  plot(res3, res3b); abline(a=0, b=1, col="red")
  
  res4 <- -4*lambda*sum(diag(aux))*Lx%*%Lx%*%alpha #right !
  
  res <- res1+res2+res3+res4
  
  return(as.numeric(res))
}


library(numDeriv)

 n <- 3
 m <- 4
 alpha0 <- matrix(seq(m*n),n,m)
 Lx <- matrix(rnorm(n*n),n,n)
 Lx <- Lx+t(Lx)
 
 phiy <- matrix(runif(m*n),n,m)
 lambda <- 0.1
 
 lossKCDC(as.numeric(alpha0), Lx, phiy, lambda, n)
 length(gradKCDC(alphaVec=as.numeric(alpha0), Lx, phiy, lambda, n))
 
 
 pm <- proc.time()
 Jnum <- jacobian(lossKCDC, as.numeric(alpha0), Lx=Lx, phiy=phiy, lambda=lambda, n=n)
 proc.time() - pm # 2 mins
 Jact <- gradKCDC(alphaVec=as.numeric(alpha0), Lx=Lx, phiy=phiy, lambda=lambda, n=n)
 length(Jact); dim(Jnum)
 plot(Jact, Jnum); abline(a=0, b=1, col="red")
 