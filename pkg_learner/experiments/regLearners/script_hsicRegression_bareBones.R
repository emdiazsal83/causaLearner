library(kernlab) # for rbfdot, kernelMatrix


remove(list=ls())

set.seed(3)
n <- 100
x <- runif(n, 0, 2*pi)
ny <- rnorm(n, 0, 0.1)
y <- sin(x) + ny
y <- y - mean(y)

plot(x, y)

# implement hsic regression in forward direction

# hsic regression 
# fucion de perdida  tr(Kxb %*% H %*% Krg %*%) + lambda* alpha %*% Kxs %*% alpha
# donde Krg = exp(-gamma ||r||^2)  r = y - kxs %*% alpha

sigma0 <- 1/median(as.numeric(dist(x)^2))

lambda <- 0.1


kernelXs <- rbfdot(sigma=sigma0)
kernelXb <- rbfdot(sigma=sigma0)

# calcular un modelo no parametrico (mars en este caso) para tener un estimación de la mediana de los residuos
library(earth)
model  <- earth(x=x[1:50], y=y[1:50], nfold=5)
yh0 <- as.numeric(predict(model, x[51:100]))
plot(x[51:100], yh0)
res0 <- y[51:100]-as.numeric(predict(model, x[51:100]))
gamma0r <- 1/median(as.numeric(dist(res0)^2))

kernelRg <- rbfdot(sigma=gamma0r)

# kernel correspondiente a yhat = Kxs %*% alpha
Kxs <- kernelMatrix(kernelXs, x)
# kernel de hsic correspondiente a los inputs
Kxb <- kernelMatrix(kernelXb, x)
N  <- nrow(Kxs)

# calcula la funcion de perdida con (rbf kernels) para una alfa dada
.hsicRegLoss <- function(alpha, y, Kxs, Kxb, kernelRg, lambda){
  
  rs <- y-Kxs%*%alpha
  # kernel de hsic correspondient a los residuos
  Krg <- kernelMatrix(kernelRg, rs)
  N = nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- Kxb%*%H
  Krgc <- Krg%*%H
  res <- sum(diag(Kxbc%*%Krgc)) + lambda*t(alpha)%*%Kxs%*%alpha  
  return(res)
}

# calcula el gradiente de la funcion de hsic regression loss (con rbf kernel) para una da alfa
.hsicRegLossGrad <- function(alpha, y, Kxs, Kxb, kernelRg, lambda){
  rs <- y-Kxs%*%alpha
  Krg <- kernelMatrix(kernelRg, rs)
  N <- nrow(Kxs)
  H <- diag(N)-matrix(1/N,N,N)
  Kxbc <- H%*%Kxb%*%H
  ones <- matrix(1,N,1)
  aux <- Kxbc*Krg*(ones%*%t(rs)-rs%*%t(ones))
  gamma <- kernelRg@kpar$sigma
  
  grad <- apply(diag(N),2, function(ek)t(ones)%*%(aux*(ones%*%t(ek)%*%Kxs-Kxs%*%ek%*%t(ones)))%*%ones)
  
  
  
  grad <- 2*gamma*grad  + 2*lambda*Kxs%*%alpha
  return(grad)
}


set.seed(1)
alpha0 <- rep(0,N)#rnorm(N)

library(Matrix) #Matrix function
alpha_krr <-  as.numeric(solve(Matrix(Kxs + diag(lambda, N))) %*% y)

# checar solución krr para funcion de perdida hsic y norma del gradiente
.hsicRegLoss(alpha_krr, y, Kxs, Kxb, kernelRg, lambda)
sqrt(sum((.hsicRegLossGrad(alpha_krr, y, Kxs, Kxb, kernelRg, lambda))^2))

library(lbfgs) #lbfgs function
pm <- proc.time()
res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha_krr, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=0, max_iterations=3000)
proc.time() - pm # 4.3 mins para 3000 iteraciones


alpha <- res$par

.hsicRegLoss(alpha, y, Kxs, Kxb, kernelRg, lambda)
sqrt(sum((.hsicRegLossGrad(alpha, y, Kxs, Kxb, kernelRg, lambda))^2))



xx <- x
yy <- Kxs%*%alpha
yy <- yy - mean(yy)

yy.krr <- Kxs%*%alpha_krr

o <- order(x)
xx <- xx[o]
yy <- yy[o]
yy.krr <- yy.krr[o]


plot(x, y, ylim=range(y,yy.krr, yy))
lines(xx, yy.krr, col="red")
lines(xx, yy, col="green")

plot(alpha, alpha_krr) # solución algo parecida a la de alpha_krr, normal porque inicializamos en esa solución



#varias (m) inicializaciones en paralelo y escoger la mejor
m <- 1000
set.seed(2)
ALPHA0 <- matrix(rnorm(m*N), N, m) 
# que tengan misma norma 2 que la de alpha_krr las m inicializaciones
ALPHA0 <- apply(ALPHA0, 2, function(col){
    res <- col/as.numeric(sqrt(t(col)%*%col))*as.numeric(sqrt(t(alpha_krr)%*%alpha_krr))
})
ALPHA0 <- as.data.frame(ALPHA0)
apply(ALPHA0, 2, function(col) as.numeric(sqrt(t(col)%*%col)))
sqrt(t(alpha_krr)%*%alpha_krr)

library(parallel) #mcmapply function
 
pm <- proc.time()
fxs <- mcmapply(FUN=function(alpha0){
  res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=1, max_iterations=10)
  return(res$value)
}, alpha0=as.data.frame(ALPHA0), mc.cores=5)
proc.time() - pm # 6 mins para 1000 inicializaciones con 10 iteraciones cada una paralelizando con 5 cores 

summary(fxs)


# escoger la mejor inicialización

indx <- which.min(fxs)
alpha0 <- ALPHA0[,indx]

# optimizar con esa inicialización pero mas iteraciones
pm <- proc.time()
res <- lbfgs(.hsicRegLoss, .hsicRegLossGrad, vars=alpha0, y=y, Kx=Kxs, Kxb=Kxb, kernelRg=kernelRg, lambda=lambda, invisible=0, max_iterations=1000)
proc.time() - pm # 1.4 mins para 1000 iteraciones 
res$value
alpha <- res$par
.hsicRegLoss(alpha, y, Kxs, Kxb, kernelRg, lambda)
sqrt(sum((.hsicRegLossGrad(alpha, y, Kxs, Kxb, kernelRg, lambda))^2))

# graficar solución obtenida
xx <- x
yy <- Kxs%*%alpha
yy <- yy - mean(yy)
yy.krr <- Kxs%*%alpha_krr
o <- order(x)
xx <- xx[o]
yy <- yy[o]
yy.krr <- yy.krr[o]

plot(x, y, ylim=range(y,yy.krr, yy))
lines(xx, yy.krr, col="red")
lines(xx, yy, col="green")

plot(alpha_krr, alpha) #muy distintas soluciones!!!




