# generate some time series to try latent noise krr 
# in order to deal with auto correlation

hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
setwd(repos)
print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

q <- 100

ts <- 1:151
w1 <- rnorm(1, sd=2)
w2 <- rnorm(1, sd=2)
z <- (sin((w1*ts/max(ts))*2*pi) + cos((w2*ts/max(ts))*2*pi)+2)/4; summary(z)
#z <- (sin((ts/max(ts))*2*pi)+1)/2; summary(z)
plot(ts, z)

m <- 100
#set.seed(123)
wx <- runif(m,-2,2)
wz <- runif(m,-2,2)
wy <- runif(m,-2,2)
alpha_xy <- rnorm(m)#runif(m, 0, 1/m)
alpha_yx <- rnorm(m)#runif(m, 0, 1/m)

prop_xy <- 0.99
prop_yx <- 1-prop_xy
prop_zy <- 0.99
prop_zx <- 1-prop_zy
prop_self <- 0.99
sigma_self <- prop_self/(1-prop_self)  
sigma_cross_xy <- prop_xy/(1-prop_xy)
sigma_cross_yx <- prop_yx/(1-prop_yx)
sigma_hidd_xy <- prop_zy/(1-prop_zy)
sigma_hidd_yx <- prop_zx/(1-prop_zx)

x0 <- 0.1
y0 <- 0.2
z0 <- 0.3
alpha <- alpha_xy
sigx <- sigma_cross_xy
sigy <- sigma_self
sigz <- sigma_hidd_xy
funcKernel <- function(x, y, z, wx, wy, wz, sigx, sigy, sigz, alpha){
  kx <- kern_rbf(matrix(x,length(x),1), matrix(wx), sigma=sigx)
  ky <- kern_rbf(matrix(y,length(y),1), matrix(wy), sigma=sigy)
  kz <- kern_rbf(matrix(z,length(z),1), matrix(wz), sigma=sigz)
  k <- kx*ky*kz    
  res <- c(k%*%alpha)
  return(res)
}

funcKernel(x0,y0,z0,wx,wy,wz,sigx,sigy,sigz, alpha)

xt <- runif(1000)
yt <- runif(1000)
zt <- runif(1000)
# fx(x,y,z)
hist(funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_self,sigy=sigma_cross_yx, sigz=sigma_hidd_yx, alpha=alpha_yx))
plot(xt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_self,sigy=sigma_cross_yx, sigz=sigma_hidd_yx, alpha=alpha_yx))
plot(yt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_self,sigy=sigma_cross_yx, sigz=sigma_hidd_yx, alpha=alpha_yx))
plot(zt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_self,sigy=sigma_cross_yx, sigz=sigma_hidd_yx, alpha=alpha_yx))
# fy(x,y,z)
hist(funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_cross_xy,sigy=sigma_self, sigz=sigma_hidd_xy, alpha=alpha_xy))
plot(xt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_cross_xy,sigy=sigma_self, sigz=sigma_hidd_xy, alpha=alpha_xy))
plot(yt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_cross_xy,sigy=sigma_self, sigz=sigma_hidd_xy, alpha=alpha_xy))
plot(zt, funcKernel(xt,yt,zt,wx,wy,wz,sigx=sigma_cross_xy,sigy=sigma_self, sigz=sigma_hidd_xy, alpha=alpha_xy))

set.seed(123)
q <- 100
writeRepos <- "./data/timeSeries/latentNoise_xyz/"
for(i in 1:q){
  # i <- 1
  #set.seed(1234)
  print(paste("i: ",i))
  ts <- 1:151
  w1 <- rnorm(1, sd=2)
  w2 <- rnorm(1, sd=2)
  z <- (sin((w1*ts/max(ts))*2*pi) + cos((w2*ts/max(ts))*2*pi)+2)/4; summary(z)
  #z <- (sin((ts/max(ts))*2*pi)+1)/2; summary(z)
  #plot(ts, z)
  x <- runif(1,-1,1); y <- runif(1,-1,1)
  coefs_fy <- rnorm(9, sd=1)
  coefs_fx <- rnorm(9, sd=1)
  polyTerms <- poly(x=matrix(c(x,y,z[1]),1,3), degree=2, raw=T)
  yzTerms <- which(sapply(strsplit(colnames(polyTerms),"\\."), function(el) any(el[2:3] %in% c("1","2"))))
  coefs_fx[yzTerms] <- coefs_fx[yzTerms]*0.1

  for(t in ts[-1]){
    #print(paste("t: ",t))
    #x <- c(x, funcKernel(x[t-1],y[t-1],z[t-1]+rnorm(1,sd=0.1),wx,wy,wz,sigx=sigma_self,sigy=sigma_cross_yx, sigz=sigma_hidd_yx, alpha=alpha_yx)+rnorm(1,sd=0.1))
    #y <- c(y, funcKernel(x[t-1],y[t-1],z[t-1]+rnorm(1,sd=0.1),wx,wy,wz,sigx=sigma_cross_xy,sigy=sigma_self, sigz=sigma_hidd_xy, alpha=alpha_xy)+rnorm(1,sd=0.1))
    polyTerms <- poly(x=matrix(c(x[t-1],y[t-1],z[t-1]),1,3), degree=2, raw=T)
    y <- c(y, tanh(c(polyTerms %*% coefs_fy)+rnorm(1,sd=0.1)))
    x <- c(x, tanh(c(polyTerms %*% coefs_fx)+rnorm(1,sd=0.1)))
  }


  df <- data.frame(x=x[51:150], y=y[51:150], z=z[51:150])
  df <- as.ts(df)
  plot(df)
  Sys.sleep(5)
  
  modMat <- data.frame(xt=x[51:150], yt=y[51:150], x1=x[50:149], y1=y[50:149], z1=z[50:149])
  modMat <- as.matrix(modMat)
  
  write.csv(modMat, file=paste(writeRepos,"ts_","latentNoise_xyz_",i,sep=""), row.names=FALSE)
}

plot(x[1:(length(ts)-1)],y[2:length(ts)])
plot(y[1:(length(ts)-1)],x[2:length(ts)])
plot(x[1:(length(ts)-1)],x[2:length(ts)])
plot(y[1:(length(ts)-1)],y[2:length(ts)])
plot(z[1:(length(ts)-1)],y[2:length(ts)])
plot(z[1:(length(ts)-1)],x[2:length(ts)])


modMat <- data.frame(xt=x[2:length(ts)], yt=y[2:length(ts)], x1=x[1:(length(ts)-1)], y1=y[1:(length(ts)-1)], z1=z[1:(length(ts)-1)])
modMat <- as.matrix(modMat)
modMat <- apply(modMat, 2, norml)
modMat2 <- as.ts(modMat)
plot(modMat2)

# with z
# x -> y
krr2_xy <- setParams(learner=krrIndSig, trainData=constructData(x=modMat[,c("x1","y1","z1")], y=modMat[,"yt"]))
krr2_xy <- krr2_xy$learn(krr2_xy)
pred_xy <- krr2_xy$predict(krr2_xy, constructData(x=modMat[,c("x1","y1","z1")], y=modMat[,"yt"]))
plot(pred_xy$gy, pred_xy$gyh)
do.call(krr2_xy$optimizeParams$losses$rmse$func, list(learner=krr2_xy, pred=pred_xy))
# y -> x
krr2_yx <- setParams(krrIndSig, trainData=constructData(x=modMat[,c("x1","y1","z1")], y=modMat[,"xt"]))
krr2_yx <- krr2_yx$learn(krr2_yx)
pred_yx <- krr2_yx$predict(krr2_yx, constructData(x=modMat[,c("x1","y1","z1")], y=modMat[,"xt"]))
plot(pred_yx$gy, pred_yx$gyh)
do.call(krr2_yx$optimizeParams$losses$rmse$func, list(learner=krr2_yx, pred=pred_yx))
# without z
# x -> y
krr2_xy <- setParams(krrIndSig, trainData=constructData(x=modMat[,c("x1","y1")], y=modMat[,"yt"]))
krr2_xy <- krr2_xy$learn(krr2_xy)
pred_xy <- krr2_xy$predict(krr2_xy, constructData(x=modMat[,c("x1","y1")], y=modMat[,"yt"]))
plot(pred_xy$gy, pred_xy$gyh)
do.call(krr2_xy$optimizeParams$losses$rmse$func, list(learner=krr2_xy, pred=pred_xy))
# y -> x
krr2_yx <- setParams(krrIndSig, trainData=constructData(x=modMat[,c("x1","y1")], y=modMat[,"xt"]))
krr2_yx <- krr2_yx$learn(krr2_yx)
pred_yx <- krr2_yx$predict(krr2_yx, constructData(x=modMat[,c("x1","y1")], y=modMat[,"xt"]))
plot(pred_yx$gy, pred_yx$gyh)
do.call(krr2_yx$optimizeParams$losses$rmse$func, list(learner=krr2_yx, pred=pred_yx))


