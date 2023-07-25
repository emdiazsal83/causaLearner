# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)


# simulation functions
source("./pkg_causaLearner/genData/func_getData_v2.R", echo=FALSE)
# learners
source("./pkg_learner/func_learners_v3.R", echo=FALSE)
source("./pkg_learner/ob_cmf_bin_learner_v1.R", echo=FALSE)
source("./pkg_learner/ob_reg_learners.R", echo=FALSE)
# for stdrize
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R", echo=FALSE)
library(reshape)
library(ggplot2)

p <- 2
nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
nodes <- rep(list(nodesX),p)
names(nodes) <- c("x","y")
dag <- matrix(c(0,0,1,0),2,2)
colnames(dag) <- rownames(dag) <- c("x","y")

n <- 100

# random function
set.seed(4)
simTest <- simRandSEM(p, n, nodes, sigma=2, sigmaErr=0, dagMat=dag)
plot(getGraph(simTest$dag))


X <- simTest$x
X <- apply(X, 2, norml)

apply(X, 2, min)
apply(X, 2, max)

x <- X[,"x"]
y <- X[,"y"]
plot(x,y)
plot(y,x)

trainDataXY <- constructData(as.matrix(x), y)
trainDataYX <- constructData(as.matrix(y), x)

#cmfm_learner_pack_none_1
#cmfm_learner_pack_kernParsXY_1
cmfm_learner <- eval(parse(text=cmfm_learner_pack_kernParsXY_1[1]))

# train hyperparameters
pm <- proc.time()
cmfm_learner_xy <- setParams(learner=cmfm_learner, trainData=trainDataXY, plot=TRUE)
cmfm_learner_yx <- setParams(learner=cmfm_learner, trainDataYX, plot=TRUE)
proc.time() - pm # 3.8 mins
# train learn parameters
cmfm_learner_xy <- cmfm_learner$learn(learner=cmfm_learner_xy)
cmfm_learner_yx <- cmfm_learner$learn(cmfm_learner_yx)
# obtain pred on mixed true/fake data to obtain loss
pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy, trainDataXY)
pred_yx <- cmfm_learner_yx$pred(learner=cmfm_learner_yx, trainDataYX)
# obtain loss
do.call(cmfm_learner_xy$optimizeParams$losses[[1]], list(learner=cmfm_learner_xy, pred=pred_xy))
do.call(cmfm_learner_yx$optimizeParams$losses[[1]], list(learner=cmfm_learner_yx, pred=pred_yx))
# calculate measures
pm <- proc.time()
cmfm_learner_xy$calcMsrs(cmemLearner=cmfm_learner_xy)
cmfm_learner_yx$calcMsrs(cmemLearner=cmfm_learner_yx)
proc.time() - pm # 0.415 secs


# look at grids

df1 <- melt(cmfm_learner_xy$hyperParams$data$optimizable$grid)
df1$dir <- c("x2y")
df2 <- melt(cmfm_learner_yx$hyperParams$data$optimizable$grid)
df2$dir <- c("y2x")
df <- rbind(df1,df2)
df <- cast(df, trainTest+sigma.rbf.Y+dir~var, value="value")

p <- ggplot(df)
p <- p + geom_line(aes(x=log(sigma.rbf.Y,10), y=negCE, colour=dir))
p <- p + geom_point(aes(x=log(sigma.rbf.Y,10), y=negCE, colour=dir))
p <- p + facet_wrap(~trainTest, scales="free")
p

# now lets set the sigma parameter of xy and yx learners
# as the direct sum of the two optimal RKHS H_y  and compare
# measures 

newSigma <- c(cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y$val,
cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.Y$val)

cmfm_learner_xy$hyperParams$data$optimizable$sigma.rbf.Y$val <- newSigma
cmfm_learner_yx$hyperParams$data$optimizable$sigma.rbf.Y$val <- newSigma

# train learn parameters
cmfm_learner_xy <- cmfm_learner$learn(learner=cmfm_learner_xy)
cmfm_learner_yx <- cmfm_learner$learn(cmfm_learner_yx)
# obtain pred on mixed true/fake data to obtain loss
pred_xy <- cmfm_learner_xy$pred(learner=cmfm_learner_xy, trainDataXY)
pred_yx <- cmfm_learner_yx$pred(learner=cmfm_learner_yx, trainDataYX)
# obtain loss
do.call(cmfm_learner_xy$optimizeParams$losses[[1]], list(learner=cmfm_learner_xy, pred=pred_xy))
do.call(cmfm_learner_yx$optimizeParams$losses[[1]], list(learner=cmfm_learner_yx, pred=pred_yx))
# calculate measures
pm <- proc.time()
cmfm_learner_xy$calcMsrs(cmemLearner=cmfm_learner_xy)
cmfm_learner_yx$calcMsrs(cmemLearner=cmfm_learner_yx)
proc.time() - pm # 0.415 secs

