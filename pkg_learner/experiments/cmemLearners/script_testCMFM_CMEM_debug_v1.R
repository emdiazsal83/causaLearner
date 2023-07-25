# Kernel Deviance first approach

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)

hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)


# p <- 2
# nodesX <- list(dist="runif", pars=list(min=0, max=2*pi), a=1, b=1)
# nodes <- rep(list(nodesX),p)
# names(nodes) <- c("x","y")
# dag <- matrix(c(0,0,1,0),2,2)
# colnames(dag) <- rownames(dag) <- c("x","y")
# nTr <- 100
# nTe <- 110
# # random function
# set.seed(4)
# simTest_tr <- simRandSEM(p, nTr, nodes, sigma=2, sigmaErr=0, dagMat=dag)
# simTest_te <- simRandSEM(p, nTe, nodes, sigma=2, sigmaErr=0, dagMat=dag)
# #plot(getGraph(simTest_tr$dag))
# Xtr <- simTest_tr$x
# Xtr <- apply(Xtr, 2, norml)
# Xte <- simTest_te$x
# Xte <- apply(Xte, 2, norml)
# xtr <- Xtr[,"x"]
# ytr <- Xtr[,"y"]
# xte <- Xte[,"x"]
# yte <- Xte[,"y"]

#experimentName <- "dag5-ME32-Add-Mult-Cmplx_cmem-none"
#experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmem-comp-logRegDists"
#experimentName <- "dag5-ME32-Add-Mult-Cmplx_cmem_comp-4-l2yKcmc_nc"
experimentName <- "dag2-ME2-Add-Mult_cmem-comp-4-l2yKcmc_nc"
dataName <- strsplit(experimentName, "_")[[1]][1]
#dataName <- "dag2-ME2-Cmplx-SinPlus-SinTimes"

block <- 10
pm <- proc.time()
load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs
dataNm <- dat$name
dag <- dat$dag

Xtr <- dat$x
Xtr <- apply(Xtr, 2, norml)
Xte <- dat$x
Xte <- apply(Xte, 2, norml)


apply(Xtr, 2, min)
apply(Xtr, 2, max)
apply(Xte, 2, min)
apply(Xte, 2, max)


indx_X <- c(1,2)[match(paste(as.numeric(dat$dag), collapse="."),  
                      c(paste(c(0,0,1,0), collapse="."),paste(c(0,1,0,0), collapse=".")))]
indx_Y <- setdiff(c(1,2), indx_X)

xtr <- Xtr[,indx_X]
ytr <- Xtr[,indx_Y]
xte <- Xte[,indx_X]
yte <- Xte[,indx_Y]


plot(xtr,ytr)
plot(ytr, xtr)

trainDataXY <- constructData(as.matrix(xtr), ytr)
trainDataYX <- constructData(as.matrix(ytr), xtr)

testDataXY <- constructData(as.matrix(xte), yte)
testDataYX <- constructData(as.matrix(yte), xte)


# cmfm_learner_pack_none_1
# cmfm_learner_pack_lambda_1
# cmfm_learner_pack_kernParsX_1
# cmfm_learner_pack_kernParsXY_1

hs_cmem_ob_version <- "v6_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)


hyp_scorers_cmem <- c(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1,
                      #cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)
                      #cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)

hyp_scorers_cmfm <- c(
  cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1,
#cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1,
cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)
#cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)

length(hyp_scorers_cmem)
length(hyp_scorers_cmfm)

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

cmfm_learners <- lapply(hyp_scorers_cmfm, function(el){
  hypScorer <- eval(parse(text=el))
  cmfm_learner1 <- hypScorer$cmemLearner1
  cmfm_learner2 <- hypScorer$cmemLearner2
  #print("***************")
  #print(cmfm_learner1)
  #print(cmfm_learner2)
  #return(NULL)
  return(list(cmfm_learner1=cmfm_learner1, cmfm_learner2=cmfm_learner2))
})

cmem_learners1 <- sapply(cmem_learners, function(el) el$cmem_learner1)
cmem_learners2 <- sapply(cmem_learners, function(el) el$cmem_learner2)


cmfm_learners1 <- sapply(cmfm_learners, function(el) el$cmfm_learner1)
cmfm_learners2 <- sapply(cmfm_learners, function(el) el$cmfm_learner2)

indx_lrn <- 9
cmem_learner1 <- cmem_learners1[indx_lrn]
cmem_learner2 <- cmem_learners2[indx_lrn]

eval(parse(text=cmem_learner1))$hyperParams$data$non_optimizable$kernelY

hypArray <- oracleME(trueDAG=dag)

expType <- toupper(strsplit(strsplit(experimentName, "_")[[1]][2], "-")[[1]][1])
#dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, sep=""))
#experimentName <- "dag5-ME32-Add-Mult-Cmplx_cmem_comp-4-l2yKcmc_nc"
#experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmem-comp-4-l2yKcmc_nc"
experimentName <- "dag2-ME2-Add-Mult_cmem-comp-4-l2yKcmc_nc"

#experimentName <- "dummy"
folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, "/",sep="")

dir(folderLearners)


which(apply(hypArray, 3, function(mat) all(dat$dag==mat)))

plotCurveUAI_master(x=Xtr, hypArray, cmemLearner1=cmem_learner1, cmemLearner2 =cmem_learner2 ,
                    ppTab=ppTabCMEM1, plot=FALSE, dataNm=dataNm, folderSave=folderLearners
                    , x0=-2, yN=1.3, xN=-0.3, y0=1.03)
#x=Xtr; cmemLearner1=cmem_learner1;cmemLearner2 =cmem_learner2;ppTab=ppTabCMEM1;plot=FALSE;folderSave=folderLearners;noiseLearner=NULL;augmentData=FALSE
scores <- cmem_hypScorer_comp(x=Xtr, hypArray, cmemLearner1=cmem_learner1, cmemLearner2 =cmem_learner2 , 
                              ppTab=ppTabCMEM1, plot=FALSE, dataNm=dataNm, folderSave=folderLearners)

# 128 for cmem (27 grid points), 278 for cmfm (27 grid points)
#  for cmem (108 grid points),  20.4 mins for cmfm (27*4=108 grid points)


learners1 <- c(cmem_learners1) #cmem_learners1, 
learners2 <- c(cmem_learners2) #cmem_learners2, 

count <<- 0
res <- mapply( function(cmfm_learner1, cmfm_learner2){
  print("*********************")
  # i <- 4; cmfm_learner1 <- learners1[i]; cmfm_learner2 <- learners2[i] 
  count <<- count + 1
  print("*************************")
  print(count)
  print(cmfm_learner1)
  print(cmfm_learner2)
  pm <- proc.time()
  scores <- cmem_hypScorer_comp(x=Xtr, hypArray, cmemLearner1=cmfm_learner1, cmemLearner2 =cmfm_learner2 , 
                              ppTab=ppTabCMFM1, plot=FALSE, dataNm=dataNm, folderSave=folderLearners)
  pm <- proc.time() - pm 
  return(list(pm=pm, scores=scores))
}, 
cmfm_learner1=learners1[(count+1):length(learners1)], 
cmfm_learner2=learners2[(count+1):length(learners2)], SIMPLIFY=FALSE)

t(sapply(res, function(el) el$pm))

t(sapply(res, function(el) el$scores))

dir("./pkg_learner/experiments/cmemLearners/")
reposSave <- "./pkg_learner/experiments/cmemLearners/"
file <- "times2by2by2Grid.RData"
save(pm, file=paste(reposSave, file, sep=""))
