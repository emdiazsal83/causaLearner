# demo of causality package

remove(list=ls())
#server <- "optimus.uv.es"
server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag5-ME32-Add-Mult-Cmplx_cmem-l2yKcmc_nc"
# folderSave <- paste("/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/CMEM/results/learners/", experimentName, "/",sep="")
# dir(folderSave)
dataName <- strsplit(experimentName, "_")[[1]][1]
hs_cmem_ob_version <- "v5_comp_ord"
hs_cmfm_ob_version <- "v5_comp"

#repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
getwd()
#repos <- paste("/media/disk/users/emiliano/causaLearner/", sep="")
#folderSave <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner/pkg_causaLearner/experiments/blocks/learners/"
repos <- paste("/home/emiliano/causaLearner/", sep="")
folderSave <- "/home/emiliano/causaLearner/pkg_causaLearner/experiments/blocks/learners/"


dir(repos)
setwd(repos)
print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

library(FNN) #  (in learn.qhsic, myGpCreate)
# FNN sometimes has a problem: cant find C functions perhaps because of some
# unloading of dll libraries: see https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16446



#########################################################################################################################################3
# Load/simulate data on which we will evaluate performance of causal learners
#########################################################################################################################################

count.graphs(type="all-dags", 5)
allDags5 <- genAllDAGS(5)
# pm <- proc.time()
# sizeME <- apply(allDags5, 3, function(dagMat) dim(getMarkovEquivClass(dagMat))[3])
# proc.time() - pm # 5 mins
# save(sizeME, file="./pkg_causaLearner/experiments/CMEM/results/sizeME_allDags5.RData")
load(file="./pkg_causaLearner/experiments/CMEM/results/sizeME_allDags5.RData")
table(sizeME)
indxME32 <- which(sizeME == 32)[1]

q <- 300
set.seed(4)
p <- 5
ps <- rep(p, q)
ns <- rep(100, q)
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)
dagMat <- allDags5[,,indxME32]
rownames(dagMat) <- colnames(dagMat) <- names(nodes)
dagMatME <- getMarkovEquivClass(dagMat)
dim(dagMatME)

apply(dagMatME, 3, sum)
apply(dagMatME, c(2,3), sum)

apply(dagMatME, c(2,3), function(from) paste(which(from==1), collapse="."))


dagMatME[,,1]
dagMatME[,,2]

permss <- perms(c(1,2,3,4,5))
dagMatME2 <- sapply(1:nrow(permss), function(i){
  # i <- 4
  res <- dagMatME[,,1]
  res <- res[permss[i,],]
  res <- res[,permss[i,]]
  return(res)
}, simplify="array")

# defAttrs <- getDefaultAttrs()
# defAttrs$node$fontsize <- "8"
# plot(as(dagMat, "graphNEL"), attrs=defAttrs)
# defAttrs$node$fontsize <- "6"
# par(mfrow=c(6,6))
# for(i in 1:dim(dagMatME)[3]) plot(as(dagMatME[,,i], "graphNEL"), attrs=defAttrs)
# plot(as(dagMatME[,,2], "graphNEL"))


 seed <- 5
 # set.seed(seed)
 # dataList1 <- simRandSEMs(q, ps, ns, nodess, sigma=1, dagMat=dagMat, markovEquiv=T)
 # pairs(dataList1$xs[[3]])
 # ggpairs(as.data.frame(dataList1$xs[[13]]),aes(alpha = 0.4))
 # dataList2 <- simRandSEMs(q, ps, ns, nodess, sigma=1, dagMat=dagMat, markovEquiv=T, geU=function(y, nois, scale, constant) y*nois)
 # pairs(dataList2$xs[[3]])
 # ggpairs(as.data.frame(dataList2$xs[[13]]),aes(alpha = 0.4))
 # dataList3 <- simRandSEMs(q, ps, ns, nodess, sigma=1, sigmaErr=1, dagMat=dagMat, markovEquiv=T, geU=function(y, nois, scale, constant) y)
 # pairs(dataList3$xs[[3]])
 # ggpairs(as.data.frame(dataList3$xs[[13]]),aes(alpha = 0.4))
 # names(dataList1); names(dataList2); names(dataList3)
 # dataList1$names
 # dataListList <- list(Add=dataList1, Mult=dataList2, Cmplx=dataList3)
 # dataList <- dataJoin(dataListList)
 # save(dataList, file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))



block <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#block <- 1

pm <- proc.time()
load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))
dat <- list(dag=dataList$dags[[block]], x=dataList$xs[[block]], nois=dataList$noiss[[block]], name=dataList$names[[block]])
proc.time() - pm # 0.058 secs

indxReady <- c(1:30, 101:104,106:113,115:116,118,120:130, 201:250)
length(indxReady)
names(dataList)
dataList$dags <- dataList$dags[indxReady]
dataList$xs <- dataList$xs[indxReady]
dataList$noiss <- dataList$noiss[indxReady]
dataList$names <- dataList$names[indxReady]

# weights for each data set
ws <- rep(1, length(dataList$names))
ws <- ws/sum(ws)
names(ws) <- dataList$names


#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm") #, "stdr"
approxDagSetMethodList <- c("oracleME1")
#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")


hyp_scorers_cmem <- c(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1,
                      #cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1,
                      cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)
#cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)


hypScorerList <- hyp_scorers_cmem

hypScorerList <- cmem_hypScorer_pack_none_1[2]
recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)

recipe <- recipe[-c(4,7),]

recipe <- recipe[c(5,8),]

dim(recipe)



#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 
scores <- applyLearnRecipesOne(recipe, data=dat, numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderSave)
proc.time()-pm 

# 

# seems like its best to parallelize on applyLearnRecipes level and on parameters but not on folsds

save("scores", file=paste("./pkg_causaLearner/experiments/blocks/", experimentName,"_",hs_cmem_ob_version, "_", block, ".RData", sep=""))

print("finish")
q("no")

experimentName2 <- "dag5-ME32-Add-Mult-Cmplx_cmem"
blockFiles <- dir("./pkg_causaLearner/experiments/blocks/")
indx <- grep(experimentName2, blockFiles)
blockFiles <- blockFiles[indx]
blockNumber <- as.numeric(sapply(strsplit(sapply(strsplit(blockFiles, "_"), function(el) el[[6]]), "\\."), function(el) el[[1]]))
(blockFiles <- blockFiles[order(blockNumber)])

scores <- lapply(blockFiles, function(el){
  # el <- blockFiles[1]
  load(file=paste("./pkg_causaLearner/experiments/blocks/", el, sep=""))
  res <- scores
})
names(scores) <- dataList$names
#scores <- scores[-c(4,5)]

scores[c(4,5)]

# For running all on my pc

pm <- proc.time()
scores <- applyLearnRecipesMany(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresHypSc=1, plot=FALSE, folderSave=folderSave)
proc.time() - pm 

# 12.4 mins for none_1      with 3-1-2 cores on my pc with 1-1 in CV
# for kernParsX_1 with 1-1-6 cores on my pc with 1-1 in CV


version <- "v_REG" #"v_none"#"vNCE_CLASS"
save("scores", file=paste("./pkg_causaLearner/experiments/CMEM/results/", 
                          experimentName, "_", hs_cmem_ob_version, "_",
                          version  ,".RData", sep=""))


unlink(paste("./pkg_causaLearner/experiments/blocks/",blockFiles, sep=""))





#########################################################################################################################################*
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################*
lambda <- log(eval(parse(text=eval(parse(text=recipe$hypScorer[1]))$cmemLearner))$hyperParams$data$non_optimizable$lambda$val,10)
load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName,"_",hs_cmem_ob_version,"_lambda_1e",lambda,"_seed_",seed, ".RData", sep=""))
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)


i <- 2
names(scores[[i]][["norm"]][["oracleME1"]])
names(scores[[i]][["norm"]][["oracleME1"]][-26])
summary(scores[[i]][["norm"]][["oracleME1"]][["hs_cmem_kern_vanANDkern_van_cc_L2_none"]])

aux <- sapply(1:25, function(j) scores[[i]][["norm"]][["oracleME1"]][[j]], simplify="array")

aux <- sapply(1:100, function(i) sapply(1:25, function(j) scores[[i]][["norm"]][["oracleME1"]][[j]], simplify="array"), simplify="array")
dim(aux)
summary(as.numeric(aux[,"KCMC",,]))


# clean all scores
# first lets see how many Inf, -Inf, NA or NULL or NaN
aux <- deliverFuncToScore_recipeTree(scores, func ="printErrors", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))

# actually we only have a problem when there is a NA, nan, or both all hypothesis are non-finite with the same sign 
# lets check first nas and nans

aux <- deliverFuncToScore_recipeTree(scores, func ="printErrors2", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))



# lets see which scores are the main problem as far as na's and nan's

aux <- deliverFuncToScore_recipeTree(scores, func ="tableErrors2")

errors <- unwrapScoreFunc(scoresFunc=aux, dataList=dataList, ws)

cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)/length(x)*100/2)

# what about non decision problem -> all hypothesis scores equal

aux <- deliverFuncToScore_recipeTree(scores, func ="printNoDecision", getDagPars = function(dag, dataName) return(list(dataName=dataName)), 
                                     getDataLevelPars = function(scoreMat, dag, reci) return(list(reci=reci)))

aux <- deliverFuncToScore_recipeTree(scores, func ="tableNoDecision")

noDecisions <- unwrapScoreFunc(scoresFunc=aux, dataList=dataList, ws)


cast(noDecisions, score~. ,value="value", fun.aggregate=function(x) sum(x, na.rm=T)/2)



#########################################################################################################################################3
# Obtain and visualize performance measures 
#########################################################################################################################################


###############################################################
# Obtain dag-classification measures
###############################################################
# Mult5- gaussloglik-L2&KCMC
# Mult32-negCE-KCMC

tmp <- scores[[105]]

tmp[[1]][[1]][c(3,5)]
tmp[[1]][[1]][[6]]
tmp[[1]][[1]][[3]] <- tmp[[1]][[1]][[6]]
tmp[[1]][[1]][[5]] <- tmp[[1]][[1]][[6]]
scores[[105]] <- tmp 

msrs <- deliverFuncToScore_recipeTree(scores, func="measuresDagDistEdge", dags=dataList$dags, ppTab=1, hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1],
                                      getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                      getDataLevelPars=getHypots)


# analyze these measures

measures <- unwrapScoreFunc(scoresFunc=msrs, dataList=dataList, ws)

version <- "v_REG2" #"v_none2"#"vNCE_CLASS_2"
save(list=c("scores","measures"), file=paste("./pkg_causaLearner/experiments/CMEM/results/", 
                          experimentName, "_", hs_cmem_ob_version, "_",
                          version  ,".RData", sep=""))


measures$dataset2 <- substr(measures$dataset,1,1)
table(measures$dataset2)

table(measures$measure)
# Total data sets per recipe and p
msrss <- c("ccr") #"edgeD","fnr","fpr","msr","nonEdgeD","npp", "npp", "ppp","sens", "spec", "totEdgeD", "tss"
indx <- with(measures, which(measure %in% msrss  ))
cast(measures[indx,], recipeFull+dataset2~measure, value="value", margins="grand_col")
# mean stats
tab <- cast(measures[indx,], dataset2+recipe~score, value="value", fun.aggregate="mean", margins="grand_col")
tab

writeRepos <- "./pkg_causaLearner/experiments/CMEM/results/"
write.csv(tab, file=paste(writeRepos, experimentName, "_edgeCCR.csv", sep=""))


msrss <- c("SHD") #"edgeD","fnr","fpr","msr","nonEdgeD","npp", "npp", "ppp","sens", "spec", "totEdgeD", "tss"
indx <- with(measures, which(measure %in% msrss  ))
cast(measures[indx,], recipeFull+dataset2~measure, value="value", margins="grand_col")
# mean stats
tab <- cast(measures[indx,], dataset2+recipe~score, value="value", fun.aggregate="mean", margins="grand_col")
tab

writeRepos <- "./pkg_causaLearner/experiments/CMEM/results/"
write.csv(tab, file=paste(writeRepos, experimentName, "_SHD.csv", sep=""))



msrss <- c("SHD")
indx <- with(measures, which(measure %in% msrss ))
p <- ggplot(measures[indx,])
p <- p + geom_boxplot(aes(x=recipeFull, y=value, color=recipeFull))
p <- p + facet_wrap(~measure, scales="free")
p


hypScs <- c("score_rf","score_lasso")
indx2 <- with(measures, which(measure %in% msrss  & dt %in%dts & hypSc %in% hypScs))
cast(measures[indx2,], recipeFull~measure, value="value", margins="grand_col")
# mean stats
tab <- cast(measures[indx2,], recipeFull~measure, value="value", fun.aggregate="mean", margins="grand_col")
tab


#############################################
# first choose parameters for aggregating 
# scores of sets of dags into two classes
#############################################


groupFuncs <- c("agg_truevsRest")
groupFuncPars <- list(list())

#############################################
# Aggregate scores according to previous 
# parameters
#############################################

# first aggregate by min score then rank by rank2func associated to the score

aggScores <- deliverFuncToScore_recipeTree(scores, func="aggregateScores", dags=dataList$dags, 
                                           getDagPars=function(dag, nm) return(list(trueDag=dag)), 
                                           getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")

rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1,
                                         hypSc_nameProxy=names(scores[[1]][["norm"]][["oracleME1"]])[1])





#############################################
# Obtain measures and analyze
#############################################

scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList, ws)
rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=dataList$names, value="value")



scoreDB2 <- unwrapScoreFunc(scoresFunc=aggScores, dataList, ws)
scoreDB2$posNeg <- as.numeric(scoreDB2$hypothesis)

tab <- cast(scoreDB2, dataset+posNeg~recipeFull, fun.aggregate="length")
tab[1:3, 1:5]

regDB <- cast(scoreDB2, dataset+posNeg~recipeFull, value="value")
dim(regDB)
colnames(regDB)[3:ncol(regDB)] <- sapply(strsplit(colnames(regDB)[3:ncol(regDB)], "\\."), 
                                         function(el) paste(paste( strsplit(strsplit(el[3],"_")[[1]][4], "AND")[[1]][1], 
                                                                  strsplit(el[3],"_")[[1]][5], sep="_"), el[4], sep="_"))
regDB$posNeg <- as.factor(regDB$posNeg)
regDB[1:3, 1:5]
summary(regDB)[,1:10]

library(glmnet)
lambda_seq <- 10^seq(2, -8, by = -.1)
X <- as.matrix(as.matrix(regDB[,3:ncol(regDB)]))
colnames(X) <- colnames(regDB)[3:ncol(regDB)]
cv_output <- cv.glmnet(x=X, y=regDB$posNeg, alpha = 1, lambda = lambda_seq, family="binomial")
best_lam <- cv_output$lambda.min
lasso_best <- glmnet(x=X, y=regDB$posNeg, alpha = 1, lambda = best_lam, family="binomial")

names(lasso_best)
as.matrix(lasso_best$beta)
rownames(as.matrix(lasso_best$beta))[which(lasso_best$beta>1e-3)]
predict(lasso_best, newx=X, type="response")

table((predict(lasso_best, newx=X, type="response")[,1]>0.5)*1, regDB$posNeg)

save(lasso_best, file="./pkg_causaLearner/experiments/CMEM/results/lassoModel.RData")


aux <- strsplit(colnames(regDB), "_")
aux1 <- sapply(aux, function(el) el[3])
aux2 <- sapply(aux, function(el) el[3])
indxVars <- which(aux %in% c("KCMC","KCRDC"))

library(randomForest)
rf <- randomForest(x=regDB[,3:ncol(regDB)], y=regDB$posNeg, ntree=50000, mtry=floor(sqrt(length(3:(ncol(regDB)-2)))), nodesize=1)
print(rf)
importance(rf)
varImpPlot(rf)

predict(rf, newdata=regDB, type="prob")[,1] # coz lower score is best

save(rf, file="./pkg_causaLearner/experiments/CMEM/results/rfModel.RData")

n.boot <- 500

statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
statsDB <- statsArrayToDB(statsArr)


msrss <- c("ccr")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value", fun.aggregate="length")
tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")

tab[order(tab$KCDC, decreasing=T),]
tab[order(tab$KCMC, decreasing=T),]
tab[order(tab$KCNSC, decreasing=T),]
tab[order(tab$KCRDC, decreasing=T),]
tab[order(tab$KCSC, decreasing=T),]
tab[order(tab$majority_score, decreasing=T),]
tab[order(tab$rf_score, decreasing=T),]
tab[order(tab$lasso_score, decreasing=T),]


msrss <- c("auc")
stats <- c("orig")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats))
tab <- cast(statsDB[indx,], dt+hypSc~cmplxFunc, value="value")
tab

tab[order(tab$KCDC, decreasing=T),]
tab[order(tab$KCMC, decreasing=T),]
tab[order(tab$KCNSC, decreasing=T),]
tab[order(tab$KCRDC, decreasing=T),]
tab[order(tab$KCSC, decreasing=T),]

msrss <- c("ccr")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("tss")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(dt~hypSc, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p




