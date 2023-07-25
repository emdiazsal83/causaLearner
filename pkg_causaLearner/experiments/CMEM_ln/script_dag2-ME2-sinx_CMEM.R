# demo of causality package

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2-ME2-Add-sinx_cmem-none"
dataName <- strsplit(experimentName, "_")[[1]][1]

repos <- "/run/user/1000/gvfs/sftp:host=optimus.uv.es/media/disk/erc/papers/CAUSALITY/causaLearner/"
#repos <- paste("/run/user/1000/gvfs/sftp:host=", server, ",user=", user, "/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
#repos <- paste("/media/disk/erc/papers/CAUSALITY/causaLearner/", sep="")
dir(repos)
setwd(repos)


print("loading causal learners functions")
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)



#########################################################################################################################################3
# Load/simulate data on which we will evaluate performance of causal learners
#########################################################################################################################################


q <- 100


edgL <- vector("list", length=2)
names(edgL) <- c("x","y")
edgL[["x"]] <- list(edges=c("y"))
XtoY <- graphNEL(nodes=c("x","y"), edgeL=edgL, edgemode="directed")
plot(XtoY)

# define function

funcs <- list(fx=function(n) n, fy=function(x, n) sin(x)*n)

# define noise
# gaussian
nodesGauss <- list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                   y=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1))
# uniform
nodesUnif <-  list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                   y=list(dist="runif", pars=list(min=-1, max=1), a=1, b=1))
# exponential
nodesExp <-  list(x=list(dist="rnorm", pars=list(mean=0, sd=1), a=1, b=1), 
                  y=list(dist="rexp", pars=list(rate=1), a=1, b=1))

# define  SEMs
for(nois in c("Gauss", "Unif", "Exp")){
      # nois <- "Gauss"
      sem <- list(dag=XtoY, funcs=funcs, simPars=list(n=100, nodes=eval(parse(text=paste("nodes", nois, sep="")))))
      assign(paste("sem_", nois, sep=""), sem)
}

# Perform all 3 type_funcs x 3 funcs per type_func x 3 noise per func = 27 simulations
# Now Perform all 3 type_funcs x 4 funcs per type_func x 3 noise per func = 36 simulations

# set.seed(23)
# plts <- list()
# nms <- character()
# 
# for(nois in c("Gauss", "Unif", "Exp")){
#       # nois <- "Gauss"
#       sim <- simSEMs(q, eval(parse(text=paste("sem_", nois, sep=""))))
#       name <- paste("sim_", nois, sep="")
#       print(paste("name: ", name))
#       assign(name, sim)
#       p <- plotPairsList(sim)
#       print(p)
#       plts <- c(plts, list(p))
#       nms <- c(nms, name)
# }
# 
# names(plts) <- nms
# 
# indx1 <- substr(ls(),1,4)=="sim_"
# sims <- ls()[which(indx1)]

# save(list=c(sims, "plts", "sims"), file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))

load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))

length(plts)
print(plts[["sim_Exp"]])




i <- 2
dim(do.call("rbind", eval(parse(text=sims[i]))$xs))
summary(do.call("rbind", eval(parse(text=sims[i]))$xs))
summary(do.call("rbind", eval(parse(text=sims[i]))$xs)[,2])
sum(is.na(do.call("rbind", eval(parse(text=sims[i]))$xs)[,2]))
sapply(sims, function(sim) sum(is.na(do.call("rbind", eval(parse(text=sim))$xs)[,2])))
sum(sapply(sims, function(sim) sum(is.na(do.call("rbind", eval(parse(text=sim))$xs)[,2]))))



# weights for each data set
ws <- rep(1, q)
ws <- ws/sum(ws)
names(ws) <- sim_Gauss$names

#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm", "stdr")
approxDagSetMethodList <- c("oracleME1")

#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")
#hypScorerList <- cmem_hypScorer_pack_none_1
hypScorerList <- cmem_hypScorer_pack_none_1
#hypScorerList <- cmem_hypScorer_pack_lambda_1[6]
#hypScorerList <- cmem_hypScorer_pack_kernParsX_1
#hypScorerList <- cmem_hypScorer_pack_kernParsXY_1




recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)




#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 


for(nois in c("Gauss", "Unif", "Exp")){
      # nois <- "Exp"
      print(paste("noise type: ", nois))
      scores <- applyLearnRecipesMany(recipe, dataList=eval(parse(text=paste("sim_", nois, sep=""))), numCoresScores=2, numCoresDt=2, numCoresHypSc=2, plot=TRUE)
      assign(paste("scores", nois, sep=""), scores)
}

proc.time()-pm 

#  22 mins for cmem_hypScorer_pack_none_1 with 2-2-2 core config
indx <- which(substr(ls(),1,6)=="scores" & nchar(ls())>6)
length(indx)


lapply(ls()[indx], function(el) class(unlist(eval(parse(text=el)))))
indx <- which(sapply(ls()[indx], function(el) class(unlist(eval(parse(text=el)))))=="character")
indx

lapply(scoresExp, function(el) class(unlist(el)))
indx <- which(sapply(scoresExp, function(el) class(unlist(el)))=="character")
indx 



# 2.83 hours with 3-2-1 cores, 2 data treatments, 9 kernels, 3 model types, 3 examples, 3 noise types
# (total 27 SEMs), 100 pairs per SEM, 2 nodes, 100 data points per data set

scores <- as.character()
for(nois in c("Gauss", "Unif", "Exp")){
      # nois <- "Gauss"
      print(paste("noise type: ", nois))
      scores <- c(scores, paste("scores", nois,sep=""))
}
  

all(scores %in% ls())






#########################################################################################################################################3
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################


# lets see which scores are the main problem as far as na's and nan's


for(nois in c("Gauss", "Unif", "Exp")){
      # nois <- "Gauss"
      print(paste("noise type: ", nois))
      aux <- deliverFuncToScore_recipeTree( eval(parse(text=paste("scores",nois, sep=""))), func ="tableErrors2")
      errors <- unwrapScoreFunc(scoresFunc=aux, dataList=eval(parse(text=paste("sim_",nois, sep=""))), ws)
      print(cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)))
}


# what about non decision problem -> all hypothesis scores equal

for(nois in c("Gauss", "Unif", "Exp")){
      #  nois <- "Gauss"
      print(paste("noise type: ", nois))
      aux <- deliverFuncToScore_recipeTree( eval(parse(text=paste("scores",nois, sep=""))), func ="tableNoDecision")
      errors <- unwrapScoreFunc(scoresFunc=aux, dataList=eval(parse(text=paste("sim_",nois, sep=""))), ws)
      print(cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)))
}



#########################################################################################################################################3
# Obtain and visualize performance measures 
#########################################################################################################################################


###############################################################
# Obtain dag-classification measures
###############################################################

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
n.boot <- 500

statsDB <- data.frame()


for(nois in c("Gauss", "Unif", "Exp")){
      # nois <- "Gauss"
      print(paste("noise type: ", nois))
      
      # first aggregate by min score then rank by rank2func associated to the score
      aggScores <- deliverFuncToScore_recipeTree(eval(parse(text=paste("scores",nois, sep=""))), func="aggregateScores", dags=eval(parse(text=paste("sim_", nois,"$dags", sep=""))), getDagPars=function(dag, nm) return(list(trueDag=dag)), getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")
      rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1)
      # Obtain measures and analyze
      scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=eval(parse(text=paste("sim_", nois, sep=""))), ws)
      rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=eval(parse(text=paste("sim_", nois, "$names", sep=""))), value="value")
      statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
      statsDBaux <- statsArrayToDB(statsArr)
      statsDBaux$noiseType <- nois
      statsDB <- rbind(statsDB, statsDBaux)

}

#save(list=ls(), file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, ".RData", sep=""))

load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, ".RData", sep=""))

msrss <- c("ccr")
stats <- c("orig")
dts <- c("norm")
cmplxFuncs <- c("KCDC", "KCDCrel") #, "KCDCrel"

indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats & dt %in% dts & cmplxFunc %in% cmplxFuncs))
tab <- cast(statsDB[indx,], cmplxFunc+hypSc~dt + noiseType, value="value")
tab

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats & dt %in% dts & cmplxFunc %in% cmplxFuncs))
tab <- cast(statsDB[indx,], cmplxFunc+hypSc~dt + noiseType, value="value")
tab

msrss <- c("ccr")
cmplxFuncs <- c("KCDC")
indx <- with(statsDB, which(measure %in% msrss & cmplxFunc %in% cmplxFuncs))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=hypSc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(.~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("tss")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(.~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(.~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p




