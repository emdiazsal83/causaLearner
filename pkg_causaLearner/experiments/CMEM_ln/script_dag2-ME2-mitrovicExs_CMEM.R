# demo of causality package

remove(list=ls())

server <- "optimus.uv.es"
#server <- "erc.uv.es"
user <- "emiliano"

experimentName <- "dag2-ME2-mitrovicExps_cmem-none"
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

# define functions
# additive
funcsAdd_fA <- list(fx=function(n) n, fy=function(x, n) x^3 + x + n)
funcsAdd_fB <- list(fx=function(n) n, fy=function(x, n) log(x+10) + x^6 + n)
funcsAdd_fC <- list(fx=function(n) n, fy=function(x, n) sin(10*x)+exp(3*x) + n)
funcsAdd_fD <- list(fx=function(n) n, fy=function(x, n) sin(pi*x) + n)

# multiplicative
funcsMult_fA <- list(fx=function(n) n, fy=function(x, n) (x^3 + x)*n)
funcsMult_fB <- list(fx=function(n) n, fy=function(x, n) (sin(10*x)+exp(3*x))*exp(n))
funcsMult_fC <- list(fx=function(n) n, fy=function(x, n) log(x+10+x^6)*exp(n))
funcsMult_fD <- list(fx=function(n) n, fy=function(x, n) x*n)

# more complex
funcsCmplx_fA <- list(fx=function(n) n, fy=function(x, n) (log(x+10)+x^2)^n)
funcsCmplx_fB <- list(fx=function(n) n, fy=function(x, n) log(x+10)+abs(x)^(2*abs(n))) #this one allowed roots of negative values so changed it
funcsCmplx_fC <- list(fx=function(n) n, fy=function(x, n) log(x^6+5)+x^5-sin((x^2)*abs(n)))
funcsCmplx_fD <- list(fx=function(n) n, fy=function(x, n) abs(cos(2*pi*x))^(sin(2*pi*n)))

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
for(typeModel in c("Add", "Mult","Cmplx")){
  # typeModel <- "Add"
  for(exampleModel in c("A","B","C","D")){
    # exampleModel <- "D"
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Cmplx"; exampleModel <- "A"; nois <- "Gauss"
      sem <- list(dag=XtoY, funcs=eval(parse(text=paste("funcs", typeModel, "_f", exampleModel, sep=""))), simPars=list(n=100, nodes=eval(parse(text=paste("nodes", nois, sep="")))))
      assign(paste("sem", typeModel,"_f", exampleModel,"_", nois, sep=""), sem)
    }
  }
}

# Perform all 3 type_funcs x 3 funcs per type_func x 3 noise per func = 27 simulations
# Now Perform all 3 type_funcs x 4 funcs per type_func x 3 noise per func = 36 simulations

# set.seed(23)
# plts <- list()
# nms <- character()
# for(typeModel in c("Add", "Mult","Cmplx")){
#   # typeModel <- "Cmplx"
#   for(exampleModel in c("A","B","C","D")){
#     # exampleModel <- "D"
#     for(nois in c("Gauss", "Unif", "Exp")){
#       # typeModel <- "Mult"; exampleModel <- "B"; nois <- "Gauss"
#       sim <- simSEMs(q, eval(parse(text=paste("sem", typeModel, "_f", exampleModel, "_", nois, sep=""))))
#       name <- paste("sim", typeModel,"_f", exampleModel,"_", nois, sep="")
#       print(paste("name: ", name))
#       assign(name, sim)
#       p <- plotPairsList(sim)
#       print(p)
#       plts <- c(plts, list(p))
#       nms <- c(nms, name)
#     }
#   }
# }
# names(plts) <- nms
# 
# indx1 <- substr(ls(),1,3)=="sim"
# indx2 <- substr(ls(),4,6) %in% c("Add","Mul","Cmp")
# length(indx1); length(indx2)
# sum(indx1 & indx2)
# sims <- ls()[which(indx1 & indx2)]

# save(list=c(sims, "plts", "sims"), file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))

load(file=paste("./data/Simulations/", dataName, "_sims.RData", sep=""))

3*4*3
length(plts)
print(plts[["simAdd_fD_Unif"]])




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
names(ws) <- simAdd_fA_Gauss$names

#########################################################################################################################################3
# Construct recipes of causal learners to compare
#########################################################################################################################################

dataTreatmentList <- c("norm", "stdr")
approxDagSetMethodList <- c("oracleME1")

#hypScorerList <- cmem_hypScorer_pack1 #c("cmem1")
#hypScorerList <- cmem_hypScorer_pack_none_1
hypScorerList <- cmem_hypScorer_pack_none_1[c(5,6)]
#hypScorerList <- cmem_hypScorer_pack_lambda_1[6]
#hypScorerList <- cmem_hypScorer_pack_kernParsX_1
#hypScorerList <- cmem_hypScorer_pack_kernParsXY_1




recipe <- genLearnRecipes(dataTreatmentList, approxDagSetMethodList, hypScorerList)




#########################################################################################################################################3
# Obtain scores 
#########################################################################################################################################


pm <- proc.time() 

for(typeModel in c("Add", "Mult","Cmplx")){
  # typeModel <- "Cmplx"
  for(exampleModel in c("A","B","C","D")){
    # exampleModel <- "B"
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Cmplx"; exampleModel <- "B"; nois <- "Exp"
      flav <- paste(typeModel, "_f", exampleModel, "_", nois, sep="")
      print(paste("type model: ", typeModel, ", example model: ", exampleModel, " noise type: ", nois))
      print(flav)
      scores <- applyLearnRecipesMany(recipe, dataList=eval(parse(text=paste("sim", flav, sep=""))), numCoresScores=2, numCoresDt=2, numCoresHypSc=2, plot=TRUE)
      assign(paste("scores", flav, sep=""), scores)
    }
  }
}
proc.time()-pm 
# 5.2 hours or  11.6 mins per 100 pairs with 3-2-1 cores, 2  recipes, 3 model types, 3 examples, 3 noise types
# (total 27 SEMs), 100 pairs per SEM, 2 nodes, 100 data points per data set

# 2.7 hours or  6 mins per 100 pairs with 3-2-1 cores, 2  recipes, 3 model types, 3 examples, 3 noise types
# (total 27 SEMs), 100 pairs per SEM, 2 nodes, 100 data points per data set for no CV (2 data-treatments, 9 hypothesis scorers, 1 approx method)

# 1.7 hours or  3.75 mins per 100 pairs with 3-1-1 cores, 1  recipes, 3 model types, 3 examples, 3 noise types
# (total 27 SEMs), 100 pairs per SEM, 2 nodes, 100 data points per data set for no CV (1 data-treatments, 1 hypothesis scorers, 1 approx method)

# 58 mins for 100 pairs for 3x3x4 data models, with 2-2-2, 2 data types and 2 hypsCoreers for a totoal of 4 recipes

indx <- which(substr(ls(),1,6)=="scores" & nchar(ls())>6)
length(indx)


lapply(ls()[indx], function(el) class(unlist(eval(parse(text=el)))))
indx <- which(sapply(ls()[indx], function(el) class(unlist(eval(parse(text=el)))))=="character")
indx

lapply(scoresCmplx_fB_Exp, function(el) class(unlist(el)))
indx <- which(sapply(scoresCmplx_fB_Exp, function(el) class(unlist(el)))=="character")
indx 



# 2.83 hours with 3-2-1 cores, 2 data treatments, 9 kernels, 3 model types, 3 examples, 3 noise types
# (total 27 SEMs), 100 pairs per SEM, 2 nodes, 100 data points per data set

scores <- as.character()
for(typeModel in c("Add", "Mult","Cmplx")){
  for(exampleModel in c("A","B","C","D")){
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Add"; exampleModel <- "A"; nois <- "Gauss"
      flav <- paste(typeModel, "_f", exampleModel, "_", nois, sep="")
      print(paste("type model: ", typeModel, ", example model: ", exampleModel, " noise type: ", nois))
      print(flav)
      scores <- c(scores, paste("scores",flav,sep=""))
    }
  }
}

all(scores %in% ls())






#########################################################################################################################################3
# Pre- processing of scores (inc. creation of new scores by combining old ones) 
#########################################################################################################################################


# lets see which scores are the main problem as far as na's and nan's

for(typeModel in c("Add", "Mult","Cmplx")){
  for(exampleModel in c("A","B","C","D")){
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Add"; exampleModel <- "A"; nois <- "Gauss"
      flav <- paste(typeModel, "_f", exampleModel, "_", nois, sep="")
      print(paste("type model: ", typeModel, ", example model: ", exampleModel, " noise type: ", nois))
      print(flav)
      aux <- deliverFuncToScore_recipeTree( eval(parse(text=paste("scores",flav, sep=""))), func ="tableErrors2")
      errors <- unwrapScoreFunc(scoresFunc=aux, dataList=eval(parse(text=paste("sim",flav, sep=""))), ws)
      print(cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)))
    }
  }
}



# what about non decision problem -> all hypothesis scores equal

for(typeModel in c("Add", "Mult","Cmplx")){
  for(exampleModel in c("A","B","C","D")){
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Add"; exampleModel <- "A"; nois <- "Gauss"
      flav <- paste(typeModel, "_f", exampleModel, "_", nois, sep="")
      print(paste("type model: ", typeModel, ", example model: ", exampleModel, " noise type: ", nois))
      print(flav)
      aux <- deliverFuncToScore_recipeTree( eval(parse(text=paste("scores",flav, sep=""))), func ="tableNoDecision")
      errors <- unwrapScoreFunc(scoresFunc=aux, dataList=eval(parse(text=paste("sim",flav, sep=""))), ws)
      print(cast(errors, score~. ,value="value", fun.aggregate=function(x) sum(x)))
    }
  }
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

for(typeModel in c("Add", "Mult","Cmplx")){
  for(exampleModel in c("A","B","C","D")){
    for(nois in c("Gauss", "Unif", "Exp")){
      # typeModel <- "Add"; exampleModel <- "A"; nois <- "Gauss"
      flav <- paste(typeModel, "_f", exampleModel, "_", nois, sep="")
      print(paste("type model: ", typeModel, ", example model: ", exampleModel, " noise type: ", nois))
      print(flav)
      
      # first aggregate by min score then rank by rank2func associated to the score
      aggScores <- deliverFuncToScore_recipeTree(eval(parse(text=paste("scores",flav, sep=""))), func="aggregateScores", dags=eval(parse(text=paste("sim",flav,"$dags", sep=""))), getDagPars=function(dag, nm) return(list(trueDag=dag)), getDataLevelPars=getHypots, groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc="aggMin")
      rnkDecs <- deliverFuncToScore_recipeTree(scores=aggScores, func="rankedDecisions", ppTab=1)
      # Obtain measures and analyze
      scoreDB <- unwrapScoreFunc(scoresFunc=rnkDecs, dataList=eval(parse(text=paste("sim",flav, sep=""))), ws)
      rnkDecsArr <- scoreDBToArray(scoreDB, names_datasets=eval(parse(text=paste("sim",flav, "$names", sep=""))), value="value")
      statsArr <- measures2DagClass(rankDecArray=rnkDecsArr, ws, n.boot)
      statsDBaux <- statsArrayToDB(statsArr)
      statsDBaux$typeModel <- typeModel
      statsDBaux$exampleModel <- exampleModel
      statsDBaux$noiseType <- nois
      statsDB <- rbind(statsDB, statsDBaux)
    }
  }
}

#save(list=ls(), file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, ".RData", sep=""))

load(file=paste("./pkg_causaLearner/experiments/CMEM/results/", experimentName, ".RData", sep=""))

msrss <- c("ccr")
stats <- c("orig")
dts <- c("norm")
cmplxFuncs <- c("KCDC", "KCDCrel") #, "KCDCrel"

for(typeMdl in c("Add","Mult","Cmplx")){
  # typeMdl <- "Add"
  
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats & dt %in% dts & cmplxFunc %in% cmplxFuncs & typeModel %in% typeMdl ))
  print(length(indx))
  tab <- cast(statsDB[indx,], typeModel+exampleModel+noiseType~hypSc+cmplxFunc, value="value")
  print(tab)
}

msrss <- c("auc")
for(typeMdl in c("Add","Mult","Cmplx")){
  # typeMdl <- "Add"
  
  indx <- with(statsDB, which(measure %in% msrss & statistic %in% stats & dt %in% dts & cmplxFunc %in% cmplxFuncs & typeModel %in% typeMdl ))
  print(length(indx))
  tab <- cast(statsDB[indx,], typeModel+exampleModel+noiseType~hypSc+cmplxFunc, value="value")
  print(tab)
}


msrss <- c("ccr")
cmplxFuncs <- c("KCDC")
indx <- with(statsDB, which(measure %in% msrss & cmplxFunc %in% cmplxFuncs))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=hypSc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(typeModel+exampleModel~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("tss")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(typeModel+exampleModel~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

msrss <- c("auc")
indx <- with(statsDB, which(measure %in% msrss))
p <- ggplot(statsDB[indx,])
p <- p + geom_jitter(aes(x=cmplxFunc, y=value, color=statistic), width=0.1)
p <- p + facet_grid(typeModel+exampleModel~noiseType, scales="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p




