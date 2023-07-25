# remove(list=ls())
# repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
# setwd(repos)
# hs_cmem_ob_version <- "v3_comp"
# hs_cmfm_ob_version <- "v5_comp"

library(graph) # graphNEL (bioconductor)
library(bnlearn) # count.graphs
library(gtools) # permutations, combinations


library(abind) #bind arrays
library(ggplot2)
library(GGally) #ggpairs
library(reshape)
library(parallel) #mcmapply, mlapply



library(pROC) # roc, auc
library(WeightedROC) # WeightedROC, WeightedAUC
library(boot) # boot

library(oro.dicom) # dec2base in ob_cme_learner_CV.R

print("loading data generating functions")
# Functions to load or generate data on which to test causality package
source("./pkg_causaLearner/genData/func_getData_v2.R")
# Utility dag functions
print("loading dag utility functions")
source("./pkg_causaLearner/utilities/func_dagStuff.R")
# Functions for contingency table based classification measures
print("loading contingency table based classification measures")
source("./pkg_causaLearner/utilities/func_classifMsrsPackage.R")


# Causal Learners

# A causal learner takes in multivariate-data and gives back a set of scored dags

# A causal learner consists of a data-treatment, an approximation dag set method and a hypothesis scorer
# 1. data treatment transforms data (eg. standardize)
# 2. approximation dag set method takes data and gives back a set of possible dags. (eg, PC algo, minimal dag)
# 3. hypothesis scorer takes in a set of dags and gives back a score or scores for each one

##########################################################################################################*
# Hypothesis Scorers
##########################################################################################################*


print("loading ANM hypothesis scorer functions")
source("./pkg_causaLearner/hypothesisScorers/pkg_ANM_hypothesisScorer/func_ANM_hypothesisScorer_v2.R")
print("loading CMEM hypothesis scorer functions")
source("./pkg_causaLearner/hypothesisScorers/cmem_hypothesisScorer/func_CMEM_hypothesisScorer_v2.R")

print("loading hypothesis scorer objects")
source("./pkg_causaLearner/hypothesisScorers/ob_hypothesisScorers_v2.R")

# random hypothesis scorer for evaluating approx dag set methods on their own

rnd_hypScorer <- function(x, hypArray, ppTab=NULL, plot=FALSE, dataNm,  folderSave){
  m <- dim(hypArray)[3]
  set.seed(34)
  scores <- matrix(rnorm(m), m, 1)
  #dagNms <- getHypID(hypArray)$id 
  dagNms <- dimnames(hypArray)$dag
  dimnames(scores) <- list(dag=dagNms, score="rnd")
  return(scores)
}

ppTabRnd <- data.frame(id=c("rnd"))
# define rank2Func transf: how t
# rnd               =                     in (-Inf,Inf) -> differ
ppTabRnd$rank2Funcs <- "differ"
ppTabRnd$rankFuncs <- "correctScoreToAdd"
ppTabRnd$probFuncs <- "scoreToProb"
ppTabRnd$argTypes <- "cmes"




##########################################################################################################*
# Data treatment functions
##########################################################################################################*

print("loading data-treatment functions")
source("./pkg_causaLearner/dataTreatments/func_dataTreatments.R")

print("loading data-treatment objects")
source("./pkg_causaLearner/dataTreatments/ob_dataTreatments.R")

###########################################################################################################*
# approx dagset methods
###########################################################################################################*

print("loading approximate dag set functions")
source("./pkg_causaLearner/approxDagSetMethods/func_approxDagSetMethods_v1.R")

print("loading approximate dag set objects")
source("./pkg_causaLearner/approxDagSetMethods/ob_approxPacks.R")

################################################################################################################*
# Functions for applying causal learners (data dreatment - approx dag set method - hypothesis scoring method)
# to lists of data 
################################################################################################################*
genLearnRecipes <- function(dataTreatmentList, approxDagSetMethodList, hypScorerList){
  # obtain N different recipes
  
  recipe <- expand.grid(dataTreatment=dataTreatmentList, approxDagSetMethod=approxDagSetMethodList,hypScorer=hypScorerList ,stringsAsFactors=F)
  recipe$id <- 1:nrow(recipe)
  
  return(recipe)
  
}

setdiffRecipes <- function(recipe, recipeSubract){
  srce <- c(rep(1, nrow(recipe)), rep(2, nrow(recipeSubtract)))
  res <- rbind(recipe, recipeSubtract)
  d1 <- duplicated(res) 
  d2 <- rev(duplicated(res[nrow(res):1,]))
  del <- which(d1 | d2 ) #| srce==2
  res <- res[-del,]
  return(res)
}

applyLearnRecipesMany <- function(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresAdsm=1, numCoresHypSc=1, plot=FALSE, folderSave){
  

  # data treatment tree - to store normalized data - same structure per data set
  # 1 level - 1 branch per unique data treatment type in recipe
  
  dts <- recipe$dataTreatment
  dts_uni <- unique(dts)
  # which data treatments will we effectively run
  dts_eff <- dts_uni
  
  dt_tree <- lapply(dts_eff, function(dt) c(dt=dt))
  names(dt_tree) <- dts_eff
  
  # approx dag set method  tree - to store hypothesis sets - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x approx dag set method in recipe
  
  adsm_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm) c(dt, adsm=adsm))
    names(res) <- adsms_uni
    return(res)
    })
  
  
  # hypothesis scorer tree - to store scores - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x hypothesis scorer in recipe
  
  hypSc_tree <- lapply(dt_tree, function(dt){ 
    hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt)]
    hypScs_uni <- unique(hypScs)
    res <- lapply(hypScs_uni, function(hypSc) c(dt, hypSc=hypSc))
    names(res) <- hypScs_uni
    return(res)
  })
  
  # recipe tree -   to store scores - same structure per data set
  # 3 levels- 1 branch per recipe (i.e. unique data treatment type x approx dag set method x hypothesis scorer)
  
  recipe_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm){
      hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt & recipe$approxDagSetMethod==adsm)]
      hypScs_uni <- unique(hypScs)
      res <- lapply(hypScs, function(hypSc) c(dt, adsm=adsm, hypSc=hypSc))
      names(res) <- hypScs_uni
      return(res)
    })
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # the adsm we run one for each leaf of adsm_tree on the data corresponding to the 1st level branch of that tree
  
  # to get the union of hypothesis, for each leaf of the hypSc_tree we look for all the leafs in the recipe tree
  # with that hypSc and take the adsm corresponding to that branch
  
  # eg - how to make union of hypothesis
  hypSc_tree.unions <- lapply(hypSc_tree, function(dt1) lapply(dt1, function(hypSc1){
    # i <- 1; j <- 1; dt1 <- hypSc_tree[[i]]; hypSc1 <- dt1[[j]]  
    # find hypSc1 in recipe tree and give back the adsm corresponding to that branch
    #print("****************") 
    #print(hypSc1)
    indx <- which(sapply(recipe_tree[[hypSc1["dt"]]], function(adsm) any(sapply(adsm, function(hypSc2){
      # k <-1; l <- 1; adsm <- recipe_tree[[hypSc1[["dt"]]]][[k]]; hypSc2 <- adsm[[l]] 
      hypSc2["hypSc"]==hypSc1["hypSc"]
      
      }))))
    res <- names(recipe_tree[[hypSc1["dt"]]])[indx]
    return(res)    
  }))
  
  # Say we have scored each hypothesis and stored in hte hypothesis scorer tree, how do put into
  # recipe tree appropriately? We traverse the recipe tree, stop at the adsm level, take the
  # hypothesis set that we want there, go up to the hypothesis level and then look for the hypothesis set
  # we just collected in the leaf of the hypothesis scorer tree corresponding to the dt-hypSc of the leaf we are in
  
  firstDag <- 1
  lastDag <- length(dataList$dags)
  
  # for each data set ...
  scores <- mcmapply(FUN=function(dag, x, n, nm) {
    # i <- 1; dag <- dataList$dags[[i]]; x <- dataList$xs[[i]]; n <- dataList$noiss[[i]]; nm <- dataList$names[i]  ;plot(getGraph(dag))
    print("*****************************************")
    print(paste("data set:  ", nm))
    if(is.character(x)){
      x <- read.csv(x, sep="", header=F)
      
      # order cause first, effect second
      x <- x[, n]
      
      colnames(x) <- colnames(dag)
      x <- as.matrix(x)
    }
    print(paste("number of samples: ", nrow(x)))
    print(head(x))
    print("true dag")
    print(dag)
    
    
    # peform data treatments
    print("performing data treatments")
    
    dt_tree.x <- lapply(dt_tree, function(dt){
      # i <- 1; dt <- dt_tree[[i]]
      #print(paste("data treatment:", dt))
      dataTreatment <- eval(parse(text=dt))
      
      if(nrow(x) > dataTreatment$maxPoints) x <- x[sample(1:nrow(x), dataTreatment$maxPoints),]
      if(dataTreatment$perm) x <- permData(x, dataTreatment$perm)
      
      x <- apply(x, 2, dataTreatment$scalingFunc)
    })
    
    # lapply(dt_tree.x, function(el) apply(el, 2, function(col) c(mu=mean(col), sigma=sd(col), low=min(col), high=max(col))))
    
    
    # get the approximate hypothesis set for each approx method x data setting : the approx method, with the specific data used
    # must be independent of the learner used. This also prevents approx methods with certain randomnmess (ex permutation based hsic used)
    # being run more than once (for diff learners) and resulting in different hypothesis sets despite having exact same parameters.
    
    print("run approximate dag set methods for each datasetting")
    
    adsm_tree.hyps <- mcmapply(FUN=function(dt){
      # i <- 1; dt <- adsm_tree[[i]]
      DAGsets <- mcmapply(FUN=function(adsm){
        # j <- 2; adsm <- dt[[j]]
        x <- dt_tree.x[[adsm["dt"]]]
      print("************************************")
      print(paste("apprx dag set method: ", adsm["adsm"]))
      apprxPack <- eval(parse(text=adsm["adsm"]))
      apprxPack$pars$dataNm <- nm
      apprxPack$pars$folderSave <- folderSave
      
      pars <- list(data=x, trueDAG=dag, pars=apprxPack$pars)
      # data <- pars$data; trueDAG <- pars$trueDAG; pars <- pars$pars
      # names(pars)
      # minDAG_cmem(data, trueDAG, pars)
      # minDAG(data, trueDAG, pars)
      # pairwise_cmem(data, trueDAG, pars)
      DAGset <- do.call(apprxPack$func, pars)
      print(paste("number of dags in set: ", dim(DAGset)[3]))
      return(DAGset)
    } , adsm=dt, mc.cores=numCoresAdsm, SIMPLIFY=FALSE)
      }, dt=adsm_tree, mc.cores=numCoresDt, SIMPLIFY=FALSE)
    
    # ... get union of hypothesis for each learner-dataSetting
    print("get union of hypotheses:")
    
    hypSc_tree.unionsHyps <- lapply(hypSc_tree, function(dt) lapply(dt, function(hypSc){
      # i <- 1; j <- 2; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
      #print(hypSc)
      dt_char <- hypSc["dt"]
      hypSc_char <- hypSc["hypSc"]
      adsms_char <- hypSc_tree.unions[[dt_char]][[hypSc_char]]
      hypsList <- adsm_tree.hyps[[dt_char]][adsms_char]
      hypsListID <- lapply(hypsList, getHypID)
      hypArray <- do.call("abind", hypsList)
      # apply union to eliminate duplicates
      hypArray <- unique(hypArray, MARGIN=3)
      dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
      print(paste("number of hypotheses to be scored for dataTreatment - hypothesis scorer ", paste(dt_char, hypSc_char, sep="-"), " is ",dim(hypArray)[3]))
      return(list(hypsListID=hypsListID, hypArray=hypArray))
    }))
    
    
  
    # get scores for each union of hypothesis
    hypSc_tree.scores  <- mcmapply(FUN=function(dt) mcmapply(FUN=function(hypSc){
      # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
      
      dt_char <- hypSc["dt"]
      hypSc_char <- hypSc["hypSc"]
      hypScorer <- eval(parse(text=hypSc_char))
      
      
      #print(paste("learner: ", hypScorer$learner, " data-regime: ", hypScorer$dataReg))
      
      pars <- hypScorer[-1]
     
      parsEff <- pars[-which(names(pars) %in% c("ppTab","plot"))] 
      #print("*****************************************")
      print(paste("hypScorer: ", hypSc_char))
      print("parameters: ")
      print(paste(names(parsEff), parsEff, sep=": "))
      
      pars$x <- dt_tree.x[[dt_char]] 
      pars$hypArray <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
      pars$plot <- plot # pars$plot <- TRUE
      pars$dataNm <- nm  
      pars$folderSave <- folderSave
      
      # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, ppTab=pars$ppTab, plot=TRUE)
      # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, plot=TRUE)
      # cmem_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
      # cmemJoint_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
      
      scores <- try(do.call(hypScorer$func, pars))
      
      return(scores)
      
    }, 
    hypSc=dt, mc.cores=numCoresHypSc, SIMPLIFY=FALSE), 
    dt=hypSc_tree, mc.cores = numCoresDt, SIMPLIFY=FALSE)
    
   
    
    # arrange scores back according to the approx dag set method - i.e. one per recipe
    recipe_tree.scores <- lapply(recipe_tree, function(dt) lapply(dt, function(adsm) lapply(adsm, function(hypSc){
      # i <- 1; j <- 1; k <-1; dt <- recipe_tree[[i]]; adsm <- dt[[j]]; hypSc <- adsm[[k]]
      #print(hypSc)
      dt_char <- hypSc["dt"]
      adsm_char <- hypSc["adsm"]
      hypSc_char <- hypSc["hypSc"]
      
      hypsIDRecipe <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypsListID[[adsm_char]]
      hypsArrRecipe <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
      
      indxInUnions <- match(hypsIDRecipe$id, as.numeric(dimnames(hypsArrRecipe)$dag))  
      
      scores <- hypSc_tree.scores[[dt_char]][[hypSc_char]] 
      if(!class(scores)=="try-error"){
        scores <- scores[indxInUnions,,drop=FALSE]
      } else{
        print("Error, recipe didnt work")
        scores <- NA
      }
      
      return(scores)
    })))
    
    
    
    return(recipe_tree.scores)
    
  }, 
  dag=dataList$dags[firstDag:lastDag], x=dataList$xs[firstDag:lastDag], 
  n=dataList$noiss[firstDag:lastDag], nm=dataList$names[firstDag:lastDag], SIMPLIFY=FALSE, mc.cores= numCoresScores)
  
  names(scores) <- dataList$names[firstDag:lastDag]
  
  return(scores)
  
}

applyLearnRecipesMany_temp <- function(recipe, dataList, numCoresScores=1, numCoresDt=1, numCoresAdsm=1, numCoresHypSc=1, plot=FALSE, folderSave){
  
  
  # data treatment tree - to store normalized data - same structure per data set
  # 1 level - 1 branch per unique data treatment type in recipe
  
  dts <- recipe$dataTreatment
  dts_uni <- unique(dts)
  # which data treatments will we effectively run
  dts_eff <- dts_uni
  
  dt_tree <- lapply(dts_eff, function(dt) c(dt=dt))
  names(dt_tree) <- dts_eff
  
  # approx dag set method  tree - to store hypothesis sets - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x approx dag set method in recipe
  
  adsm_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm) c(dt, adsm=adsm))
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # hypothesis scorer tree - to store scores - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x hypothesis scorer in recipe
  
  hypSc_tree <- lapply(dt_tree, function(dt){ 
    hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt)]
    hypScs_uni <- unique(hypScs)
    res <- lapply(hypScs_uni, function(hypSc) c(dt, hypSc=hypSc))
    names(res) <- hypScs_uni
    return(res)
  })
  
  # recipe tree -   to store scores - same structure per data set
  # 3 levels- 1 branch per recipe (i.e. unique data treatment type x approx dag set method x hypothesis scorer)
  
  recipe_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm){
      hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt & recipe$approxDagSetMethod==adsm)]
      hypScs_uni <- unique(hypScs)
      res <- lapply(hypScs, function(hypSc) c(dt, adsm=adsm, hypSc=hypSc))
      names(res) <- hypScs_uni
      return(res)
    })
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # the adsm we run one for each leaf of adsm_tree on the data corresponding to the 1st level branch of that tree
  
  # to get the union of hypothesis, for each leaf of the hypSc_tree we look for all the leafs in the recipe tree
  # with that hypSc and take the adsm corresponding to that branch
  
  # eg - how to make union of hypothesis
  hypSc_tree.unions <- lapply(hypSc_tree, function(dt1) lapply(dt1, function(hypSc1){
    # i <- 1; j <- 1; dt1 <- hypSc_tree[[i]]; hypSc1 <- dt1[[j]]  
    # find hypSc1 in recipe tree and give back the adsm corresponding to that branch
    #print("****************") 
    #print(hypSc1)
    indx <- which(sapply(recipe_tree[[hypSc1["dt"]]], function(adsm) any(sapply(adsm, function(hypSc2){
      # k <-1; l <- 1; adsm <- recipe_tree[[hypSc1[["dt"]]]][[k]]; hypSc2 <- adsm[[l]] 
      hypSc2["hypSc"]==hypSc1["hypSc"]
      
    }))))
    res <- names(recipe_tree[[hypSc1["dt"]]])[indx]
    return(res)    
  }))
  
  # Say we have scored each hypothesis and stored in hte hypothesis scorer tree, how do put into
  # recipe tree appropriately? We traverse the recipe tree, stop at the adsm level, take the
  # hypothesis set that we want there, go up to the hypothesis level and then look for the hypothesis set
  # we just collected in the leaf of the hypothesis scorer tree corresponding to the dt-hypSc of the leaf we are in
  
  firstDag <- 1
  lastDag <- length(dataList$dags)
  
  # for each data set ...
  scores <- mcmapply(FUN=function(dag, x, n, nm) {
    # i <- 1; dag <- dataList$dags[[i]]; x <- dataList$xs[[i]]; n <- dataList$noiss[[i]]; nm <- dataList$names[i]  ;plot(getGraph(dag))
    print("*****************************************")
    print(paste("data set:  ", nm))
    if(is.character(x)){
      x <- read.csv(x, sep="", header=F)
      
      # order cause first, effect second
      x <- x[, n]
      
      colnames(x) <- colnames(dag)
      x <- as.matrix(x)
    }
    print(paste("number of samples: ", nrow(x)))
    print(head(x))
    print("true dag")
    print(dag)
    
    
    # peform data treatments
    print("performing data treatments")
    
    dt_tree.x <- lapply(dt_tree, function(dt){
      # i <- 1; dt <- dt_tree[[i]]
      #print(paste("data treatment:", dt))
      dataTreatment <- eval(parse(text=dt))
      
      if(nrow(x) > dataTreatment$maxPoints) x <- x[sample(1:nrow(x), dataTreatment$maxPoints),]
      if(dataTreatment$perm) x <- permData(x, dataTreatment$perm)
      
      x <- apply(x, 2, dataTreatment$scalingFunc)
    })
    
    # lapply(dt_tree.x, function(el) apply(el, 2, function(col) c(mu=mean(col), sigma=sd(col), low=min(col), high=max(col))))
    
    
    # get the approximate hypothesis set for each approx method x data setting : the approx method, with the specific data used
    # must be independent of the learner used. This also prevents approx methods with certain randomnmess (ex permutation based hsic used)
    # being run more than once (for diff learners) and resulting in different hypothesis sets despite having exact same parameters.
    
    print("run approximate dag set methods for each datasetting")
    
    adsm_tree.hyps <- mcmapply(FUN=function(dt){
      # i <- 1; dt <- adsm_tree[[i]]
      DAGsets <- mcmapply(FUN=function(adsm){
        # j <- 2; adsm <- dt[[j]]
        x <- dt_tree.x[[adsm["dt"]]]
        print("************************************")
        print(paste("apprx dag set method: ", adsm["adsm"]))
        apprxPack <- eval(parse(text=adsm["adsm"]))
        apprxPack$pars$dataNm <- nm
        apprxPack$pars$folderSave <- folderSave
        
        pars <- list(data=x, trueDAG=dag, pars=apprxPack$pars)
        # data <- pars$data; trueDAG <- pars$trueDAG; pars <- pars$pars
        # names(pars)
        # minDAG_cmem(data, trueDAG, pars)
        # minDAG(data, trueDAG, pars)
        # pairwise_cmem(data, trueDAG, pars)
        DAGset <- do.call(apprxPack$func, pars)
        print(paste("number of dags in set: ", dim(DAGset)[3]))
        return(DAGset)
      } , adsm=dt, mc.cores=numCoresAdsm, SIMPLIFY=FALSE)
    }, dt=adsm_tree, mc.cores=numCoresDt, SIMPLIFY=FALSE)
    
    # ... get union of hypothesis for each learner-dataSetting
    print("get union of hypotheses:")
    
    hypSc_tree.unionsHyps <- lapply(hypSc_tree, function(dt) lapply(dt, function(hypSc){
      # i <- 1; j <- 2; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
      #print(hypSc)
      dt_char <- hypSc["dt"]
      hypSc_char <- hypSc["hypSc"]
      adsms_char <- hypSc_tree.unions[[dt_char]][[hypSc_char]]
      hypsList <- adsm_tree.hyps[[dt_char]][adsms_char]
      hypsListID <- lapply(hypsList, getHypID)
      hypArray <- do.call("abind", hypsList)
      # apply union to eliminate duplicates
      hypArray <- unique(hypArray, MARGIN=3)
      dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
      print(paste("number of hypotheses to be scored for dataTreatment - hypothesis scorer ", paste(dt_char, hypSc_char, sep="-"), " is ",dim(hypArray)[3]))
      return(list(hypsListID=hypsListID, hypArray=hypArray))
    }))
    
    
    
    # get scores for each union of hypothesis
    hypSc_tree.scores  <- mcmapply(FUN=function(dt) mcmapply(FUN=function(hypSc){
      # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
      
      dt_char <- hypSc["dt"]
      hypSc_char <- hypSc["hypSc"]
      hypScorer <- eval(parse(text=hypSc_char))
      
      
      #print(paste("learner: ", hypScorer$learner, " data-regime: ", hypScorer$dataReg))
      
      pars <- hypScorer[-1]
      
      parsEff <- pars[-which(names(pars) %in% c("ppTab","plot"))] 
      #print("*****************************************")
      print(paste("hypScorer: ", hypSc_char))
      print("parameters: ")
      print(paste(names(parsEff), parsEff, sep=": "))
      
      pars$x <- dt_tree.x[[dt_char]] 
      pars$hypArray <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
      pars$plot <- plot # pars$plot <- TRUE
      pars$dataNm <- nm  
      pars$folderSave <- folderSave
      
      # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, ppTab=pars$ppTab, plot=TRUE)
      # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, plot=TRUE)
      # cmem_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
      # cmemJoint_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
      
      scores <- try(do.call(hypScorer$func, pars))
      
      return(scores)
      
    }, 
    hypSc=dt, mc.cores=numCoresHypSc, SIMPLIFY=FALSE), 
    dt=hypSc_tree, mc.cores = numCoresDt, SIMPLIFY=FALSE)
    
    
    
    
    
    
    return(hypSc_tree.scores)
    
  }, 
  dag=dataList$dags[firstDag:lastDag], x=dataList$xs[firstDag:lastDag], 
  n=dataList$noiss[firstDag:lastDag], nm=dataList$names[firstDag:lastDag], SIMPLIFY=FALSE, mc.cores= numCoresScores)
  
  names(scores) <- dataList$names[firstDag:lastDag]
  
  return(scores)
  
}


applyLearnRecipesOne <- function(recipe, data,  numCoresDt=1, numCoresHypSc=2, plot=FALSE, folderSave){
  
  
  # data treatment tree - to store normalized data - same structure per data set
  # 1 level - 1 branch per unique data treatment type in recipe
  
  {
  dts <- recipe$dataTreatment
  dts_uni <- unique(dts)
  # which data treatments will we effectively run
  dts_eff <- dts_uni
  
  dt_tree <- lapply(dts_eff, function(dt) c(dt=dt))
  names(dt_tree) <- dts_eff
  
  # approx dag set method  tree - to store hypothesis sets - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x approx dag set method in recipe
  
  adsm_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm) c(dt, adsm=adsm))
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # hypothesis scorer tree - to store scores - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x hypothesis scorer in recipe
  
  hypSc_tree <- lapply(dt_tree, function(dt){ 
    hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt)]
    hypScs_uni <- unique(hypScs)
    res <- lapply(hypScs_uni, function(hypSc) c(dt, hypSc=hypSc))
    names(res) <- hypScs_uni
    return(res)
  })
  
  # recipe tree -   to store scores - same structure per data set
  # 3 levels- 1 branch per recipe (i.e. unique data treatment type x approx dag set method x hypothesis scorer)
  
  recipe_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm){
      hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt & recipe$approxDagSetMethod==adsm)]
      hypScs_uni <- unique(hypScs)
      res <- lapply(hypScs, function(hypSc) c(dt, adsm=adsm, hypSc=hypSc))
      names(res) <- hypScs_uni
      return(res)
    })
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # the adsm we run one for each leaf of adsm_tree on the data corresponding to the 1st level branch of that tree
  
  # to get the union of hypothesis, for each leaf of the hypSc_tree we look for all the leafs in the recipe tree
  # with that hypSc and take the adsm corresponding to that branch
  
  # eg - how to make union of hypothesis
  hypSc_tree.unions <- lapply(hypSc_tree, function(dt1) lapply(dt1, function(hypSc1){
    # i <- 1; j <- 1; dt1 <- hypSc_tree[[i]]; hypSc1 <- dt1[[j]]  
    # find hypSc1 in recipe tree and give back the adsm corresponding to that branch
    #print("****************") 
    #print(hypSc1)
    indx <- which(sapply(recipe_tree[[hypSc1["dt"]]], function(adsm) any(sapply(adsm, function(hypSc2){
      # k <-1; l <- 1; adsm <- recipe_tree[[hypSc1[["dt"]]]][[k]]; hypSc2 <- adsm[[l]] 
      hypSc2["hypSc"]==hypSc1["hypSc"]
      
    }))))
    res <- names(recipe_tree[[hypSc1["dt"]]])[indx]
    return(res)    
  }))
  
  # Say we have scored each hypothesis and stored in hte hypothesis scorer tree, how do put into
  # recipe tree appropriately? We traverse the recipe tree, stop at the adsm level, take the
  # hypothesis set that we want there, go up to the hypothesis level and then look for the hypothesis set
  # we just collected in the leaf of the hypothesis scorer tree corresponding to the dt-hypSc of the leaf we are in
  

  dag <- data$dag
  x <- data$x
  n <- data$nois
  nm <- data$name  
  #print("*****************************************")
  #print(paste("data set:  ", nm))
  if(is.character(x)){
    x <- read.csv(x, sep="", header=F)
    # order cause first, effect second
    x <- x[, n]
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
  }
  #print(paste("number of samples: ", nrow(x)))
  #print(head(x))
  #print("true dag")
  #print(dag)
    
  # peform data treatments
  #print("performing data treatments")
    
  dt_tree.x <- lapply(dt_tree, function(dt){
      # i <- 1; dt <- dt_tree[[i]]
      #print(paste("data treatment:", dt))
      dataTreatment <- eval(parse(text=dt))
      
      if(nrow(x) > dataTreatment$maxPoints) x <- x[sample(1:nrow(x), dataTreatment$maxPoints),]
      if(dataTreatment$perm) x <- permData(x, dataTreatment$perm)
      
      x <- apply(x, 2, dataTreatment$scalingFunc)
    })
    
  # lapply(dt_tree.x, function(el) apply(el, 2, function(col) c(mu=mean(col), sigma=sd(col), low=min(col), high=max(col))))
    
  # get the approximate hypothesis set for each approx method x data setting : the approx method, with the specific data used
  # must be independent of the learner used. This also prevents approx methods with certain randomnmess (ex permutation based hsic used)
  # being run more than once (for diff learners) and resulting in different hypothesis sets despite having exact same parameters.
    
  #print("run approximate dag set methods for each datasetting")
    
  adsm_tree.hyps <- lapply(adsm_tree, function(dt){
    # i <- 1; dt <- adsm_tree[[i]]
    DAGsets <- lapply(dt, function(adsm){
      # j <- 5; adsm <- dt[[j]]
      x <- dt_tree.x[[adsm["dt"]]]
      print("************************************")
      print(paste("apprx dag set method: ", adsm["adsm"]))
      apprxPack <- eval(parse(text=adsm["adsm"]))
      apprxPack$pars$dataNm <- nm
      apprxPack$pars$folderSave <- folderSave
      
      pars <- list(data=x, trueDAG=dag, pars=apprxPack$pars)
      # data <- pars$data; trueDAG <- pars$trueDAG; pars <- pars$pars
      # names(pars)
      # minDAG_cmem(data, trueDAG, pars)
      # minDAG(data, trueDAG, pars)
      # pairwise_cmemComp(data, trueDAG, pars)
      DAGset <- do.call(apprxPack$func, pars)
      #print(paste("number of dags in set: ", dim(DAGset)[3]))
      return(DAGset)
    } )
  })
    
  # ... get union of hypothesis for each learner-dataSetting
  #print("get union of hypotheses:")
    
  
  
  hypSc_tree.unionsHyps <- lapply(hypSc_tree, function(dt) lapply(dt, function(hypSc){
    # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
    #print(hypSc)
    dt_char <- hypSc["dt"]
    hypSc_char <- hypSc["hypSc"]
    adsms_char <- hypSc_tree.unions[[dt_char]][[hypSc_char]]
    hypsList <- adsm_tree.hyps[[dt_char]][adsms_char]
    hypsListID <- lapply(hypsList, getHypID)
    hypArray <- do.call("abind", hypsList)
    # apply union to eliminate duplicates
    hypArray <- unique(hypArray, MARGIN=3)
    #names(dimnames(hypArray)) <- c("from","to","dag")
    dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
    #print(paste("number of hypotheses to be scored for dataTreatment - hypothesis scorer ", paste(dt_char, hypSc_char, sep="-"), " is ",dim(hypArray)[3]))
    return(list(hypsListID=hypsListID, hypArray=hypArray))
  }))
  
  }
  
  # get scores for each union of hypothesis
  hypSc_tree.scores  <- mcmapply(FUN=function(dt) mcmapply(FUN=function(hypSc){
    # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
      
    dt_char <- hypSc["dt"]
    hypSc_char <- hypSc["hypSc"]
    hypScorer <- eval(parse(text=hypSc_char))
      
      
    
      
    pars <- hypScorer[-1]
      
    parsEff <- pars[-which(names(pars) %in% c("ppTab","plot"))] 
    
     #print(paste("learner: ", hypScorer$learner, " data-regime: ", hypScorer$dataReg))
     print("*****************************************")
     print(paste("hypScorer: ", hypSc_char))
     print("parameters: ")
     print(paste(names(parsEff), parsEff, sep=": "))
      
    pars$x <- dt_tree.x[[dt_char]] 
    pars$hypArray <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
    pars$plot <- plot # pars$plot <- TRUE
    pars$dataNm <- nm  
    pars$folderSave <- folderSave  
    pars$hypSc_char <- hypSc_char
    
    # x=pars$x; hypArray=pars$hypArray; dataNm=pars$dataNm; folderSave=pars$folderSave
    # rnd_hypScorer(x, hypArray, ppTab=NULL, plot=FALSE, dataNm,  folderSave)
    # x=pars$x; hypArray=pars$hypArray; dataReg=pars$dataReg; learner=pars$learner; complexityPack=pars$complexityPack; ppTab=pars$ppTab; plot=TRUE; dataNm=pars$dataNm; folderSave=pars$folderSave; hypSc_char= pars$hypSc_char
    # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, ppTab=pars$ppTab, plot=TRUE, dataNm=pars$dataNm, folderSave=pars$folderSave, hypSc_char= pars$hypSc_char)
    # x=pars$x; hypArray=pars$hypArray; dataReg=pars$dataReg; learner=pars$learner; complexityPack=pars$complexityPack;numBoots=pars$numBoot; numPerBoot=pars$numPerBoot; ppTab=pars$ppTab; plot=TRUE; dataNm=pars$dataNm; folderSave=pars$folderSave; hypSc_char= pars$hypSc_char
    # anm_hypScorer_boot(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack,numBoots=pars$numBoot, numPerBoot=pars$numPerBoot, ppTab=pars$ppTab, plot=TRUE, dataNm=pars$dataNm, folderSave=pars$folderSave, hypSc_char= pars$hypSc_char)
    # x=pars$x; hypArray=pars$hypArray; jointFeats=pars$jointFeats; smoothFeats=pars$smoothFeats; ppTab=pars$ppTab; plot=pars$plot; dataNm=pars$dataNm; folderSave=pars$folderSave
    # boot_cmem_hypScorer_eqSig(x=pars$x, hypArray=pars$hypArray, jointFeats=pars$jointFeats,  ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm, folderSave=pars$folderSave)
    # x=pars$x; hypArray=pars$hypArray; cmemLearner=pars$cmemLearner; ppTab=pars$ppTab; plot=pars$plot; dataNm=pars$dataNm; folderSave=pars$folderSave
    # cmem_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm, folderSave=pars$folderSave)
    # cmemJoint_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
    # x=pars$x; hypArray=pars$hypArray; cmemLearner1=pars$cmemLearner1; cmemLearner2=pars$cmemLearner2; noiseLearner=pars$noiseLearner; augmentData=pars$augmentData;ppTab=pars$ppTab; plot=pars$plot; dataNm=pars$dataNm;  folderSave=pars$folderSave
    # cmem_hypScorer_comp(x=pars$x, hypArray=pars$hypArray, cmemLearner1=pars$cmemLearner1, cmemLearner2=pars$cmemLearner2, noiseLearner=pars$noiseLearner, augmentData=pars$augmentData,ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm,  folderSave=pars$folderSave)
    # x=pars$x; hypArray=pars$hypArray; cmemLearner1=pars$cmemLearner1; cmemLearner2=pars$cmemLearner2; numBoots=pars$numBoots; numPerBoot=pars$numPerBoot; noiseLearner=pars$noiseLearner; augmentData=pars$augmentData;ppTab=pars$ppTab; plot=pars$plot; dataNm=pars$dataNm;  folderSave=pars$folderSave
    # cmem_hypScorer_comp_boot(x=pars$x, hypArray=pars$hypArray, cmemLearner1=pars$cmemLearner1, cmemLearner2=pars$cmemLearner2, numBoots=pars$numBoots, numPerBoot=pars$numPerBoot,noiseLearner=pars$noiseLearner, augmentData=pars$augmentData,ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm,  folderSave=pars$folderSave)
    # cmem_hypScorer_comp_nn(x=pars$x, hypArray=pars$hypArray, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm,  folderSave=pars$folderSave)
    scores <- try(do.call(hypScorer$func, pars))
      
    return(scores)
      
    }, 
  hypSc=dt, mc.cores=numCoresHypSc, SIMPLIFY=FALSE),  
  dt=hypSc_tree, mc.cores = numCoresDt, SIMPLIFY=FALSE)
  
  # arrange scores back according to the approx dag set method - i.e. one per recipe
  recipe_tree.scores <- lapply(recipe_tree, function(dt) lapply(dt, function(adsm) lapply(adsm, function(hypSc){
    # i <- 1; j <- 1; k <-1; dt <- recipe_tree[[i]]; adsm <- dt[[j]]; hypSc <- adsm[[k]]
    #print(hypSc)
    dt_char <- hypSc["dt"]
    adsm_char <- hypSc["adsm"]
    hypSc_char <- hypSc["hypSc"]
      
    hypsIDRecipe <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypsListID[[adsm_char]]
    hypsArrRecipe <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
    
    #indxInUnions <- match(hypsIDRecipe$id, as.numeric(dimnames(hypsArrRecipe)$dag))   
    indxInUnions <- match(as.character(hypsIDRecipe$id), as.character(dimnames(hypsArrRecipe)$dag))  
      
    scores <- hypSc_tree.scores[[dt_char]][[hypSc_char]]
    if(!class(scores)=="try-error"){
      scores <- scores[indxInUnions,,drop=FALSE]
    } else{
      print("Error, recipe didnt work")
      scores <- NA
    }
    return(scores)
  })))
    
    

  return(recipe_tree.scores)
  
}

# THIS is a hack I made for the confounding case implementing Cieran lees paper
applyLearnRecipesOne_temp <- function(recipe, data,  numCoresDt=1, numCoresHypSc=2, plot=FALSE, folderSave){
  
  
  # data treatment tree - to store normalized data - same structure per data set
  # 1 level - 1 branch per unique data treatment type in recipe
  
  {
  dts <- recipe$dataTreatment
  dts_uni <- unique(dts)
  # which data treatments will we effectively run
  dts_eff <- dts_uni
  
  dt_tree <- lapply(dts_eff, function(dt) c(dt=dt))
  names(dt_tree) <- dts_eff
  
  # approx dag set method  tree - to store hypothesis sets - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x approx dag set method in recipe
  
  adsm_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm) c(dt, adsm=adsm))
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # hypothesis scorer tree - to store scores - same structure per data set
  # 2 levels - 1 branch per unique data treatment type x hypothesis scorer in recipe
  
  hypSc_tree <- lapply(dt_tree, function(dt){ 
    hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt)]
    hypScs_uni <- unique(hypScs)
    res <- lapply(hypScs_uni, function(hypSc) c(dt, hypSc=hypSc))
    names(res) <- hypScs_uni
    return(res)
  })
  
  # recipe tree -   to store scores - same structure per data set
  # 3 levels- 1 branch per recipe (i.e. unique data treatment type x approx dag set method x hypothesis scorer)
  
  recipe_tree <- lapply(dt_tree, function(dt){ 
    adsms <- recipe$approxDagSetMethod[which(recipe$dataTreatment==dt)]
    adsms_uni <- unique(adsms)
    res <- lapply(adsms_uni, function(adsm){
      hypScs <- recipe$hypScorer[which(recipe$dataTreatment==dt & recipe$approxDagSetMethod==adsm)]
      hypScs_uni <- unique(hypScs)
      res <- lapply(hypScs, function(hypSc) c(dt, adsm=adsm, hypSc=hypSc))
      names(res) <- hypScs_uni
      return(res)
    })
    names(res) <- adsms_uni
    return(res)
  })
  
  
  # the adsm we run one for each leaf of adsm_tree on the data corresponding to the 1st level branch of that tree
  
  # to get the union of hypothesis, for each leaf of the hypSc_tree we look for all the leafs in the recipe tree
  # with that hypSc and take the adsm corresponding to that branch
  
  # eg - how to make union of hypothesis
  hypSc_tree.unions <- lapply(hypSc_tree, function(dt1) lapply(dt1, function(hypSc1){
    # i <- 1; j <- 1; dt1 <- hypSc_tree[[i]]; hypSc1 <- dt1[[j]]  
    # find hypSc1 in recipe tree and give back the adsm corresponding to that branch
    #print("****************") 
    #print(hypSc1)
    indx <- which(sapply(recipe_tree[[hypSc1["dt"]]], function(adsm) any(sapply(adsm, function(hypSc2){
      # k <-1; l <- 1; adsm <- recipe_tree[[hypSc1[["dt"]]]][[k]]; hypSc2 <- adsm[[l]] 
      hypSc2["hypSc"]==hypSc1["hypSc"]
      
    }))))
    res <- names(recipe_tree[[hypSc1["dt"]]])[indx]
    return(res)    
  }))
  
  # Say we have scored each hypothesis and stored in hte hypothesis scorer tree, how do put into
  # recipe tree appropriately? We traverse the recipe tree, stop at the adsm level, take the
  # hypothesis set that we want there, go up to the hypothesis level and then look for the hypothesis set
  # we just collected in the leaf of the hypothesis scorer tree corresponding to the dt-hypSc of the leaf we are in
  
  
  dag <- data$dag
  x <- data$x
  n <- data$nois
  nm <- data$name  
  #print("*****************************************")
  #print(paste("data set:  ", nm))
  if(is.character(x)){
    x <- read.csv(x, sep="", header=F)
    # order cause first, effect second
    x <- x[, n]
    colnames(x) <- colnames(dag)
    x <- as.matrix(x)
  }
  #print(paste("number of samples: ", nrow(x)))
  #print(head(x))
  #print("true dag")
  #print(dag)
  
  # peform data treatments
  #print("performing data treatments")
  
  dt_tree.x <- lapply(dt_tree, function(dt){
    # i <- 1; dt <- dt_tree[[i]]
    #print(paste("data treatment:", dt))
    dataTreatment <- eval(parse(text=dt))
    
    if(nrow(x) > dataTreatment$maxPoints) x <- x[sample(1:nrow(x), dataTreatment$maxPoints),]
    if(dataTreatment$perm) x <- permData(x, dataTreatment$perm)
    
    x <- apply(x, 2, dataTreatment$scalingFunc)
  })
  
  # lapply(dt_tree.x, function(el) apply(el, 2, function(col) c(mu=mean(col), sigma=sd(col), low=min(col), high=max(col))))
  
  # get the approximate hypothesis set for each approx method x data setting : the approx method, with the specific data used
  # must be independent of the learner used. This also prevents approx methods with certain randomnmess (ex permutation based hsic used)
  # being run more than once (for diff learners) and resulting in different hypothesis sets despite having exact same parameters.
  
  #print("run approximate dag set methods for each datasetting")
  
  adsm_tree.hyps <- lapply(adsm_tree, function(dt){
    # i <- 1; dt <- adsm_tree[[i]]
    DAGsets <- lapply(dt, function(adsm){
      # j <- 5; adsm <- dt[[j]]
      x <- dt_tree.x[[adsm["dt"]]]
      print("************************************")
      print(paste("apprx dag set method: ", adsm["adsm"]))
      apprxPack <- eval(parse(text=adsm["adsm"]))
      apprxPack$pars$dataNm <- nm
      apprxPack$pars$folderSave <- folderSave
      
      pars <- list(data=x, trueDAG=dag, pars=apprxPack$pars)
      # data <- pars$data; trueDAG <- pars$trueDAG; pars <- pars$pars
      # names(pars)
      # minDAG_cmem(data, trueDAG, pars)
      # minDAG(data, trueDAG, pars)
      # pairwise_cmemComp(data, trueDAG, pars)
      DAGset <- do.call(apprxPack$func, pars)
      #print(paste("number of dags in set: ", dim(DAGset)[3]))
      return(DAGset)
    } )
  })
  
  
  
  
  # ... get union of hypothesis for each learner-dataSetting
  #print("get union of hypotheses:")
  
  hypSc_tree.unionsHyps <- lapply(hypSc_tree, function(dt) lapply(dt, function(hypSc){
    # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
    #print(hypSc)
    dt_char <- hypSc["dt"]
    hypSc_char <- hypSc["hypSc"]
    adsms_char <- hypSc_tree.unions[[dt_char]][[hypSc_char]]
    hypsList <- adsm_tree.hyps[[dt_char]][adsms_char]
    hypsListID <- lapply(hypsList, getHypID)
    hypArray <- do.call("abind", hypsList)
    # apply union to eliminate duplicates
    hypArray <- unique(hypArray, MARGIN=3)
    dimnames(hypArray) <- list(from=colnames(dag), to=colnames(dag), dag=getHypID(hypArray)$id)
    #print(paste("number of hypotheses to be scored for dataTreatment - hypothesis scorer ", paste(dt_char, hypSc_char, sep="-"), " is ",dim(hypArray)[3]))
    return(list(hypsListID=hypsListID, hypArray=hypArray))
  }))
  
  }
  
  # get scores for each union of hypothesis
  hypSc_tree.scores  <- mcmapply(FUN=function(dt) mcmapply(FUN=function(hypSc){
    # i <- 1; j <- 1; dt <- hypSc_tree[[i]]; hypSc <- dt[[j]]
    
    dt_char <- hypSc["dt"]
    hypSc_char <- hypSc["hypSc"]
    hypScorer <- eval(parse(text=hypSc_char))
    
    
    
    
    pars <- hypScorer[-1]
    
    parsEff <- pars[-which(names(pars) %in% c("ppTab","plot"))] 
    
    #print(paste("learner: ", hypScorer$learner, " data-regime: ", hypScorer$dataReg))
    print("*****************************************")
    print(paste("hypScorer: ", hypSc_char))
    print("parameters: ")
    print(paste(names(parsEff), parsEff, sep=": "))
    
    pars$x <- dt_tree.x[[dt_char]] 
    pars$hypArray <- hypSc_tree.unionsHyps[[dt_char]][[hypSc_char]]$hypArray
    pars$plot <- plot # pars$plot <- TRUE
    pars$dataNm <- nm  
    pars$folderSave <- folderSave  
    
    # anm_hypScorer(x=pars$x, hypArray=pars$hypArray, dataReg=pars$dataReg, learner=pars$learner, complexityPack=pars$complexityPack, ppTab=pars$ppTab, plot=TRUE, dataNm=pars$dataNm, folderSave=pars$folderSave)
    # x=pars$x; hypArray=pars$hypArray; cmemLearner=pars$cmemLearner; ppTab=pars$ppTab; plot=pars$plot
    # cmem_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm, folderSave=pars$folderSave)
    # cmemJoint_hypScorer(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot)
    # cmem_hypScorer_comp(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm,  folderSave=pars$folderSave)
    # cmem_hypScorer_comp_nn(x=pars$x, hypArray=pars$hypArray, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm,  folderSave=pars$folderSave)
    # x=pars$x; hypArray=pars$hypArray; cmemLearner=pars$cmemLearner; ppTab=pars$ppTab; plot=pars$plot; dataNm=pars$dataNm; folderSave=pars$folderSave
    # cmem_hypScorer_confounder_isomap(x=pars$x, hypArray=pars$hypArray, cmemLearner=pars$cmemLearner, ppTab=pars$ppTab, plot=pars$plot, dataNm=pars$dataNm, folderSave=pars$folderSave)
    scores <- try(do.call(hypScorer$func, pars))
    
    return(scores)
    
  }, 
  hypSc=dt, mc.cores=numCoresHypSc, SIMPLIFY=FALSE),  
  dt=hypSc_tree, mc.cores = numCoresDt, SIMPLIFY=FALSE)
  
  return(hypSc_tree.scores)
  
}


getData <- function(dataList, nm){
  dag <- dataList$dags[[nm]]
  x <- dataList$xs[[nm]]
  nois <- dataList$noiss[[nm]]
  name <- nm
  res <- list(dag=dag, x=x, nois=nois, name=name)
  return(res)  
}

getSubsetData <- function(dataList, indxs){
  dags <- dataList$dags[indxs]
  xs <- dataList$xs[indxs]
  noiss <- dataList$noiss[indxs]
  names <- dataList$names[indxs]
  res <- list(dags=dags, xs=xs, noiss=noiss, names=names)
  return(res)  
}

#Functions to retrieve lambda, gamma, sigma.rbf.X and sigma.rbfY
# once an experiment is done

# this one is for when i calculated the grid of each
# loss function (cmem_L2, gll, cmem_rel, NCE) in separate learner
# have to run with appropriate cmem_ob
getExpParams1 <- function(recipe, dataList, folderSave){
  
  cmem_learners <- lapply(recipe$hypScorer, function(el){
    hypScorer <- eval(parse(text=el))
    cmem_learner1 <- hypScorer$cmemLearner1
    cmem_learner2 <- hypScorer$cmemLearner2
    return(list(cmem_learner1=cmem_learner1, cmem_learner2=cmem_learner2))
  })
  
  cmem_learners1 <- sapply(cmem_learners, function(el) el$cmem_learner1)
  cmem_learners2 <- sapply(cmem_learners, function(el) el$cmem_learner2)
  
  params <- lapply(dataList$names, function(dataNm){
    #dataNm <- dataList$names[1]
    print(paste("dataNm: ", dataNm))
    paramsData <- mcmapply(function(cmem_learner1, cmem_learner2){
      #indx_lrn <- 4; cmem_learner1 <- cmem_learners1[indx_lrn]; cmem_learner2 <- cmem_learners2[indx_lrn]
      
      cmemLearner1 <- cmem_learner1
      aux <- strsplit(cmemLearner1, "_")[[1]]
      indxDEL1 <- which(aux=="DEL1")
      indxDEL2 <- which(aux=="DEL2")
      cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
      nm1 <- paste(dataNm, cmemLearner1_save, "1",sep="_")
      nm <- nm1
      paramNodeTo_1 <- lapply(colnames(dataList$xs[[dataNm]]), function(nodeTo){
        # nodeTo <- colnames(dataList$xs[[dataNm]])[1]
        regressors <- setdiff(colnames(dataList$xs[[dataNm]]), nodeTo)
        regressorsChar <- paste(regressors, collapse="-")
        regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
        fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
        fileSaveCmemLearner %in% dir(folderSave)
        
        load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
        lambda <- getHyperPar(cmemLearnerAux,"lambda")
        sigma.rbf.X <- getHyperPar(cmemLearnerAux,"sigma.rbf.X")
        sigma.rbf.Y <- getHyperPar(cmemLearnerAux,"sigma.rbf.Y")
        pars <- list(lambda=lambda, sigma.rbf.X=sigma.rbf.X, sigma.rbf.Y=sigma.rbf.Y)
        return(pars)
      })
      names(paramNodeTo_1) <- colnames(dataList$xs[[dataNm]])
      nm1 <- paste(dataNm, cmem_learner1, "1",sep="_")
      nm2 <- paste(cmem_learner2, "2",sep="_")
      nm <- paste("nm1",nm1,"nm2",nm2, sep="_")
      paramNodeTo_2 <- lapply(colnames(dataList$xs[[dataNm]]), function(nodeTo){
        # nodeTo <- colnames(dataList$xs[[dataNm]])[1]
        regressors <- setdiff(colnames(dataList$xs[[dataNm]]), nodeTo)
        regressorsChar <- paste(regressors, collapse="-")
        regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
        fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
        fileSaveCmemLearner %in% dir(folderSave)
        if(fileSaveCmemLearner %in% dir(folderSave)){
          load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
          lambda <- getHyperPar(cmemLearnerAux,"lambda")
          sigma.rbf.X <- getHyperPar(cmemLearnerAux,"sigma.rbf.X")
          pars <- list(lambda=lambda, sigma.rbf.X=sigma.rbf.X)
        } else{
          pars <- list(lambda=NA, sigma.rbf.X=NA)
        }
        return(pars)
      })
      names(paramNodeTo_2) <- colnames(dataList$xs[[dataNm]])
      paramNodeTo <- list(round1=paramNodeTo_1, round2=paramNodeTo_2)
      #paramNodeTo <- unlist(paramNodeTo)
      
      return(paramNodeTo)   
    }, cmem_learner1=cmem_learners1, cmem_learner2=cmem_learners2, SIMPLIFY=FALSE)
    names(paramsData) <- recipe$hypScorer
    return(paramsData)
  })
  names(params) <- dataList$names
  return(params)
}
# this one is for when i calculated the grid of each
# loss function in one learner so it involves re-defining the 
# mainLoss and optFunc(grid)
getExpParams2 <- function(recipe, dataList, folderSave){
  
  cmem_learners <- lapply(recipe$hypScorer, function(el){
    hypScorer <- eval(parse(text=el))
    cmem_learner1 <- hypScorer$cmemLearner1
    cmem_learner2 <- hypScorer$cmemLearner2
    return(list(cmem_learner1=cmem_learner1, cmem_learner2=cmem_learner2))
  })
  
  cmem_learners1 <- sapply(cmem_learners, function(el) el$cmem_learner1)
  cmem_learners2 <- sapply(cmem_learners, function(el) el$cmem_learner2)
  
  params <- lapply(dataList$names, function(dataNm){
    #dataNm <- dataList$names[98]
    print(paste("dataNm: ", dataNm))
    paramsData <- mcmapply(function(cmem_learner1, cmem_learner2){
      #indx_lrn <- 4; cmem_learner1 <- cmem_learners1[indx_lrn]; cmem_learner2 <- cmem_learners2[indx_lrn]
      
      cmemLearner1 <- cmem_learner1
      aux <- strsplit(cmemLearner1, "_")[[1]]
      indxDEL1 <- which(aux=="DEL1")
      indxDEL2 <- which(aux=="DEL2")
      cmemLearner1_save <- paste(c(aux[1:(indxDEL1-1)],aux[(indxDEL2+1):length(aux)]), collapse="_")
      nm1 <- paste(dataNm, cmemLearner1_save, "1",sep="_")
      nm <- nm1
      paramNodeTo_1 <- lapply(colnames(dataList$xs[[dataNm]]), function(nodeTo){
        # nodeTo <- colnames(dataList$xs[[dataNm]])[1]
        regressors <- setdiff(colnames(dataList$xs[[dataNm]]), nodeTo)
        regressorsChar <- paste(regressors, collapse="-")
        regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
        fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
        fileSaveCmemLearner %in% dir(folderSave)
        cmemLearnerAux <- eval(parse(text=cmemLearner1))
        cmemLearnerAuxOrig <- cmemLearnerAux
        load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
        cmemLearnerAux$optimizeParams$mainLoss <- cmemLearnerAuxOrig$optimizeParams$mainLoss
        cmemLearnerAux$optimizeParams$optLossFunc <- cmemLearnerAuxOrig$optimizeParams$optLossFunc
        print("re-sets optimal params based on ready made grid and right mainLoss")
        pm <- proc.time()
        cmemLearnerAux <- setParams(learner=cmemLearnerAux, trainData=cmemLearnerAux$hyperParams$trainData, mc_cores=1, numCoresFold=1)
        pm <- proc.time() - pm
        lambda <- getHyperPar(cmemLearnerAux,"lambda")
        sigma.rbf.X <- getHyperPar(cmemLearnerAux,"sigma.rbf.X")
        sigma.rbf.Y <- getHyperPar(cmemLearnerAux,"sigma.rbf.Y")
        pars <- list(lambda=lambda, sigma.rbf.X=sigma.rbf.X, sigma.rbf.Y=sigma.rbf.Y)
        return(pars)
      })
      names(paramNodeTo_1) <- colnames(dataList$xs[[dataNm]])
      nm1 <- paste(dataNm, cmem_learner1, "1",sep="_")
      nm2 <- paste(cmem_learner2, "2",sep="_")
      nm <- paste("nm1",nm1,"nm2",nm2, sep="_")
      paramNodeTo_2 <- lapply(colnames(dataList$xs[[dataNm]]), function(nodeTo){
        # nodeTo <- colnames(dataList$xs[[dataNm]])[1]
        regressors <- setdiff(colnames(dataList$xs[[dataNm]]), nodeTo)
        regressorsChar <- paste(regressors, collapse="-")
        regressionChar <- paste(nodeTo, "on", regressorsChar, sep="")
        fileSaveCmemLearner <- paste(nm, "_", regressionChar, ".RData", sep="")
        fileSaveCmemLearner %in% dir(folderSave)
        if(fileSaveCmemLearner %in% dir(folderSave)){
          load(file=paste(folderSave, fileSaveCmemLearner, sep=""))
          lambda <- getHyperPar(cmemLearnerAux,"lambda")
          sigma.rbf.X <- getHyperPar(cmemLearnerAux,"sigma.rbf.X")
          sigma.rbf.Y <- getHyperPar(cmemLearnerAux,"sigma.rbf.Y")
          if(length(sigma.rbf.Y)==1) sigma.rbf.Y <- rep(sigma.rbf.Y,2)
          names(sigma.rbf.Y) <- c("sigma.rbf.Y.1","sigma.rbf.Y.2")
          pars <- list(lambda=lambda, sigma.rbf.X=sigma.rbf.X, sigma.rbf.Y.1=sigma.rbf.Y[1], sigma.rbf.Y.2=sigma.rbf.Y[2])
        } else{
          pars <- list(lambda=NA, sigma.rbf.X=NA, sigma.rbf.Y.1=NA, sigma.rbf.Y.2=NA)
        }
        return(pars)
      })
      names(paramNodeTo_2) <- colnames(dataList$xs[[dataNm]])
      paramNodeTo <- list(round1=paramNodeTo_1, round2=paramNodeTo_2)
      #paramNodeTo <- unlist(paramNodeTo)
      
      return(paramNodeTo)   
    }, cmem_learner1=cmem_learners1, cmem_learner2=cmem_learners2, SIMPLIFY=FALSE)
    names(paramsData) <- recipe$hypScorer
    return(paramsData)
  })
  names(params) <- dataList$names
  return(params)
}


################################################################################################################*
# Functions for process scores of different causal learners and datasets 
################################################################################################################*


# to apply a function at the data-recipe level you need, the scores, the function to apply,
# you might also need:
# the true dag of each dag, the complexity pack used in order to transform scores, 
# a function at the dataset level to get info about the data or dag,
# a function ath the dataset_recipe level to get info about the scores (which hypotheses, which scores, etc)

deliverFuncToScore_recipeTree <- function(scores, func, dags=rep(NA, length(scores)), ppTab=NULL, hypSc_nameProxy=NULL,
                               getDagPars=function(dag, dataName) return(list()), 
                               getDataLevelPars=function(scoreMat, dag, reci) return(list()), ...){
  
  
  pars <- list(...)
  # pars <- list()
  # for aggregateScores:
  # pars <- list(groupFuncs=groupFuncs, groupFuncPars=groupFuncPars, aggFunc=aggFunc)
  # pars <- list(rank2Func="probRD")
  #print("length(pars)")
  #print(length(pars))
  #print("names(pars)")
  #print(names(pars))
  
  
  # pars <- list()
  
  start <- 1
  finish <- length(scores)
  res <- mapply(FUN=function(score_tree, nm, dag){
    # i <- 1; score_tree <- scores[[i]]; dag <- dags[[i]]; nm <- names(scores)[i]
    print("***********")
    print(paste("score name: ", nm))
    # print("true Dag")
    # print(dag)
    # plot(as(dags[[3]], "graphNEL"))
    
    parsDag <- getDagPars(dag, nm) 
    res <- mapply(FUN=function(dt, dt_name) mapply(FUN=function(adsm, adsm_name) mapply(FUN=function(hypSc, hypSc_name){
      # j <- 1; k <- 1; l <- 1; dt <- score_tree[[j]]; dt_name <- names(score_tree)[j]; adsm <- dt[[k]]; adsm_name <- names(dt)[k]; hypSc <- adsm[[l]]; hypSc_name <- names(adsm)[l]
       print(paste("recipe: ", dt_name, " - ", adsm_name, " - ", hypSc_name ))
      
      # in case of mixed recipes we take the complexity pack corresponding to the first hypothesis scorer
      #hypSc_nameProxy <- strsplit(hypSc_name, "__")[[1]][1]
      if(!is.null(ppTab)) ppTab <- eval(parse(text=hypSc_nameProxy))$ppTab
      pars <- c(pars, ppTab=ppTab)
      reci <- paste(dt_name, adsm_name, hypSc_name, sep=".")
      # scoreMat <- hypSc
      parsDataScoreRecipe <- getDataLevelPars(hypSc, dag, reci)
      # names(parsDataScoreRecipe)
      pars2 <- c(pars, parsDag, parsDataScoreRecipe)
      pars2$scoreMat <- hypSc
      
      names(pars2)
      
      # ppTab <- pars2$ppTab; scoreMat <- pars2$scoreMat
      # rankedDecisions(scoreMat, ppTab)
      # rank2Func <- pars2$rank2Func; scoreMat <- pars2$scoreMat  
      # groupFuncs <- pars2$groupFuncs; groupFuncPars <- pars2$groupFuncPars; aggFunc <- pars2$aggFunc; trueDag <- pars2$trueDag; dagHyps <- pars2$dagHyps; scoreMat <- pars2$scoreMat
      # aggregateScores(scoreMat, dagHyps, trueDag, groupFuncs, groupFuncPars, aggFunc)
      # ppTab <- pars2$ppTab; trueDag <- pars2$trueDag; dagHyps <- pars2$dagHyps; scoreMat <- pars2$scoreMat
      # 
      # measuresDagDistEdge
      # names(pars2)
      # scoreMat <- pars2$scoreMat; dagHyps <- pars2$dagHyps; trueDag <- pars2$trueDag; ppTab <- pars2$ppTab
      # measuresDagDistEdge(scoreMat, dagHyps, trueDag, ppTab)
      res <- do.call(func, pars2)
      
      return(res)
    }, hypSc=adsm, hypSc_name=names(adsm), SIMPLIFY=FALSE), adsm=dt, adsm_name=names(dt),SIMPLIFY=FALSE), dt=score_tree, dt_name=names(score_tree), SIMPLIFY=FALSE)
    return(res)
  }, score_tree=scores[start:finish], nm=names(scores)[start:finish], dag=dags[start:finish], SIMPLIFY=FALSE)
  return(res)
}


# usually the getDataLevel Pars will be "getHypots" which obtains hypotheses in id form which have been scored 
# transforms to matrix form
getHypots <- function(scoreMat, dag, reci){
  p <- nrow(dag)
  ids <- as.numeric(rownames(scoreMat))
  dagHyps <- sapply(ids, function(id) dagIDtoMatrix(dagID=id,p), simplify="array")
  return(list(dagHyps=dagHyps))
}


# functions corresponding to "func" argument of deliverFuncToScore function
#a print infinite, not-a-numbers or nas
printErrors <- function(scoreMat, dataName, reci){
  if(any(is.infinite(scoreMat) | is.nan(scoreMat) | is.na(scoreMat))){
    print(paste("for data ", dataName, " and recipe ", reci, " there are ", sum(is.infinite(scoreMat)) + sum(is.nan(scoreMat)) + sum(is.na(scoreMat)),
                " non valid scores out of ", nrow(scoreMat), " hypothesis x ", ncol(scoreMat), " scores = ", length(scoreMat)))
  } 
  return(NULL)
}

printErrors2 <- function(scoreMat, dataName, reci){
  if(any( is.nan(scoreMat) | is.na(scoreMat))){
    print(paste("for data ", dataName, " and recipe ", reci, " there are ",  sum(is.nan(scoreMat)) + sum(is.na(scoreMat)),
                " non valid scores out of ", nrow(scoreMat), " hypothesis x ", ncol(scoreMat), " scores = ", length(scoreMat)))
  } 
  return(NULL)
}

printNoDecision <- function(scoreMat, dataName, reci){
  chk <- apply(scoreMat, 2, function(col) (all(col==col[1]) & length(col)>1))
  #print(paste("chk: ", chk))
  if(any(chk, na.rm=T)){
    print(paste("for data ", dataName, " and recipe ", reci, " there are ",  sum(chk, na.rm=T) ,
                " non valid scores out of ", nrow(scoreMat), " hypothesis x ", ncol(scoreMat), " scores = ", length(scoreMat)))
  } 
  return(NULL)
}


# make a table of errors
tableErrors <- function(scoreMat){
  res <- apply(scoreMat, 2, function(col) sum(is.infinite(col) | is.nan(col) | is.na(col)))
  nms <- names(res)
  res <- matrix(res, 1, length(res))
  colnames(res) <- nms
  rownames(res) <- "countErrors"
  names(dimnames(res)) <- c("countErrors","score")
  return(res)
}

tableErrors2 <- function(scoreMat){
  res <- apply(scoreMat, 2, function(col) sum(is.nan(col) | is.na(col)))
  nms <- names(res)
  res <- matrix(res, 1, length(res))
  colnames(res) <- nms
  rownames(res) <- "countErrors"
  names(dimnames(res)) <- c("countErrors","score")
  return(res)
}

tableNoDecision <- function(scoreMat){
  res <- apply(scoreMat, 2, function(col) (all(col == col[1]) & length(col)>1))
  nms <- names(res)
  res <- matrix(res, 1, length(res))
  colnames(res) <- nms
  rownames(res) <- "countNoDecisions"
  names(dimnames(res)) <- c("countNoDecisions","score")
  return(res)
}


#b correct not-a-numbers or nas
correctNaNaNs <- function(scoreMat, funcReplace){
  
  scoreMat2 <- apply(scoreMat, 2, function(col){
    indxError <- which(is.nan(col) | is.na(col))
    if(length(indxError)>0){
      print("indxError")
      print(indxError)
      col[indxError] <- do.call(funcReplace, list(col[-indxError]))
    }
    return(col)
  })
  dim(scoreMat2) <- dim(scoreMat)
  dimnames(scoreMat2) <- dimnames(scoreMat)
  #rownames(scoreMat2) <- rownames(scoreMat)
  #colnames(scoreMat2) <- colnames(scoreMat)
  return(scoreMat2)
}

correctNoDecisions <- function(scoreMat){
  
  scoreMat2 <- apply(scoreMat, 2, function(col){
    error <- (all(col==col[1]) & length(col)>1)
    if(error){
      d1 <- duplicated(col)
      d2 <- rev(duplicated(rev(col)))
      dup <- which(d1 | d2)
      print("no decision error")
      if(is.infinite(col[1])){
        
        repl <- rep(Inf, length(dup))
        repl[1:floor(length(dup)/2)] <- -Inf
        col[dup] <- repl
      } else{
        col[dup] <- col[dup] + rnorm(length(dup), sd=10^-8)
      }
      
    }
    return(col)
  })
  dim(scoreMat2) <- dim(scoreMat)
  dimnames(scoreMat2) <- dimnames(scoreMat)
  #rownames(scoreMat2) <- rownames(scoreMat)
  #colnames(scoreMat2) <- colnames(scoreMat)
  return(scoreMat2)
}


#c normalize scores to prob
correctToProb <- function(scoreMat, dagHyps, trueDag){
  # I converted to prob first and corrected for similar hypothesis later so as to deal
  # with non-finite scores
   
  scoreT <- apply(scoreMat, 2, scoreToProb)
  #scoreT <- apply(scoreMat, 2, correctScoreToAdd, hyps=dagHyps)
  dim(scoreT) <- dim(scoreMat)
  scoreT <- apply(scoreT, 2, correctScoreToAdd, hyps=dagHyps)
  #scoreT <- apply(scoreT, 2, scoreToProb)
  dim(scoreT) <- dim(scoreMat)
  dimnames(scoreT) <- dimnames(scoreMat)
  return(scoreT)
}

#d meta score functions

score_pUnifdHSIC <- function(scoreMat){
  # scoreMat <- scores[[1]][[1]]
  newVar <- "pUnifpdHSIC"; expr <- "pvalUnifPart * pdHSIC"
  newScore <- as.matrix(within(as.data.frame(scoreMat), eval(parse(text=paste(newVar, expr, sep=" <- ")))))
  
  dimnms <- dimnames(scoreMat)
  dimnms$score <- c(dimnms$score, newVar)
  
  dimnames(newScore) <- dimnms
  
  return(newScore)
}

# see also aggregateScores and rankedDecisions below

# unwrap scores list structure (data list, then recipe list, each with score matrix) into a database long format
unwrapScoreFunc <- function(scoresFunc, dataList, ws){
  innerDims <- names(dimnames(scoresFunc[[1]][[1]][[1]][[1]]))
  scoresFuncs <- mapply(FUN=function(scrFuncs, dta, w, nms){
    # i <- 1; scrFuncs <- scoresFunc[[i]]; nms <- names(scoresFunc)[i]; dta <- dataList$xs[[i]]; w <- ws[i]
    # print(paste("dataset: ", nms))
    
    res1 <- lapply(names(scrFuncs) ,function(dt){
      # dt <- names(scrFuncs)[1]
      res2 <- do.call(rbind, lapply(names(scrFuncs[[dt]]), function(adsm){ 
        # adsm <- names(scrFuncs[[dt]])[1]
        res3 <- do.call(rbind, lapply(names(scrFuncs[[dt]][[adsm]]), function(hypSc){
          # j <- 1; k <- 1; l <- 2; dt <- names(scrFuncs)[j]; adsm <- names(scrFuncs[[dt]])[k]; hypSc <- names(scrFuncs[[dt]][[adsm]])[l]
          # dt; adsm; hypSc
          #print(paste("hypSc: ", hypSc))
          res4 <- melt(scrFuncs[[dt]][[adsm]][[hypSc]])
        
          res4$dt <- dt
          res4$adsm <- adsm
          res4$hypSc <- hypSc
          return(res4)
        }))
      }))
    })
    res <- do.call(rbind, res1)
    
    if(is.character(dta)) dta <- read.csv(dta, sep="", header=F)
    
    res$dataset <- nms
    res$p <- ncol(dta)
    res$n <- nrow(dta)
    res$w <- w
    res$weighted_value <- res$w*res$value
    
    colnames(res) <- c(innerDims,"value","dt","adsm","hypSc","dataset","p","n", "w","weighted_value")
    res <- res[,c("dataset", "p", "n", "w", "dt","adsm","hypSc", rev(innerDims),"value", "weighted_value")]
    return(res)
  }, scrFuncs=scoresFunc, dta=dataList$xs, w=ws, nms=names(scoresFunc), SIMPLIFY=FALSE)
  
  scoresFuncs <- do.call(rbind, scoresFuncs)
  rownames(scoresFuncs) <- 1:nrow(scoresFuncs)
  
  scoresFuncs$recipe <- apply(scoresFuncs[,c("dt","adsm","hypSc")], 1, paste, collapse=".")
  scoresFuncs$recipeFull <- apply(scoresFuncs[,c("recipe",innerDims[2:length(innerDims)])], 1, paste, collapse=".")
  
  scoresFuncs <- scoresFuncs[,c("dataset", "p", "n", "w","recipeFull", "recipe", "dt", "adsm", "hypSc", rev(innerDims),"value", "weighted_value")]
  return(scoresFuncs)
}

# uwrap scores with lis structure and with same number of hypothesis into an array format
scoreDBToArray <- function(scoreDB, names_datasets, value=c("value", "weighted_value")){
  
  value <- match.arg(value)
  
  # count <- 1
  scoreArr <- sapply(unique(scoreDB$recipeFull), function(reci){
    # reci <- unique(scoreDB$recipeFull)[16]
    # reci <- "stdr.oracleME1.krr2_re_full.SME_spacing_LL.agg_truevsRest"
    # print(paste("count: ", count)); count <<- count + 1
    # print(paste("recipe: ", reci))
    
    indxReci <- which(scoreDB$recipeFull==reci)
    mat <- scoreDB[indxReci,]
    mat <- cast(mat, dataset~hypothesis, value="value")
    
    if(nrow(mat) < length(names_datasets)){
      print(paste("for recipe ", reci, " there are only scores for", nrow(mat), " datasets" ))
    }
    
    hyps <- colnames(mat)[c(2,3)]
    # meta learners are not available for all datasets
    mat <- mat[match(as.character(names_datasets), mat$dataset),]
    mat <- as.matrix(mat[, c("hyp0pos","hyp1pos")])
    
    matFin <- mat<Inf & mat>-Inf & !is.na(mat)
    maxFini <- max(mat[matFin])
    minFini <- min(mat[matFin])
    medFini <- median(mat[matFin])
    if(sum(matFin)==0) {
      mat <- sign(mat)
      print(paste("error for recipe ", reci, " all scores NA"))
      #stop("error!")
    } else{
    
      # our AUC functions cant deal with infinite values so we simply give them
      # the next best or worst rank-score
      mat[mat==-Inf] <- minFini
      mat[mat==Inf] <- maxFini  
      
      # NA's correspond to groups of hypothesis which we were not able to score so
      # we assign the median score
      mat[is.na(mat)] <- medFini
    }
        
    
    
    
    
    dimnames(mat) <- list(dataset=names_datasets, hypothesis=hyps)
    return(mat)
  }, simplify="array")
  
  names(dimnames(scoreArr))[3] <- "recipeFull"
  
  return(scoreArr)
}




#This one only combines recipes where the approx method (and parameters!) and data setting are the same 
# and only the learner changes: guaranteeing that the same hypotheses will be scored by all learners
# with score tree structure we should only combine branches from same dt and adsm

addMetaLearner <- function(scores, func, score_name=NULL, ...){
  pars <- list(...)
  # pars <- list()
  # pars <- list(scrNms=scrNms, lrnNms=NULL)
  # pars <- list(mod=mod)
  # pars <- list(lambdas=lambdas, crit=crit)
  #print("names(pars)")
  #print(names(pars))
  #print(pars)
  
  score_tree <- scores[[1]]
  
 
  anyMixing <- sapply(names(score_tree), function(dt) sapply(names(score_tree[[dt]]), 
                                                       function(adsm){ 
                                                         # dt <- names(score_tree)[1]; adsm <- names(score_tree[[dt]])[1]
                                                         if(length(score_tree[[dt]][[adsm]])>1){
                                                           anyMixing <- TRUE
                                                           print(paste("for data treatment ", dt, 
                                                                       " and approx dag method- ", adsm, " mixing following recipes (", 
                                                                       length(score_tree[[dt]][[adsm]]), ", total)"))
                                                           sapply(names(score_tree[[dt]][[adsm]]), 
                                                                              function(hypSc){ 
                                                                                print(paste(dt, adsm, hypSc, sep="-"))
                                                                                })}
                                                         else{
                                                           # anyMixing <- FALSE
                                                           anyMixing <- TRUE # made it true to see if i can also use to combine scores for one same learner
                                                         }
                                                         return(anyMixing)
                                                         }))
  
  numMixing <- sum(unlist(anyMixing))
  print(paste("mixing ", numMixing, " recipes per data:"))
  
  
  
  
  if(numMixing>0){
    newScores <- mapply(FUN=function(score_tree, nm){
      # i <- 81; score_tree <- scores[[i]]; nm <- names(scores)[i]
      
      #print("*****************************************")
      #print(paste("data set:  ", nm))
      
      # we then form a list of extra-scores, one for each list of recipes to merge and its corresponding resulting
      # hypothesis set
      newScores <- lapply(score_tree, function(dt) lapply(dt, function(adsm){
        # i <- 1; j <- 1; dt <- score_tree[[i]]; adsm <- dt[[j]]
        
        if(length(adsm)>0){ #>1
          # we form an array with the recipes to merge and the corresponding hypotheses
          #scoreArr <- sapply(adsm, function(hypSc) hypSc, simplify="array")
          #names(dimnames(scoreArr))[3] <- "hypScorer"
          pars2 <- pars
          pars2$scoreList <- adsm
          # names(pars2)
          # score_majority(scoreList=pars2$scoreList, scrNms=pars2$scrNms, lrnNms=pars2$lrnNms)
          # score_rf(scoreList=pars2$scoreList, mod=pars2$mod)
          # score_bestTrainingError(scoreList=pars2$scoreList, crit=pars2$crit)
          extra_score <- do.call(func, pars2)
          if(!is.list(extra_score)){
            extra_score <- list(extra_score)
            if(is.null(score_name)){
              nam <- paste(names(adsm), collapse="__")
            } else{
              nam <- score_name
            }
            names(extra_score) <- nam
          }
          newScore <- c(adsm, extra_score)
        } else{
          newScore <- adsm
        }    
      return(newScore)
      }))  
    }, score_tree=scores, nm=names(scores), SIMPLIFY=FALSE)
    names(newScores) <- names(scores)
  } else{
    newScores <- scores
  }
  
  return(newScores)
  
}

# functions to addMetaLearner should recieve a score array where:
# 1d: is the dag/hypothesis scored
# 2d: is the score
# 3d: is the learner

score_bestpUnifpdHSIC <- function(scoreArr){
  if(dim(scoreArr)[1]>1){
    diffPvals <- apply(scoreArr[,"pvalUnifPart",], "hypScorer", function(x) min(x[-which.min(x)])-min(x))
  } else{
    diffPvals <- scoreArr[1,"pvalUnifPart",]
  }
  
  res <- scoreArr[,"pdHSIC",which.max(diffPvals), drop=FALSE]    
  dimnames(res)$score <- "bestpUnifpdHSIC"
  res <- adrop(res, drop=3)
  return(res)
}

score_majority <- function(scoreList, scrNms=NULL, lrnNms=NULL){
  
  if(is.null(lrnNms)) lrnNms <- names(scoreList)
  if(is.null(scrNms)) scrNms <- unique(unlist(sapply(scoreList, function(el) colnames(el))))
  
  rankList <- lapply(scoreList[lrnNms], function(scr) scr[,which(colnames(scr)%in% scrNms),drop=FALSE])
  #sum(sapply(scoreList[lrnNms], function(el) ncol(el)))
  
  rankMat <- do.call(cbind, rankList)
  aux <- unlist(lapply(scoreList[lrnNms], function(el) colnames(el)[which(colnames(el) %in% scrNms)]))
  colnames(rankMat) <- paste(names(aux), aux, sep=".")
  rankMat <- apply(rankMat, 2, function(scr) rank(scr))
  
  
  res1 <- apply(rankMat, 1, sum)
  
  res2 <- apply(rankMat, 1, function(mat) sum(mat>1))
  
  res <- cbind(res1, res2)
  
  colnames(res) <- c("rankMajority_score","majority_score")
  rownames(res) <- rownames(rankMat)
  names(dimnames(res)) <- c("dag","score")
  
  return(res)
}

score_rf <- function(scoreList, mod){
  
  regDB <- melt(scoreList)
  colnames(regDB)[4] <- "hypScorer"
  regDB$recipeShort <- paste( sapply(strsplit(as.character(regDB$hypScorer), "_"), function(el) strsplit(el[4],"AND")[[1]][1]), 
                              sapply(strsplit(as.character(regDB$hypScorer), "_"), function(el) el[5]),regDB$score,sep="_")
  
  length(unique(regDB$recipeShort))
  

  #tab <- cast(regDB, dag~recipeShort, fun.aggregate="length")
  #tab[,1:4]
  
  regDB <- cast(regDB, dag~recipeShort, value="value")
  regDB[1:3, 1:5]
  
  regDB <- regDB[match(rownames(scoreList[[1]]), regDB$dag),]
  
  data.frame(regDB$dag, rownames(scoreList[[1]]))
  
  res <- predict(mod, newdata=regDB, type="prob")[,1]
  
  
  res <- matrix(res)
  colnames(res) <- c("rf_score")
  rownames(res) <- rownames(scoreList[[1]])
  names(dimnames(res)) <- c("dag","score")
  return(res)
}

score_lasso <- function(scoreList, mod){
  
  regDB <- melt(scoreList)
  colnames(regDB)[4] <- "hypScorer"
  regDB$recipeShort <- paste( sapply(strsplit(as.character(regDB$hypScorer), "_"), function(el) strsplit(el[4],"AND")[[1]][1]), 
                              sapply(strsplit(as.character(regDB$hypScorer), "_"), function(el) el[5]),regDB$score,sep="_")
  
  length(unique(regDB$recipeShort))
  
  #tab <- cast(regDB, dag~recipeShort, fun.aggregate="length")
  #tab[1:4,1:4]
  
  regDB <- cast(regDB, dag~recipeShort, value="value")
  regDB[1:3, 1:5]
  
  regDB <- regDB[match(rownames(scoreList[[1]]), regDB$dag),]
  
  data.frame(regDB$dag, rownames(scoreList[[1]]))
  
  X <- as.matrix(as.matrix(regDB[,2:ncol(regDB)]))
  colnames(X) <- colnames(regDB)[2:ncol(regDB)]
  
  indx <- which(colnames(X) %in% rownames(mod$beta))
  X <- X[,indx]
  
  X <- X[,match(colnames(X), rownames(mod$beta))]
  
  data.frame(colnames(X), rownames(mod$beta))
  
  # its 1-prob coz if score is low it means it is better
  res <- 1-predict(mod, newx=X, type="response") 
  
  
  res <- matrix(res)
  colnames(res) <- c("lasso_score")
  rownames(res) <- rownames(scoreList[[1]])
  names(dimnames(res)) <- c("dag","score")
  return(res)
}

# compare accross kernels and kernel parameters a
score_bestTrainingError <- function(scoreList,  crit){
  
  lambdaScores <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[9]))
  kernelXScores <- sapply(strsplit(names(scoreList), "_"), function(el) strsplit(el[4], "AND")[[1]][1])
  featureYScores <- sapply(strsplit(names(scoreList), "_"), function(el) el[5])
  kernX <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[10]))
  featY <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[11]))
  regDB <- melt(scoreList)
  regDB$kernelX <- sapply(strsplit(regDB$L1, "_"), function(el) strsplit(el[4],"AND")[[1]][1])
  regDB$featureY <- sapply(strsplit(regDB$L1, "_"), function(el) el[5])
  regDB$lambda <- as.numeric(sapply(strsplit(regDB$L1, "_"), function(el) el[9]))
  regDB$kernX <- sapply(strsplit(regDB$L1, "_"), function(el) el[10])
  regDB$featY <- sapply(strsplit(regDB$L1, "_"), function(el) el[11])
  tab <- cast(regDB[which(regDB$score==crit),], featureY~featY~dag~kernelX~kernX~lambda, value="value")
  #tab <- abind(tab, abs(tab[,,1,,,,drop=F]-tab[,,2,,,,drop=F]), along=3)
  
  
  #apply(tab, c(1,4), function(mat) print(dim(mat)))
  parsOpt <- apply(tab, c(1,3,4), function(mat){
    #mat <- tab[1,,1,1,,]
    indxMat <- mat==min(mat)
    res <- sapply(1:length(dim(indxMat)), function(k) which(apply(indxMat, k, sum)==1))
    #indxMat[matrix(res,nrow=1,byrow=T)]
    #mat[matrix(res,nrow=1,byrow=T)]
    res <- as.numeric(names(res))
    names(res) <- names(dimnames(indxMat))
    return(res)
    })
  
  
  
  
  indxScores <- sapply(dimnames(parsOpt)$kernelX, 
                                            function(kX) sapply(dimnames(parsOpt)$featureY, 
                                                                function(fY){
   #kX <- dimnames(parsOpt)$kernelX[1]; fY <- dimnames(parsOpt)$featureY[1]
                                                                  res <- sapply(rev(dimnames(parsOpt)$dag), 
                                                                                    function(dag) {
    #dag <- dimnames(parsOpt)$dag[1]; 
    indx <- which(kernelXScores==kX & featureYScores==fY & lambdaScores==parsOpt["lambda",fY,dag,kX] 
                  & featY==parsOpt["featY",fY,dag,kX] & kernX==parsOpt["kernX",fY,dag,kX])
    
    res <- scoreList[[indx]][dag,]
    
    return(res)
  }, simplify="array")
      res <- t(res)
      dimnames(res) <- dimnames(scoreList[[1]])
      
      return(res)
 }, simplify="array"), simplify="array")
  
  dim(indxScores)
  dimnames(indxScores)[3:4] <- list(dimnames(parsOpt)$featureY, dimnames(parsOpt)$kernelX)
  names(dimnames(indxScores))[3:4] <- c("featureY","kernelX")
  
  nmsScores <- sapply(strsplit(names(scoreList), "_"), function(el) paste(c(el[1:8], rep("opt",3)), collapse="_"))
  nmsScores <- unique(nmsScores)
  
  res <- lapply(nmsScores, function(el){
    # el <- nmsScores[1]
    kernX <- strsplit(strsplit(el, "_")[[1]][4],"AND")[[1]][1]
    featY <- strsplit(el, "_")[[1]][5]
    dimnames(indxScores)
    res <- indxScores[,,kernX,featY]
    return(res)
  })
  
  names(res) <- nmsScores
  return(res)
}

score_phaseTrajectory <- function(scoreList){
  
  lambdaScores <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[9]))
  kernelXScores <- sapply(strsplit(names(scoreList), "_"), function(el) strsplit(el[4], "AND")[[1]][1])
  featureYScores <- sapply(strsplit(names(scoreList), "_"), function(el) el[5])
  kernX <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[10]))
  featY <- as.numeric(sapply(strsplit(names(scoreList), "_"), function(el) el[11]))
  regDB <- melt(scoreList)
  regDB$kernelX <- sapply(strsplit(regDB$L1, "_"), function(el) strsplit(el[4],"AND")[[1]][1])
  regDB$featureY <- sapply(strsplit(regDB$L1, "_"), function(el) el[5])
  regDB$lambda <- as.numeric(sapply(strsplit(regDB$L1, "_"), function(el) el[9]))
  regDB$kernX <- sapply(strsplit(regDB$L1, "_"), function(el) el[10])
  regDB$featY <- sapply(strsplit(regDB$L1, "_"), function(el) el[11])
  
  #indx <- which(regDB$score %in% c("EDML","KCDC"))
  cast(regDB, dag+kernX+featY+lambda~score, value="value", fun.aggregate=length)
  regDB2 <- cast(regDB, dag+kernX+featY+lambda~score, value="value")
  
  # Test CV error
  p <- ggplot(regDB2)
  p <- p + geom_point(aes(x=rev(lambda), y=log(CVTeRNE,10), colour=factor(dag)))
  p <- p + facet_grid(kernX~featY, scales="free")
  p
  
  # KCDC
  p <- ggplot(regDB2)
  p <- p + geom_point(aes(x=rev(lambda), y=KCDC, colour=factor(dag)))
  p <- p + geom_line(aes(x=rev(lambda), y=KCDC, colour=factor(dag)))
  p <- p + facet_grid(kernX~featY, scales="free")
  p
  
  tab <- cast(regDB[which(regDB$score==crit),], featureY~featY~dag~kernelX~kernX~lambda, value="value")
  #tab <- abind(tab, abs(tab[,,1,,,,drop=F]-tab[,,2,,,,drop=F]), along=3)
  
  
  #apply(tab, c(1,4), function(mat) print(dim(mat)))
  parsOpt <- apply(tab, c(1,3,4), function(mat){
    #mat <- tab[1,,1,1,,]
    indxMat <- mat==min(mat)
    res <- sapply(1:length(dim(indxMat)), function(k) which(apply(indxMat, k, sum)==1))
    #indxMat[matrix(res,nrow=1,byrow=T)]
    #mat[matrix(res,nrow=1,byrow=T)]
    res <- as.numeric(names(res))
    names(res) <- names(dimnames(indxMat))
    return(res)
  })
  
  
  
  
  indxScores <- sapply(dimnames(parsOpt)$kernelX, 
                       function(kX) sapply(dimnames(parsOpt)$featureY, 
                                           function(fY){
                                             #kX <- dimnames(parsOpt)$kernelX[1]; fY <- dimnames(parsOpt)$featureY[1]
                                             res <- sapply(rev(dimnames(parsOpt)$dag), 
                                                           function(dag) {
                                                             #dag <- dimnames(parsOpt)$dag[1]; 
                                                             indx <- which(kernelXScores==kX & featureYScores==fY & lambdaScores==parsOpt["lambda",fY,dag,kX] 
                                                                           & featY==parsOpt["featY",fY,dag,kX] & kernX==parsOpt["kernX",fY,dag,kX])
                                                             
                                                             res <- scoreList[[indx]][dag,]
                                                             
                                                             return(res)
                                                           }, simplify="array")
                                             res <- t(res)
                                             dimnames(res) <- dimnames(scoreList[[1]])
                                             
                                             return(res)
                                           }, simplify="array"), simplify="array")
  
  dim(indxScores)
  dimnames(indxScores)[3:4] <- list(dimnames(parsOpt)$featureY, dimnames(parsOpt)$kernelX)
  names(dimnames(indxScores))[3:4] <- c("featureY","kernelX")
  
  nmsScores <- sapply(strsplit(names(scoreList), "_"), function(el) paste(c(el[1:8], rep("opt",3)), collapse="_"))
  nmsScores <- unique(nmsScores)
  
  res <- lapply(nmsScores, function(el){
    # el <- nmsScores[1]
    kernX <- strsplit(strsplit(el, "_")[[1]][4],"AND")[[1]][1]
    featY <- strsplit(el, "_")[[1]][5]
    dimnames(indxScores)
    res <- indxScores[,,kernX,featY]
    return(res)
  })
  
  names(res) <- nmsScores
  return(res)
}


########################################################################################*
# Performance functions for causal learning methods
########################################################################################*
# overview

# Performance evaluation

# combine accross learners (but asame approx method and data setting to keep same hypothesis space) or within learners (combine different scores
# for same learner-aprroxDagset-recycle)

measuresDagDistEdge <- function(scoreMat, dagHyps, trueDag, ppTab){
  
  p <- dim(trueDag)[1]
  offDiag <- !diag(p)
  
  # Transform scores to weights according to complexity pack chosen
  tab <- eval(parse(text=ppTab))
  scoreIds <- sapply(strsplit(colnames(scoreMat), "_"), function(el) el[[1]]) 
  indx <- match(scoreIds, sapply(strsplit(as.character(tab$id), "_"), function(el) el[1]))
  rankFunc <- tab$rankFunc[indx]
  probFunc <- tab$probFunc[indx]
  
  
  ws <- sapply(1:ncol(scoreMat), function(i){
    #print(i)
    # i <- 1
    #print("******************")
    #print(paste("score :", colnames(scoreMat)[i]))
    #print(scoreMat[,i])
    #res <- do.call(rankFunc[i], list(x=scoreMat[,i], hyps=dagHyps))
    #res <- do.call(probFunc[i], list(res, dagHyps))
    
    # if we get the prob first we can handle Infinite values and then correct for similarity of hypothesis
    res <- do.call(probFunc[i], list(x=scoreMat[,i], hyps=dagHyps))
    res <- do.call(rankFunc[i], list(x=res, hyps=dagHyps))
    
    
  }, simplify="matrix")
  dim(ws) <- dim(scoreMat)
  
  
  # Best dag per score
  indxWinner <- apply(scoreMat, 2, which.min)
  chosenDag <- dagHyps[,,indxWinner, drop=FALSE]
  
  # plot(as(chosenDag[,,1], "graphNEL"))
  # plot(dag2cpdag(as(chosenDag[,,1], "graphNEL")))
  chosenCpdag <- sapply(1:(dim(chosenDag)[3]), function(i) dag2cpdag(chosenDag[,,i])*1, simplify="array")
  
  # plot(as(dag, "graphNEL"))
  # plot(dag2cpdag(as(dag, "graphNEL")))
  # dag; dag2cpdag(dag)*1
  trueCpdag <- dag2cpdag(trueDag)*1
  
  cnst <- choose(p,2)
  
  # edge distances best score
  SHDs <- apply(chosenDag, 3, SHD, dagTrue=trueDag)/cnst
  edgeDists <- apply(chosenDag, 3, edgeDist, dagTrue=trueDag)/cnst
  nonEdgeDists <- apply(chosenDag, 3, nonEdgeDist, dagTrue=trueDag)/cnst
  totalEdgeDists <- (edgeDists + nonEdgeDists)/2
  #apply(chosenDag, 3, totalEdgeDist, dagTrue=trueDag)/(2*cnst) 
  
  #weighted edgeDistances
  edgeDistsW <- apply(dagHyps, 3, edgeDist, dagTrue=trueDag)/cnst
  nonEdgeDistsW <- apply(dagHyps, 3, nonEdgeDist, dagTrue=trueDag)/cnst
  #totalEdgeDistsW <- apply(dagHyps, 3, totalEdgeDist, dagTrue=trueDag)/(2*cnst)
  edgeDistsW <- apply(ws, 2, function(col) sum(col*edgeDistsW))
  nonEdgeDistsW <- apply(ws, 2, function(col) sum(col*nonEdgeDistsW))
  #totalEdgeDistsW <- apply(ws, 2, function(col) sum(col*totalEdgeDistsW))
  totalEdgeDistsW <- (edgeDistsW + nonEdgeDistsW)/2
  
  resEdgeDists <- rbind(SHDs, edgeDists, nonEdgeDists, totalEdgeDists, edgeDistsW, nonEdgeDistsW, totalEdgeDistsW)
  rownames(resEdgeDists) <- c("SHD","edgeD", "nonEdgeD", "totEdgeD","WedgeD", "WnonEdgeD", "WtotEdgeD")
  colnames(resEdgeDists) <- colnames(scoreMat)
  
  
  # edge distances CPDAG best score
  SHDs <- apply(chosenCpdag, 3, SHD, dagTrue=trueCpdag)/cnst
  edgeDists <- apply(chosenCpdag, 3, edgeDist, dagTrue=trueCpdag)/cnst
  nonEdgeDists <- apply(chosenCpdag, 3, nonEdgeDist, dagTrue=trueCpdag)/cnst
  totalEdgeDists <- (edgeDists + nonEdgeDists)/2
  #apply(chosenCpdag, 3, totalEdgeDist, dagTrue=trueCpdag)/(2*cnst) 
  
  #weighted edgeDistances
  edgeDistsW <- apply(dagHyps, 3, edgeDist, dagTrue=trueCpdag)/cnst
  nonEdgeDistsW <- apply(dagHyps, 3, nonEdgeDist, dagTrue=trueCpdag)/cnst
  #totalEdgeDistsW <- apply(dagHyps, 3, totalEdgeDist, dagTrue=trueCpdag)/(2*cnst)
  edgeDistsW <- apply(ws, 2, function(col) sum(col*edgeDistsW))
  nonEdgeDistsW <- apply(ws, 2, function(col) sum(col*nonEdgeDistsW))
  #totalEdgeDistsW <- apply(ws, 2, function(col) sum(col*totalEdgeDistsW))
  totalEdgeDistsW <- (edgeDistsW + nonEdgeDistsW)/2
  
  resCpdagEdgeDists <- rbind(SHDs, edgeDists, nonEdgeDists, totalEdgeDists, edgeDistsW, nonEdgeDistsW, totalEdgeDistsW)
  rownames(resCpdagEdgeDists) <- paste(c("SHD","edgeD", "nonEdgeD", "totEdgeD","WedgeD", "WnonEdgeD", "WtotEdgeD"), "cpdag", sep="_")
  colnames(resCpdagEdgeDists) <- colnames(scoreMat)
  
  # edgewise contingency table based measures
  conTabs <- apply(chosenDag, 3, function(pred) contingencyTable(pred[offDiag], trueDag[offDiag]))
  ccrs  <- apply(conTabs, 2, correctCR)
  msrs  <- apply(conTabs, 2, misCR)
  ppps  <- apply(conTabs, 2, posPP) 
  npps  <- apply(conTabs, 2, negPP)
  senss <- apply(conTabs, 2, sensitivity)
  specs <- apply(conTabs, 2, specificity) 
  fprs  <- apply(conTabs, 2, fpr)
  fnrs  <- apply(conTabs, 2, fnr)
  tsss  <- apply(conTabs, 2, tss) 
  
  # weighted edgewise contingency table based measures
  conTabs <- apply(dagHyps, 3, function(pred) contingencyTable(pred[offDiag], trueDag[offDiag]))
  ccrsW  <- apply(conTabs, 2, correctCR)
  msrsW  <- apply(conTabs, 2, misCR)
  pppsW  <- apply(conTabs, 2, posPP) 
  nppsW  <- apply(conTabs, 2, negPP)
  senssW <- apply(conTabs, 2, sensitivity)
  specsW <- apply(conTabs, 2, specificity) 
  fprsW  <- apply(conTabs, 2, fpr)
  fnrsW  <- apply(conTabs, 2, fnr)
  tsssW  <- apply(conTabs, 2, tss) 
  ccrsW  <- apply(ws, 2, function(col) sum(col*ccrsW))
  msrsW  <- apply(ws, 2, function(col) sum(col*msrsW))
  pppsW  <- apply(ws, 2, function(col) sum(col*pppsW)) 
  nppsW  <- apply(ws, 2, function(col) sum(col*nppsW))
  senssW <- apply(ws, 2, function(col) sum(col*senssW))
  specsW <- apply(ws, 2, function(col) sum(col*specsW)) 
  fprsW  <- apply(ws, 2, function(col) sum(col*fprsW))
  fnrsW  <- apply(ws, 2, function(col) sum(col*fnrsW))
  tsssW  <- apply(ws, 2, function(col) sum(col*tsssW))
  
  #AUC
  preds <- apply(dagHyps, 3, function(mat) as.numeric(mat[offDiag])) %*% ws
  
  obs <- trueDag[offDiag]
  numTrials <- length(obs)
  numCasesGoal <- numTrials/2
  numCasesNow <- sum(obs)
  numCasesAdd <- numCasesGoal - numCasesNow
  indx0 <- which(obs==0)
  indx1 <- which(obs==1)
  smpl <- as.numeric()
  
  #print("numCasesAdd")
  #print(numCasesAdd)
  
  if(numCasesAdd>0) smpl <- sample(indx0, size=numCasesAdd)
  if(numCasesAdd<0) smpl <- sample(indx1, size=-numCasesAdd)
  
  obsNew <- obs
  obsNew[smpl] <- (!obs[smpl])*1
  predsNew <- preds
  predsNew[smpl,] <- 1-preds[smpl,] 
  
  
  aucs <- apply(predsNew, 2, function(pred) as.numeric(roc(obsNew, pred)$auc))
  
  
  resEdgeCont <- rbind(ccrs,  msrs,  ppps,  npps,  senss,  specs,  fprs,  fnrs, tsss, 
                       ccrsW, msrsW, pppsW, nppsW, senssW, specsW, fprsW, fnrsW, tsssW, aucs)
  rownames(resEdgeCont) <- c("ccr",  "msr",  "ppp", "npp",  "sens",  "spec", "fpr", "fnr", "tss",
                             "Wccr", "Wmsr", "Wppp","Wnpp", "Wsens", "Wspec","Wfpr","Wfnr","Wtss", "edgeAUC")
  colnames(resEdgeCont) <- colnames(scoreMat)
  
  
  # edgewise CPDAG contingency table based measures
  conTabs <- apply(chosenCpdag, 3, function(pred) contingencyTable(pred[offDiag], trueCpdag[offDiag]))
  ccrs  <- apply(conTabs, 2, correctCR)
  msrs  <- apply(conTabs, 2, misCR)
  ppps  <- apply(conTabs, 2, posPP) 
  npps  <- apply(conTabs, 2, negPP)
  senss <- apply(conTabs, 2, sensitivity)
  specs <- apply(conTabs, 2, specificity) 
  fprs  <- apply(conTabs, 2, fpr)
  fnrs  <- apply(conTabs, 2, fnr)
  tsss  <- apply(conTabs, 2, tss) 
  
  # weighted CPDAG edgewise contingency table based measures
  conTabs <- apply(dagHyps, 3, function(pred) contingencyTable(pred[offDiag], trueCpdag[offDiag]))
  ccrsW  <- apply(conTabs, 2, correctCR)
  msrsW  <- apply(conTabs, 2, misCR)
  pppsW  <- apply(conTabs, 2, posPP) 
  nppsW  <- apply(conTabs, 2, negPP)
  senssW <- apply(conTabs, 2, sensitivity)
  specsW <- apply(conTabs, 2, specificity) 
  fprsW  <- apply(conTabs, 2, fpr)
  fnrsW  <- apply(conTabs, 2, fnr)
  tsssW  <- apply(conTabs, 2, tss) 
  ccrsW  <- apply(ws, 2, function(col) sum(col*ccrsW))
  msrsW  <- apply(ws, 2, function(col) sum(col*msrsW))
  pppsW  <- apply(ws, 2, function(col) sum(col*pppsW)) 
  nppsW  <- apply(ws, 2, function(col) sum(col*nppsW))
  senssW <- apply(ws, 2, function(col) sum(col*senssW))
  specsW <- apply(ws, 2, function(col) sum(col*specsW)) 
  fprsW  <- apply(ws, 2, function(col) sum(col*fprsW))
  fnrsW  <- apply(ws, 2, function(col) sum(col*fnrsW))
  tsssW  <- apply(ws, 2, function(col) sum(col*tsssW))
  
  #AUC
  preds <- apply(dagHyps, 3, function(mat) as.numeric(mat[offDiag])) %*% ws
  
  obs <- trueCpdag[offDiag]
  numTrials <- length(obs)
  numCasesGoal <- numTrials/2
  numCasesNow <- sum(obs)
  numCasesAdd <- numCasesGoal - numCasesNow
  indx0 <- which(obs==0)
  indx1 <- which(obs==1)
  smpl <- as.numeric()
  
  #print("numCasesAdd")
  #print(numCasesAdd)
  
  if(numCasesAdd>0) smpl <- sample(indx0, size=numCasesAdd)
  if(numCasesAdd<0) smpl <- sample(indx1, size=-numCasesAdd)
  
  obsNew <- obs
  obsNew[smpl] <- (!obs[smpl])*1
  predsNew <- preds
  predsNew[smpl,] <- 1-preds[smpl,] 
  
  
  aucs <- apply(predsNew, 2, function(pred) as.numeric(roc(obsNew, pred)$auc))
  
  
  resCpdagEdgeCont <- rbind(ccrs,  msrs,  ppps,  npps,  senss,  specs,  fprs,  fnrs, tsss, 
                       ccrsW, msrsW, pppsW, nppsW, senssW, specsW, fprsW, fnrsW, tsssW, aucs)
  rownames(resCpdagEdgeCont) <- paste(c("ccr",  "msr",  "ppp", "npp",  "sens",  "spec", "fpr", "fnr", "tss",
                             "Wccr", "Wmsr", "Wppp","Wnpp", "Wsens", "Wspec","Wfpr","Wfnr","Wtss", "edgeAUC"), "cpdag", sep="_")
  colnames(resCpdagEdgeCont) <- colnames(scoreMat)
  
  
  res <- rbind(resEdgeDists, resCpdagEdgeDists, resEdgeCont, resCpdagEdgeCont)
  
  names(dimnames(res)) <- c("measure","score")
  
  return(res)
}


# I. Dag/Edge wise evaluation: Distance from true hyp to chosen hyp 

#   1) number of edges, non-edges and total diff between true and chosen hyp

#   2) get distancest to all hyp then weight by prob/conf

# II. Edge wise evaluation

#  1) compare true hyp to a single chosen hyp
#    a) accuracy, ROC/AUC, TSS, recall, precison, etc 
#  2) compare true hyp to all hyps and then weight measures by prob
#    b) accuracy, ROC/AUC, TSS, recall, precison, etc

# III. Dag-set wise evaluation

# Ranking of two hypotheses as in Mooij et al, pg 32-33

# 1. divide hypothesis space in two - part with true hyp and other

aggregateScores <- function(scoreMat, dagHyps, trueDag, groupFuncs, groupFuncPars, aggFunc){
  
  # construct names for the different bipartition methods groupFuncs and their parameters groupFuncPars
  nms <- mapply(function(func, pars){
    # i <- 2; func <- groupFuncs[i]; pars <- groupFuncPars[[i]]
    aux <- unlist(pars)
    res <- func
    if(!is.null(aux)){
      # ED- edge distance, NED-non-edge distance, TED- total edge distance
      bit1 <- c("ED","NED","TED")[match(aux["distFun"], c("edgeDist","nonEdgeDist","totalEdgeDist"))]
      bit2 <- paste("k", aux["k"], sep="")
      res <- paste(res, bit1, bit2, sep="")
    }
    return(res)
  }, func=groupFuncs, pars=groupFuncPars, USE.NAMES = F)
  
  
  # Aggregate score matrix
  #print("aggregating score matrix")
  aggScoreMat <- mapply(function(grpFunc, grpFuncPars){
    # i <- 1; grpFunc <- groupFuncs[i]; grpFuncPars <- groupFuncPars[[i]]
    #print("*************")
    #print(paste("group function: ", grpFunc))
    pars <- grpFuncPars
    pars$dagHyps <- dagHyps
    pars$trueDag <- trueDag
    # agg_MECvsRest(dagHyps, trueDag)
    grp <- do.call(grpFunc, pars)
    indx0 <- which(grp==0)
    indx1 <- which(grp==1)
    # x <- scoreMat[,1]
    res <- apply(scoreMat, 2, aggFunc, indx1=indx0, indx2=indx1)
    return(res)
  }, grpFunc=groupFuncs, grpFuncPars=groupFuncPars, SIMPLIFY="array")
  
  names(dimnames(aggScoreMat))[c(1,3)] <- c("hypothesis", "bipartCrit")
  dimnames(aggScoreMat)$bipartCrit <- nms
  
  return(aggScoreMat)
  
}


# a ) true hyp vs else

agg_truevsRest <- function(dagHyps, trueDag){
  grp <- apply(dagHyps, 3, function(dag) any(duplicated(abind(dag, trueDag, along=3), MARGIN=3)))*1
  return(grp)
}

# b )  +- k edges 

agg_distKvsRest <- function(dagHyps, trueDag, distFun, k){
  # is this function the same as agg_truevsRest if I pass k=0??
  distss <- apply(dagHyps, 3, distFun, dagTrue=trueDag)
  grp <- (distss <= k)*1  
  return(grp)
}

#     i) true hyp +- k edges vs else
#     ii) true hyp + k edges vs else (less independencies)
#     iii) true hyp - k edges vs else (more independencies)
# c)   Independencies
#    i) Markov equivalence class vs else

agg_MECvsRest <- function(dagHyps, trueDag){
  
  MEdagTrue <- getMarkovEquivClass(trueDag)
  #p <- nrow(trueDag)
  #dag <- as(trueDag, "graphNEL")
  #cpdag <- dag2cpdag(dag)
  #cpdagMat <- as(cpdag, "matrix")
  #MEdagTrue <- pdag2allDags(cpdagMat)$dags
  #MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
  grp <- (duplicated(abind(MEdagTrue, dagHyps ,along=3), MARGIN=3)[(dim(MEdagTrue)[3]+1):(dim(dagHyps)[3]+dim(MEdagTrue)[3])])*1
  return(grp)
}



#    ii) Markov equivalence class +- k edges

agg_distKMECvsRest <- function(dagHyps, trueDag, distFun, k){
  # is this the same as agg_distKvsRest if I pass k=0??
  MEdagTrue <- getMarkovEquivClass(trueDag)
  #p <- nrow(trueDag)
  #MEdagTrue <- pdag2allDags(pcalg:::dag2cpdag(trueDag))$dags
  #MEdagTrue <- sapply(1:nrow(MEdagTrue),  function(i) matrix(MEdagTrue[i,], p, p, byrow=T), simplify="array")
  distss <- apply(MEdagTrue, 3, function(dag) apply(dagHyps, 3, distFun, dagTrue=dag))
  dim(distss) <- c(dim(dagHyps)[3], dim(MEdagTrue)[3])
  grp <- apply(distss,1, function(row) any(row <= k)*1)
  return(grp)
}

#     A) true hyp +- k edges vs else
#     B) true hyp + k edges vs else (less independencies)
#     C) true hyp - k edges vs else (more independencies)        




# 2: aggregate scores by hypothesis group (max?, mean?, dircetly on score and on prob)

aggMin <- function(x, indx1, indx2){
  
    if(length(indx1)==0){
      # if a certain hypothesis group wasn't scored then it should have very bad "Inf" score
      res1 <- Inf 
    } else{
      if(sum(!is.na(x[indx1]) & !is.nan(x[indx1]))==0){
        # if a certain hypothesis group was scored but produced only NA or NaN errors then
        # we will assign NA aggregate score and later give it the median score
        res1 <- NA
      } else{
        res1 <- min(x[indx1], na.rm=T)
      }
    }
  
    if(length(indx2)==0){
      res2 <- Inf #max(x, na.rm=T)+100 # Inf
    } else{
      if(sum(!is.na(x[indx2]) & !is.nan(x[indx2]))==0){
        res2 <- NA
      } else{
        res2 <- min(x[indx2], na.rm=T)
      }
    }
  
  
  res <- c(res1, res2)
  names(res) <- c("0","1")
  return(res)
  
}

aggSum <- function(x, indx1, indx2){
  
  if(sum(!is.na(x) & ! is.nan(x))==0){
    res1 <- 2
    res2 <- 2
  } else{
    
    if(sum(!is.na(x[indx1]) & ! is.nan(x[indx1]))==0){
      res1 <- Inf 
    } else{
      res1 <- sum(x[indx1], na.rm=T)
    }
  
    if(sum(!is.na(x[indx2]) & ! is.nan(x[indx2]))==0){
      res2 <- Inf
    } else{
      res2 <- sum(x[indx2], na.rm=T)
    }
  }
  
  res <- c(res1, res2)
  names(res) <- c("0","1")
  return(res)
  
}


# 3. Add metascore -> combine scores from same recipe to obtain new score

# 4. Add metalearner -> combine scores from different learner-dataSetting 
# (but same approxDagSet method to keep hypothesis space the same) 
# eg. use  uniformity p-value score accross learners to decide which to use

# 5.  Rank transformation on score (This is only for 2 hypotheses so must be done after dividing hypothesis space)

rankedDecisionsKeep <- function(scoreMat, ppTab, rank2Func=NULL){
  # scoreMat <- aggScores[[1]][[1]][[1]][[1]]
  # For Ranked decisions
  
  # Transform scores  according to complexity pack chosen
  #print("obtaining rank decision functions")
  if(! is.null(rank2Func)){
    rank2Func <- rep(rank2Func, ncol(scoreMat))
  } else{
    tab <- eval(parse(text=ppTab))
    scoreIds <- sapply(strsplit(colnames(scoreMat), "_"), function(el) el[[1]]) 
    indx <- match(scoreIds, tab$id)
    rank2Func <- tab$rank2Func[indx]
  }
  
  # we add both ways because later we will use one or the other to balance "positive" cases
  #print("applying rank decision functions")
  # for each score we apply rank function
  rnkDec <- sapply(1:dim(scoreMat)[2], function(i){
    # i <- 4
    res <- apply( adrop(scoreMat[,i,, drop=FALSE],2), 2, function(hyps){
      print(hyps)
      res1 <- do.call(rank2Func[i], list(hyps[1], hyps[2]))
      res2 <- do.call(rank2Func[i], list(hyps[2], hyps[1]))
      res <- c(res1, res2)
      names(res) <- c("hyp1pos","hyp0pos")
      return(res)
    })
  }, simplify="array")
  
  rnkDec <- aperm(rnkDec, c(1,3,2))
  dimnames(rnkDec)[[2]] <- dimnames(scoreMat)[[2]]
  names(dimnames(rnkDec))[2] <- names(dimnames(scoreMat))[2]
  #dimnames(rnkDec)[[3]] <- dimnames(scoreMat)[[3]]
  #names(dimnames(rnkDec))[3] <- names(dimnames(scoreMat))[3]
  
  return(rnkDec)
}


rankedDecisions <- function(scoreMat, ppTab, rank2Func=NULL){
  # scoreMat <- aggScores[[1]][[1]][[1]][[1]]
  # For Ranked decisions
  
  # Transform scores  according to complexity pack chosen
  #print("obtaining rank decision functions")
  if(! is.null(rank2Func)){
    rank2Func <- rep(rank2Func, ncol(scoreMat))
  } else{
    tab <- eval(parse(text=ppTab))
    scoreIds <- sapply(strsplit(colnames(scoreMat), "_"), function(el) el[[1]]) 
    indx <- sapply(scoreIds, function(id) which(regexpr(id, as.character(tab$id))==1)[1])#match(scoreIds, tab$id)
    rank2Func <- tab$rank2Func[indx]
  }
  
  # we add both ways because later we will use one or the other to balance "positive" cases
  #print("applying rank decision functions")
  # for each score we apply rank function
  rnkDec <- sapply(1:dim(scoreMat)[2], function(i){
    # i <- 4
    res <- apply( adrop(scoreMat[,i,, drop=FALSE],2), 2, function(hyps){
      #print(hyps)
      res1 <- do.call(rank2Func[i], list(hyps[1], hyps[2]))
      res2 <- do.call(rank2Func[i], list(hyps[2], hyps[1]))
      res <- c(res1, res2)
      names(res) <- c("hyp1pos","hyp0pos")
      return(res)
    })
  }, simplify="array")
  
  rnkDec <- aperm(rnkDec, c(1,3,2))
  dimnames(rnkDec)[[2]] <- dimnames(scoreMat)[[2]]
  names(dimnames(rnkDec))[2] <- names(dimnames(scoreMat))[2]
  #dimnames(rnkDec)[[3]] <- dimnames(scoreMat)[[3]]
  #names(dimnames(rnkDec))[3] <- names(dimnames(scoreMat))[3]
  
  return(rnkDec)
}


# 6: Calculate 2-class performance measures: accuracy, ROC/AUC, TSS, recall, precison, etc

statsArrayToDB <- function(statsArr){
  statss <- melt(statsArr)
  
  aux <- strsplit(as.character(statss$fullRecipe), "\\.")
  statss$dt <- sapply(aux, function(el) el[1])
  statss$adsm <- sapply(aux, function(el) el[2])
  statss$hypSc <- sapply(aux, function(el) el[3])
  statss$cmplxFunc <- sapply(aux, function(el) el[4])
  statss$bipartCrit <- sapply(aux, function(el) el[5])
  statss <- statss[,c("fullRecipe", "dt", "adsm", "hypSc", "cmplxFunc","bipartCrit","measure","statistic","value")]
  return(statss)  
  
}

measures2DagClassDeprecated <- function(rankDecArray, posNeg, ws, n.boot){
  n <- dim(rankDecArray)[1]
  indxMat <- cbind(1:n, posNeg+1)
  pm <- proc.time()
  statss <- apply(rankDecArray, 3, function(mat){
    # mat <- rnkDecsArr[,,1]
    predHard <- (apply(mat,1, which.max)-1)
    predHard[which(!posNeg)] <- (!predHard[which(!posNeg)])*1 
    predSoft <- mat[indxMat]
    contBased <- getStats(posNeg, predHard, ws, n.boot, stat="cntgncyStats")
    rankBased <- getStats(posNeg, predSoft, ws, n.boot, stat="aucStat")
    stats <- cbind(contBased, rankBased)
    return(stats)
  })
  proc.time()-pm # 33 mins for 2000 bootstrap samples for 50 data sets
  
  nmsFullRecipe <- dimnames(statss)[2]
  
  # inner stats matrix which we can't access from outer environment has 7 x 9:
  # 7 stats  5 quantiles, the mean and the original 
  # 10 measures: ccr, mcr, ppp, npp, fpr, fnr, sens, spec, tss, auc
  
  mat <- rnkDecsArr[,,1]
  # higher scores better sow e take the max
  predHard <- (apply(mat,1, which.max)-1)
  predHard[which(!posNeg)] <- (!predHard[which(!posNeg)])*1 
  predSoft <- mat[indxMat]
  contBased <- getStats(posNeg, predHard, ws, n.boot, stat="cntgncyStats")
  rankBased <- getStats(posNeg, predSoft, ws, n.boot, stat="aucStat")
  stats <- cbind(contBased, rankBased)
  
  dim(statss) <- c(dim(stats), dim(statss)[2])
  dimnames(statss) <- c(dimnames(stats), nmsFullRecipe)
  names(dimnames(statss)) <- c("statistic", "measure", "fullRecipe")
  
  return(statss)
}

getStatsBoot <- function(obs, pred, ws, n.boot, stat){
  
  
  data <- cbind(obs=obs, pred=pred, ws=ws)
  boots <- boot(data, statistic=eval(parse(text=stat)), R=n.boot, sim="ordinary", stype="i")
  qs <- apply(boots$t, 2 , quantile, probs=c(0.025, 0.18, 0.5, 0.84, 0.975))
  mu <- apply(boots$t,2, mean)
  dim(mu) <- c(1, ncol(qs))
  stats <- rbind(qs[1:2,,drop=F], mu, qs[3,,drop=F], boots$t0, qs[4:5,,drop=F])
  rownames(stats) <- c("q025", "q18", "mu", "q50", "orig", "q84", "q975")
  return(stats)
}

measures2DagClass <- function(rankDecArray, ws, n.boot){
  n <- dim(rankDecArray)[1]
  
  pm <- proc.time()
  statss <- apply(rankDecArray, 3, function(mat){
    # mat <- rnkDecsArr[,,1]
    
    stats <- getStatsDoubleBoot(mat, ws, n.boot)
    
    return(stats)
  })
  proc.time()-pm # 33 mins for 2000 bootstrap samples for 50 data sets
  
  nmsFullRecipe <- dimnames(statss)[2]
  
  # inner stats matrix which we can't access from outer environment has 7 x 9:
  # 7 stats  5 quantiles, the mean and the original 
  # 10 measures: ccr, mcr, ppp, npp, fpr, fnr, sens, spec, tss, auc
  
  # Run Stats on first full recipe to get dimensions along which we shd reshape stats
  i <- 1
  mat <- rankDecArray[,,i]
  stats <- getStatsDoubleBoot(mat, ws, n.boot)
  
  
  # reshape statss
  
  dim(statss) <- c(dim(stats), dim(statss)[2])
  dimnames(statss) <- c(dimnames(stats), nmsFullRecipe)
  names(dimnames(statss)) <- c("statistic", "measure", "fullRecipe")
  
  #test that re-shaping done appropriately
  #all(statss[,,i]==stats)
  
  return(statss)
}

getStatsDoubleBoot <- function(mat, ws, n.boot){
  
  n <- nrow(mat)
  numPos <- ceil(n/2)
  numNeg <- n - numPos 
  posNegAux <- c(rep(0, numNeg), rep(1,numPos))
  
  n.boot1 <- floor(sqrt(n.boot))
  n.boot2 <- ceiling(n.boot/n.boot1)
  
  set.seed(1234)
  posNeg <- sapply(1:n.boot1, function(i){
    smpl <- sample(n)
    res <- posNegAux[smpl]  
    return(res)
  })
  
  predHardAux <- (apply(mat,1, which.max)-1) 
  
  boots <- sapply(as.data.frame(posNeg), function(obs){
    # obs <- posNeg[,1]
    predHard <- predHardAux
    predHard[which(!obs)] <- (!predHardAux[which(!obs)])*1 
    indxMat <- cbind(1:n, obs+1)
    predSoft <- mat[indxMat]
    
    dataHard <- cbind(obs=obs, pred=predHard, ws=ws)
    dataSoft <- cbind(obs=obs, pred=predSoft, ws=ws)
    contBased <- boot(dataHard, statistic=cntgncyStats, R=n.boot2, sim="ordinary", stype="i")
    rankBased <- boot(dataSoft, statistic=aucStat     , R=n.boot2, sim="ordinary", stype="i")
    
    contBased <- rbind(contBased$t0, contBased$t)
    rankBased <- rbind(rankBased$t0, rankBased$t)
    
    boot <- cbind(contBased, rankBased)

    return(boot)
  }, simplify="array")
  
  
  boots <- aperm(boots, c(1,3,2))
  
  t0 <- boots[1,,] 
  t0 <- apply(t0, 2, mean)
  nmsMsrs <- names(t0)
  dim(t0) <- c(1, length(t0))
  
  boots <- boots[-1,,]
  
  dim(boots) <- c(dim(boots)[1]*dim(boots)[2], dim(boots)[3])
  
  
  
  qs <- apply(boots, 2 , quantile, probs=c(0.025, 0.18, 0.5, 0.84, 0.975))
  mu <- apply(boots, 2, mean)
  dim(mu) <- c(1, ncol(qs))
  stats <- rbind(qs[1:2,,drop=F], mu, qs[3,,drop=F], t0, qs[4:5,,drop=F])
  rownames(stats) <- c("q025", "q18", "mu", "q50", "orig", "q84", "q975")
  colnames(stats) <- nmsMsrs
  return(stats)
}


cntgncyStats <- function(data, indx=1:nrow(data)){
  
  pred <- data[indx,"pred"]
  obs <- data[indx, "obs"]
  ws <- data[indx, "ws"]
  
  
  # don't name the same as something outside otherwise it screws up (in this case)
  tab <- contingencyTable(pred, obs, ws)
  
  ccr <- correctCR(tab)
  mcr <- misCR(tab) 
  ppp <- posPP(tab)
  npp <- negPP(tab)
  sens <- sensitivity(tab)
  spec <- specificity(tab)
  fpr <- fpr(tab)
  fnr <- fnr(tab)
  tss <- tss(tab)
  res <- c(ccr, mcr, ppp, npp, sens, spec, fpr, fnr, tss)
  names(res) <- c("ccr", "mcr", "ppp", "npp", "sens", "spec", "fpr", "fnr", "tss")
  return(res)
} 

aucStat <- function(data, indx=1:nrow(data)){
  res <- WeightedAUC(WeightedROC(data[indx,"pred"], data[indx,"obs"], data[indx,"ws"]))
  names(res) <- "auc"
  return(res)
} 

getAucsDoubleBoot <- function(mat, ws, n.boot){
  
  n <- nrow(mat)
  numPos <- ceil(n/2)
  numNeg <- n - numPos 
  posNegAux <- c(rep(0, numNeg), rep(1,numPos))
  
  n.boot1 <- floor(sqrt(n.boot))
  n.boot2 <- ceiling(n.boot/n.boot1)
  n.boot; n.boot1*n.boot2
  
  set.seed(1234)
  posNeg <- sapply(1:n.boot1, function(i){
    smpl <- sample(n)
    res <- posNegAux[smpl]  
    return(res)
  })
  
  boots <- sapply(as.data.frame(posNeg), function(obs){
    # obs <- posNeg[,1]
    indxMat <- cbind(1:n, obs+1)
    predSoft <- mat[indxMat]
    dataSoft <- cbind(obs=obs, pred=predSoft, ws=ws)
    rankBased <- boot(dataSoft, statistic=aucStat     , R=n.boot2, sim="ordinary", stype="i")
    
    boot <-  c(rankBased$t)
    
    return(boot)
  }, simplify="array")
  
  boots <- c(boots)
  
  return(boots)
}

aucDist <- function(rankDecArray, ws, n.boot){
  n <- dim(rankDecArray)[1]
  
  pm <- proc.time()
  aucss <- apply(rankDecArray, 3, function(mat){
    # mat <- rankDecArray[,,1]
    
    aucs_boots <- getAucsDoubleBoot(mat, ws, n.boot)
    
    return(aucs_boots)
  })
  proc.time()-pm # 33 mins for 2000 bootstrap samples for 50 data sets
  
  dimnames(aucss)[[1]] <- 1:dim(aucss)[1]
  names(dimnames(aucss))[1] <- "boot"
  
  
  return(aucss)
}


#################################################################################################################*
# Miscellaneous

# for pairs plots, to plot smoother and calculate 1-pval_hsic as a dependence strength
panel.hsic <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- 1-dhsic.test(list(x,y))$p.value
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) 
    cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r) 
}


# for pairs plots, to calculate 1-pval_hsic as a dependence strength
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x,y)
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) 
    cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r) 
}

# for pairs plots, to smooth linearly
panel.lin <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.line="red") {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(list(x=x, y=predict(lm(y~x), data=x)), col = col.line)
}

# for pairs plots to add x=y line
panel.xy <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.line="red") {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(list(x=x, y=x), col = col.line)
}
