
# for cmfm.kcmc check what gammas were fitted to see if 
# the extra regularizer is helping

remove(list=ls())

repos <- "/home/emiliano/Documents/ISP/proyectos/causality/causaLearner_R_pkg/causaLearner"
dir(repos)
setwd(repos)

hs_cmem_ob_version <- "v5_comp"
hs_cmfm_ob_version <- "v5_comp"
source("./pkg_causaLearner/func_causaLearners_v2.R", echo=FALSE)

# KCMC
experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-kcmc-4"
expType <- toupper(strsplit(strsplit(experimentName, "_")[[1]][2], "-")[[1]][1])
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", sep=""))
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, sep=""))
folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, "/",sep="")
files <- dir(folderLearners)
indx <- grep("KCMC", files)
length(indx)
300*4*2 # 300 data sets, 4 learners, 2 causal hypothesis
files <- files[indx]
gammas <- sapply(files, function(file){
  #file <- files[1]
  #print(file)
  load(file=paste(folderLearners, file, sep=""))
  gamma <- cmemLearnerAux$hyperParams$data$optimizable$gamma$val
  return(gamma)
})
table(log(gammas, 10))/length(gammas)*100

# KCSC
experimentName <- "dag2-ME2-Cmplx-SinPlus-SinTimes_cmfm-comp-kcsc-4"
expType <- toupper(strsplit(strsplit(experimentName, "_")[[1]][2], "-")[[1]][1])
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", sep=""))
dir(paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, sep=""))
folderLearners <- paste(repos, "/pkg_causaLearner/experiments/", expType, "/results/learners/", experimentName, "/",sep="")
files <- dir(folderLearners)
indx <- grep("KCSC", files)
length(indx)
300*4*2 # 300 data sets, 4 learners, 2 causal hypothesis
files <- files[indx]
gammas <- sapply(files, function(file){
  #file <- files[1]
  #print(file)
  load(file=paste(folderLearners, file, sep=""))
  gamma <- cmemLearnerAux$hyperParams$data$optimizable$gamma$val
  return(gamma)
})
table(log(gammas, 10))/length(gammas)*100
