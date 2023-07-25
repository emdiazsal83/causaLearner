# construct hypothesis scorers
print("ob_hypothesisScorers.R")

# they have: data_reg, other_pars and function parameters

# other_pars in case of model based hoyer-style are: learner, complexity_pack


# A hypothesis scorerer object should have

# 1. a hypothesis scoring function taking in a hypothesis (dag) list, data and giving back a score matrix
# one row per hyothesis one per (different flavour) score
# 2. a post-processing table (ppTab) giving a unique id to the hypothesis scorer (one per col of score matrix)
# and also indicates how post-processsing should be done: aggregating, ranking, other transformations
# 3. other parameters needed by hypothesis scoring function
#   a. ANMs: dataRegiment, learner, residualComplexity functions
#   b. LVMs: ???
#   c. CMEMs: dep-vars (x) kernel, indep-vars (y) kernel, measures (KCDC, KCSC, etc)
#   d. DMs: ???

# ANMs: additive noise model (based on Mooij et al 2016), also includes PNL-style modeling

print("loading learner objects")
source("./pkg_learner/ob_reg_learners.R")
print("loading complexity objects")
source("./pkg_cmplxtyScores/ob_complexPacks.R")
        
# random

rnd1 <- list(func="rnd_hypScorer", ppTab="ppTabRnd")

# Additive: y = f(x) + err
krr1_re <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="krr1", complexityPack="cmplxScorePack1")

krr2_re_red <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="krr2", complexityPack="cmplxScorePackRed")

krr2_ho_red <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="holdout", learner="krr2", complexityPack="cmplxScorePackRed")

krr2_re <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="krr2", complexityPack="cmplxScorePackMooij")

krr2_ho <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="holdout", learner="krr2", complexityPack="cmplxScorePackMooij")

krr2_re_full <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="krr2", complexityPack="cmplxScorePack1")

krr2_ho_full <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="holdout", learner="krr2", complexityPack="cmplxScorePack1")

qhsic_re <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="qhsic", complexityPack="cmplxScorePack1")

qhsic_ho <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="holdout", learner="qhsic", complexityPack="cmplxScorePack1")

gptk_re <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="recycle", learner="gptk2", complexityPack="cmplxScorePack1")

gptk_ho <- list(func="anm_hypScorer", ppTab="ppTab1", dataReg="holdout", learner="gptk2", complexityPack="cmplxScorePack1")

gptk_re2 <- list(func="anm_hypScorer", ppTab="ppTabMooij", dataReg="recycle", learner="gptk2", complexityPack="cmplxScorePackMooij")

gptk_ho2 <- list(func="anm_hypScorer", ppTab="ppTabMooij", dataReg="holdout", learner="gptk2", complexityPack="cmplxScorePackMooij")

# Non additive: g(y) = f(x) + err

# PNL: post non linear = y = g(f(x) + err) (g invertible)

# LVMs: latent variable models (based on Mooij et al 2010): y = f(x, e)
mooij2010 <- list(func="lvm_hypScorer", ppTab="", dataReg="", priors="", optim="")

# CMEM: conditional mean embedding measures (based on Mitrovic et al 2018)

print("loading CMEM and CMFM objects")
source(paste("./pkg_learner/ob_cme_learner_", hs_cmem_ob_version, ".R", sep=""))

#cmem1 <- list(func="cmem_hypScorer", ppTab="ppTabCMEM1", kernelPack="kernelPack1", measurePack="measurePack1")
# cmem1 <- list(func="cmem_hypScorer", ppTab="ppTabCMEM1", cmemLearner="cmem_rbf1")

cmem_hypScorer_pack_none_1       <- paste("hs", cmem_learner_pack_none_1, sep="_")
cmem_hypScorer_pack_lambda_1     <- paste("hs", cmem_learner_pack_lambda_1, sep="_") 
cmem_hypScorer_pack_kernParsX_1  <- paste("hs", cmem_learner_pack_kernParsX_1, sep="_")
cmem_hypScorer_pack_kernParsXY_1 <- paste("hs", cmem_learner_pack_kernParsXY_1, sep="_")

if(length(cmem_learner_pack_none_1)>0) for(i in 1:length(cmem_learner_pack_none_1)) assign(cmem_hypScorer_pack_none_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_none_1[i]))
if(length(cmem_learner_pack_lambda_1)>0) for(i in 1:length(cmem_learner_pack_lambda_1)) assign(cmem_hypScorer_pack_lambda_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_lambda_1[i]))
if(length(cmem_learner_pack_kernParsX_1)>0) for(i in 1:length(cmem_learner_pack_kernParsX_1)) assign(cmem_hypScorer_pack_kernParsX_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_kernParsX_1[i]))
if(length(cmem_learner_pack_kernParsXY_1)>0) for(i in 1:length(cmem_learner_pack_kernParsXY_1)) assign(cmem_hypScorer_pack_kernParsXY_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_kernParsXY_1[i]))

# comparison cmem hyothesis scorers
# lambda_kernParsX_lambda_kernParsX_L2 - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))

cmem_hypScorer_comp_nn1 <- list(func="cmem_hypScorer_comp_nn", ppTab="ppTabCMEM1") 

# lambda_kernParsXY_lambda_kernParsX_L2 - choose y-param with suvervision in 1st round

lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsX_lambdaGamma_kernParsX_KCMC  - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCMC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsXY_lambdaGamma_kernParsX_KCMC - choose y-param with suvervision in 1st round

lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCMC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsX_lambdaGamma_kernParsX_KCSC  - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCSC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))

# lambda_kernParsXY_lambdaGamma_kernParsX_KCSC - choose y-param with suvervision in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCSC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)>0) 
  for(i in 1:length(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)) 
    assign(cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMEM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))



cmemJ_hypScorer_pack_none_1       <- paste("hs", sub("cmem","cmemJ",cmem_learner_pack_none_1), sep="_")
cmemJ_hypScorer_pack_lambda_1     <- paste("hs", sub("cmem","cmemJ",cmem_learner_pack_lambda_1), sep="_") 
cmemJ_hypScorer_pack_kernParsX_1  <- paste("hs", sub("cmem","cmemJ",cmem_learner_pack_kernParsX_1), sep="_")
cmemJ_hypScorer_pack_kernParsXY_1 <- paste("hs", sub("cmem","cmemJ",cmem_learner_pack_kernParsXY_1), sep="_")

if(length(cmem_learner_pack_none_1)>0) for(i in 1:length(cmem_learner_pack_none_1)) assign(cmemJ_hypScorer_pack_none_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_none_1[i]))
if(length(cmem_learner_pack_lambda_1)>0) for(i in 1:length(cmem_learner_pack_lambda_1)) assign(cmemJ_hypScorer_pack_lambda_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_lambda_1[i]))
if(length(cmem_learner_pack_kernParsX_1)>0) for(i in 1:length(cmem_learner_pack_kernParsX_1)) assign(cmemJ_hypScorer_pack_kernParsX_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_kernParsX_1[i]))
if(length(cmem_learner_pack_kernParsXY_1)>0) for(i in 1:length(cmem_learner_pack_kernParsXY_1)) assign(cmemJ_hypScorer_pack_kernParsXY_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMEM1", cmemLearner=cmem_learner_pack_kernParsXY_1[i]))


# CMFM: conditional mean feature measures (based on Mitrovic et al 2018)

print("loading CMFM objects")
print(hs_cmfm_ob_version)
print(getwd())
source(paste("./pkg_learner/ob_cmf_learner_", hs_cmfm_ob_version, ".R", sep=""))
print(cmfm_learner_pack_none_1)
cmfm_hypScorer_pack_none_1       <- paste("hs", cmfm_learner_pack_none_1, sep="_")
cmfm_hypScorer_pack_lambda_1     <- paste("hs", cmfm_learner_pack_lambda_1, sep="_") 
cmfm_hypScorer_pack_kernParsX_1  <- paste("hs", cmfm_learner_pack_kernParsX_1, sep="_")
cmfm_hypScorer_pack_kernParsXY_1 <- paste("hs", cmfm_learner_pack_kernParsXY_1, sep="_")

if(length(cmfm_learner_pack_none_1)>0) for(i in 1:length(cmfm_learner_pack_none_1)) assign(cmfm_hypScorer_pack_none_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_none_1[i]))
if(length(cmfm_learner_pack_lambda_1)>0) for(i in 1:length(cmfm_learner_pack_lambda_1)) assign(cmfm_hypScorer_pack_lambda_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_lambda_1[i]))
if(length(cmfm_learner_pack_kernParsX_1)>0) for(i in 1:length(cmfm_learner_pack_kernParsX_1)) assign(cmfm_hypScorer_pack_kernParsX_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_kernParsX_1[i]))
if(length(cmfm_learner_pack_kernParsXY_1)>0) for(i in 1:length(cmfm_learner_pack_kernParsXY_1)) assign(cmfm_hypScorer_pack_kernParsXY_1[i], list(func="cmem_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_kernParsXY_1[i]))

# comparison cmfm hyothesis scorers
# lambda_kernParsX_lambda_kernParsX_L2 - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambda_kernParsX_L2_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))

# lambda_kernParsXY_lambda_kernParsX_L2 - choose y-param with suvervision in 1st round

lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambda_kernParsX_L2_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsX_lambdaGamma_kernParsX_KCMC  - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCMC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCMC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsXY_lambdaGamma_kernParsX_KCMC - choose y-param with suvervision in 1st round

lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCMC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))


# lambda_kernParsX_lambdaGamma_kernParsX_KCSC  - choose y-param with heuristic in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss == "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCSC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2 <- strsplit(lrn2, "_")
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsX_L2_lambdaGamma_kernParsX_KCSC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))

# lambda_kernParsXY_lambdaGamma_kernParsX_KCSC - choose y-param with suvervision in 1st round
lrn1 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsXY" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_L2")]
lrn2 <- kernelTab1$learnerName[which(kernelTab1$optParms=="lambda_kernParsX" & kernelTab1$mainLoss != "cmem_L2_f" & kernelTab1$betaLearns == "learnBlambda_KCSC")]
aux1 <- strsplit(lrn1, "_")
aux1_2 <- sapply(aux1, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
aux1 <- sapply(aux1, function(el) paste(el[1:(length(el)-3)], collapse="_"))
aux2_2 <- sapply(aux2, function(el) paste(el[(length(el)-2):length(el)], collapse="_"))
cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1 <- paste("hs", aux1, aux1_2, aux2_2, sep="_")
if(length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)>0) 
  for(i in 1:length(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1)) 
    assign(cmfm_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCSC_1[i], 
           list(func="cmem_hypScorer_comp", ppTab="ppTabCMFM1", 
                cmemLearner1=lrn1[i], cmemLearner2=lrn2[i]))



# Joint hypScorer

cmfmJ_hypScorer_pack_none_1       <- paste("hs", sub("cmfm", "cmfmJ", cmfm_learner_pack_none_1), sep="_")
cmfmJ_hypScorer_pack_lambda_1     <- paste("hs", sub("cmfm", "cmfmJ", cmfm_learner_pack_lambda_1), sep="_") 
cmfmJ_hypScorer_pack_kernParsX_1  <- paste("hs", sub("cmfm", "cmfmJ", cmfm_learner_pack_kernParsX_1), sep="_")
cmfmJ_hypScorer_pack_kernParsXY_1 <- paste("hs", sub("cmfm", "cmfmJ", cmfm_learner_pack_kernParsXY_1), sep="_")

if(length(cmfm_learner_pack_none_1)>0) for(i in 1:length(cmfm_learner_pack_none_1)) assign(cmfmJ_hypScorer_pack_none_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_none_1[i]))
if(length(cmfm_learner_pack_lambda_1)>0) for(i in 1:length(cmfm_learner_pack_lambda_1)) assign(cmfmJ_hypScorer_pack_lambda_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_lambda_1[i]))
if(length(cmfm_learner_pack_kernParsX_1)>0) for(i in 1:length(cmfm_learner_pack_kernParsX_1)) assign(cmfmJ_hypScorer_pack_kernParsX_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_kernParsX_1[i]))
if(length(cmfm_learner_pack_kernParsXY_1)>0) for(i in 1:length(cmfm_learner_pack_kernParsXY_1)) assign(cmfmJ_hypScorer_pack_kernParsXY_1[i], list(func="cmemJoint_hypScorer", ppTab="ppTabCMFM1", cmemLearner=cmfm_learner_pack_kernParsXY_1[i]))


# DM: distribution classifier (based on David Lopez-Paz et al 2015)

lopezPaz2015 <- list(func="dm_hypScorer", ppTab="", dataReg="", distClassif="", motherDist="")

