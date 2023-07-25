# complexity function packs
print("in ob_complexPacks.R")


##################################################################################################################
# A general complexity and post-processing dictionary from which to build specific ones
##################################################################################################################

entropyFuncs <- c("Shannon_1sp", "Shannon_Gauss","Shannon_kNN_k", "Shannon_spacing_V", "Shannon_spacing_Vb",
                  "Shannon_spacing_Vpconst","Shannon_spacing_Vplin", "Shannon_spacing_Vplin2", "Shannon_spacing_VKDE",
                  "Shannon_spacing_LL", "Shannon_PSD_SzegoT", "Shannon_Edgeworth", "Shannon_MaxEnt1",
                  "Shannon_MaxEnt2", "Shannon_expF", "Shannon_KDP", "Shannon_vME")


entropyFuncsMV <- c("Shannon_kNN_k","Shannon_Edgeworth","Shannon_expF", "Shannon_KDP", "Shannon_vME")

cmplxFuncs <- c("score_sumLogLik","score_pvalHSIC", "score_pHSIC", "score_HSIC", "score_HSIC_fix", 
                "score_pvaldHSIC","score_pdHSIC","score_pdHSICgMin","score_pdHSICgMean" ,"score_dHSIC", "score_dHSIC_fix",
                "score_pvalUnifPart", "score_pvalUnifBoot",
                rep("score_sumMarginalEntropies", length(entropyFuncs)),
                rep("score_entropy", length(entropyFuncsMV)))


cmplxPars <- c(list(list()),lapply(1:5, function(i) list(method="gamma")),
               lapply(1:3, function(i) list(method="gamma")), #permutation
               lapply(1:2, function(i) list(method="gamma")),
               list(list(method="gamma", numParts=10)),
               list(list(method="gamma", numSmpls=100)),
               lapply(entropyFuncs[1:2], function(i) list(type=i)),
               list(list(type="Shannon_kNN_k", k=3)),
               lapply(entropyFuncs[4:length(entropyFuncs)], function(i) list(type=i)),
               list(list(type="Shannon_kNN_k", k=3)),
               lapply(entropyFuncsMV[2:length(entropyFuncsMV)], function(i) list(type=i)))

ids <- c("logLik","pvalHSIC", "pHSIC", "HSIC", "HSIC_fix", 
         "pvaldHSIC","pdHSIC", "pdHSICgMin","pdHSICgMean","dHSIC", "dHSIC_fix",
         "pvalUnifPart", "pvalUnifBoot",
         paste("SME", sapply(strsplit(entropyFuncs[1:2], "_"), function(el) paste(el[2:length(el)], collapse="_")), sep="_"),
         "SME_3NN",
         paste("SME", sapply(strsplit(entropyFuncs[4:length(entropyFuncs)], "_"), function(el) paste(el[2:length(el)], collapse="_")), sep="_"),
         "ent_3NN",
         paste("ent", sapply(strsplit(entropyFuncsMV[2:length(entropyFuncsMV)], "_"), function(el) paste(el[2:length(el)], collapse="_")), sep="_"))

args <- c(list("pred"), lapply(1:6, function(i) "resid"), lapply(1:2, function(i) c("resid","grp")), 
                        lapply(1:(length(ids)-9), function(i) "resid"))

cmplxScores <- mapply(FUN=function(a, b, c, d){
  # i <- 8; a <- ids[i]; b <- cmplxFuncs[i]; c <- cmplxPars[[i]]; d <- args[[i]]
  res <- list(id=a)
  res$func <- b
  res$pars <- c
  res$arg <- d
  return(res)
}, a=ids, b=cmplxFuncs , c=cmplxPars, d=args, SIMPLIFY=FALSE)


names(cmplxScores) <- ids

# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

ppTab <- data.frame(id=unique(sapply(strsplit(ids, "_"), function(el) el[[1]])))

# lets add pUnifHSIC and bestpUnifpHSIC which are functions of scores added after

ppTab <- rbind(ppTab, data.frame(id=c("pUnifpdHSIC","bestpUnifpdHSIC","timeLearn")))



# define rank2Func transf: how t
# logLik          = -sum(log(liks))                  in (0,Inf) -> differ
# pvalHSIC        = (1-pvalHSIC)                     in (0,1)   -> addQuot
# pHSIC           = -log(pvalHSIC)                  in (0,Inf)  -> quot (as in Mooij et al 2016)
# HSIC            = HSIC                            in (0,Inf)  -> addQuot (as in Mooij et al 2016)
# pvaldHSIC       = (1-pvaldHSIC)                   in (0,1),   -> addQuot
# pdHSIC          = -log(pvaldHSIC)                 in (0,Inf)  -> quot (as in Mooij et al 2016)
# dHSIC           = dHSIC                           in (0,Inf)  -> addQuot (as in Mooij et al 2016)
# pvalUnifPart    = -log(pvaldHSIC or pvalKS_dHSIC) in (0,Inf)  -> quot
# pvalUnifBoot    = -log(1-KS_dHSIC)                in (0,Inf)  -> quot 
# SME             = sum(ents)                       in (0,Inf)  -> differ
# ent             = ent                             in (0,Inf)  -> differ
# pUnifpdHSIC     = pUnifPart*pdHSIC                 in (0,Inf) -> quot
# bestpUnifpdHSIC = pHSIC of max diff between 2nd 
#                  lowest and lowest pvalUnifPart 
#                  among diff recipes              in (0,Inf)  -> quot


ppTab$rank2Funcs <- c("quot","addQuot","differ")[c(3,2,1,2,2,1,3,3,2,1,1,3,3,1,1,1)]
ppTab$rankFuncs <- rep("correctScoreToAdd", nrow(ppTab))
ppTab$probFuncs <- rep("scoreToProb", nrow(ppTab))
ppTab$argTypes <- c("negLogLiks", rep("residuals", 5), rep("residuals_grp",2),rep("residuals", 5),"recipe_score_mat", "recipeList_score_array","time")



##################################################################################################################
# Mooij Complexity Pack and post-processing table
##################################################################################################################



#ones from Mooij et al 2016 (which I have available)
names(cmplxScores)
indx <- c(3, 4, 5, 14, 17, 29, 24, 16, 25, 28, 26, 27)
cmplxScorePackMooij <- cmplxScores[indx]
names(cmplxScorePackMooij)
ppTabMooij <- ppTab


# a reduced version
indx <- c(29, 26) #5
cmplxScorePackRed <- cmplxScores[indx]
names(cmplxScorePackRed)
ppTabRed <- ppTab


# a reduced version for kernel quantile reg 
indx <- c(1, 2, 7, 8, 9, 29, 26) #5
cmplxScorePackRedq <- cmplxScores[indx]
names(cmplxScorePackRedq)
ppTabRedq <- ppTab


##################################################################################################################
# General Complexity Pack and post-processing table for flex script
##################################################################################################################

names(cmplxScores)
indx <- c(30, 35)
names(cmplxScores)[indx]
cmplxScorePack1 <- cmplxScores[-indx]
names(cmplxScorePack1)

ppTab1 <- ppTab
