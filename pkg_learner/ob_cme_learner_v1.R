# Construct Conditional Mean Embedding Measures



cmemPack1 <- list()

measurePack1 <- list(KCDC=list(func="KCDC", aggFunc="sum",pars=list(lambda=10^-5)),
                     KCDCrel=list(func="KCDCrel", aggFunc="prod", pars=list(lambda=10^-5, numPerms=100)))


kernelPack1 <- list(kernelsX=list(rbfMed=list(func="kern_rbf", setParams="medSig"), 
                                  log=list(func="kern_log", setParams="deg1"),
                                  quad=list(func="kern_quad", setParams="off1")), 
                    kernelsY=list(rbfMed=list(func="kern_rbf", setParams="medSig"), 
                                  log=list(func="kern_log", setParams="deg1"),
                                  quad=list(func="kern_quad", setParams="off1")))

# Lets build post-processing table indicating what the post-processing (aggregation, ranking, etc) should be for each
# type of score function. The first part of the ID (up to the first underscore) should be enough to determine this

ppTab <- data.frame(id=c("KCDC","KCDCrel"))




# define rank2Func transf: how t
# KCDC              =                     in (0,Inf)   -> addQuot
# KCDCrel           =                     in (0,Inf)   -> addQuot


ppTab$rank2Funcs <- c("quot","addQuot","differ")[c(2,2)]
ppTab$rankFuncs <- rep("correctScoreToAdd", 2)
ppTab$probFuncs <- rep("scoreToProb", 2)
ppTab$argTypes <- c(rep("cmes", 2))


ppTabCMEM1 <- ppTab