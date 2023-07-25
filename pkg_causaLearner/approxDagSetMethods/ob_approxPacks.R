# DAG approximation packs


# approx methods should take data x and pars 




exhaustive1 <- list(func="exhaustive", pars=list())

oracleME1 <- list(func="oracleME", pars=list())

timeSeriesXYZ_ts1 <- list(func="timeSeriesXYZ_ts", pars=list())
timeSeriesXYZ_cmpr <- list(func="timeSeriesXYZ_cmpr", pars=list())

pcSkeleton1 <- list(func="pcSkeleton", pars=list(ic.method="hsic.gamma", indepTest="kernelCItest", alpha=0.02))

pcME1 <- list(func="pcMarkovEquiv", pars=list(ic.method="hsic.gamma", indepTest="kernelCItest", alpha=0.01))

minDAG_gptk_ae <- list(func="minDAG", pars=list(learner="gptk2", dataRegime="recycle", method="permutation", resolution="addEdges"))

minDAG_gptk_dag <- list(func="minDAG", pars=list(learner="gptk2", dataRegime="recycle", method="permutation", resolution="dag"))

minDAG_gptk_co <- list(func="minDAG", pars=list(learner="gptk2", dataRegime="recycle", method="permutation", resolution="co"))

minDAG_krr2_dag <- list(func="minDAG", pars=list(learner="krr2", dataRegime="recycle", method="permutation", resolution="dag"))

minDAG_krr2_co <- list(func="minDAG", pars=list(learner="krr2", dataRegime="recycle", method="permutation", resolution="co"))


# msrCO - causal ordering
# msrPD - parent deletion
# ae - add edges
# dag - take off edges 
# co - causal ordering dag

# minDAG -cmem

minDAG_KCDC_HSIC_ae <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                      msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="addEdges"))


minDAG_KCDC_HSIC_dag <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="dag"))

minDAG_KCDC_HSIC_co <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="co"))


minDAG_KCMC_HSIC_ae <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="addEdges"))


minDAG_KCMC_HSIC_dag <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="dag"))

minDAG_KCMC_HSIC_co <- list(func="minDAG_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="co"))


# pairwise -cmem


pairwise_KCDC_HSIC_max <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=FALSE))

pairwise_KCDC_HSIC_min <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="min", withinMEC=FALSE))


pairwise_KCMC_HSIC_max <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                           msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=FALSE))

pairwise_KCMC_HSIC_min <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                          msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="min", withinMEC=FALSE))


pairwise_cmemComp_KCDC_HSIC_min <- list(func="pairwise_cmemComp", pars=list(hypScorer=cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[3],
                                                                   ppTab=ppTabCMEM1, msrCO="KCDC_reglr_comp",  msrPD="HSIC_cmem", numPerms=100, resolution="min"))
pairwise_cmemComp_KCDC_HSIC_max <- list(func="pairwise_cmemComp", pars=list(hypScorer=cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[3],
                                                                            ppTab=ppTabCMEM1, msrCO="KCDC_reglr_comp",  msrPD="HSIC_cmem", numPerms=100, resolution="max"))

pairwise_cmemComp_KCMC_HSIC_min <- list(func="pairwise_cmemComp", pars=list(hypScorer=cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[3],
                                                                            ppTab=ppTabCMEM1, msrCO="KCMC_reglr_comp",  msrPD="HSIC_cmem", numPerms=100, resolution="min"))
pairwise_cmemComp_KCMC_HSIC_max <- list(func="pairwise_cmemComp", pars=list(hypScorer=cmem_hypScorer_pack_compCombos_lambda_kernParsXY_L2_lambdaGamma_kernParsX_KCMC_1[3],
                                                                            ppTab=ppTabCMEM1, msrCO="KCMC_reglr_comp",  msrPD="HSIC_cmem", numPerms=100, resolution="max"))




pairwise_KCDC_HSIC_max_MEC <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                               msrCO="KCDC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=TRUE))


pairwise_KCMC_HSIC_max_MEC <- list(func="pairwise_cmem", pars=list(cmemLearner="cmem_kern_quadANDkern_quad_nn_L2_none",
                                                               msrCO="KCMC", msrPD="HSIC_cmem", numPerms=100, resolution="max", withinMEC=TRUE))


