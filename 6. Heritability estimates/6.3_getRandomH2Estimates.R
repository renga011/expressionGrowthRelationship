#!/usr/bin/env Rscript --vanilla
setwd("~/ExpressionPhenotypeProject/")

#get estimate of proportion variance explained for any 102 markers in the genome

library(lme4qtl)
library(parallel)
library(parallelly)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "randomMarkersForRandomHeritabilityEstimate.rda"))

#get the argument from cmd line ----
i = as.numeric(commandArgs(trailingOnly = TRUE))
i = i-1

#get the heritability estimate for random 102 marker set
start = i*100 + 1
end = (i+1)*100

randomMarkers = randomMarkers[start:end]

randomH2Estimate = mclapply(1:20, function(n){
  print(n)
  #make GRM
  grm_random = cov(t(dplyr::select(genotypesCommonSegregants, randomMarkers[[n]])))
  colnames(grm_random) = 1:979
  row.names(grm_random) = 1:979 
  
  #proportion of variance explained
  getPropVarExplained = lapply(colnames(traitCommonSegregants_std), function(condition_i){
    print(condition_i)
    
    data1 = data.frame(trait = unlist(traitCommonSegregants_std[condition_i]),
                       R = 1:979)
    
    #get proportion of variance explained
    
    m_random = relmatLmer(trait ~ (1|R), data1, relmat = list(R = grm_random))
    
    #get proportion of variance explained
    R = VarProp(m_random)
    
    return(R$prop[which(R$grp == "R")])  
    
  })
  
  return(unlist(getPropVarExplained))
  
}, mc.cores = parallelly::availableCores())

randomH2Estimates_allConditions = do.call(rbind, randomH2Estimate)

ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "randomH2Estimates"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "randomH2Estimates"))), FALSE)

save(randomH2Estimates_allConditions, file = paste0(results_dir, RObj_dir, "randomH2Estimates/", "h2Estimates_randomMarkers_allConditions_", i, ".rda"))



