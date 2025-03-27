#!/usr/bin/env Rscript --vanilla
setwd("home/albertf/renga011/ExpressionPhenotypeProject/")

#get the h2 explained by hotspots for each of the 46 conditions 

devtools::install_github("variani/lme4qtl")

library(lme4qtl)
library(parallel)
library(parallelly)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#get the argument from cmd line ----
i = as.numeric(commandArgs(trailingOnly = TRUE))
i = i-1

start = i*100 + 1
end = (i+1)*100

#load data
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

#function to make GRMs and get h2 estimates ------

#Function to get the estimates of h2
geth2 = function(GRM, condition, traitCommonSegregants_std_i){
  data1 = data.frame(trait = unlist(traitCommonSegregants_std_i[condition]), X = 1:979)
  model = relmatLmer(trait ~ (1|X), data1, relmat = list(X = GRM))
  #proportion of variance explained
  
  X = VarProp(model)
  
  return(X$prop[which(X$grp == "X")])
}

getH2Estimates = function(genotypeCommonSegregants_i, traitCommonSegregants_std_i){
  #grm for hotspot markers
  grm_h = cov(t(dplyr::select(genotypeCommonSegregants_i, hotspotData_Albert2018$hotspotMarker)))
  colnames(grm_h) = 1:979
  row.names(grm_h) = 1:979
  
  #grm all markers
  grm = cov(t(genotypeCommonSegregants_i))
  colnames(grm)= 1:979
  row.names(grm)= 1:979
  
  #grm all trans-eQTL markers
  eQTL_transMarkers = unique(eQTL_Albert2018$pmarker[which(eQTL_Albert2018$cis=="FALSE")])
  grm_trans = cov(t(dplyr::select(genotypeCommonSegregants_i, eQTL_transMarkers)))
  colnames(grm_trans) = 1:979
  row.names(grm_trans) = 1:979
  
  #grm all cis-eQTL markers
  eQTL_cisMarkers = unique(eQTL_Albert2018$pmarker[which(eQTL_Albert2018$cis=="TRUE")])
  grm_cis = cov(t(dplyr::select(genotypeCommonSegregants_i, eQTL_cisMarkers)))
  colnames(grm_cis) = 1:979
  row.names(grm_cis) = 1:979
  
  #MAKE MIXED MODEL AND GET PROPORTION OF VAR EXPLAINED BY THE GRMs -----
  getPropVarExplained = parallel::mclapply(colnames(traitCommonSegregants_std_i), 
                                           function(condition_i){
                                             print(condition_i)
                                             
                                             #all genotype markers ---
                                             A = geth2(grm, condition_i, traitCommonSegregants_std_i)
                                             
                                             #all hotspots ---
                                             H = geth2(grm_h, condition_i, traitCommonSegregants_std_i)
                                             
                                             #all gQTLs ---
                                             ##make a grm
                                             gQTLs = dplyr::filter(BloomQTL, Trait == condition_i)
                                             genotypeMarkers = stringr::str_split(colnames(genotypeCommonSegregants_i), pattern = "_", simplify = TRUE)[,1]
                                             col_i = which(genotypeMarkers %in% gQTLs$growthQTL)
                                             genotypeMatrix = genotypeCommonSegregants_i[col_i]
                                             
                                             grm_g = cov(t(genotypeMatrix))
                                             colnames(grm_g) = 1:979
                                             row.names(grm_g) = 1:979
                                             
                                             ##mixed model
                                             G = geth2(grm_g, condition_i, traitCommonSegregants_std_i)
                                             
                                             ##top 3 hotspots ---
                                             ##make grm
                                             grm_top3Hotspots = cov(t(dplyr::select(genotypeCommonSegregants_i,
                                                                                    `chrXIV:466588_T/G`, #MKT1
                                                                                    `chrXII:657022_T/C`, #HAP1
                                                                                    `chrXV:171150_T/C`, #IRA2
                                             )))
                                             colnames(grm_top3Hotspots) = 1:979
                                             row.names(grm_top3Hotspots) = 1:979
                                             
                                             ##mixed model
                                             T3 = geth2(grm_top3Hotspots, condition_i, traitCommonSegregants_std_i)
                                             
                                             ## all trans eQTLs ---
                                             allTrans = geth2(grm_trans, condition_i, traitCommonSegregants_std_i)
                                             
                                             ## all cis eQTLs ---
                                             allCis = geth2(grm_cis, condition_i, traitCommonSegregants_std_i)
                                             
                                             return(c("A" = A, 
                                                      "H" = H, 
                                                      "allTrans" = allTrans, 
                                                      "allCis" = allCis, 
                                                      "G" = G, 
                                                      "T3" = T3))    
                                             
                                           }, mc.cores = parallelly::availableCores())
  propVarExplained_allConditions = as.data.frame(do.call(rbind, getPropVarExplained))
  propVarExplained_allConditions$condition  = colnames(traitCommonSegregants_std_i)
  return(propVarExplained_allConditions)
}

#get h2 estimates for actual genotype and trait matrix
if(i == 0){ #so that this operation is done only once.
propVarExplained_allConditions_actual = getH2Estimates(genotypesCommonSegregants, traitCommonSegregants_std)

#save 
save(propVarExplained_allConditions_actual, 
     file = paste0(results_dir, RObj_dir, "h2Estimates_allVsHotspot_allConditions.rda"))
}

#make bootstrapped matrices
propVarExplained_allConditions_boot = mclapply(start:end, function(i){
  set.seed(i)
  print(i)
  boot_indices = sample(1:nrow(genotypesCommonSegregants), replace = TRUE, 
                        size = nrow(genotypesCommonSegregants))
  genotypes_i = genotypesCommonSegregants[boot_indices,]
  traits_i = traitCommonSegregants_std[boot_indices,]
  propVarExplained_allConditions_i = getH2Estimates(genotypes_i, traits_i)
  return(propVarExplained_allConditions_i)
})

#save and plot ---
save(propVarExplained_allConditions_boot,
     file = paste0(results_dir, RObj_dir, "h2Estimates_boot_", i, ".rda"))
