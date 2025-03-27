#!/usr/bin/env Rscript --vanilla
setwd("~/ExpressionPhenotypeProject/")

##compute nulls for local eQtl, trans eQTL and hotspot regress out of GC analysis

library(dplyr)
library(ggplot2)
library(reshape2)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

#get the argument from cmd line ----
i = as.numeric(commandArgs(trailingOnly = TRUE))
condition_i = colnames(traitCommonSegregants_std)[i]

#functions -----
regressOutQTLEffects = function(df_i){
  model_growth = lm(growth~.-expression, data = df_i)
  growth_res = residuals(model_growth)
  
  model_expression = lm(expression~.-growth, data= df_i)
  expression_res = residuals(model_expression)
  
  geneticR_regressedOut = cor.test(x = expression_res, growth_res, method = "pearson")
  r_out = as.numeric(geneticR_regressedOut$estimate)
  p_out = geneticR_regressedOut$p.value
  
  return(c("r_out" = r_out, "p_out" = p_out))
} #df_i is data frame with the genotype vectors, growth and expression vectors

sampleRandomMarkerSetOfGivenSize = function(n){
  marker_set = sample(colnames(genotypesCommonSegregants), size = n)
  return(marker_set)
}

regressOutQTLEffectsAndComputePhenoR_nulls = function(gene_i, condition_i){
  print(gene_i)
  localEQTLs = dplyr::filter(eQTL_Albert2018, gene== gene_i & cis == "TRUE")
  transEQTLs = dplyr::filter(eQTL_Albert2018, gene== gene_i & cis == "FALSE")
  
  growth = unlist(traitCommonSegregants_std[condition_i])
  expression = unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])
  
  #regress out local-eQTL effects and compute genetic r
  if(nrow(localEQTLs) > 0){
    ##local null: 500 set
    localOut_null_500 = lapply(1:500, 
                               function(i){
                                 set.seed(i)
                                 randomMarkerSet_i = sampleRandomMarkerSetOfGivenSize(nrow(localEQTLs))
                                 random_QTLs = dplyr::select(genotypesCommonSegregants, randomMarkerSet_i)
                                 random_QTLs$growth = growth
                                 random_QTLs$expression = expression
                                 random_QTLs = random_QTLs[complete.cases(random_QTLs),]
                                 
                                 randomOut = regressOutQTLEffects(random_QTLs)
                               })
    localOut_null_500 = as.data.frame(do.call(rbind,localOut_null_500))
  } else{
    localOut_null_500 = as.data.frame(matrix(NA, nrow = 500, ncol = 2))
  }
  
  colnames(localOut_null_500) = c("local_r_out", "local_p_out")
  
  #regress out trans-eQTL effects and compute genetic r
  if(nrow(transEQTLs) > 0){
    ##trans null: 500 set
    transOut_null_500 = lapply(1:500, 
                               function(i){
                                 set.seed(i)
                                 randomMarkerSet_i = sampleRandomMarkerSetOfGivenSize(nrow(transEQTLs))
                                 random_QTLs = dplyr::select(genotypesCommonSegregants, randomMarkerSet_i)
                                 random_QTLs$growth = growth
                                 random_QTLs$expression = expression
                                 random_QTLs = random_QTLs[complete.cases(random_QTLs),]
                                 
                                 randomOut = regressOutQTLEffects(random_QTLs)
                               })
    transOut_null_500 = as.data.frame(do.call(rbind,transOut_null_500))
  } else {
    transOut_null_500 = as.data.frame(matrix(NA, nrow = 500, ncol = 2))
  }
  
  colnames(transOut_null_500) = c("trans_r_out", "trans_p_out")
  
  ##hotspot null: 500 set
  hotspotOut_null_500 = lapply(1:500, 
                               function(i){
                                 set.seed(i)
                                 randomMarkerSet_i = sampleRandomMarkerSetOfGivenSize(102)
                                 random_QTLs = dplyr::select(genotypesCommonSegregants, randomMarkerSet_i)
                                 random_QTLs$growth = growth
                                 random_QTLs$expression = expression
                                 random_QTLs = random_QTLs[complete.cases(random_QTLs),]
                                 
                                 randomOut = regressOutQTLEffects(random_QTLs)
                               })
  hotspotOut_null_500 = as.data.frame(do.call(rbind,hotspotOut_null_500))
  colnames(hotspotOut_null_500) = c("hotspot_r_out", "hotspot_p_out")
  
  nulls_df = cbind(localOut_null_500, transOut_null_500, hotspotOut_null_500)
  
  return(nulls_df)
}


#compute regressed out nulls ------

print(condition_i)
sigGenes = dplyr::filter(phenotypicCorrelations[[condition_i]], p < 0.05 & q < 0.05)$gene

effects_i_nulls = parallel::mclapply(sigGenes,
                                     FUN = regressOutQTLEffectsAndComputePhenoR_nulls,
                                     condition_i = condition_i,
                                     mc.cores = parallelly::availableCores())

names(effects_i_nulls) = sigGenes

if(!dir.exists(paste0(results_dir, RObj_dir, "nullsForLocalVsTransVsHotspot/"))){
  dir.create(paste0(results_dir, RObj_dir, "nullsForLocalVsTransVsHotspot/"))
} else {FALSE}

#save nulls for condition_i

save(effects_i_nulls, file = paste0(results_dir, RObj_dir, 
                                    "nullsForLocalVsTransVsHotspot/",
                                    paste0("nulls_", condition_i, ".rda")))