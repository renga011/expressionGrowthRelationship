#!/usr/bin/env Rscript --vanilla
setwd("~/ExpressionPhenotypeProject/")

##mediation for the candidate hotspot

library(mediation)
library(readr)
library(parallel)
library(parallelly)
library(dplyr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
otherFiles_dir = "otherFiles_101522/"

##load data ------
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx")) #Albert et al. 2018 - Source data 8
  hotspotData_Albert2018$ID = paste0("h", 1:nrow(hotspotData_Albert2018))

hotspotEffects_expression = read_csv(paste0(results_dir, otherFiles_dir, "hotspotEffectsOnExpression.csv")) #Albert et al. 2018 - Source data 9: the b matrix

##
hotspot = "chrXV:171150_T/C"
hotspotID = "IRA2"
condition = "Hydrogen_Peroxide"
g = unlist(genotypesCommonSegregants[hotspot])
he = unlist(hotspotEffects_expression[hotspot])

##
getMediationEstimatesForGene = function(gene){
  #gene = hotspotEffects_expression$gene[which(hotspotEffects_expression$h88 == mVal)]
  dat$mediator = unlist(expressionCommonSegregants_batchODCorrected_std[gene])
  dat= na.omit(dat)
  
  #x = as.formula(paste0("mediator~", hotspot))
  fit.mediator = lm(mediator~genotype, data = dat)
  
  fit.dv = lm(trait~genotype + mediator, data = dat)
  
  results = mediate(fit.mediator, fit.dv, treat = "genotype", mediator = 'mediator', boot = T)
  mediation_summary = summary(results)
  #print(mediation_summary)
  
  
  return(c("mediated_effect" = mediation_summary$d.avg,
           "mediated_p" = mediation_summary$d.avg.p,
           "direct_effect" = mediation_summary$z.avg,
           "direct_p" = mediation_summary$z.avg.p,
           "total_effect" = mediation_summary$tau.coef,
           "total_p" = mediation_summary$tau.p,
           "propMediated" = mediation_summary$n.avg, 
           "p" = mediation_summary$n.avg.p))
}

trait_i = unlist(traitCommonSegregants_std[condition])
hotspotTargets_i = hotspotEffects_expression$gene[which(he!=0)]
dat = data.frame(genotype = g,
                trait= trait_i)
mediationEstimatesForAllHotspotTargets = mclapply(hotspotTargets_i, 
                                                FUN = getMediationEstimatesForGene,
                                                mc.cores = availableCores())
mediationEstimatesForAllHotspotTargets = as.data.frame(do.call(rbind, mediationEstimatesForAllHotspotTargets))

mediationEstimatesForAllHotspotTargets$gene = hotspotTargets_i

mediationEstimatesForAllHotspotTargets$fdr = p.adjust(mediationEstimatesForAllHotspotTargets$p, method = "fdr")

##save
ifelse(!dir.exists(paste0(results_dir, RObj_dir, "mediationForHotspotTargets_examples/")),
       dir.create(paste0(results_dir, RObj_dir, "mediationForHotspotTargets_examples/")),
       FALSE)

save(mediationEstimatesForAllHotspotTargets, file = paste0(results_dir, RObj_dir, hotspotID, "__", condition, ".rda"))
