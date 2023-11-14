#!/usr/bin/env Rscript --vanilla
#commands required for MSI 
setwd("~/ExpressionPhenotypeProject/")

library(dplyr)
library(qlcMatrix)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load the random marker set ----
load(file = paste0(results_dir, RObj_dir, "randomMarkersForRandomHeritabilityEstimate.rda"))
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

#input
index = as.numeric(commandArgs(trailingOnly = TRUE))
start_index = (index-1)*100 + 1
end_index = index*100
randomMarkers_sub = randomMarkers[start_index:end_index]

#compute eQTL effects and gQTL effects for a random marker set ----
getGrowthEffectsForRandomMarkers = function(randomMarkerSet_i){
  growthCorrected_allConditions = lapply(colnames(traitCommonSegregants_std),
                                         FUN = function(condition_i){
                      growth_corrected = unlist(traitCommonSegregants_std[condition_i])
                      completeCases = complete.cases(growth_corrected)
                      growth_corrected = growth_corrected[completeCases]
                      
                      growthEffects_i = rep(0, length(randomMarkerSet_i))
                      names(growthEffects_i) = randomMarkerSet_i
                      
                      for(i in 1:length(randomMarkerSet_i)){
                        pMarker = randomMarkerSet_i[i]
                        
                        #growth-correlation
                        genotypeAtpMarker = unlist(genotypesCommonSegregants[pMarker])
                        genotypeAtpMarker = genotypeAtpMarker[completeCases]
                        computeRGrowth = cor.test(growth_corrected, genotypeAtpMarker)
                        r_growth = computeRGrowth$estimate
                        p_growth = computeRGrowth$p.value
                        
                        growthEffects_i[pMarker] = ifelse(p_growth < 0.05, r_growth, 0)
                        
                        growthModel =lm(growth_corrected ~ genotypeAtpMarker)
                        growth_corrected = residuals(growthModel)
                      }
                
                        return(growthEffects_i)
                                         })
  names(growthCorrected_allConditions) = colnames(traitCommonSegregants_std)
  growthCorrected_allConditions = as.data.frame(growthCorrected_allConditions)
  return(growthCorrected_allConditions)
}

getExpressionEffectsForRandomMarkers = function(randomMarkerSet_i){
  expressionCorrected_allConditions = parallel::mclapply(colnames(expressionCommonSegregants_batchODCorrected_std),
                                         FUN = function(gene_i){
                                           print(gene_i)
                                           expression_corrected = unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])
                                           completeCases = complete.cases(expression_corrected)
                                           expression_corrected = expression_corrected[completeCases]
                                           
                                           expressionEffects_i = rep(0, length(randomMarkerSet_i))
                                           names(expressionEffects_i) = randomMarkerSet_i
                                           
                                           for(i in 1:length(randomMarkerSet_i)){
                                             pMarker = randomMarkerSet_i[i]
                                             
                                             #growth-correlation
                                             genotypeAtpMarker = unlist(genotypesCommonSegregants[pMarker])
                                             genotypeAtpMarker = genotypeAtpMarker[completeCases]
                                             computeRExpression = cor.test(expression_corrected, genotypeAtpMarker)
                                             r_expression = computeRExpression$estimate
                                             p_expression = computeRExpression$p.value
                                             
                                             expressionEffects_i[pMarker] = ifelse(p_expression < 0.05, r_expression, 0)
                                             
                                             expressionModel =lm(expression_corrected ~ genotypeAtpMarker)
                                             expression_corrected = residuals(expressionModel)
                                           }
                                           
                                           return(expressionEffects_i)
                                         },
                                         mc.cores = parallelly::availableCores())
  names(expressionCorrected_allConditions) = colnames(expressionCommonSegregants_batchODCorrected_std)
  expressionCorrected_allConditions = as.data.frame(expressionCorrected_allConditions)
  return(expressionCorrected_allConditions)
}

#compute "hotspot" effect r for random marker set
hotspotEffectCorrelations = lapply(randomMarkers_sub,
                                   FUN = function(rs_i){
                                     growthEffects_i = getGrowthEffectsForRandomMarkers(rs_i)
                                     expressionEffects_i = getExpressionEffectsForRandomMarkers(rs_i)
                                     effectsCorrelation_i = corSparse(as.matrix(growthEffects_i), as.matrix(expressionEffects_i))
                                     row.names(effectsCorrelation_i) = colnames(traitCommonSegregants_std)
                                     colnames(effectsCorrelation_i) = colnames(expressionCommonSegregants_batchODCorrected_std)
                                     return(effectsCorrelation_i)
                                   })

#save random "hotspot" effect table
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "randomHotspotEffectCorrelationEstimates/"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "randomHotspotEffectCorrelationEstimates/"))), FALSE)
 
save(hotspotEffectCorrelations, file = paste0(results_dir, RObj_dir, "randomHotspotEffectCorrelationEstimates/", "randomHotspotEffectCor_", index, ".rda"))

### compute the expression and growth effects for the actual hotspot markers again and then compare the correlations once again...

hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

rs_i = hotspotData_Albert2018$hotspotMarker

growthEffects_i = getGrowthEffectsForRandomMarkers(rs_i)
expressionEffects_i = getExpressionEffectsForRandomMarkers(rs_i)
effectsCorrelation_i = corSparse(as.matrix(growthEffects_i), as.matrix(expressionEffects_i))
row.names(effectsCorrelation_i) = colnames(traitCommonSegregants_std)
colnames(effectsCorrelation_i) = colnames(expressionCommonSegregants_batchODCorrected_std)

#Hotspot effect correlation
hotspotEffectCorrelations = lapply(colnames(growthEffects_i),
FUN = function(condition_i){
  print(condition_i)
  growth_i = unlist(growthEffects_i[condition_i])
  r_allGenes = parallel::mclapply(colnames(expressionEffects_i),FUN = function(gene_i){expression_i = unlist(expressionEffects_i[gene_i])
  corTest = cor.test(growth_i, expression_i, method = "pearson")
  return(c("r" = corTest$estimate,"p" = corTest$p.value))
  }, mc.cores = parallelly::availableCores())
  r_allGenes = as.data.frame(do.call(rbind, r_allGenes))
  r_allGenes$gene = colnames(expressionEffects_i)
  r_allGenes$q = p.adjust(r_allGenes$p, method = "fdr")
  return(r_allGenes)
})

names(hotspotEffectCorrelations) = colnames(traitCommonSegregants_std)

##save the recomputed hotspot effect correlations, growth and expression effects ---
save(effectsCorrelation_i, expressionEffects_i, growthEffects_i, file = paste0(results_dir, RObj_dir, "hotspotEffectCorrelations_recomputed_onlyMyScanningMethod.rda"))

save(hotspotEffectCorrelations, file = paste0(results_dir, RObj_dir, "hotspotEffectCorrelations.rda"))

##############################

############ ANALYSE HEC ---
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations_recomputed_onlyMyScanningMethod.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#
nGenes_HEC = sapply(colnames(traitCommonSegregants_std),
                    FUN = function(condition_i){
                      HEC_i = hotspotEffectCorrelations[[condition_i]]
                      colnames(HEC_i) = c("r", "p", "gene", "q")
                      nGenes = dplyr::filter(HEC_i, q < 0.05 & p < 0.05)$gene
                      return(length(nGenes))
                    })

##agreement with GC

#nGenes --
ft_GCVsHEC = sapply(colnames(traitCommonSegregants_std),
                    function(condition_i){
                      GC = phenotypicCorrelations[[condition_i]]
                      genes_GC = dplyr::filter(GC, q< 0.05 & p < 0.05)$gene
                      
                      HEC = hotspotEffectCorrelations[[condition_i]]
                      genes_HEC = dplyr::filter(HEC, q<0.05 & p< 0.05)$gene
                      
                      t_t = length(intersect(genes_GC, genes_HEC))
                      t_f = length(setdiff(genes_GC, genes_HEC))
                      f_t = length(setdiff(genes_HEC, genes_GC))
                      f_f = length(colnames(expressionCommonSegregants_batchODCorrected_std)) - length(union(genes_GC, genes_HEC))
                      
                      twoByTwo = matrix(c(t_t, t_f, f_t, f_f), nrow = 2, ncol = 2)
                      ft = fisher.test(twoByTwo)
                      
                      return(c("oddsRatio" = ft$estimate,
                               "p" = ft$p.value))
                    })

ft_GCVsHEC = as.data.frame(t(ft_GCVsHEC))


##genetic vs hotspot effect R
r_GCVsHEC_recompute = sapply(colnames(traitCommonSegregants_std),
                          function(condition_i){
                            GC_i = phenotypicCorrelations[[condition_i]]
                            HC_i = as.data.frame(t(effectsCorrelation_i))
                            HC_i = HC_i[condition_i]
                            HC_i$gene = row.names(HC_i)
                            GC_i$r_HEC = NA
                            GC_i$r_HEC = HC_i[match(GC_i$gene, HC_i$gene), 1]
                            
                            r = cor.test(GC_i$r, GC_i$r_HEC, method = "pearson")$estimate
                            return(r)
                          })

r_GCVsHEC_recompute = as.data.frame(r_GCVsHEC_recompute)
r_GCVsHEC_recompute$condition = colnames(traitCommonSegregants_std)

##


