#!/usr/bin/env Rscript --vanilla

## Comment 10: The rho between GC and HEC could be inflated -----
#split the data in half. Estimate GCs with one half and HECs with the other

#change home directory ----
setwd("/home/albertf/renga011/ExpressionPhenotypeProject/")

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))


#input
i = as.numeric(commandArgs(trailingOnly = TRUE))

# split the data in half ----
set.seed(i)
index_split = sample(1:nrow(genotypesCommonSegregants), size = 0.5*nrow(genotypesCommonSegregants))
                     
# Data for GC ----

trait_GC = traitCommonSegregants_std[index_split,]
expression_GC = expressionCommonSegregants_batchODCorrected_std[index_split,]

# Data for HEC ----

genotypes_HEC = genotypesCommonSegregants[-index_split,]
trait_HEC = traitCommonSegregants_std[-index_split,]
expression_HEC = expressionCommonSegregants_batchODCorrected_std[-index_split,]

#functions to compute GC and HEC ----

computeGC = function(trait_data, expression_data){
  r_expressionGrowth = function(gene, condition){
    traitValues = unlist(trait_data[condition])
    expressionForThisGene = unlist(expression_data[,gene])
    cor_expressionGrowth = cor.test(expressionForThisGene, traitValues)
    return(list(r = cor_expressionGrowth$estimate,
                p = cor_expressionGrowth$p.value))
  }
  
  phenotypicCorrelations = parallel::mclapply(colnames(trait_data),
                                              FUN = function(x){
                                                print(x)
                                                df = parallel::mclapply(unique(eQTL_Albert2018$gene),
                                                                        FUN = r_expressionGrowth,
                                                                        condition = x, 
                                                                        mc.cores = parallelly::availableCores())
                                                df = as.data.frame(do.call(rbind, df))
                                                df = as.data.frame(apply(df, 2, FUN = unlist))
                                                df = as.data.frame(apply(df, 2, FUN = as.numeric))
                                                df$gene = unique(eQTL_Albert2018$gene)
                                                df$q = qvalue::qvalue(df$p)$qvalue
                                                return(df)
                                              }, 
                                              mc.cores = parallelly::availableCores())
  
  names(phenotypicCorrelations) = colnames(trait_data)
  
  phenotypicCorrelations_df = lapply(colnames(trait_data),
                                     FUN = function(x){
                                       df = phenotypicCorrelations[[x]]
                                       df$condition = x
                                       return(df)
                                     })
  
  phenotypicCorrelations_df = as.data.frame(do.call(rbind, phenotypicCorrelations_df))
  
  return(phenotypicCorrelations_df)
}

computeHEC = function(trait_data, expression_data, genotype_data){
  getGrowthEffectsForRandomMarkers = function(randomMarkerSet_i){
    growthCorrected_allConditions = lapply(colnames(trait_data),
                                           FUN = function(condition_i){
                                             growth_corrected = unlist(trait_data[condition_i])
                                             completeCases = complete.cases(growth_corrected)
                                             growth_corrected = growth_corrected[completeCases]
                                             
                                             growthEffects_i = rep(0, length(randomMarkerSet_i))
                                             names(growthEffects_i) = randomMarkerSet_i
                                             
                                             for(i in 1:length(randomMarkerSet_i)){
                                               pMarker = randomMarkerSet_i[i]
                                               
                                               #growth-correlation
                                               genotypeAtpMarker = unlist(genotype_data[pMarker])
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
    names(growthCorrected_allConditions) = colnames(trait_data)
    growthCorrected_allConditions = as.data.frame(growthCorrected_allConditions)
    return(growthCorrected_allConditions)
  }
  
  getExpressionEffectsForRandomMarkers = function(randomMarkerSet_i){
    expressionCorrected_allConditions = parallel::mclapply(colnames(expression_data),
                                                           FUN = function(gene_i){
                                                             print(gene_i)
                                                             expression_corrected = unlist(expression_data[gene_i])
                                                             completeCases = complete.cases(expression_corrected)
                                                             expression_corrected = expression_corrected[completeCases]
                                                             
                                                             expressionEffects_i = rep(0, length(randomMarkerSet_i))
                                                             names(expressionEffects_i) = randomMarkerSet_i
                                                             
                                                             for(i in 1:length(randomMarkerSet_i)){
                                                               pMarker = randomMarkerSet_i[i]
                                                               
                                                               #growth-correlation
                                                               genotypeAtpMarker = unlist(genotype_data[pMarker])
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
    names(expressionCorrected_allConditions) = colnames(expression_data)
    expressionCorrected_allConditions = as.data.frame(expressionCorrected_allConditions)
    return(expressionCorrected_allConditions)
  }
  
  rs_i = hotspotData_Albert2018$hotspotMarker
  
  growthEffects_i = getGrowthEffectsForRandomMarkers(rs_i)
  expressionEffects_i = getExpressionEffectsForRandomMarkers(rs_i)
  
  #Hotspot effect correlation
  hotspotEffectCorrelations = parallel::mclapply(colnames(growthEffects_i),
                                                 FUN = function(condition_i){
                                                   print(condition_i)
                                                   growth_i = unlist(growthEffects_i[condition_i])
                                                   r_allGenes = parallel::mclapply(colnames(expressionEffects_i), FUN = function(gene_i){expression_i = unlist(expressionEffects_i[gene_i])
                                                   corTest = cor.test(growth_i, expression_i, method = "pearson")
                                                   return(c("r" = corTest$estimate,"p" = corTest$p.value))
                                                   }, mc.cores = parallelly::availableCores())
                                                   r_allGenes = as.data.frame(do.call(rbind, r_allGenes))
                                                   r_allGenes$gene = colnames(expressionEffects_i)
                                                   r_allGenes$q = p.adjust(r_allGenes$p, method = "fdr")
                                                   return(r_allGenes)
                                                 },
                                                 mc.cores = parallelly::availableCores())
  
  names(hotspotEffectCorrelations) = colnames(trait_data)
  
  hotspotEffectCorrelations_df = lapply(colnames(trait_data),
                                        FUN = function(x){
                                          df = hotspotEffectCorrelations[[x]]
                                          df$condition = x
                                          return(df)
                                        })
  hotspotEffectCorrelations_df = as.data.frame(do.call(rbind, hotspotEffectCorrelations_df))
  return(hotspotEffectCorrelations_df)
}

#compute GC ----
GC_i = computeGC(trait_GC, expression_GC)

#compute HEC ----
HEC_i = computeHEC(trait_HEC, expression_HEC, genotypes_HEC)

# collate into list and save
GCHEC = list("GC" = GC_i, "HEC" = HEC_i)

# save
save(GCHEC, file = paste0(results_dir, RObj_dir, "GC_HEC_", i, ".rda"))












