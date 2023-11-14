library(dplyr)
library(ggpubr)
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

regressOutQTLEffectsAndComputePhenoR = function(gene_i, condition_i){
  print(gene_i)
  localEQTLs = dplyr::filter(eQTL_Albert2018, gene== gene_i & cis == "TRUE")
  transEQTLs = dplyr::filter(eQTL_Albert2018, gene== gene_i & cis == "FALSE")
  
  QTL_pmarkerList = c(localEQTLs$pmarker, transEQTLs$pmarker, hotspotData_Albert2018$hotspotMarker)

  geno_QTLs = dplyr::select(genotypesCommonSegregants, QTL_pmarkerList)
  geno_QTLs$growth = unlist(traitCommonSegregants_std[condition_i])
  geno_QTLs$expression = unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])
  
  #genetic r 
  geneticR_old = cor.test(x = geno_QTLs$expression, y = geno_QTLs$growth)
  r_old = as.numeric(geneticR_old$estimate)
  p_old = geneticR_old$p.value
  
  #regress out local-eQTL effects and compute genetic r
  if(nrow(localEQTLs) > 0){
  df_local = dplyr::select(geno_QTLs, localEQTLs$pmarker, growth, expression)
  df_local = df_local[complete.cases(df_local),] #remove NA entries
  
  localOut = regressOutQTLEffects(df_local)
  } else { 
    localOut = c("r_out" = r_old, "p_out" = p_old)
  }
  
  #regress out trans-eQTL effects and compute genetic r
  if(nrow(transEQTLs) > 0){
    df_trans = dplyr::select(geno_QTLs, transEQTLs$pmarker, growth, expression)
    df_trans = df_trans[complete.cases(df_trans),] #remove NA entries
    
    transOut = regressOutQTLEffects(df_trans)
  } else { 
    transOut = c("r_out" = r_old, "p_out" = p_old)
  }
  
  #regress out hotspot effects and compute genetic r
    df_hotspot = dplyr::select(geno_QTLs, hotspotData_Albert2018$hotspotMarker, growth, expression)
    df_hotspot = df_hotspot[complete.cases(df_hotspot),] #remove NA entries
    hotspotOut = regressOutQTLEffects(df_hotspot)
  
  return(c("r_old" = r_old, "p_old" = p_old,
           "r_localOut" = localOut['r_out'], "p_localOut" = localOut['p_out'],
           "r_transOut" = transOut['r_out'], "p_transOut" = transOut['p_out'],
           "r_hotspotOut" = hotspotOut['r_out'], "p_hotspotOut" = hotspotOut['p_out']))
}

geneticR_oldLocalTransHotspotRegressOut = lapply(colnames(traitCommonSegregants_std),
       FUN = function(condition_i){
         print(condition_i)
         sigGenes = dplyr::filter(phenotypicCorrelations[[condition_i]], p < 0.05 & q < 0.05)$gene
         effects_i = parallel::mclapply(sigGenes,
                            FUN = regressOutQTLEffectsAndComputePhenoR,
                            condition_i = condition_i,
                            mc.cores = parallelly::availableCores())
         effects_i = as.data.frame(do.call(rbind, effects_i))
         colnames(effects_i) = c("r", "p", "r_localOut", "p_localOut", "r_transOut", "p_transOut",
                                 "r_hotspotOut", "p_hotspotOut")
         effects_i$gene = sigGenes
         return(effects_i)
       })
names(geneticR_oldLocalTransHotspotRegressOut) = colnames(traitCommonSegregants_std)

save(geneticR_oldLocalTransHotspotRegressOut, file = paste0(results_dir, RObj_dir, "geneticR_afterRegressingOutLocalTransHotspotEffects.rda"))

####PLOT -----
load(paste0(results_dir, RObj_dir, "geneticR_afterRegressingOutLocalTransHotspotEffects.rda"))
for(condition_i in names(geneticR_oldLocalTransHotspotRegressOut)){
  geneticR_oldLocalTransHotspotRegressOut[[condition_i]]$condition = condition_i
}

df_all = as.data.frame(do.call(rbind, geneticR_oldLocalTransHotspotRegressOut))
df_all = dplyr::filter(df_all, r != r_localOut & r != r_transOut)

nGenes = df_all %>% group_by(condition) %>% summarize(n = sum(p < 0.05), n_localOut = sum(p_localOut < 0.05), n_transOut = sum(p_transOut < 0.05), n_hotspotOut = sum(p_hotspotOut < 0.05))

df_melt = df_all %>% select(r, contains("r_"), condition) %>% melt(id.vars = "condition")
df_melt$value_abs = abs(df_melt$value)

nGenes_melt = melt(nGenes, id.vars= "condition")

pdf(paste0(results_dir, plotting_dir, "geneticR_afterRegressingOutLocalTransHotspotEffects_conditionWise.pdf"),
    width = 15, height = 20)

ggboxplot(df_melt,
          x = "variable",
          y = "value_abs",
          facet.by = "condition",
          outlier.shape = NA) + stat_compare_means(ref.group = "r", label = "p.signif", method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab("abs(r)") + scale_y_continuous(limits= c(0, 0.25))


ggbarplot(nGenes_melt,
          x = "variable",
          y = "value",
          facet.by = "condition") + #stat_compare_means(ref.group = "r", "p.signif") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab("nSiginficantGenes")

dev.off()

pdf(paste0(results_dir, plotting_dir, "geneticR_afterRegressingOutLocalTransHotspotEffects.pdf"))

ggboxplot(df_melt,
          x = "variable",
          y = "value_abs",
          outlier.shape = NA) + stat_compare_means(ref.group = "r", label = "p.format", method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab("abs(r)") + scale_y_continuous(limits= c(0, 0.25))

dev.off()

#compute drop in median between different categories
(median(abs(df_all$r_localOut))-median(abs(df_all$r)))/median(abs(df_all$r))
#0.1%
(median(abs(df_all$r_transOut))-median(abs(df_all$r)))/median(abs(df_all$r))
#0.52
(median(abs(df_all$r_hotspotOut))-median(abs(df_all$r)))/median(abs(df_all$r))
#0.32


wilcox.test(abs(df_all$r), abs(df_all$r_localOut), paired = TRUE)
# p = 1e-07

wilcox.test(abs(df_all$r), abs(df_all$r_transOut), paired = TRUE)
# p < 2e-16

wilcox.test(abs(df_all$r), abs(df_all$r_hotspotOut), paired = TRUE)$p.value
# p < 2e-16



