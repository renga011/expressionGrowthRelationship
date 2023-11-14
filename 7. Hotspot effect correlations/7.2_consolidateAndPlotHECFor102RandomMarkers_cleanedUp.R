library(dplyr)
library(ggpubr)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"

#load the random marker set ----
load(file = paste0(results_dir, RObj_dir, "randomMarkersForRandomHeritabilityEstimate.rda"))
load(file = paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(file= paste0(results_dir, RObj_dir, "hotspotEffectCorrelations.rda"))
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations_recomputed_onlyMyScanningMethod.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#append all random hotspot r into one list of 1000, each entry corresponding to one random marker
listOfAllRandomHotspotR = list()

for(index in 1:10){
  print(index)
  load(file = paste0(results_dir, RObj_dir, "randomHotspotEffectCorrelationEstimates/", "randomHotspotEffectCor_", index, ".rda"))
  
  listOfAllRandomHotspotR = c(listOfAllRandomHotspotR, hotspotEffectCorrelations)
}

#get median hotspot effect correlation per condition for each random marker set
medianHotspotEffect_i = lapply(listOfAllRandomHotspotR, 
                               FUN = function(set_i){
                                 set_i = t(set_i) #change to dim 5720 X 46, where 46 is number of conditions
                                 a = apply(set_i, 2, median, na.rm = TRUE)
                                 names(a) = colnames(set_i)
                                 return(a)
                               })

medianHotspotEffect_i = do.call(rbind, medianHotspotEffect_i)

#plot the median hotspot effect - actual vs random hotspots like in heritability plot
medianHotspotEffect_actual = lapply(hotspotEffectCorrelations, 
                                    FUN = function(condition_set){
                                      colnames(condition_set) = c("r", "p", "gene", "q")
                                      return(median(condition_set$r, na.rm = TRUE))
                                    })
medianHotspotEffect_actual = as.data.frame(do.call(rbind, medianHotspotEffect_actual))
medianHotspotEffect_actual$condition = row.names(medianHotspotEffect_actual)
colnames(medianHotspotEffect_actual) = c("median_hotspotEffectR", "condition")


#
randomMedianHotspotCor = reshape2::melt(medianHotspotEffect_i)
colnames(randomMedianHotspotCor) = c("id", "condition", "median_hotspotEffectR")

pdf(file = paste0(results_dir, plotting_dir, "medianHotspotEffectCorrelations_withRandomHotspotEstimates.pdf"), width = 10, height = 8)

ggplot() +
  geom_boxplot(randomMedianHotspotCor,
               mapping = aes(x = condition, y = median_hotspotEffectR), color = "darkgrey", 
               alpha = 0.5, outlier.alpha = 0) +
  geom_point(medianHotspotEffect_actual,
             mapping = aes(x = condition, y = median_hotspotEffectR), color = "darkred") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Median Hotspot effect correlation")

dev.off()

#genetic vs hotspot effect correlation
##get correlation of genetic vs hotspot effect correlation for random marker sets - condition wise

corGeneticVsHotspot = lapply(colnames(traitCommonSegregants_std),
                             FUN = function(condition_i){
                               print(condition_i)
                               GC_i_df = phenotypicCorrelations[[condition_i]]
                               HCVsGC = sapply(1:length(listOfAllRandomHotspotR),
                                               function(i){
                                                 #print(i)
                                                 rs_i = listOfAllRandomHotspotR[[i]]
                                                 HC_random_i = rs_i[,GC_i_df$gene] #common genes
                                                 HC_random_i = t(HC_random_i)
                                                 HC_i = unlist(HC_random_i[,condition_i])
                                                 GC_i = unlist(GC_i_df$r)
                                                 names(GC_i) = GC_i_df$gene
                                                 if(length(which(is.finite(HC_i) == "TRUE"))!=0){
                                                 r = cor.test(HC_i, GC_i, 
                                                              method = "pearson")$estimate
                                                 } else { r = NA}
                                                 return(r)
                                               })
                               return(HCVsGC)
                             })
names(corGeneticVsHotspot) = colnames(traitCommonSegregants_std)
corGeneticVsHotspot = as.data.frame(corGeneticVsHotspot)
colnames(corGeneticVsHotspot) = colnames(traitCommonSegregants_std)

##plot
melt_GCVsHEC_random = reshape2::melt(corGeneticVsHotspot)
colnames(melt_GCVsHEC_random) = c("condition", "GCVsHEC")

r_GCVsHEC_actual = sapply(colnames(traitCommonSegregants_std),
                          function(condition_i){
                            GC_i = phenotypicCorrelations[[condition_i]]
                            HC_i = hotspotEffectCorrelations[[condition_i]]
                            GC_i$r_HEC = NA
                            GC_i$r_HEC = HC_i$r.cor[match(GC_i$gene, HC_i$gene)]
                            
                            r = cor.test(GC_i$r, GC_i$r_HEC, method = "pearson")$estimate
                            return(r)
                          })
r_GCVsHEC_actual = as.data.frame(r_GCVsHEC_actual)
r_GCVsHEC_actual$condition = colnames(traitCommonSegregants_std)


#get p-value
for(i in 1:nrow(r_GCVsHEC_actual)){
  r_i = r_GCVsHEC_actual$r_GCVsHEC_actual[i]
  condition_i = r_GCVsHEC_actual$condition[i]
  
  n = length(which(unlist(corGeneticVsHotspot[condition_i]) >= r_i))
  r_GCVsHEC_actual$p[i] = n/nrow(corGeneticVsHotspot)
  
  max_i = max(unlist(corGeneticVsHotspot[condition_i]))
  r_GCVsHEC_actual$maxRandom[i] = max_i
}

#
pdf(file = paste0(results_dir, plotting_dir, "GeneticVsHotspotEffectCorrelation_withRandomEstimates_recomputed.pdf"), width = 10, height = 8)

ggplot() +
  geom_boxplot(melt_GCVsHEC_random,
               mapping = aes(x = condition, y = GCVsHEC), color = "darkgrey", 
               alpha = 0.5, outlier.alpha = 0) +
  geom_point(r_GCVsHEC_actual,
             mapping = aes(x = condition, y = r_GCVsHEC_actual), color = "darkred") +
  geom_text(r_GCVsHEC_actual,
            mapping = aes(x = condition, y = maxRandom + 0.05, label = p), size = 3, angle = 90) +
  geom_text(r_GCVsHEC_actual,
            mapping = aes(x = condition, y = maxRandom + 0.05, 
                          label = ifelse(p < 0.05, "*", "")), vjust = -0.75, size = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Genetic vs Hotspot effect correlation")

dev.off()