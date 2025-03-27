## see how much heritability is captured by top X hotspots

#devtools::install_github("variani/lme4qtl")

library(lme4qtl)
library(parallel)
library(parallelly)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

#
hotspotData_Albert2018 = hotspotData_Albert2018[order(hotspotData_Albert2018$numberNonzeroEffects,
                                                      decreasing = TRUE),]

h2_topNHotspots = lapply(colnames(traitCommonSegregants_std),
                         FUN = function(condition_i){
                           print(condition_i)
                           trait_i = unlist(traitCommonSegregants_std[condition_i])
                           h2_topNHotspots_i = sapply(2:102, FUN = function(n){
                             print(paste0("n = ", n))
  data1 = data.frame(trait = trait_i,
                     Tn = 1:979)
  
  hotspotList = hotspotData_Albert2018$hotspotMarker[1:n]
  grm_topNHotspots = cov(t(dplyr::select(genotypesCommonSegregants, hotspotList)
  ))
  colnames(grm_topNHotspots) = 1:979
  row.names(grm_topNHotspots) = 1:979
  
  ##mixed model
  m_topNHotspots = relmatLmer(trait ~ (1|Tn), data1, relmat = list(Tn = grm_topNHotspots))
  
  ##proportion variance explained
  Tn = VarProp(m_topNHotspots)
  
  return(Tn$prop[which(Tn$grp == "Tn")])
})
                           df_topN_i = data.frame(n = 2:102,
                                                  h2 = h2_topNHotspots_i)
                           return(df_topN_i)
                         })

names(h2_topNHotspots) = colnames(traitCommonSegregants_std)

## plot the h2 trails for all conditions

h2_topNHotspots = lapply(colnames(traitCommonSegregants_std),
                         FUN = function(condition_i){
                           df = h2_topNHotspots[[condition_i]]
                           df$condition = condition_i
                           return(df)
                         })
names(h2_topNHotspots) = colnames(traitCommonSegregants_std)
  
h2_df = as.data.frame(do.call(rbind, h2_topNHotspots))

save(h2_topNHotspots, file = paste0(results_dir, RObj_dir, "propVarExplained_topNHotspots.rda"))



# load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations_recomputed_onlyMyScanningMethod.rda"))
# 
# hotspotEffectsOnExpression_boolean = mutate(as.data.frame(t(expressionEffects_i)), across(where(is.numeric), ~ +as.logical(.x)))
# 
# nGenes_nonZeroEffects_nHotspots = sapply(2:102, function(n){
#   hotspotList = hotspotData_Albert2018$hotspotMarker[1:n]
#   df_hotspotEffects = dplyr::select(hotspotEffectsOnExpression_boolean, hotspotList)
#   df_hotspotEffects$geneEffectsVector_boolean = apply(df_hotspotEffects, 1, FUN = sum)
#   nGenes_nonZeroEffects_boolean = length(which(df_hotspotEffects$geneEffectsVector_boolean != 0))
#   
#   #df_hotspotEffects$geneEffectsVector = apply(df_hotspotEffects, 1, FUN = sum)
#   #nGenes_nonZeroEffects = length(which(df_hotspotEffects$geneEffectsVector != 0))
#   
#   #df = data.frame(nGenes_boolean = nGenes_nonZeroEffects_boolean,
#   #                nGenes = nGenes_nonZeroEffects)
#   return(nGenes_nonZeroEffects_boolean)
# })
# 
# propGenes_nonZeroEffects_nHotspots = nGenes_nonZeroEffects_nHotspots/5720
# df_propGenes = data.frame(n = 2:102,
#                           propGenesWithEffects = propGenes_nonZeroEffects_nHotspots)



###### plotting ####

h2_halfMax = h2_df %>% group_by(condition) %>% summarise(halfMax = 0.5 * max(h2), eightyMax = 0.8 * max(h2))

load(paste0(results_dir, RObj_dir, "propVarExplained_topNHotspots.rda"))

pdf(paste0(results_dir, plotting_dir, "proportionVarianceExplainedByHotspots.pdf"),
    width = 12, height = 12)
ggscatter(data = h2_df,
          x = "n",
          y = "h2",
          facet.by = "condition",
          size = 0.5,
          xlab = "Number of hotspots",
          ylab = "Proportion of variance explained") + 
  geom_hline(data = h2_halfMax,
             mapping = aes(yintercept = halfMax), color = "deeppink3", alpha = 0.5, width = 3) +
  geom_hline(data = h2_halfMax, aes(yintercept = eightyMax), color = "deeppink3", alpha = 0.8, width = 3) +
  #geom_point(data = df_propGenes, aes(x = n, y = propGenesWithEffects), color = "darkblue", size = 0.5) + 
  theme_minimal()

dev.off()


# nHotspots_q2 = sapply(h2_topNHotspots, 
#                       FUN = function(df_i){
#                         max_h2 = max(df_i$h2)
#                         #range_h2 = range(df_i$h2)[2] - range(df_i$h2)[1]
#                         halfMax = max_h2*0.50
#                         max80 = max_h2*0.80
#                         n_q2_50 = min(df_i$n[which(df_i$h2 >= halfMax)])
#                         n_q2_80 = min(df_i$n[which(df_i$h2 >= max80)])
#                         
#                         return(c("50% increase" = n_q2_50, "80% increase" = n_q2_80))
#                       })
# 
# nHotspots_q2 = as.data.frame(t(nHotspots_q2))
# 
#   #some summary statistics
#   summary(nHotspots_q2)
# 
# pdf(paste0(results_dir, plotting_dir, "nHotspotsForHalfMaxH2_histogram.pdf"))
# hist(nHotspots_q2)
# dev.off()



## correlate proportion of h2 contributed by hotspot and number of genes affected by it

load(paste0(results_dir, RObj_dir, "propVarExplained_topNHotspots.rda"))

hotspots = hotspotData_Albert2018[order(hotspotData_Albert2018$numberNonzeroEffects,
                                                      decreasing = TRUE),]
hotspots = select(hotspots, hotspotMarker, numberNonzeroEffects)

h2_diffBetweenHotspots = lapply(colnames(traitCommonSegregants_std),
                         function(condition_i){
                           df = h2_topNHotspots[[condition_i]]
                           diff_hotspots = diff(df$h2) #difference between consecutive rows
                           diff_hotspots = c(NA, NA, diff_hotspots) #adding 2 NAs for 1st & 2nd hotspot (descending order)
                           return(diff_hotspots)
                         })
names(h2_diffBetweenHotspots) = colnames(traitCommonSegregants_std)
h2_diffBetweenHotspots = as.data.frame(h2_diffBetweenHotspots)
h2_diffBetweenHotspots_mean = as.data.frame(apply(h2_diffBetweenHotspots, 1, mean))
colnames(h2_diffBetweenHotspots_mean) = "meanH2Contribution"

h2_df = cbind(hotspots, h2_diffBetweenHotspots_mean)

pdf(paste0(results_dir, plotting_dir, "h2ContributionByHotspotsVsNumberOfGenesAffected.pdf"))

ggplot(h2_df, aes(x = log10(numberNonzeroEffects), y = log10(meanH2Contribution)))+
  geom_point(shape = 21) +
  xlab("log10(Number of genes affected by hotspot)") +
  ylab("log10(Average h2 explained by hotspot)") +
  stat_cor(method = "spearman") +
  theme_classic()

dev.off()

