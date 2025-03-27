#consolidate nulls

library(dplyr)
library(ggpubr)
library(reshape2)

#load data -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "geneticR_afterRegressingOutLocalTransHotspotEffects.rda"))

#genetic correlation after regressing out local, trans, hotspot effects out
for(condition_i in names(geneticR_oldLocalTransHotspotRegressOut)){
  geneticR_oldLocalTransHotspotRegressOut[[condition_i]]$condition = condition_i
}

df_all = as.data.frame(do.call(rbind, geneticR_oldLocalTransHotspotRegressOut))
df_all = dplyr::filter(df_all, r != r_localOut & r != r_transOut)

nGenes = df_all %>% group_by(condition) %>% summarize(n = sum(p < 0.05), n_localOut = sum(p_localOut < 0.05), n_transOut = sum(p_transOut < 0.05), n_hotspotOut = sum(p_hotspotOut < 0.05))

#consolidate the null estimates for each category

consolidateNulls = lapply(colnames(traitCommonSegregants_std),
                          FUN = function(condition_i){
                            if(file.exists(paste0(results_dir, RObj_dir, 
                                                  "nullsForLocalVsTransVsHotspot/",
                                                  "nulls_", condition_i, ".rda"))){
                              load(paste0(results_dir, RObj_dir, 
                                          "nullsForLocalVsTransVsHotspot/",
                                          "nulls_", condition_i, ".rda"))
                              
                              effects_i_nulls = purrr::compact(effects_i_nulls)
                              
                              avg_list = lapply(effects_i_nulls,
                                                FUN = function(df_i){
                                                  r = dplyr::select(df_i, contains("_r_"))
                                                  avg_r = apply(r, 2, mean, na.rm = TRUE)
                                                  return(avg_r)
                                                })
                              
                              avg_df = as.data.frame(do.call(rbind, avg_list))
                              avg_df$condition = condition_i
                              return(avg_df)
                            } else{
                              NULL
                            }
                          })

consolidateNulls_all = as.data.frame(do.call(rbind, consolidateNulls))

colnames(consolidateNulls_all) = c(paste0("r_", c("localOut", "transOut", "hotspotOut"), "_nulls"), "condition")

#melt and prep for plotting

  #actual estimates
df_melt_actuals = df_all %>% select(r, contains("r_"), condition) %>% melt(id.vars = "condition")
df_melt_actuals$value_abs = abs(df_melt_actuals$value)

nGenes_melt = melt(nGenes, id.vars= "condition")

  #null estimates
df_melt_null = consolidateNulls_all %>% melt(id.vars = "condition")
df_melt_null$value_abs = abs(df_melt_null$value)

  #combine and plot
df_melt = rbind(df_melt_actuals, df_melt_null)
df_melt$variable = factor(df_melt$variable, levels = c("r", "r_localOut", "r_localOut_nulls",
                                                       "r_transOut", "r_transOut_nulls",
                                                       "r_hotspotOut", "r_hotspotOut_nulls"))

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


#
consolidateNulls = lapply(colnames(traitCommonSegregants_std),
                          FUN = function(condition_i){
                            if(file.exists(paste0(results_dir, RObj_dir, 
                                                      "nullsForLocalVsTransVsHotspot/",
                                                      "nulls_", condition_i, ".rda"))){
                            load(paste0(results_dir, RObj_dir, 
                                        "nullsForLocalVsTransVsHotspot/",
                                        "nulls_", condition_i, ".rda"))
                              
                            effects_i_nulls = purrr::compact(effects_i_nulls)
                              
                            avg_list = lapply(effects_i_nulls,
                                              FUN = function(df_i){
                                                r = dplyr::select(df_i, contains("_r_"))
                                                avg_r = apply(r, 2, mean, na.rm = TRUE)
                                                return(avg_r)
                                              })
                            
                            avg_df = as.data.frame(do.call(rbind, avg_list))
                            return(avg_df)
                            } else{
                              NULL
                            }
                          })
consolidateNulls_all = as.data.frame(do.call(rbind, consolidateNulls))
summary(abs(consolidateNulls_all))

#
wilcox.test(df_all$r_localOut, consolidateNulls_all$r_localOut_nulls)$p.value
  # p = 0.73

wilcox.test(df_all$r_transOut, consolidateNulls_all$r_transOut_nulls)
  # p = 0

wilcox.test(df_all$r_hotspotOut, consolidateNulls_all$r_hotspotOut_nulls)
  # p = 0
