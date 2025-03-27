library(dplyr)
library(ggpubr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load datasets----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "h2Estimates_allVsHotspot_allConditions.rda"))
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))
load(paste0(results_dir, RObj_dir, "h2_se.rda"))

#consolidate the random h2 estimates across conditions ----
randomH2Estimates = data.frame()
for(i in 1:10){
  load(paste0(results_dir, RObj_dir, "randomH2Estimates/", "h2Estimates_randomMarkers_allConditions_", i, ".rda"))
  randomH2Estimates = rbind(randomH2Estimates, randomH2Estimates_allConditions)
}

colnames(randomH2Estimates) = propVarExplained_allConditions_actual$condition

#compute p-value ----
for(i in 1:nrow(propVarExplained_allConditions_actual)){
  h_i = propVarExplained_allConditions_actual$H[i]
  condition_i = propVarExplained_allConditions_actual$condition[i]
  
  n = length(which(unlist(randomH2Estimates[condition_i]) >= h_i))
  propVarExplained_allConditions_actual$p[i] = n/1000
}

#add bootstrapped se to the same table as h2 estimates ----
row.names(se_h2) = se_h2$condition
se_h2 = select(se_h2, -condition)
colnames(se_h2) = paste0(colnames(se_h2), "_se")
propVarExplained_allConditions_actual = cbind(propVarExplained_allConditions_actual, se_h2)

#plot----
  #prepare data
randomH2Estimates_molten = reshape2::melt(randomH2Estimates)
colnames(randomH2Estimates_molten) = c("condition", "h2")

propVarExplained_allConditions_actual = dplyr::arrange(propVarExplained_allConditions_actual, desc(A))

pdf(file = paste0(results_dir, plotting_dir, "h2_allVsHotspots_allConditions_withRandomH2Estimates.pdf"), width = 10, height = 8)
ggplot() +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = reorder(condition, -A),y = A, color = "darkred")) +
   geom_errorbar(propVarExplained_allConditions_actual,
                 mapping = aes(x = reorder(condition, -A), ymin= A - A_se , ymax=A + A_se), 
                 color = "darkred", alpha = 0.5) +
  geom_boxplot(randomH2Estimates_molten,
               mapping = aes(x = condition, y = h2), outlier.shape = NA, color = "grey") +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = condition, y = H, color = "darkblue")) +
   geom_errorbar(propVarExplained_allConditions_actual,
                 mapping = aes(x = condition, ymin= H - H_se , ymax= H + H_se),
                 color = "darkblue", alpha = 0.5) +
  # geom_point(propVarExplained_allConditions_actual,
  #            mapping = aes(x = condition, y = allTrans, color = "darkcyan")) +
  #  geom_errorbar(propVarExplained_allConditions_actual,
  #                mapping = aes(x = condition, ymin= allTrans-allTrans_se ,ymax= allTrans + allTrans_se),
  #                color = "darkcyan", alpha = 0.5) +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = condition, y = G, color = "gold3")) +
   geom_errorbar(propVarExplained_allConditions_actual,
                 mapping = aes(x = condition, ymin= G - G_se ,ymax= G + G_se),
                color = "gold3", alpha = 0.5) +
  geom_text(propVarExplained_allConditions_actual,
            mapping = aes(x = condition, y = A + 0.05, label = p), size = 3, angle = 90) +
  geom_text(propVarExplained_allConditions_actual,
            mapping = aes(x = condition, y = A + 0.05, 
                          label = ifelse(p < 0.05, "*", "")), vjust = -2, size = 5) +
  ylab("Additive heritability (h2)") +
  scale_color_identity(name = "h2 explained by:",
                       breaks = c("darkred", "darkblue", "darkcyan", "gold3"),
                      labels = c("All genetic markers", "trans-eQTL hotspots", "all trans-eQTLs", "all gQTLs of condition"),
                      guide = "legend") +
  guides(color = guide_legend(title = "h2 explained by:", 
                              title.position = "left")) +
  theme_bw() +
  theme_textProperties +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1))
dev.off()

## reviewer figure --- 5A with all trans
  
pdf(file = paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/h2_allVsHotspots_allConditions_withRandomH2Estimates.pdf"), width = 10, height = 8)
ggplot() +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = reorder(condition, -A),y = A, color = "darkred")) +
  geom_errorbar(propVarExplained_allConditions_actual,
                mapping = aes(x = reorder(condition, -A), ymin= A - A_se , ymax=A + A_se), 
                color = "darkred", alpha = 0.5) +
  geom_boxplot(randomH2Estimates_molten,
               mapping = aes(x = condition, y = h2), outlier.shape = NA, color = "grey") +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = condition, y = H, color = "darkblue")) +
  geom_errorbar(propVarExplained_allConditions_actual,
                mapping = aes(x = condition, ymin= H - H_se , ymax= H + H_se),
                color = "darkblue", alpha = 0.5) +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = condition, y = allTrans, color = "darkcyan")) +
  geom_errorbar(propVarExplained_allConditions_actual,
                mapping = aes(x = condition, ymin= allTrans-allTrans_se ,ymax= allTrans + allTrans_se),
                color = "darkcyan", alpha = 0.5) +
  geom_point(propVarExplained_allConditions_actual,
             mapping = aes(x = condition, y = G, color = "gold3")) +
  geom_errorbar(propVarExplained_allConditions_actual,
                mapping = aes(x = condition, ymin= G - G_se ,ymax= G + G_se),
                color = "gold3", alpha = 0.5) +
  geom_text(propVarExplained_allConditions_actual,
            mapping = aes(x = condition, y = A + 0.05, label = p), size = 3, angle = 90) +
  geom_text(propVarExplained_allConditions_actual,
            mapping = aes(x = condition, y = A + 0.05, 
                          label = ifelse(p < 0.05, "*", "")), vjust = -2, size = 5) +
  ylab("Additive heritability (h2)") +
  scale_color_identity(name = "h2 explained by:",
                       breaks = c("darkred", "darkblue", "darkcyan", "gold3"),
                       labels = c("All genetic markers", "trans-eQTL hotspots", "all trans-eQTLs", "all gQTLs of condition"),
                       guide = "legend") +
  guides(color = guide_legend(title = "h2 explained by:", 
                              title.position = "left")) +
  theme_bw() +
  theme_textProperties +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1))
dev.off()
