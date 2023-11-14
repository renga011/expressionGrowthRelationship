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

#consolidate the random h2 estimates across conditions ----
randomH2Estimates = data.frame()
for(i in 1:10){
  load(paste0(results_dir, RObj_dir, "randomH2Estimates/", "h2Estimates_randomMarkers_allConditions_", i, ".rda"))
  randomH2Estimates = rbind(randomH2Estimates, randomH2Estimates_allConditions)
}

colnames(randomH2Estimates) = propVarExplained_allConditions$condition

#compute p-value ----
for(i in 1:nrow(propVarExplained_allConditions)){
  h_i = propVarExplained_allConditions$H[i]
  condition_i = propVarExplained_allConditions$condition[i]
  
  n = length(which(unlist(randomH2Estimates[condition_i]) >= h_i))
  propVarExplained_allConditions$p[i] = n/1000
}

#plot----
randomH2Estimates_molten = reshape2::melt(randomH2Estimates)
colnames(randomH2Estimates_molten) = c("condition", "h2")

propVarExplained_allConditions = dplyr::arrange(propVarExplained_allConditions, desc(A))

pdf(file = paste0(results_dir, plotting_dir, "h2_allVsHotspots_allConditions_withRandomH2Estimates.pdf"), width = 10, height = 8)
ggplot() +
  geom_point(propVarExplained_allConditions,
             mapping = aes(x = reorder(condition, -A), y = A, color = "darkred")) +
  geom_boxplot(randomH2Estimates_molten,
               mapping = aes(x = condition, y = h2), outlier.shape = NA, color = "grey") +
  geom_point(propVarExplained_allConditions,
             mapping = aes(x = condition, y = H, color = "darkblue")) +
  geom_point(propVarExplained_allConditions,
             mapping = aes(x = condition, y = G, color = "gold3")) +
  geom_text(propVarExplained_allConditions,
            mapping = aes(x = condition, y = A + 0.05, label = p), size = 3, angle = 90) +
  geom_text(propVarExplained_allConditions,
            mapping = aes(x = condition, y = A + 0.05, 
                          label = ifelse(p < 0.05, "*", "")), vjust = -2, size = 5) +
  ylab("Additive heritability (h2)") +
  scale_color_identity(name = "h2 explained by:",
                       breaks = c("darkred", "darkblue", "gold3"),
                      labels = c("All genetic markers", "trans-eQTL hotspots", "all gQTLs of condition"),
                      guide = "legend") +
  guides(color = guide_legend(title = "h2 explained by:", 
                              title.position = "left")) +
  theme_bw() +
  theme_textProperties +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1))
dev.off()
  