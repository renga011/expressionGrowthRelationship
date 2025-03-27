## Reviewer 2, Comment 1: Batch-OD covariate 

library(ggpubr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load files ----

covariates_2018 = readr::read_csv(paste0(results_dir, otherFiles_dir,"covariates_2018.csv")) #AB18 - source data 2

#

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "batchVsOD.pdf"))
ggplot(data = covariates_2018) +
  geom_boxplot(aes(x = as.factor(batch), y = OD_covariate)) +
  ylab("Growth covariate") +
  xlab("Batch") +
  theme_bw()

ggplot(data = covariates_2018, aes(x = batch, y = OD_covariate)) +
  geom_point(shape = 21, alpha = 0.6) +
  ylab("Growth covariate") +
  xlab("Batch") +
  stat_cor(method = "spearman") +
  scale_x_continuous(breaks = seq(1,13, by = 1)) +
  theme_bw()

dev.off()

