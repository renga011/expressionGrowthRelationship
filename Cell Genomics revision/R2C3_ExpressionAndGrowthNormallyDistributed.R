# Reviewer 2, Comment 3: are growth and expression distributions normally distributed? 
## Response: Show the distributions of growth values for all 46 conditions and expression values for 50 genes 

library(ggpubr)
# library(car)
# library(cowplot)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

# load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

# reshape data ----
growth_molten = reshape2::melt(traitCommonSegregants_std)
colnames(growth_molten) = c("condition", "value")

set.seed(1234)
colNos = sample(1:ncol(expressionCommonSegregants_batchODCorrected_std),size = 100)
expression_molten = reshape2::melt(expressionCommonSegregants_batchODCorrected_std[, colNos])
colnames(expression_molten) = c("gene", "value")

#plot distributions

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/growthValuesDistributions_qqPlot.pdf"), width = 10, height = 10)

ggplot(data = growth_molten,
       aes(sample = value)) +
  stat_qq(shape = 21, alpha = 0.5) + facet_wrap(~condition, ncol = 5) +
  geom_abline(slope = 1, color = "darkred") +
  ylab("colony size growth (standardized normal values)") +
  theme_bw()

dev.off()

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/expressionValuesDistributions_qqPlot.pdf"), width = 20, height = 20)

ggplot(data = expression_molten, 
       aes(sample = value)) +
  stat_qq(shape = 21, alpha = 0.5) + facet_wrap(~gene, ncol = 10) +
  geom_abline(slope = 1, color = "darkred") +
  ylab("mRNA expression values (standardized normal values)") +
  theme_bw()

dev.off()

