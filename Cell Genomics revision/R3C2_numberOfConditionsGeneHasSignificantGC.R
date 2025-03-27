## COMMENT R3C2: number of conditions where a gene has significant GC
library(ggpubr)
library(dplyr)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#
phenotypicCorrelations = lapply(colnames(traitCommonSegregants_std),
                                FUN = function(condition_i){
                                  df = phenotypicCorrelations[[condition_i]]
                                  df$condition = condition_i
                                  return(df)
                                })
names(phenotypicCorrelations) = colnames(traitCommonSegregants_std)
GC_all = rbind(do.call(rbind, phenotypicCorrelations))

table(GC_all$gene)

  #filter significant entriess only at FDR 5% and then count how many times each gene appears
nSignificantConditions_geneWise = GC_all %>% filter(p < 0.05 & q< 0.05) %>% group_by(gene) %>% summarise(nConditions = n())

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "histogramOfNConditionsGeneHasSignificantGC.pdf"))

gghistogram(nSignificantConditions_geneWise,
            x = "nConditions") + theme_classic() + 
  xlab("number of growth conditions gene has significant genetic correlation") +
  ylab("number of genes")

dev.off()

write.csv(nSignificantConditions_geneWise,
          file = paste0(results_dir, otherFiles_dir, "supplementaryTables/",
                        "revision_tables_CellGenomics/",
                        "R3C2_numberOfConditionsGeneHasSignificantGC.csv"), 
          quote = FALSE, row.names = FALSE)
