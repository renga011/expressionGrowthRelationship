## Reviewer 2, Comment 4: is there a relationship between expression level and detecting a cis eQTL
## get the batch and OD corrected expression levels and compare for genes with at least one cis-eQTL and those without even 1 cis eQTL using a wilcoxon test.

library(dplyr)
library(ggpubr)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#get expression data which is not standardized but batch and OD corrected because the standardized measures are both positive and negative which starts to become weird for the wilcoxon test.

expression1000segregants = readr::read_delim(paste0(results_dir, otherFiles_dir, "SI_Data_01_expressionValues.txt"), "\t", escape_double = FALSE, trim_ws = TRUE) # AB18 - source data 1
colnames(expression1000segregants)[1] = "wellID"

covariates_2018 = readr::read_csv(paste0(results_dir, otherFiles_dir,"covariates_2018.csv")) #AB18 - source data 2
colnames(covariates_2018)[1] = "wellID"

#batch and OD correct expression of 1012 segregants ---
batch = covariates_2018$batch
OD = covariates_2018$OD_covariate
expression1000segregants = select(expression1000segregants, -wellID)

expression1000Segregants_batchODCorrected = sapply(expression1000segregants,
                                                   FUN = function(x){
                                                     model = lm(x ~0 + batch + OD) #0 ensures that the model doesn't correct for the intercept - values stay positive.
                                                     x_batchODCorrected = residuals(model)
                                                     return(x_batchODCorrected)
                                                   })
expression1000Segregants_batchODCorrected = as.data.frame(expression1000Segregants_batchODCorrected)

genes_cis = eQTL_Albert2018 %>% group_by(gene) %>% summarize(nCis = sum(cis == "TRUE"))
genes_cis$atleastOneCis = ifelse(genes_cis$nCis !=0, TRUE, FALSE)
for(i in 1:nrow(genes_cis)){
  gene_i = genes_cis$gene[i]
  genes_cis$correctedExpression_avg[i] = median(unlist(expression1000Segregants_batchODCorrected[gene_i]), na.rm = TRUE)
  genes_cis$expression_avg[i] = median(unlist(expression1000segregants[gene_i]), na.rm = TRUE)
}

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/presenceOfCisQTLVsExpressionLevel.pdf"))

ggviolin(genes_cis,
         x = "atleastOneCis",
         y = "correctedExpression_avg") + theme_bw() + stat_compare_means(method = "wilcox.test") + xlab("has at least one detectable local eQTL in Albert et al. (2018)") + ylab("median corrected mRNA expression (across 1012 segregants)")

ggviolin(genes_cis,
         x = "atleastOneCis",
         y = "expression_avg") + theme_bw() + stat_compare_means(method = "wilcox.test") + xlab("has at least one detectable local eQTL in Albert et al. (2018)") + ylab("median uncorrected mRNA expression (across 1012 segregants)")

dev.off()

#get medians of both distributions
genes_cis %>% group_by(atleastOneCis) %>% summarize(median_correctedExpression = median(correctedExpression_avg), median_uncorrectedExpression = median(expression_avg))

