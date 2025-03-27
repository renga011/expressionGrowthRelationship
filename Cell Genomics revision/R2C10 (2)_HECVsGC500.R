# analyze HEC vs GC data -----
library(dplyr)
library(ggpubr)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations.rda"))

#
GCHEC_all = list()

for(i in 1:100){
    load(paste0(results_dir, RObj_dir, "HECVsGC_5050/", "GC_HEC_", i, ".rda"))
    GCHEC_all[[i]] = GCHEC
}

#
compareGCHEC = function(GC, HEC){
  commonGenes = intersect(GC$gene, HEC$gene)
  GC = dplyr::filter(GC, gene %in% commonGenes)
  GC$HEC = HEC$r.cor[match(paste0(GC$gene,GC$condition), paste0(HEC$gene, HEC$condition))]
  ct = cor.test(GC$r, GC$HEC, method = "spearman")
  return(c("r" = as.numeric(ct$estimate), "p" = ct$p.value))
}

GC_all = lapply(colnames(traitCommonSegregants_std),
                FUN = function(condition_i){
                  df = phenotypicCorrelations[[condition_i]]
                  df$condition = condition_i
                  return(df)
                })
GC_all = as.data.frame(do.call(rbind, GC_all))

HEC_all = lapply(colnames(traitCommonSegregants_std),
                FUN = function(condition_i){
                  df = hotspotEffectCorrelations[[condition_i]]
                  df$condition = condition_i
                  return(df)
                })
HEC_all = as.data.frame(do.call(rbind, HEC_all))

#compare for actual data ----
GCHEC_actual = compareGCHEC(GC_all, HEC_all)

# compare for 50% split data ----
GCHEC_5050 = sapply(GCHEC_all,
                    function(list_i){
                      GC_i = list_i$GC
                      HEC_i = list_i$HEC
                      x = compareGCHEC(GC_i, HEC_i)
                      return(x)
                    })
GCHEC_5050 = t(GCHEC_5050)

#plot histogram ----
pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "hist_GCVsHEC_50pcVsall.pdf"))

gghistogram(as.data.frame(GCHEC_5050),
            x = "r") + 
  geom_vline(xintercept = GCHEC_actual["r"], color = "darkred") +
  xlab("rho") +
  theme_classic()

dev.off()
