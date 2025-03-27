# analyze QEC vs GC data -----
library(ggpubr)
library(dplyr)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----

load(file = paste0(results_dir, RObj_dir, "QTLEffectsCorrelation_atEQTLs.rda"))
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#
GCQEC_all = list()

for(i in 1:100){
  load(paste0(results_dir, RObj_dir, "QECVsGC_5050/", "GC_QEC_", i, ".rda"))
  GCQEC_all[[i]] = GCQEC
}


#
compareGCQEC = function(GC, QEC){
  commonGenes = intersect(GC$gene, QEC$gene)
  GC = dplyr::filter(GC, gene %in% commonGenes)
  GC$QEC = QEC$r[match(paste0(GC$gene,GC$condition), paste0(QEC$gene, QEC$condition))]
  ct = cor.test(GC$r, GC$QEC, method = "spearman")
  return(c("r" = as.numeric(ct$estimate), "p" = ct$p.value))
}

GC_all = lapply(colnames(traitCommonSegregants_std),
                FUN = function(condition_i){
                  df = phenotypicCorrelations[[condition_i]]
                  df$condition = condition_i
                  return(df)
                })
GC_all = as.data.frame(do.call(rbind, GC_all))

QEC_all = lapply(colnames(traitCommonSegregants_std),
                 FUN = function(condition_i){
                   df = QEC[[condition_i]]
                   df$condition = condition_i
                   return(df)
                 })
QEC_all = as.data.frame(do.call(rbind, QEC_all))

#compare for actual data ----
GCQEC_actual = compareGCQEC(GC_all, QEC_all)

# compare for 50% split data ----
GCQEC_5050 = sapply(GCQEC_all,
                    function(list_i){
                      GC_i = list_i$GC
                      QEC_i = list_i$QEC
                      x = compareGCQEC(GC_i, QEC_i)
                      return(x)
                    })
GCQEC_5050 = t(GCQEC_5050)

#plot histogram ----
pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "hist_GCVsQEC_50pcVsall.pdf"))

gghistogram(as.data.frame(GCQEC_5050),
            x = "r") + 
  geom_vline(xintercept = GCQEC_actual["r"], color = "darkred") +
  xlab("rho") +
  theme_classic()

dev.off()
