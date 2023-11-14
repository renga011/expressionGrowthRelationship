#plotting and analyses for phenotypic correlations----
library(dplyr)
library(ggpubr)
library(ggpattern)

#global variables ---

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
plotting_dir = "plots_092522/"
RObj_dir = "RObjects_092522/"

#create directories ---
ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "phenotypicCorrelations/"))), 
       dir.create(file.path(paste0(results_dir, plotting_dir, "phenotypicCorrelations/"))), FALSE)


#functions -----
fisherTestTable = function(true_true, true_false, false_true, false_false){
  x = t(matrix(c(true_true, true_false, false_true, false_false), nrow = 2, ncol = 2))
  return(x)
}
getFisherP = function(fisherTable){return(fisher.test(fisherTable)$p.value)}

#data -----

load(paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))
load(paste0(results_dir, RObj_dir, "colorPaletteForConditions.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))


#plot the p-value histograms ----
ifelse(!dir.exists(file.path(results_dir)), 
       dir.create(file.path(results_dir)), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir))), 
       dir.create(file.path(paste0(results_dir, plotting_dir))), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "phenotypicCorrelations/"))), 
       dir.create(file.path(paste0(results_dir, plotting_dir, "phenotypicCorrelations/"))), FALSE)

pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations/", "pValueHistograms.pdf"))

lapply(colnames(traitCommonSegregants_std), FUN = function(condition){
  plot = gghistogram(phenotypicCorrelations[[condition]],
                     x = "p",
                     bins = 50,
                     title = condition,
                     xlab = "p") + geom_vline(xintercept = 0.05, color = "maroon") + theme_light() +
    theme_textProperties
})

dev.off()

#phenotypic r vs p across all conditions ----

pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations/", "phenoRVsP.pdf"))  
lapply(colnames(traitCommonSegregants_std), 
       FUN = function(condition_i){
         pheno_i = phenotypicCorrelations[[condition_i]]
         pheno_i$significance_FDR = ifelse(pheno_i$p < 0.05 & pheno_i$q < 0.05, "TRUE", "FALSE")
         pheno_i$logP = -log10(pheno_i$p)
         plot1 = ggscatter(pheno_i,
                           x = "r",
                           y = "logP",
                           shape = 21,
                           color = "significance_FDR",
                           alpha = 0.6,
                           title = condition_i) + 
           geom_hline(yintercept = -log10(0.05), color = "deep pink") +
           ylab("-log10(p)") + scale_color_manual(values = c("TRUE" = "dark blue", "FALSE" = "black")) + theme_light() + theme_textProperties + guides(color = guide_legend(title = "Significance @ 5% FDR"))
         
       })

dev.off()

#get genes at nominal bonferonni significance and 5% FDR ----

##bonferonni significance ---
phenotypicCorrelations_bonferonniSignificant = lapply(colnames(traitCommonSegregants_std),
                                                      FUN = function(condition){
                                                        df = phenotypicCorrelations[[condition]]
                                                        bonfFilter = which(df$p < 0.05/(5643 * 46)) #total number of tests across conditions is number of genes * number of conditions
                                                        df$condition = condition
                                                        return(df[bonfFilter,])
                                                      })

names(phenotypicCorrelations_bonferonniSignificant)= colnames(traitCommonSegregants_std)

df_bonf = do.call(rbind, phenotypicCorrelations_bonferonniSignificant)
df_bonf$significanceLevel = "bonferonni"

## 5% FDR ---

phenotypicCorrelations_significantAtFDR5 = lapply(colnames(traitCommonSegregants_std), FUN= function(condition){
  df = dplyr::filter(phenotypicCorrelations[[condition]], p < 0.05 & q < 0.05)
  df$condition = condition
  return(df)
})
names(phenotypicCorrelations_significantAtFDR5)= colnames(traitCommonSegregants_std)

df_FDR5 = do.call(rbind, phenotypicCorrelations_significantAtFDR5)
df_FDR5$significanceLevel = "5% FDR"


##nominal significance ---


df_significant = do.call(rbind, lapply(colnames(traitCommonSegregants_std), 
                                       FUN = function(condition){
                                         df = phenotypicCorrelations[[condition]]
                                         df$condition = condition
                                         return(dplyr::filter(df, p< 0.05))
                                       }))

df_significant$significanceLevel = "nominal"

##consolidate number of genes at different significance thresholds ---

nGenes_summary = rbind(df_significant, df_FDR5, df_bonf) %>% group_by(condition, significanceLevel) %>% summarise(nSignificantGenes = n(), nSignificantGenes_positive = sum(r >=0), nSignificantGenes_negative = sum(r < 0), proportionPositiveAssociations = ifelse(nSignificantGenes_positive> nSignificantGenes_negative, nSignificantGenes_positive/nSignificantGenes, -nSignificantGenes_positive/nSignificantGenes))

#plot ---

pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations/", "nGenes_phenotypicCorrelations_rectangularBarPlot.pdf" ),
    width = 15, height = 10)

nGenes_summary$significanceLevel = factor(nGenes_summary$significanceLevel, levels = c("nominal", "5% FDR", "bonferonni"))

ggplot() +
  geom_col(data = dplyr::filter(nGenes_summary, significanceLevel == "5% FDR"),
           aes(x = condition, y = nSignificantGenes_positive, fill = "blue",
               alpha = 0.3),                  
            width = 0.5) +
  geom_col(data = dplyr::filter(nGenes_summary, significanceLevel == "5% FDR"),
           aes(x = condition, y = -nSignificantGenes_negative, fill = "red",
               alpha = 0.3),
           width = 0.5) +
  geom_col(data = dplyr::filter(nGenes_summary, significanceLevel == "bonferonni"),
           aes(x = condition, y = nSignificantGenes_positive, fill = "blue", alpha = 1),
           width = 0.5) +
  geom_col(data = dplyr::filter(nGenes_summary, significanceLevel == "bonferonni"),
           aes(x = condition, y = -nSignificantGenes_negative, fill = "red", alpha = 1),
           width = 0.5) +
  theme_classic() +
  ylab("Number of genes") + 
  theme_textProperties +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust= 1, vjust = 0.5)) +
  scale_alpha_identity(name =  "Significance Level",
                      breaks = c(0.3, 1),
                      labels = c("5% FDR", "Bonferonni"),
                      guide = "legend") +
  scale_fill_identity(name =  "Direction of correlation:",
                         breaks = c("blue", "red"),
                         labels = c("r > 0", "r < 0"),
                         guide = "legend") +
  guides(alpha = guide_legend(title = "Significance Level:",
                             title.position = "left",
                             title.vjust = 1),
         fill = guide_legend(title = "Direction of correlation:",
                              title.position = "left",
                              title.vjust = 1)) +
  scale_y_continuous(limits = c(-3000, 4000), 
                     breaks = c(-3000, -1500, 0, 1500, 3000),
                     labels = c(3000, 1500, 0, 1500, 3000))



