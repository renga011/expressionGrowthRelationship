#plot the genetic correlation of HSP12-growth in copper
#plot the eQTL effects on expression and growth on a chromosome plot
#plot all GC vs QEC
#plot all hotspot effects on expression and growth on a chromosome plot
#plot all HEC vs GC

library(dplyr)
library(ggplot2)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

# load data
load(paste0(results_dir, "chromosomeSizes.rda"))
    chromosomeSizes$plotting_locus= 0.5* (chromosomeSizes$chromosomeStart_plotting + chromosomeSizes$chromosomeEnd_plotting)
    chromosomeSizes$label= as.character(as.roman(chromosomeSizes$chromosome))
    chromosomeSizes = dplyr::filter(chromosomeSizes, chromosome != 0)

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(file = paste0(results_dir, RObj_dir, "QTLEffectsCorrelation_atEQTLs.rda"))
load(file = paste0(results_dir, RObj_dir, "expressionAndGrowthEffectsAtEQTLs.rda")) #r for expression and growth at different eQTLs
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations_recomputed_onlyMyScanningMethod.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations.rda"))
                                                                                                #functions ---
addPlottingCoords=  function(QTLs){
  for(i in 1:nrow(QTLs)){
    chr_i = QTLs$chromosome[i]
    additiveFactor = chromosomeSizes$chromosomeStart_plotting[which(chromosomeSizes$chromosome == chr_i)]
    CILeft_i = QTLs$CILeft[i]
    CIRight_i = QTLs$CIRight[i]
    QTLs$plotting_start[i] = CILeft_i + additiveFactor
    QTLs$plotting_end[i] = CIRight_i + additiveFactor
  }
  return(QTLs)
}

#genetic correlation
condition_i = "Copper"
df_GC = data.frame(expression = unlist(expressionCommonSegregants_batchODCorrected_std["YFL014W"]),
                growth = unlist(traitCommonSegregants_std[condition_i]))

pdf(paste0(results_dir, plotting_dir, "HSP12_GC.pdf"))

plot1 = ggplot(df_GC, aes(x = expression, y = growth)) + 
  #geom_point(alpha = 0.5, shape = 21) +
  #geom_smooth(method = "lm") +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c() + 
  theme_bw() +
  xlab("Expression") + ylab ("Growth") + stat_cor(method= "spearman")

print(plot1)

dev.off()

  #chromosome plot ---

#QEC- expression vs growth for HSP12, 
condition_i = "Copper"
df = dplyr::filter(computeRGrowth[[condition_i]], gene == "YFL014W")
  for(i in 1:nrow(df)){
    expression_i = unlist(expressionCommonSegregants_batchODCorrected_std["YFL014W"])
    pmarker = df$pmarker[i]
    genotype_i = unlist(genotypesCommonSegregants[pmarker])
    
    growth_i = unlist(traitCommonSegregants_std[condition_i])
    
    ct_expression = cor.test(x = genotype_i, y = expression_i)
    df$r_expression[i] = ct_expression$estimate
    df$r_expression_CILeft[i] = ct_expression$conf.int[1]
    df$r_expression_CIRight[i] = ct_expression$conf.int[2]
    
    ct_growth = cor.test(x = genotype_i, y = growth_i)
    df$r_growth_CILeft[i] = ct_growth$conf.int[1]
    df$r_growth_CIRight[i] = ct_growth$conf.int[2]
  }

pdf(paste0(results_dir, plotting_dir, "HSP12_QEC_fullGenomePlot.pdf"), width = 10, height = 8)


  plot1 = ggplot() +
    geom_point(df, mapping = aes(x = position_plotting, y = r_expression), color = "maroon",
               shape = 21, size = 7) +
    geom_errorbar(df, mapping = aes(x = position_plotting, ymin = r_expression_CILeft,
                                    ymax = r_expression_CIRight), color = "maroon",alpha = 0.6) +
    geom_point(df, mapping = aes(x = position_plotting,  y = r_growth), color = "blue4",
               shape = 21, size = 7) +
    geom_errorbar(df, mapping = aes(x = position_plotting, ymin = r_growth_CILeft,
                                    ymax = r_growth_CIRight), color = "blue4", alpha = 0.6) +
    theme_light() +
    geom_vline(data = chromosomeSizes,
               aes(xintercept = chromosomeEnd_plotting),
               color = "gray", linetype = "dotted") +
    geom_vline(xintercept = 0,
               color = "gray", linetype = "dotted") +
    scale_x_continuous(breaks = chromosomeSizes$plotting_locus, 
                       labels = chromosomeSizes$label,
                       name = "chromosome") + 
    ylab("r") +
    scale_color_identity(name =  "Trait:",
                         breaks = c("maroon", "blue4"),
                         labels = c("Expression", "Growth"),
                         guide = "legend") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  print(plot1)

dev.off()

pdf(paste0(results_dir, plotting_dir, "HSP12_QEC_ExpressionVsGrowth.pdf"))
  
  plot2 = ggplot(df) +
    geom_point(mapping = aes(x = r_expression, y = r_growth, 
                             size = (ci_expression*ci_growth)^-1)) +
    geom_errorbar(mapping = aes(x = r_expression, y = r_growth,
                                xmin = r_expression_CILeft, 
                                xmax = r_expression_CIRight),
                  alpha = 0.6) +
    geom_errorbar(mapping = aes(x = r_expression, y = r_growth,
                                ymin = r_growth_CILeft, 
                                ymax = r_growth_CIRight),
                  alpha = 0.6) +
    theme_bw() +
    theme(legend.position = "top") +
    geom_smooth(method = "lm", 
                mapping = aes(x = r_expression, y = r_growth,
                              weight = (ci_expression*ci_growth)^-1))
  
  print(plot2)
dev.off()

#correlation between expression and growth HSP12 - put as label on graph in illustrator
  ct_GCVsQEC_HSP12 = cor.test(x = df$r_expression, y = df$r_growth,
                            method = "spearman")

##all GC vs QEC
GC_all = do.call(rbind, phenotypicCorrelations)
QEC_all = do.call(rbind, QEC)
df_plot = data.frame(GC = GC_all$r,
                    GC_p = GC_all$p,
                    GC_q = GC_all$q,
                    QEC= QEC_all$r_QEC,
                    QEC_p = QEC_all$p_QEC,
                    QEC_q = QEC_all$q_QEC,
                    gene = GC_all$gene)
plot3 = ggplot(df_plot, aes(x = GC, y = QEC)) + 
  #geom_point(alpha = 0.1) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c() + 
    theme_bw() +
  xlab("Genetic r") + ylab ("QTL-effects r")

pdf(paste0(results_dir, plotting_dir, "GCVsQEC.pdf"), width = 15, height = 10)
print(plot3)
dev.off()

#correlation between expression and growth- all genes - GC vs QEC - put as label on graph in illustrator
  ct_GCVsQEC_all = cor.test(x = df_plot$GC, y = df_plot$QEC,
                            method = "spearman")
# GC vs HEC
#make hotspot table
df_i = data.frame(hotspot_marker = row.names(expressionEffects_i),
                  expression = unlist(expressionEffects_i["YFL014W"]),
                  growth = unlist(growthEffects_i["Copper"]))
df_i$peakMarkerPosition = tidyr::extract_numeric(df_i$hotspot_marker)
df_i$chr = stringr::str_remove(hotspotData_Albert2018$chromosome, pattern = "chr")
df_i$chr = as.numeric(as.roman(df_i$chr))
  #add plotting position
for(i in 1:nrow(df_i)){
  chr_i = df_i$chr[i]
  additiveFactor = chromosomeSizes$chromosomeStart_plotting[which(chromosomeSizes$chromosome == chr_i)]
  df_i$markerPosition_plotting[i] = df_i$peakMarkerPosition[i] + additiveFactor
}

  #chromosome plot
pdf(paste0(results_dir, plotting_dir, "HSP12_HEC_fullGenomePlot.pdf"), width = 10, height = 8)

plot3 = ggplot() +
  geom_point(df_i, 
             mapping = aes(x = markerPosition_plotting, y = expression,
                           alpha = ifelse(expression == 0, 0.2, 1)), 
             shape = 21, size = 6, color= "maroon") +
  geom_point(df_i, mapping = aes(x = markerPosition_plotting,  y = growth,
                                 alpha = ifelse(growth == 0, 0.2, 1)),
             shape = 21, size = 6, color= "blue4") +
  theme_light() +
  geom_vline(data = chromosomeSizes,
             aes(xintercept = chromosomeEnd_plotting),
             color = "gray", linetype = "dotted") +
  geom_vline(xintercept = 0,
             color = "gray", linetype = "dotted") +
  scale_x_continuous(breaks = chromosomeSizes$plotting_locus, 
                     labels = chromosomeSizes$label,
                     name = "chromosome") + 
  ylab("r") +
  scale_color_identity(name =  "Trait:",
                       breaks = c("maroon", "blue4"),
                       labels = c("Expression", "Growth"),
                       guide = "legend") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

print(plot3)

dev.off()

pdf(paste0(results_dir, plotting_dir, "HSP12_HEC_ExpressionVsGrowth.pdf"))
plot4 = ggplot(df_i) +
  geom_point(mapping = aes(x = expression, y = growth), size = 6, alpha = 1, shape = 21) +
  theme_bw() +
  geom_smooth(method = "lm", 
              mapping = aes(x = expression, y = growth))

print(plot4)
dev.off()

#correlation between expression and growth HSP12 - put as label on graph in illustrator
ct_GCVsHEC_HSP12 = cor.test(x = df_i$expression, y = df_i$growth,
                          method = "spearman")
##all GC vs HEC
phenotypicCorrelations = lapply(colnames(traitCommonSegregants_std),
                                function(condition_i){
                                  x = phenotypicCorrelations[[condition_i]]
                                  x$condition = condition_i
                                  return(x)
                                })
names(phenotypicCorrelations) = colnames(traitCommonSegregants_std)

hotspotEffectCorrelations = lapply(colnames(traitCommonSegregants_std),
                                function(condition_i){
                                  x = hotspotEffectCorrelations[[condition_i]]
                                  x$condition = condition_i
                                  return(x)
                                })
names(hotspotEffectCorrelations) = colnames(traitCommonSegregants_std)

GEC_all = do.call(rbind, phenotypicCorrelations)
HEC_all = do.call(rbind, hotspotEffectCorrelations)

for(i in 1:nrow(GEC_all)){
  print(i)
  gene_i = GEC_all$gene[i]
  condition_i = GEC_all$condition[i]
  row_i = which(HEC_all$gene==gene_i & HEC_all$condition==condition_i)
  if(!(rlang::is_empty(row_i))){
  GEC_all$r_HEC[i] = HEC_all$r.cor[row_i]
  } else{
    GEC_all$r_HEC[i] = NA
  }
}

plot6 = ggplot(GEC_all, aes(x = r, y = r_HEC)) + 
  #geom_point(alpha = 0.1) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c() + 
  theme_bw() +
  xlab("Genetic r") + ylab ("Hotspot effects r")

pdf(paste0(results_dir, plotting_dir, "GCVsHEC.pdf"), width = 15, height = 10)
print(plot6)
dev.off()

  #correlation between expression and growth - all genes, GC vs HEC - put as label on graph in illustrator
ct_GCVsHEC_all = cor.test(x = GEC_all$r, y = GEC_all$r_HEC,
                          method = "spearman")



###
df_i = dplyr::filter(df_plot, QEC_p < 0.05 & QEC_q < 0.05)
cor.test(df_i$GC, df_i$QEC)
