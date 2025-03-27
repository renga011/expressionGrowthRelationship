#Reviewer 1 comment 3: Causal hotspot locus

#load libraries
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

## load files -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "propVarExplained_topNHotspots.rda")) # from the R script h2ExplainedByDifferentHotspots.r
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

##

hotspots = hotspotData_Albert2018[order(hotspotData_Albert2018$numberNonzeroEffects,
                                        decreasing = TRUE),]
hotspots = select(hotspots, hotspotMarker, numberNonzeroEffects)

h2_diffBetweenHotspots = lapply(colnames(traitCommonSegregants_std),
                                function(condition_i){
                                  df = h2_topNHotspots[[condition_i]]
                                  diff_hotspots = diff(df$h2) #difference between consecutive rows
                                  diff_hotspots = c(NA, NA, diff_hotspots) #adding 2 NAs for 1st & 2nd hotspot (descending order)
                                  return(diff_hotspots)
                                })

names(h2_diffBetweenHotspots) = colnames(traitCommonSegregants_std)
h2_diffBetweenHotspots = as.data.frame(h2_diffBetweenHotspots)
row.names(h2_diffBetweenHotspots) = stringr::str_split(hotspotData_Albert2018$hotspotMarker, "_", simplify = TRUE)[,1]

# plot ----
#prep data to remove NAs
h2_diffBetweenHotspots = h2_diffBetweenHotspots[complete.cases(h2_diffBetweenHotspots),]

# Define the custom color palette
# Map -0.3 to blue, 0 to white, and 0.3 to red
my_colors = colorRampPalette(c("blue", "white", "red"))(100)

# Define the breaks for the heatmap
# Ensure that -0.3, 0, and 0.3 are nicely aligned with the colors
breaks = seq(-0.2, 0.2, length.out = 101)

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "propOfh2ByEachHotspotX46Conditions.pdf"), width = 15, height = 15)
pheatmap(as.matrix(t(h2_diffBetweenHotspots)), color = my_colors,
         breaks = breaks, cluster_cols = FALSE)
dev.off()

# make corresponding supplementary table -----
write.csv(h2_diffBetweenHotspots, file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", "propOfGrowthVarianceByEachHotspotX46Conditions.csv"), quote = FALSE)

