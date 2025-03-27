library(ggpubr)
library(dendextend)
library(dplyr)
library(Hmisc)
library(ComplexHeatmap)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)

#
traitCorrelationMatrix = rcorr(as.matrix(traitCommonSegregants_std),type = "spearman")

#square heatmap
my_colors = circlize::colorRamp2(c(-1, 0, 1), c("navyblue", "white", "maroon"))

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "squareHeatmapOfTraitCorrelations.pdf"), width = 10, height = 10)
ht_traitCorrelation = Heatmap(traitCorrelationMatrix$r, col = my_colors,
                use_raster = FALSE, show_row_names = TRUE,
                #width = unit(40, "cm"),
                #height = unit(50, "cm"),
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                #column_dend_gp = gpar(col = col_dend$annot$grp_color),
                #column_split = 2,
                #column_km = 3,
                #cluster_columns = col_dend$dendrogram,
                #row_names_gp = gpar(fontsize = 20),
                #column_names_gp = gpar(fontsize = 22),
                #top_annotation = ha2,
                #bottom_annotation = ha3,
                #cluster_rows = FALSE,
                #cluster_columns = FALSE,
                row_title_rot = 0, 
                name = "GO fold enrichment (count/expected count):", 
                row_names_side = "left",
                heatmap_legend_param = list(direction = "horizontal", side = "bottom",
                                            title = "rho:"),
                row_gap = unit(1, "mm"))
draw(ht_traitCorrelation, heatmap_legend_side = "top",
     padding = unit(1, "cm"))

dev.off()

# GO matrix in Figure 6

load(paste0(results_dir, RObj_dir, "GO_geneticCorrelations_QTLEffectsCorrelation_fdr5_HotspoteffectCorrelation_plusMinus_allP_myOwnFishersTest.rda"))

col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

df = matrix(0, nrow= nrow(GoSlimTerms), ncol = ncol(traitCommonSegregants_std))
colnames(df) = names(GO_GC)
GoSlimTerms_mat_GC = cbind(GoSlimTerms, df)

for(condition in colnames(traitCommonSegregants_std)){
  GO_GC_i = GO_GC[[condition]]
  GO_GC_plus_i = GO_GC_plus[[condition]]
  GO_GC_minus_i = GO_GC_minus[[condition]]
  
  #filter for significant terms
  go_plus = dplyr::filter(GO_GC_plus_i, p < 0.001 & log2FC > 0)
  go_minus = dplyr::filter(GO_GC_minus_i, p < 0.001 & log2FC > 0)
  
  #make log2FC + for go_plus and - for go_minus
  xPlus = go_plus$log2FC[match(row.names(GoSlimTerms), go_plus$GOID)] * 1
  names(xPlus) = GoSlimTerms$Term
  xMinus = go_minus$log2FC[match(row.names(GoSlimTerms), go_minus$GOID)] * -1
  names(xMinus) = GoSlimTerms$Term
  
  x = coalesce(xPlus, xMinus)
  
  #for GO-BP terms with nominally significant enrichment in both plus and minus, we retain the one where the fold enrichment is higher.
  rowsWithNonNA = which(!is.na(xPlus) & !is.na(xMinus))
  for(i in rowsWithNonNA){
    x[i] = ifelse(x[i] == max(abs(xPlus[i]), abs(xMinus[i])), xPlus[i], xMinus[i])
  }
  
  GoSlimTerms_mat_GC[condition] = x
  
}

row.names(GoSlimTerms_mat_GC) = GoSlimTerms_mat_GC$Term
GoSlimTerms_mat_GC = dplyr::select(GoSlimTerms_mat_GC, -GOID, -Term, -OC)
GoSlimTerms_mat_GC = as.matrix(GoSlimTerms_mat_GC)
GoSlimTerms_mat_GC[is.na(GoSlimTerms_mat_GC) | is.infinite(GoSlimTerms_mat_GC) | is.nan(GoSlimTerms_mat_GC)] = 0 #replace NA terms with 0

ht_GC = Heatmap(GoSlimTerms_mat_GC, col = col,
                use_raster = FALSE, show_row_names = TRUE,
                #width = unit(40, "cm"),
                #height = unit(50, "cm"),
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                #column_dend_gp = gpar(col = col_dend$annot$grp_color),
                #column_split = 2,
                column_km = 3,
                #cluster_columns = col_dend$dendrogram,
                #row_names_gp = gpar(fontsize = 20),
                #column_names_gp = gpar(fontsize = 22),
                #top_annotation = ha2,
                #bottom_annotation = ha3,
                #cluster_rows = FALSE,
                #cluster_columns = FALSE,
                row_title_rot = 0, 
                name = "GO fold enrichment (count/expected count):", 
                row_names_side = "left",
                heatmap_legend_param = list(direction = "horizontal", side = "bottom",
                                            title = "GO term fold enrichment:"),
                row_gap = unit(1, "mm"))



# compare dendrograms by comparing their distance matrices (using correlation based distance)
  # Step 1: Extract the pairwise distance matrices (dissimilarity matrices)
  getCorrDistance = function(matrix1){
    zero_variance_cols <- which(apply(matrix1, 2, var, na.rm = TRUE) == 0)
    matrix_clean <- matrix1[, -zero_variance_cols]
    dist1 = as.dist(1-cor(matrix_clean, use = "pairwise.complete.obs"))  # Correlation distance matrix for the first set of data
    return(dist1)
  }
  
  dist_GC = getCorrDistance(GoSlimTerms_mat_GC)
  dist_traitCorrelations = as.dist(1-cor(traitCorrelationMatrix_1, use = "pairwise.complete.obs"))

    # Correlation distance matrix for the first set of data
  
  # Step 2: Flatten the distance matrices into vectors (since they are symmetric)
  dist1_vector <- as.vector(dist_GC)
  dist2_vector <- as.vector(dist_traitCorrelations)
  
  # correlation test
  cor_test <- cor.test(dist1_vector, dist2_vector, method = "spearman")

  pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "comparisonOfDendrogramsOfTraitClusters_GCVsTraitCorrelation.pdf"))
  ggscatter(data.frame("GC_distances" = dist1_vector,
                       "traitCorrelation_distances" = dist2_vector),
            y = "GC_distances",
            x = "traitCorrelation_distances",
            shape = 21) + stat_cor(method = "spearman")
  dev.off()
  