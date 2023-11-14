#Make GO Term panel

library(ComplexHeatmap)
library(ggpubr)
library(RColorBrewer)
library(circlize)
library(dplyr)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

# load data
load(paste0(results_dir, RObj_dir, 
            "GO_geneticCorrelations_QTLEffectsCorrelation_fdr5_HotspoteffectCorrelation_plusMinus_allP_myOwnFishersTest.rda"))
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#
col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
col2 = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(256)

GoSlimTerms = dplyr::filter(GoSlimTerms, OC == "BP")

#make a dendrogram 
render_col_dendrogram = function(GOEnrichments_plus, GOEnrichments_minus, km){
  df = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
  colnames(df) = colnames(traitCommonSegregants_std)
  GoSlimTerms = cbind(GoSlimTerms, df)
  
  for(condition_i in colnames(traitCommonSegregants_std)){
    df_plus = GOEnrichments_plus[[condition_i]]
    df_plus = dplyr::filter(df_plus, p < 0.001)
    #df_plus$logFold = log2(df_plus$Count/df_plus$ExpCount)
    
    df_minus = GOEnrichments_minus[[condition_i]]
    df_minus = dplyr::filter(df_minus, p < 0.001)
    #df_minus$logFold = log2(df_minus$Count/df_minus$ExpCount)

    xPlus = df_plus$log2FC[match(row.names(GoSlimTerms), df_plus$GOID)] * 1
    names(xPlus) = GoSlimTerms$Term
    xMinus = df_minus$log2FC[match(row.names(GoSlimTerms), df_minus$GOID)] * -1
    names(xMinus) = GoSlimTerms$Term
    
    x = coalesce(xPlus, xMinus)
    
    #for GO-BP terms with nominally significant enrichment in both plus and minus, we retain the one where the fold enrichment is higher.
    rowsWithNonNA = which(!is.na(xPlus) & !is.na(xMinus))
    for(i in rowsWithNonNA){
      x[i] = ifelse(x[i] == max(abs(xPlus[i]), abs(xMinus[i])), xPlus[i], xMinus[i])
    }
    
    GoSlimTerms[condition_i] = x
  }
  
  row.names(GoSlimTerms) = GoSlimTerms$Term
  GoSlimTerms = dplyr::select(GoSlimTerms,-GOID, -Term, -OC)
  GoSlimTerms = as.matrix(GoSlimTerms)
  GoSlimTerms[is.na(GoSlimTerms) | is.infinite(GoSlimTerms) | is.nan(GoSlimTerms)] = 0 #replace NA terms with 0
  #GoSlimTerms = GoSlimTerms[which(rowSums(GoSlimTerms) != 0),]
  GoSlimTerms = t(GoSlimTerms)
  
  library(dendextend)
  col_dend = as.dendrogram(hclust(dist(as.matrix(GoSlimTerms))))
  color_colBranches = color_branches(col_dend, k = km) #3 clusters
  colAnno = data.frame(condition = labels(color_colBranches),
                       grp_color = get_leaves_branches_col(color_colBranches))
  colAnno = colAnno[order(match(colAnno$condition, row.names(GoSlimTerms))),] #sort the conditions in the same order as in allPearsonRs matrix
  row.names(colAnno) = row.names(GoSlimTerms)
  
  return(list(dendrogram = col_dend,
              annot = colAnno))
}

#####

##plot the logFC for GO-BP slim terms with p < 0.001 on heatmap
fn_plotGO = function(GOEnrichments_plus, GOEnrichments_minus, km,
                     nGenes_vector, title){
  df = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
  colnames(df) = colnames(traitCommonSegregants_std)
  GoSlimTerms = cbind(GoSlimTerms, df)
  
  for(condition_i in colnames(traitCommonSegregants_std)){
    print(condition_i)
    df_plus = GOEnrichments_plus[[condition_i]]
    df_plus = dplyr::filter(df_plus, p < 0.001)
    #df_plus$logFold = log2(df_plus$Count/df_plus$ExpCount)
    
    df_minus = GOEnrichments_minus[[condition_i]]
    df_minus = dplyr::filter(df_minus, p < 0.001)
    #df_minus$logFold = log2(df_minus$Count/df_minus$ExpCount)
    
    xPlus = df_plus$log2FC[match(row.names(GoSlimTerms), df_plus$GOID)] * 1
    names(xPlus) = GoSlimTerms$Term
    xMinus = df_minus$log2FC[match(row.names(GoSlimTerms), df_minus$GOID)] * -1
    names(xMinus) = GoSlimTerms$Term
    
    x = coalesce(xPlus, xMinus)
    
    #for GO-BP terms with nominally significant enrichment in both plus and minus, we retain the one where the fold enrichment is higher.
    rowsWithNonNA = which(!is.na(xPlus) & !is.na(xMinus))
    for(i in rowsWithNonNA){
      x[i] = ifelse(x[i] == max(abs(xPlus[i]), abs(xMinus[i])), xPlus[i], xMinus[i])
    }
    
    GoSlimTerms[condition_i] = x
  }
  
  row.names(GoSlimTerms) = GoSlimTerms$Term
  GoSlimTerms = dplyr::select(GoSlimTerms, -GOID, -Term, -OC)
  GoSlimTerms = as.matrix(GoSlimTerms)
  GoSlimTerms[is.na(GoSlimTerms) | is.infinite(GoSlimTerms) | is.nan(GoSlimTerms)] = 0 #replace NA terms with 0
  #GoSlimTerms = GoSlimTerms[which(rowSums(GoSlimTerms) != 0),]
  
  # #annotate each condition-column with number of significant genes considered for GO enrichments
   ha2 = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:km+1)),
                           nGenes = anno_barplot(nGenes_vector, gp = gpar(fill = c("red3", "blue3"))),
                           show_annotation_name = FALSE, height = unit(2, "cm"))
   
  
  #dendrogram 
  col_dend = render_col_dendrogram(GOEnrichments_plus = GO_GC_plus, 
                                   GOEnrichments_minus = GO_GC_minus,
                                   km = km)
  
  ht = Heatmap(as.matrix(GoSlimTerms), col = col,
               use_raster = FALSE, show_row_names = TRUE,
               width = unit(40, "cm"),
               height = unit(50, "cm"),
               show_row_dend = FALSE,
               show_column_dend = TRUE,
               #column_dend_gp = gpar(col = col_dend$annot$grp_color),
               #column_split = 2,
               column_km = km,
               #cluster_columns = col_dend$dendrogram,
               row_names_gp = gpar(fontsize = 20),
               column_names_gp = gpar(fontsize = 22),
               column_title = title,
               top_annotation = ha2,
               #bottom_annotation = ha3,
               #cluster_rows = FALSE,
               #cluster_columns = FALSE,
               row_title_rot = 0, 
               name = "GO fold enrichment (count/expected count):", 
               row_names_side = "left",
               heatmap_legend_param = list(direction = "horizontal", side = "bottom",
                                           title = "GO term fold enrichment:"),
               row_gap = unit(1, "mm"))
  
  #draw(ht, heatmap_legend_side = "top")
  return(ht)
}

nGenes_GC = lapply(colnames(traitCommonSegregants_std),
                   FUN = function(condition_i){
                     df = dplyr::filter(phenotypicCorrelations[[condition_i]], p < 0.05 &q < 0.05)
                     summary_df = summarise(df, plus = sum(r > 0), minus = sum(r < 0))
                     return(summary_df)
                   })
nGenes_GC = do.call(rbind, nGenes_GC)
ht_GC = fn_plotGO(GO_GC_plus, GO_GC_minus,
                  km = 3, nGenes_GC, title = "Genetic correlation")

## get the brauer and ESR heat maps now ------
load(file = paste0(results_dir, RObj_dir, "aveGeneticR.rda"))

#function to plot the heatmap
plotHt = function(df, heatmapLegend, colorPal, name){
  #colnames(df) = colnames_withoutPrefix
  ht1_genetic = Heatmap(as.matrix(df), 
                        col =colorPal, use_raster = FALSE,
                        row_order = factor(row.names(df), levels = unique(row.names(df))),
                        cluster_columns = FALSE,
                        show_heatmap_legend = heatmapLegend, 
                        height = unit(2, "cm"), 
                        heatmap_legend_param = list(direction = "horizontal", side = "right"),
                        show_row_dend = FALSE,
                        name = name,
                        row_gap = unit(3, "mm"),
                        column_gap = unit(5, "mm"),
                        #border = TRUE,
                        row_title = "Gene",
                        column_title = "Growth condition"
  )
  return(ht1_genetic)
}

#brauer
ht1_brauer = plotHt(aveGeneticR_brauerGroups,
                    TRUE, col2, "av. r")
#ESR
ht1_ESR = plotHt(aveGeneticR_ESRGroups,
                 TRUE, col2, "av. r")

##combine GO-heatmap, brauer and ESR
ht =ht_GC %v%
  ht1_brauer %v% ht1_ESR

pdf(paste0(results_dir, plotting_dir, "GOSlims_GC_p05_3cuts.pdf"), width = 25, height = 25)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right",
     padding = unit(1, "cm"))
dev.off()

#### get average brauer and ESR for each of the clusters ----
df = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
colnames(df) = colnames(traitCommonSegregants_std)
GoSlimTerms = cbind(GoSlimTerms, df)

for(condition_i in colnames(traitCommonSegregants_std)){
  print(condition_i)
  df_plus = GO_GC_plus[[condition_i]]
  df_plus = dplyr::filter(df_plus, p < 0.001)
  #df_plus$logFold = log2(df_plus$Count/df_plus$ExpCount)
  
  df_minus = GO_GC_minus[[condition_i]]
  df_minus = dplyr::filter(df_minus, p < 0.001)
  #df_minus$logFold = log2(df_minus$Count/df_minus$ExpCount)
  
  xPlus = df_plus$log2FC[match(row.names(GoSlimTerms), df_plus$GOID)] * 1
  names(xPlus) = GoSlimTerms$Term
  xMinus = df_minus$log2FC[match(row.names(GoSlimTerms), df_minus$GOID)] * -1
  names(xMinus) = GoSlimTerms$Term
  
  x = coalesce(xPlus, xMinus)
  
  #for GO-BP terms with nominally significant enrichment in both plus and minus, we retain the one where the fold enrichment is higher.
  rowsWithNonNA = which(!is.na(xPlus) & !is.na(xMinus))
  for(i in rowsWithNonNA){
    x[i] = ifelse(x[i] == max(abs(xPlus[i]), abs(xMinus[i])), xPlus[i], xMinus[i])
  }
  
  GoSlimTerms[condition_i] = x
}

row.names(GoSlimTerms) = GoSlimTerms$Term
GoSlimTerms = dplyr::select(GoSlimTerms, -GOID, -Term, -OC)
GoSlimTerms = as.matrix(GoSlimTerms)
GoSlimTerms[is.na(GoSlimTerms) | is.infinite(GoSlimTerms) | is.nan(GoSlimTerms)] = 0 #replace NA terms with 0

ht = Heatmap(as.matrix(GoSlimTerms), column_km = 3, cluster_columns = TRUE)

trait_clusters = column_order(draw(ht))
names(trait_clusters) = c("green", "blue", "yellow")

GC_avg = lapply(trait_clusters, 
                    function(cluster_i){
                      df_all = rbind(aveGeneticR_brauerGroups, aveGeneticR_ESRGroups)
                      df_all_i = df_all[,cluster_i]
                      avg_allGrps = apply(df_all_i, 1, FUN = median, na.rm = TRUE)
                      return(avg_allGrps)
                    })

