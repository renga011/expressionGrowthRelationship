## Reviewer 1, Comment 6: GO Enrichment without the plus minus division -----

library(org.Sc.sgd.db)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

# load files -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))

  # get genes belonging to each GO term -----

gn2go <- mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db), "GOALL", "ORF", multiVals = "list")
#'select()' returned 1:many mapping between keys and columns
## filter out NA mappings
gn2go <- gn2go[!sapply(gn2go, function(x) all(is.na(x)))]

#shortlist only the GO-Slim terms
gn2GO_slim = lapply(gn2go, function(gene_i){
  slim_index = which(gene_i %in% row.names(GoSlimTerms))
  return(gene_i[slim_index])
})

GO2gn = lapply(row.names(GoSlimTerms),
               FUN = function(slimTerm){
                 genes_i = names(gn2GO_slim)[sapply(1:length(gn2GO_slim),
                                                    function(x){slimTerm %in% gn2GO_slim[[x]]})]
                 return(genes_i)
               })
names(GO2gn) = row.names(GoSlimTerms)

##convert common names to systematic names
GO2gn = lapply(GO2gn, function(genes_GOTerm){
  for(i in 1: length(genes_GOTerm)){
    gene_i = genes_GOTerm[i]
    a = which(geneList_SGD$stdName == gene_i)
    if(!(rlang::is_empty(a))){
      genes_GOTerm[i] = geneList_SGD$sysName[a]
    }
    return(genes_GOTerm)
  }})

# function to do our own fisher's gene enrichment test ------
GO_FT = function(GOterm, geneList, 
                 universe = colnames(expressionCommonSegregants_batchODCorrected_std)){
  genes_GOTerm = GO2gn[[GOterm]]
  t_t = length(intersect(genes_GOTerm, geneList))
  t_f = length(setdiff(genes_GOTerm, geneList))
  f_t = length(setdiff(geneList, genes_GOTerm))
  f_f = length(setdiff(universe, union(geneList, genes_GOTerm)))
  
  twoByTwo = matrix(c(t_t, t_f, f_t, f_f), nrow =2, ncol = 2)
  
  expectedCount = (sum(twoByTwo[1,])* sum(twoByTwo[,1]))/sum(unlist(twoByTwo))
  
  ft = fisher.test(twoByTwo)
  return(c("size" = length(genes_GOTerm),
           "p" = ft$p.value, "odds_ratio" = ft$estimate, 
           "count" = t_t, "expected_count" = expectedCount, "log2FC" = log2(t_t/expectedCount)))
}

# do the GO term analysis (doing only for genetic correlation list) ----
GO_GC = lapply(colnames(traitCommonSegregants_std),
               FUN = function(condition_i){
                 print(condition_i)
                 df = phenotypicCorrelations[[condition_i]]
                 sigGenes_i = df$gene[which(df$p < 0.05 & df$q < 0.05)]
                 
                 GO_i = lapply(row.names(GoSlimTerms),
                               FUN = function(GOTerm_i){
                                 GO_FT_i = GO_FT(GOTerm_i,
                                                 sigGenes_i)
                                 return(GO_FT_i)
                               })
                 GO_i = as.data.frame(do.call(rbind, GO_i))
                 GO_i$q = p.adjust(GO_i$p, method = "fdr")
                 table = as.data.frame(cbind(GoSlimTerms, GO_i))
                 table$condition = condition_i
                 
                 return(table)
               })

names(GO_GC) = colnames(traitCommonSegregants_std)

save(GO_GC,
     file = paste0(results_dir, RObj_dir, "GO_geneticCorrelations_QTLEffectsCorrelation_fdr5_HotspoteffectCorrelation_noDirections_allP_myOwnFishersTest.rda"))


# plot the GO term analysis -----

# color palettes and GO slim terms
col = colorRamp2(c(0, 4), c("white", "darkred"))
col2 = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(256)

GoSlimTerms = dplyr::filter(GoSlimTerms, OC == "BP")

#make a dendrogram 
render_col_dendrogram = function(GOEnrichments, km){
  df = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
  colnames(df) = colnames(traitCommonSegregants_std)
  GoSlimTerms = cbind(GoSlimTerms, df)
  
  for(condition_i in colnames(traitCommonSegregants_std)){
    df = GOEnrichments[[condition_i]]
    df = dplyr::filter(df, p < 0.001 & log2FC > 0)
    
    x = df$log2FC[match(row.names(GoSlimTerms), df$GOID)] * 1
    names(x) = GoSlimTerms$Term
    
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

##plot the logFC for GO-BP slim terms with p < 0.001 on heatmap
fn_plotGO = function(GOEnrichments, km,
                     nGenes_vector, title){
  df = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
  colnames(df) = colnames(traitCommonSegregants_std)
  GoSlimTerms = cbind(GoSlimTerms, df)
  
  for(condition_i in colnames(traitCommonSegregants_std)){
    print(condition_i)
    df = GOEnrichments[[condition_i]]
    df = dplyr::filter(df, p < 0.001 & log2FC > 0)
    
    x = df$log2FC[match(row.names(GoSlimTerms), df$GOID)] * 1
    names(x) = GoSlimTerms$Term
    
    GoSlimTerms[condition_i] = x
  }
  
  row.names(GoSlimTerms) = GoSlimTerms$Term
  GoSlimTerms = dplyr::select(GoSlimTerms, -GOID, -Term, -OC)
  GoSlimTerms = as.matrix(GoSlimTerms)
  GoSlimTerms[is.na(GoSlimTerms) | is.infinite(GoSlimTerms) | is.nan(GoSlimTerms)] = 0 #replace NA terms with 0
  #GoSlimTerms = GoSlimTerms[which(rowSums(GoSlimTerms) != 0),]
  
  # #annotate each condition-column with number of significant genes considered for GO enrichments
  ha2 = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:km+1)),
                          nGenes = anno_barplot(nGenes_vector, gp = gpar(fill = c("darkgray"))),
                          show_annotation_name = FALSE, height = unit(2, "cm"))
  
  
  #dendrogram 
  col_dend = render_col_dendrogram(GOEnrichments = GO_GC,
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
                     summary_df = summarise(df, nGenes = n())
                     return(summary_df)
                   })

nGenes_GC = do.call(rbind, nGenes_GC)
ht_GC = fn_plotGO(GO_GC,
                  km = 3, nGenes_GC, title = "Genetic correlation")

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/GOSlims_GC_noDirection_p05_3cuts.pdf"), width = 25, height = 25)
draw(ht_GC, heatmap_legend_side = "top", annotation_legend_side = "right",
     padding = unit(1, "cm"))
dev.off()


## compare the GO term matrices - with and without binning -----
# GO matrix - directional binning ----
load(paste0(results_dir, RObj_dir, 
            "GO_geneticCorrelations_QTLEffectsCorrelation_fdr5_HotspoteffectCorrelation_plusMinus_allP_myOwnFishersTest.rda"))

df_1 = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
colnames(df_1) = colnames(traitCommonSegregants_std)
GoSlimTerms_0 = cbind(GoSlimTerms, df_1)

for(condition_i in colnames(traitCommonSegregants_std)){
  print(condition_i)
  df_plus = GO_GC_plus[[condition_i]]
  df_plus = dplyr::filter(df_plus, p < 0.001)
  #df_plus$logFold = log2(df_plus$Count/df_plus$ExpCount)
  
  df_minus = GO_GC_minus[[condition_i]]
  df_minus = dplyr::filter(df_minus, p < 0.001)
  #df_minus$logFold = log2(df_minus$Count/df_minus$ExpCount)
  
  xPlus = df_plus$log2FC[match(row.names(GoSlimTerms_0), df_plus$GOID)] * 1
  names(xPlus) = GoSlimTerms_0$Term
  xMinus = df_minus$log2FC[match(row.names(GoSlimTerms_0), df_minus$GOID)] * -1
  names(xMinus) = GoSlimTerms_0$Term
  
  x = coalesce(xPlus, xMinus)
  
  #for GO-BP terms with nominally significant enrichment in both plus and minus, we retain the one where the fold enrichment is higher.
  rowsWithNonNA = which(!is.na(xPlus) & !is.na(xMinus))
  for(i in rowsWithNonNA){
    x[i] = ifelse(x[i] == max(abs(xPlus[i]), abs(xMinus[i])), xPlus[i], xMinus[i])
  }
  
  GoSlimTerms_0[condition_i] = x
}

row.names(GoSlimTerms_0) = GoSlimTerms_0$Term
GoSlimTerms_0 = dplyr::select(GoSlimTerms_0, -GOID, -Term, -OC)
GoSlimTerms_0 = as.matrix(GoSlimTerms_0)
GoSlimTerms_0[is.na(GoSlimTerms_0) | is.infinite(GoSlimTerms_0) | is.nan(GoSlimTerms_0)] = 0 #replace NA terms with 0

# GO matrix - no directional binning ----
df_2 = matrix(0, nrow= nrow(GoSlimTerms), ncol = length(colnames(traitCommonSegregants_std)))
colnames(df_2) = colnames(traitCommonSegregants_std)
GoSlimTerms_1 = cbind(GoSlimTerms, df_2)

for(condition_i in colnames(traitCommonSegregants_std)){
  print(condition_i)
  df_2 = GO_GC[[condition_i]]
  df_2 = dplyr::filter(df_2, p < 0.001 & log2FC > 0)
  
  x = df_2$log2FC[match(row.names(GoSlimTerms), df_2$GOID)] * 1
  names(x) = GoSlimTerms_1$Term
  
  GoSlimTerms_1[condition_i] = x
}

row.names(GoSlimTerms_1) = GoSlimTerms_1$Term
GoSlimTerms_1 = dplyr::select(GoSlimTerms_1, -GOID, -Term, -OC)
GoSlimTerms_1 = as.matrix(GoSlimTerms_1)
GoSlimTerms_1[is.na(GoSlimTerms_1) | is.infinite(GoSlimTerms_1) | is.nan(GoSlimTerms_1)] = 0 #replace NA terms with 0


# compare both the GO matrices ----
molten_df_0 = reshape2::melt(GoSlimTerms_0)
molten_df_0$dataTreatment = "signBinned"
molten_df_1 = reshape2::melt(GoSlimTerms_1)
molten_df_1$dataTreatment = "noSignBinned"

#fisher's test
t_t = length(which(molten_df_0$value != 0 & molten_df_1$value != 0))
t_f = length(which(molten_df_0$value != 0 & molten_df_1$value == 0))
f_t = length(which(molten_df_0$value == 0 & molten_df_1$value != 0))
f_f = length(which(molten_df_0$value == 0 & molten_df_1$value == 0))

fisher.test(matrix(c(t_t, t_f, f_t, f_f), nrow = 2, ncol = 2))

#paired t test
molten_df = rbind(molten_df_0, molten_df_1)
molten_df$value = abs(molten_df$value) # the signs corresponding to positive and negative genetic correlations are absoluted as they don't mean anything

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/comparingGOTerms_signedVsNoSign.pdf.pdf"), width = 25, height = 25)
ggpaired(data = molten_df,
         x = "dataTreatment",
         y = "value",
         facet.by = "Var2",
         alpha = 0.5) + stat_compare_means(paired = TRUE) 
dev.off()




