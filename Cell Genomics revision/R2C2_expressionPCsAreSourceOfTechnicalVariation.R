## Reviewer 2, Comment 2: Expression PCs are a source of technical variation? 

# Make a scree plot for expression data 
library(psych)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

# get the expression PCs and make scree plot

  #load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))

  #do PCA ----
fit = prcomp(expressionCommonSegregants_batchODCorrected_std, center = TRUE, scale = TRUE)
Y = summary(fit)

Y_t = as.data.frame(t(Y$importance))
Y_t$eigenValue = Y_t$`Standard deviation`^2
Y_t$nPCs = 1:nrow(Y_t)

  # table showing proportion of variance explained by first 10 PCs
write.csv( Y_t[1:10,], 
           file = paste0(results_dir, otherFiles_dir, "revision_tables_cellGenomics/", "R2C2_proportionOfVarianceExplainedByFirst10PCs.csv"), 
           quote = FALSE, row.names = FALSE)

  #plot scree plot ----
pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/","PC-ScreePlot_expression.pdf"))

ggplot(data = Y_t,
       aes(x = nPCs, y = log10(eigenValue)))+
  #geom_line(linetype = 'dotted') +
  geom_point(aes(#fill = ifelse(eigenValue > 1, "black", "lightgray"),
    color = ifelse(eigenValue > 1, "deeppink4", "lightgray")),
    shape = 21, alpha = 0.75) +
  scale_fill_manual(values = c("deeppink4" = "deeppink4", "lightgray" = "lightgray")) +
  scale_color_manual(values = c("deeppink4" = "deeppink4", "lightgray" = "lightgray")) +
  ylim(c(-1,3)) +
  geom_hline(yintercept = log2(1), linetype = "dashed") +
  guides(fill = "none", color = "none") +
  theme(legend.position = "none") +
  xlab("Number of Principal components") +
  ylab("Eigen value") +
  
  theme_classic()

ggplot(data = Y_t,
       aes(x = nPCs, y = `Cumulative Proportion`))+
  #geom_line(linetype = 'dotted') +
  geom_point(aes(#fill = ifelse(eigenValue > 1, "black", "lightgray"),
    color = ifelse(eigenValue > 1, "deeppink4", "lightgray")),
    shape = 21, alpha = 0.75) +
  scale_fill_manual(values = c("deeppink4" = "deeppink4", "lightgray" = "lightgray")) +
  scale_color_manual(values = c("deeppink4" = "deeppink4", "lightgray" = "lightgray")) +
  guides(fill = "none", color = "none") +
  theme(legend.position = "none") +
  xlab("Number of Principal components") +
  ylab("Cumulative proportion of variance explained") +
  theme_classic()

dev.off()


#regress out the first 10 PCs and recompute the genetic correlations ----- DO THIS to check for technical covariates vs biological covariates -----

# First 3 PCs account for 50% of variance. Let's explore what about the first 10 PCs -----
  #function to regress out first 10 PCs from expression of gene
pcr = function(df){
  model = lm(expression~., data = df)
  expression_residualized = residuals(model)
  return(expression_residualized)
}

GC_pcOut = lapply(colnames(traitCommonSegregants_std),
       FUN = function(condition_i){
         print(condition_i)
         growth = unlist(traitCommonSegregants_std[condition_i])
         GC_pcOut = lapply(colnames(expressionCommonSegregants_batchODCorrected_std),
                           FUN = function(gene_i){
                             print(gene_i)
                             expression_i = unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])
                             
                             pcs_df = as.data.frame(fit$x)
                             pcs_df_10 = pcs_df[1:10]
                             pcs_df_10$expression = expression_i
                             
                             #get expression after regressing out PCs 
                             expression_pcOut = pcr(pcs_df_10)
                             
                             #compute GC 
                              #uncorrected
                             GC = cor.test(growth, expression_i)
                             r = as.numeric(GC$estimate)
                             p = GC$p.value
                              #corrected
                             GC_pcOut = cor.test(growth, expression_pcOut)
                             r_pcOut = as.numeric(GC_pcOut$estimate)
                             p_pcOut = GC_pcOut$p.value
                             
                             return(c("r" = r, "p" = p,
                                      "r_pcOut" = r_pcOut, "p_pcOut" = p_pcOut))
                           }) 
         
         GC_pcOut = as.data.frame(do.call(rbind, GC_pcOut))
         GC_pcOut$gene = colnames(expressionCommonSegregants_batchODCorrected_std)
         GC_pcOut$condition = condition_i
         return(GC_pcOut)
       })
GC_pcOut = lapply(GC_pcOut,
                  FUN = function(df){
                    df$q = qvalue::qvalue(df$p)$q
                    df$q_pcOut = qvalue::qvalue(df$p_pcOut)$q
                    return(df)
                  })

GC_pcOut_all = as.data.frame(do.call(rbind, GC_pcOut))
nSigGenes = GC_pcOut_all %>% group_by(condition) %>% summarize(before = sum(p< 0.05 & q < 0.05),
                                                               after = sum(p_pcOut < 0.05 & q_pcOut<0.05))

nSigGenes_melt = reshape2::melt(nSigGenes, id.vars = "condition")

pdf(paste0(results_dir, plotting_dir, 
           "revision_figures_CellGenomics/","nSignificantGenesAfter10PCsOut.pdf"))

ggbarplot(nSigGenes_melt,
          x = 'condition',
          y = 'value',
          fill = 'variable')+ theme_textProperties + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of Significant genes @ FDR < 5%")

ggpaired(nSigGenes_melt,
         x = 'variable',
         y = 'value') + stat_compare_means(method = "wilcox.test", paired = TRUE) + theme_textProperties + ylab("Number of Significant genes @ FDR < 5%")

ggscatter(GC_pcOut_all,
          x = "r",
          y = "r_pcOut",
          shape = 21,
          alpha = 0.3) + 
  xlab("Genetic correlation coefficient before regressing out first 10 PCs") +
  ylab("Genetic correlation coefficient after regressing out first 10 PCs") +
  stat_cor() + theme_bw() + theme_textProperties
  

dev.off()

# How much are the first 10 PCs correlated with the 46 growth conditions ----
PCs = Y$x[,1:10]
GC_PCs = lapply(1:10, function(i){
  PC_i = unlist(PCs[,i])
  gc_pc_i = lapply(colnames(traitCommonSegregants_std),
                   function(condition_i){
                     growth_i = unlist(traitCommonSegregants_std[condition_i])
                     ct = cor.test(PC_i, growth_i)
                     return(c("r" = as.numeric(ct$estimate), 
                              "p" = ct$p.value, 
                              "condition" = condition_i,
                              "PC" = i))
                   }) 
  gc_pc_i = as.data.frame(do.call(rbind, gc_pc_i))
  return(gc_pc_i)
})

GC_PCs = as.data.frame(do.call(rbind, GC_PCs))
GC_PCs$r = as.numeric(GC_PCs$r)
GC_PCs$p = as.numeric(GC_PCs$p)
GC_PCs$PC = as.factor(GC_PCs$PC)
GC_PCs$q = qvalue::qvalue(GC_PCs$p)$q

# plot

  #for heatmap (GC vs PC)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100) #color palette

GC_PC_mat = dplyr::filter(GC_PCs, p < 0.05 & q < 0.05)
GC_PC_mat = dplyr::select(GC_PC_mat, condition, PC, r)

GC_PC_mat = reshape2::dcast(GC_PC_mat, PC~condition, value.var = "r")
GC_PC_mat$PC = as.numeric(GC_PC_mat$PC)
row.names(GC_PC_mat) = paste0("PC", GC_PC_mat$PC)
GC_PC_mat = dplyr::select(GC_PC_mat, -PC)

GC_PC_mat = t(GC_PC_mat)
GC_PC_mat[is.na(GC_PC_mat)] = 0
GC_PC_mat[!is.finite(GC_PC_mat)] = 0


pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "PCsVsGrowth_differentConditions.pdf"), width = 10, height = 10)
#pheatmap(as.matrix(go_mat), color = my_palette, breaks = breaks, cluster_cols = FALSE)
pheatmap::pheatmap(GC_PC_mat, color = my_colors, breaks = seq(-0.5, 0.5, length.out = 101), cluster_cols = FALSE)

ggplot(data = GC_PCs) +
  geom_point(aes(x = PC, y = r, color = ifelse(p < 0.05 & q < 0.05, "deeppink3", "lightgray"))) +
  scale_color_manual(values = c("deeppink3" = "deeppink3", "lightgray" = "lightgray")) +
  facet_wrap(~condition) +
  guides(color = "none") +
  scale_x_discrete(breaks = levels(GC_PCs$PC)) +
  theme_bw()


dev.off()



# what genes are enriched for each PC? ------
#for each PC loading, you pick out the top 50% genes. Then, do a fisher's test for GO enrichment 

library(org.Sc.sgd.db)

# A data frame with gene names and their PC loadings for set of PCs
pc_loadings = as.data.frame(fit$rotation)

# load GO term data
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))

  # get genes belonging to each GO term -----

gn2go = mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db), "GOALL", "ORF", multiVals = "list")

## filter out NA mappings
gn2go = gn2go[!sapply(gn2go, function(x) all(is.na(x)))]

#shortlist only the GO-Slim terms
gn2GO_slim = lapply(gn2go, function(gene_i){
  slim_index = which(gene_i %in% row.names(GoSlimTerms))
  return(gene_i[slim_index])
})

GO2gn = lapply(row.names(GoSlimTerms), FUN = function(slimTerm){
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

  # function to perform our own GO enrichment using a fisher's exact test ----
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

  # Perform enrichment analysis for each PC ----
go_enrichment = lapply(1:10,
                       FUN = function(pc){
                         gene_names = row.names(pc_loadings)
                         
                         #split the pc loadings based on positive and negative first
                         pc_loadings_plus = pc_loadings[which(pc_loadings[pc] >= 0),]
                         pc_loadings_minus = pc_loadings[which(pc_loadings[pc] < 0),]
                         
                         #now get the top 10% genes in each of these sets
                         sorted_pc_loadings_plus = pc_loadings_plus %>%
                           arrange(desc(abs(pc_loadings_plus[pc])))
                         top_genes_plus = row.names(sorted_pc_loadings_plus)[1:(nrow(sorted_pc_loadings_plus)/10)]
                         
                         sorted_pc_loadings_minus = pc_loadings_minus %>%
                           arrange(desc(abs(pc_loadings_minus[pc])))
                         top_genes_minus = row.names(sorted_pc_loadings_minus)[1:(nrow(sorted_pc_loadings_minus)/10)]
                        
                          #do the GO enrichment for these sets now
                         
                          #plus
                         GO_i_plus = lapply(row.names(GoSlimTerms), 
                                       FUN = function(go_term){
                                                GO_FT_i = GO_FT(go_term, top_genes_plus)
                                                return(GO_FT_i)
                                              })
                         
                         GO_i_plus = as.data.frame(do.call(rbind, GO_i_plus))
                         GO_i_plus$q = p.adjust(GO_i_plus$p, method = "fdr")
                         table_plus = as.data.frame(cbind(GoSlimTerms, GO_i_plus))
                         table_plus$PC = pc 
                         
                          #minus
                         GO_i_minus = lapply(row.names(GoSlimTerms), 
                                            FUN = function(go_term){
                                              GO_FT_i = GO_FT(go_term, top_genes_minus)
                                              return(GO_FT_i)
                                            })
                         
                         GO_i_minus = as.data.frame(do.call(rbind, GO_i_minus))
                         GO_i_minus$q = p.adjust(GO_i_minus$p, method = "fdr")
                         table_minus = as.data.frame(cbind(GoSlimTerms, GO_i_minus))
                         table_minus$PC = pc 
                         
                         #table
                          table = list(table_plus, table_minus)
                          names(table) = c("plus", "minus")
                         return(table)
                       })

names(go_enrichment) = paste0("PC",1:10)

  # plot enrichment results -----
library(pheatmap)
library(ComplexHeatmap)
  ## color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-3, 3, length.out = 101)

  ##Only GO-BP terms
GoSlimTerms = dplyr::filter(GoSlimTerms, OC == "BP")

  ##prepare matrix to plot the log2FC GO enrichment
df = matrix(0, nrow= nrow(GoSlimTerms), ncol = 10)
colnames(df) = names(go_enrichment)
GoSlimTerms_mat = cbind(GoSlimTerms, df)

for(pc in 1:10){
  go_enrichment_i = go_enrichment[[paste0("PC", pc)]]
  go_plus = go_enrichment_i[["plus"]]
  go_minus = go_enrichment_i[["minus"]]
  
  #filter for significant terms
  go_plus = dplyr::filter(go_plus, p < 0.001 & log2FC > 0)
  go_minus = dplyr::filter(go_minus, p < 0.001 & log2FC > 0)
  
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
  
  GoSlimTerms_mat[paste0("PC", pc)] = x
  
}

row.names(GoSlimTerms_mat) = GoSlimTerms_mat$Term
GoSlimTerms_mat = dplyr::select(GoSlimTerms_mat, -GOID, -Term, -OC)
GoSlimTerms_mat = as.matrix(GoSlimTerms_mat)
GoSlimTerms_mat[is.na(GoSlimTerms_mat) | is.infinite(GoSlimTerms_mat) | is.nan(GoSlimTerms_mat)] = 0 #replace NA terms with 0

  #plot heatmap with the row order in current Figure 6

  ## plot top heatmap of Fig 6 again ---
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


  #get row order of above heatmap
ro1 = row_order(draw(ht_GC))

  #plot GO-PC enrichments using the above row order. we are again splitting into 3 clusters 
ht_PC = Heatmap(GoSlimTerms_mat, col = col,
             use_raster = FALSE, show_row_names = TRUE,
             #width = unit(40, "cm"),
             #height = unit(50, "cm"),
             show_row_dend = FALSE,
             show_column_dend = TRUE,
             row_order = ro1,
             #column_dend_gp = gpar(col = col_dend$annot$grp_color),
             #column_split = 2,
             #column_km = 3,
             #cluster_columns = col_dend$dendrogram,
             #row_names_gp = gpar(fontsize = 20),
             #column_names_gp = gpar(fontsize = 22),
             #top_annotation = ha2,
             #bottom_annotation = ha3,
             #cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_title_rot = 0, 
             name = "GO fold enrichment (count/expected count):", 
             row_names_side = "left",
             heatmap_legend_param = list(direction = "horizontal", side = "bottom",
                                         title = "GO term fold enrichment:"),
             row_gap = unit(1, "mm"))
  
  # get the average PC loadings for Brauer and ESR groups -----
  # get brauer and ESR genes
Brauer_S1 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"Brauer_S1.xls"),
                               col_types = c("text", "text", "text", "text", "numeric", 
                                             "numeric","numeric", "text", "text", "text",
                                             "numeric", "numeric", "numeric", "numeric", 
                                             "numeric", "text")) #supplementary table 1 from Brauer et al. 2008

up = Brauer_S1$ORF[which(Brauer_S1$`Growth Rate Response @1.5SD`=="up")]
down = Brauer_S1$ORF[which(Brauer_S1$`Growth Rate Response @1.5SD`=="down")]
unresponsive = Brauer_S1$ORF[which(Brauer_S1$`Unresponsive?` == "yes")]

up_ESR = Brauer_S1$ORF[which(Brauer_S1$ESR =="up")]
down_ESR = Brauer_S1$ORF[which(Brauer_S1$ESR =="down")]


pc_brauer = lapply(1:10,
                   FUN = function(pc){
                     gene_names = row.names(pc_loadings)
                     pc = unlist(pc_loadings[pc])
                     names(pc) = gene_names
                     
                     up_pc = mean(pc[up],na.rm = TRUE)
                     down_pc = mean(pc[down], na.rm = TRUE)
                     unresponsive_pc = mean(pc[unresponsive], na.rm = TRUE)
                     up_ESR_pc = mean(pc[up_ESR], na.rm = TRUE)
                     down_ESR_pc = mean(pc[down_ESR], na.rm = TRUE)
                     
                     return(c("up" = up_pc,
                              "down" = down_pc,
                              "unresponsive" = unresponsive_pc,
                              "up_ESR" = up_ESR_pc,
                              "down_ESR" = down_ESR_pc))
                   })

pc_brauer = do.call(rbind, pc_brauer)
row.names(pc_brauer) = paste0("PC", 1:10)
pc_brauer = t(pc_brauer)
 
  #plot of brauer ESR genes avg PC loadings
ht1_brauerESR_PC = Heatmap(pc_brauer, 
                           col =col2, use_raster = FALSE,
                           row_order = factor(row.names(pc_brauer), 
                                              levels = unique(row.names(pc_brauer))),
                           cluster_columns = FALSE,
                           show_heatmap_legend = TRUE, 
                           height = unit(2, "cm"), 
                           heatmap_legend_param = list(direction = "horizontal", side = "left"),
                           show_row_dend = FALSE,
                           name = "av. pc loading",
                           row_gap = unit(3, "mm"),
                           column_gap = unit(5, "mm"))

  #plot the GO enrichments for the expression PCs

ht_PC_GO = ht_PC %v% ht1_brauerESR_PC

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "ExpressionPCs_GOEnrichment.pdf"), width = 10, height = 15)
draw(ht_PC_GO, heatmap_legend_side = "top",padding = unit(2, "cm"))
dev.off()

                   
  # correlate this GO matrix for PCs with the GO matrix for the 46 traits to see if the PCs equal to states. Should give you another heatmap -----

 #pairwise correlation
pairwise_correlation = cor(GoSlimTerms_mat, GoSlimTerms_mat_GC, method = "spearman")

pairwise_correlation[is.na(pairwise_correlation) | is.infinite(pairwise_correlation) | is.nan(pairwise_correlation)] = 0 #replace NA terms with 0

library(pheatmap)
pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "pairwiseCorrelation_PCVsTrait_GOEnrichmentBased.pdf"), width = 10, height = 10)
ht_pairwise = Heatmap(pairwise_correlation, #col = col,
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
                cluster_rows = FALSE,
                cluster_columns = TRUE,
                row_title_rot = 0, 
                name = "GO fold enrichment (count/expected count):", 
                row_names_side = "left",
                heatmap_legend_param = list(direction = "horizontal", side = "bottom",
                                            title = "rho:"),
                row_gap = unit(1, "mm"))

draw(ht_pairwise, heatmap_legend_side = "top", padding = unit(1, "cm"))

dev.off()

# # # correlate the dendrograms of GC- GO enrichment cluster and PC - GO enrichment cluster
#  dend_GC = hclust(dist(as.matrix(t(GoSlimTerms_mat_GC)))) #clustering of traits based on GO on GC
#  dend_PC = hclust(dist(as.matrix(GoSlimTerms_mat))) # clustering of GO terms based on GO on PC
#  dend_pairwiseCorrelations = hclust(dist(as.matrix(t(pairwise_correlation)))) #clustering of traits based on correlation of GO logFC enrichments using GC and PC
#  
#  dend_list = dendlist(dend_GC, dend_pairwiseCorrelations)
#   dend_list %>%
#     untangle(method = "step1side") %>% # Find the best alignment layout
#     tanglegram()    
#   cor.dendlist(dend_list, method = "cophenetic")
#   
# # dend_list = dendlist(dend_GO, dend_PC)
# # 
# 
# # 
# # cor.dendlist(dend_list, method = "cophenetic")
# # cor.dendlist(dend_list, method = "baker")
# 
# dist1 = dist(as.matrix(GoSlimTerms_mat_GC))
# dist2 = dist(as.matrix(GoSlimTerms_mat))
# dist3 = dist(pairwise_correlation)
# 
# cor.test(dist1, dist2, method = "spearman")$p.value
