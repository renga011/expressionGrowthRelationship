## analyze mediation results

library(dplyr)
library(org.Sc.sgd.db)
library(readr)
library(ggpubr)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "IRA2__Hydrogen_Peroxide.rda"))
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))

#some stats ----
##number of significant genes at 5% FDR
significantMediators_IRA2 = mediationEstimatesForAllHotspotTargets$gene[which(mediationEstimatesForAllHotspotTargets$p < 0.05 & mediationEstimatesForAllHotspotTargets$fdr < 0.05)]
length(significantMediators_IRA2)
# 380

#GO analyses ----

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

##do your own fisher's gene enrichment test for GO term enrichment
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

##apply the above function to the list of genes with significant mediation
 GO_i = lapply(row.names(GoSlimTerms), FUN = function(GOTerm_i){
                                      GO_FT_i = GO_FT(GOTerm_i,
                                                      significantMediators_IRA2)
                                      return(GO_FT_i)
                                    })
 GO_i = as.data.frame(do.call(rbind, GO_i))
 GO_i$q = p.adjust(GO_i$p, method = "fdr")
 GO_IRA2_H2O2_significantMediators = as.data.frame(cbind(GoSlimTerms, GO_i))
 
 #write the mediation table as a supplementary table - Supp Table 6
 library(openxlsx)
 write.xlsx(GO_IRA2_H2O2_significantMediators, file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", "Supplementary Table 7 - GO_IRA2_H2O2_significantMediators.xlsx"))

 
 #MSN2 targets analysis
 write(significantMediators_IRA2, file = paste0(results_dir, otherFiles_dir, "IRA2_significantMediators.txt"))
 
 #analyze this on the yeastract database for Msn2 enrichment
 #now load the results and add details to the mediation estimates table - that becomes Supp. Table 6
 MSN2_targets_IRA2_significantMediators = read_delim(paste0(results_dir, otherFiles_dir,
                                                    "MSN2_targets_IRA2_significantMediators.csv"), delim = ";", escape_double = FALSE, trim_ws = TRUE)
 
 documentedTargets = stringr::str_split(MSN2_targets_IRA2_significantMediators$`Documented Targets`[1],pattern = " ", simplify = TRUE)[1,]
 
 documentedTargets_sysName = geneList_SGD$sysName[match(documentedTargets, geneList_SGD$stdName)]
 
 mediationEstimatesForAllHotspotTargets$gene_stdName = geneList_SGD$stdName[(match(mediationEstimatesForAllHotspotTargets$gene, geneList_SGD$sysName))]
 
 mediationEstimatesForAllHotspotTargets$MSN2_documentedTarget = ifelse(mediationEstimatesForAllHotspotTargets$gene %in% documentedTargets_sysName, TRUE, FALSE)
 
 #get average proportion mediated for MSN2 targets
 summary(abs(filter(mediationEstimatesForAllHotspotTargets, MSN2_documentedTarget == "TRUE")$propMediated))

 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.02514 0.06403 0.09593 0.12387 0.14686 0.51741 
 
 #save as Supplementary Table 7
 mediationEstimatesForAllHotspotTargets = select(mediationEstimatesForAllHotspotTargets, gene, gene_stdName, MSN2_documentedTarget, everything())
 
 write.xlsx(mediationEstimatesForAllHotspotTargets,
            file = paste0(results_dir, otherFiles_dir, "supplementaryTables/",
                          "Supplementary Table 6 - mediationResults_IRA2_H2O2.xlsx"))
 
 # plot
 df_i = dplyr::filter(mediationEstimatesForAllHotspotTargets, p < 0.05 & fdr < 0.05 | gene_stdName == "MSN2")
  
  df_i$label = ifelse(df_i$gene_stdName %in% c("MSN2", "PNC1", "TPS2", "HSP12"), TRUE, FALSE)
  
  df_i = df_i %>%
    arrange(propMediated)
  
  # df_i$cumulative_fraction = 1:nrow(df_i)
  # df_i$cumulative_frequency = cumsum(df_i$cumulative_fraction)/sum(df_i$cumulative_fraction)
  
  df_i_MSN2 = dplyr::filter(df_i, MSN2_documentedTarget == "TRUE")
  df_i_noMSN2 = dplyr::filter(df_i, MSN2_documentedTarget == "FALSE")
  
  df_i_MSN2$cumulative_fraction = 1:nrow(df_i_MSN2)
  df_i_MSN2$cumulative_frequency =cumsum(df_i_MSN2$cumulative_fraction)/sum(df_i_MSN2$cumulative_fraction)
  
  df_i_noMSN2$cumulative_fraction = 1:nrow(df_i_noMSN2)
  df_i_noMSN2$cumulative_frequency =cumsum(df_i_noMSN2$cumulative_fraction)/sum(df_i_noMSN2$cumulative_fraction)
  
  df_i = rbind(df_i_MSN2, df_i_noMSN2)
  
  pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "R3C10_fig7BRevised.pdf"))
  
  ggplot(df_i, aes(x = propMediated, y = cumulative_frequency, 
                   color = factor(MSN2_documentedTarget), 
                   size = factor(label))) +
    geom_point(shape = 21) +
   #geom_jitter(position = position_jitter(width = 0.02, height = 0.02), shape = 21)+
    scale_color_manual(values = c("TRUE" = "deeppink4", "FALSE" = "blue4")) +  # Customize colors if needed
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2)) +
    geom_vline(xintercept = 0, color = "black") +
    #facet_wrap(~MSN2_documentedTarget) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      x = "propMediated",
      y = "Cumulative Distribution",
      color = "MSN2 Documented Target",
      size = "Label"
    )
  
  dev.off()
  