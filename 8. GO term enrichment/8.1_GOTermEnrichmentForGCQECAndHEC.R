#get GO terms and genes mapping to GO term

library(org.Sc.sgd.db)
library(dplyr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(file = paste0(results_dir, RObj_dir, "QTLEffectsCorrelation_atEQTLs.rda"))
#load(paste0(results_dir, RObj_dir, "QTLOverlap_CorrelationOfQTLRs.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))
load(paste0(results_dir, RObj_dir, "hotspotEffectCorrelations.rda"))
load(paste0(results_dir, RObj_dir,"GoSlimTerms.rda"))
  GoSlimTerms$GOID = row.names(GoSlimTerms)
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))

#

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

##do your own fisher's gene enrichment test
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

GO_GC_plus = lapply(colnames(traitCommonSegregants_std),
                    FUN = function(condition_i){
                      print(condition_i)
                      df = phenotypicCorrelations[[condition_i]]
                      sigGenes_i = df$gene[which(df$p < 0.05 & df$q < 0.05 & df$r >0)]
                      
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

names(GO_GC_plus) = colnames(traitCommonSegregants_std)

GO_GC_minus = lapply(colnames(traitCommonSegregants_std),
                     FUN = function(condition_i){
                       print(condition_i)
                       df = phenotypicCorrelations[[condition_i]]
                       sigGenes_i = df$gene[which(df$p < 0.05 & df$q < 0.05 & df$r <0)]
                       
                       GO_i = lapply(row.names(GoSlimTerms),
                                     FUN = function(GOTerm_i){
                                       GO_FT_i = GO_FT(GOTerm_i,
                                                       sigGenes_i)
                                       return(GO_FT_i)
                                     })
                       GO_i = as.data.frame(do.call(rbind, GO_i))
                       table = as.data.frame(cbind(GoSlimTerms, GO_i))
                       table$condition = condition_i
                       return(table)
                     })

names(GO_GC_minus) = colnames(traitCommonSegregants_std)

GO_EQTLEC_plus = lapply(colnames(traitCommonSegregants_std),
                        FUN = function(condition_i){
                          print(condition_i)
                          df = QEC[[condition_i]]
                          sigGenes_i = df$gene[which(df$q_QEC < 0.2 & df$r_QEC > 0)]
                          
                          GO_i = lapply(row.names(GoSlimTerms),
                                        FUN = function(GOTerm_i){
                                          GO_FT_i = GO_FT(GOTerm_i,
                                                          sigGenes_i)
                                          return(GO_FT_i)
                                        })
                          GO_i = as.data.frame(do.call(rbind, GO_i))
                          table = as.data.frame(cbind(GoSlimTerms, GO_i))
                          table$condition = condition_i
                          return(table)
                        })

names(GO_EQTLEC_plus) = colnames(traitCommonSegregants_std)

GO_EQTLEC_minus = lapply(colnames(traitCommonSegregants_std),
                         FUN = function(condition_i){
                           print(condition_i)
                           df = QEC[[condition_i]]
                           sigGenes_i = df$gene[which(df$p_QEC < 0.05 & df$q_QEC < 0.2 & df$r_QEC < 0)]
                           GO_i = lapply(row.names(GoSlimTerms),
                                         FUN = function(GOTerm_i){
                                           GO_FT_i = GO_FT(GOTerm_i,
                                                           sigGenes_i)
                                           return(GO_FT_i)
                                         })
                           GO_i = as.data.frame(do.call(rbind, GO_i))
                           table = as.data.frame(cbind(GoSlimTerms, GO_i))
                           table$condition = condition_i
                           return(table)
                         })

names(GO_EQTLEC_minus) = colnames(traitCommonSegregants_std)

GO_HEC_plus = lapply(colnames(traitCommonSegregants_std),
                     FUN = function(condition_i){
                       print(condition_i)
                       df = hotspotEffectCorrelations[[condition_i]]
                       colnames(df) = c("r", "p", "gene", "q")
                       sigGenes_i = df$gene[which(df$p < 0.05 & df$q < 0.05 & df$r >0)]
                       
                       GO_i = lapply(row.names(GoSlimTerms),
                                     FUN = function(GOTerm_i){
                                       GO_FT_i = GO_FT(GOTerm_i,
                                                       sigGenes_i)
                                       return(GO_FT_i)
                                     })
                       GO_i = as.data.frame(do.call(rbind, GO_i))
                       table = as.data.frame(cbind(GoSlimTerms, GO_i))
                       table$condition = condition_i
                       return(table)
                     })

names(GO_HEC_plus) = colnames(traitCommonSegregants_std)

GO_HEC_minus = lapply(colnames(traitCommonSegregants_std),
                      FUN = function(condition_i){
                        print(condition_i)
                        df = hotspotEffectCorrelations[[condition_i]]
                        colnames(df) = c("r", "p", "gene", "q")
                        sigGenes_i = df$gene[which(df$p < 0.05 & df$q < 0.05 & df$r <0)]
                        
                        GO_i = lapply(row.names(GoSlimTerms),
                                      FUN = function(GOTerm_i){
                                        GO_FT_i = GO_FT(GOTerm_i,
                                                        sigGenes_i)
                                        return(GO_FT_i)
                                      })
                        GO_i = as.data.frame(do.call(rbind, GO_i))
                        table = as.data.frame(cbind(GoSlimTerms, GO_i))
                        table$condition = condition_i
                        return(table)
                      })

names(GO_HEC_minus) = colnames(traitCommonSegregants_std)


save(GO_GC_plus, GO_GC_minus, GO_EQTLEC_plus, GO_EQTLEC_minus, 
     GO_HEC_plus, GO_HEC_minus,
     file = paste0(results_dir, RObj_dir, "GO_geneticCorrelations_QTLEffectsCorrelation_fdr5_HotspoteffectCorrelation_plusMinus_allP_myOwnFishersTest.rda"))

