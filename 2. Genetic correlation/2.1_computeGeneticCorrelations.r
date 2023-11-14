library(dplyr)

#global variables ---

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"

#load data -----
load(paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))

#compute phenotypic correlations -----
r_expressionGrowth = function(gene, condition){
  traitValues = unlist(traitCommonSegregants_std[condition])
  expressionForThisGene = unlist(expressionCommonSegregants_batchODCorrected_std[,gene])
  cor_expressionGrowth = cor.test(expressionForThisGene, traitValues)
  return(list(r = cor_expressionGrowth$estimate,
              p = cor_expressionGrowth$p.value))
}

phenotypicCorrelations = lapply(colnames(traitCommonSegregants_std),
                                FUN = function(x){
                                  print(x)
                                  df = parallel::mclapply(unique(eQTL_Albert2018$gene),
                                                          FUN = r_expressionGrowth,
                                                          condition = x, 
                                                          mc.cores = parallelly::availableCores())
                                  df = as.data.frame(do.call(rbind, df))
                                  df = as.data.frame(apply(df, 2, FUN = unlist))
                                  df = as.data.frame(apply(df, 2, FUN = as.numeric))
                                  df$gene = unique(eQTL_Albert2018$gene)
                                  df$q = qvalue::qvalue(df$p)$qvalue
                                  return(df)
                                })

names(phenotypicCorrelations) = colnames(traitCommonSegregants_std)

#save -----

save(phenotypicCorrelations, file = paste0(results_dir, RObj_dir,
                                           "phenotypicCorrelationTable_pearson.rda"))
 
