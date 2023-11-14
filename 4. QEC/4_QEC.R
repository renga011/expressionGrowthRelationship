##code to compute the correlation between the r_expression and r_growth at the eQTLs for each gene..

library(dplyr)
library(ggpubr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

# Calculate growth effects at all eQTLs of gene
computeRGrowth = lapply(colnames(traitCommonSegregants_std),
                        FUN = function(condition_i){
                          growth_i = unlist(traitCommonSegregants_std[condition_i])
                          
                          for(i in 1:nrow(eQTL_Albert2018)){
                            pmarker_i = eQTL_Albert2018$pmarker[i]
                            genotype_i = unlist(genotypesCommonSegregants[pmarker_i])

                            gene_i = eQTL_Albert2018$gene[i]
                            expression_i = unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])

                            ct_expression = cor.test(expression_i, genotype_i)
                            eQTL_Albert2018$r_expression[i] = ct_expression$estimate
                            eQTL_Albert2018$ci_expression[i] = ct_expression$conf.int[2] - ct_expression$conf.int[1]

                            ct_i = cor.test(growth_i, genotype_i)
                            eQTL_Albert2018$r_growth[i] = ct_i$estimate
                            eQTL_Albert2018$p_growth[i] = ct_i$p.value
                            eQTL_Albert2018$r_growth_ci_left[i] = ct_i$conf.int[1]
                            eQTL_Albert2018$r_growth_ci_right[i] = ct_i$conf.int[2]
                            eQTL_Albert2018$ci_growth[i] = ct_i$conf.int[2] - ct_i$conf.int[1]
                          }
                          
                          return(eQTL_Albert2018)
                        })
names(computeRGrowth) = colnames(traitCommonSegregants_std)

# compute QEC
computeQEC = function(eQTLTable, condition_i){
  qec_i = sapply(unique(eQTLTable$gene),
                 FUN = function(gene_i){
                   df_i = dplyr::filter(eQTLTable, gene == gene_i)
                   wts = (df_i$ci_expression*df_i$ci_growth)^-1
                   w = try({weights::wtd.cor(x = df_i$r_expression, y = df_i$r_growth, 
                                             weight = wts)})
                   if(is(w, 'try-error')|is(w, 'error')|nrow(df_i)==2){
                     r = NA
                     p = NA
                   } else {
                     r = as.numeric(w[1, "correlation"])
                     p = as.numeric(w[1, "p.value"])
                   }
                   return(c("r_QEC" = r, "p_QEC" = p))
                 })
  qec_i = as.data.frame(t(qec_i))
  qec_i$gene = unique(eQTLTable$gene)
  qec_i$condition = condition_i
  
  return(qec_i)
  
}

QEC = lapply(colnames(traitCommonSegregants_std),
             FUN = function(condition_i){
               table_i = computeRGrowth[[condition_i]]
               QEC_i = computeQEC(table_i, condition_i)
               return(QEC_i)
             })
names(QEC) = colnames(traitCommonSegregants_std)

QEC = lapply(colnames(traitCommonSegregants_std),
             FUN = function(condition_i){
               table_i = QEC[[condition_i]]
               q = p.adjust(table_i$p_QEC, method = "fdr")
               table_i$q_QEC = q
               return(table_i)
             })

QEC_all = as.data.frame(do.call(rbind, QEC))

nQTL_summary = eQTL_Albert2018 %>% group_by(gene) %>% summarise(n = n())
QEC_all$nQTL = nQTL_summary$n[match(QEC_all$gene, nQTL_summary$gene)]

## save ----
write.csv(QEC_all, file = paste0(results_dir, otherFiles_dir, "tableOfQTLEffectsCorrelation_atEQTLs.csv"), quote = FALSE, row.names = FALSE)

save(QEC, file = paste0(results_dir, RObj_dir, "QTLEffectsCorrelation_atEQTLs.rda"))

