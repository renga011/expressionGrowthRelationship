#!/usr/bin/env Rscript --vanilla

## Comment 10: The rho between GC and QEC could be inflated -----
#split the data in half. Estimate GCs with one half and QECs with the other

#change home directory ----
setwd("/home/albertf/renga011/ExpressionPhenotypeProject/")

## GLOBAL VARIABLES -----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#input
i = as.numeric(commandArgs(trailingOnly = TRUE))

# split the data in half ----
set.seed(i)
index_split = sample(1:nrow(genotypesCommonSegregants), size = 0.5*nrow(genotypesCommonSegregants))

# Data for GC ----

trait_GC = traitCommonSegregants_std[index_split,]
expression_GC = expressionCommonSegregants_batchODCorrected_std[index_split,]

# Data for QEC ----

genotypes_QEC = genotypesCommonSegregants[-index_split,]
trait_QEC = traitCommonSegregants_std[-index_split,]
expression_QEC = expressionCommonSegregants_batchODCorrected_std[-index_split,]

#functions to compute GC and QEC ----

computeGC = function(trait_data, expression_data){
  r_expressionGrowth = function(gene, condition){
    traitValues = unlist(trait_data[condition])
    expressionForThisGene = unlist(expression_data[,gene])
    cor_expressionGrowth = cor.test(expressionForThisGene, traitValues)
    return(list(r = cor_expressionGrowth$estimate,
                p = cor_expressionGrowth$p.value))
  }
  
  phenotypicCorrelations = parallel::mclapply(colnames(trait_data),
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
                                              }, 
                                              mc.cores = parallelly::availableCores())
  
  names(phenotypicCorrelations) = colnames(trait_data)
  
  phenotypicCorrelations_df = lapply(colnames(trait_data),
                                     FUN = function(x){
                                       df = phenotypicCorrelations[[x]]
                                       df$condition = x
                                       return(df)
                                     })
  
  phenotypicCorrelations_df = as.data.frame(do.call(rbind, phenotypicCorrelations_df))
  
  return(phenotypicCorrelations_df)
}

computeQEC = function(trait_data, expression_data, genotype_data){
  # Calculate growth effects at all eQTLs of gene
  computeRGrowth = lapply(colnames(trait_data),
                          FUN = function(condition_i){
                            growth_i = unlist(trait_data[condition_i])
                            df_i = parallel::mclapply(1:nrow(eQTL_Albert2018),
                                                      FUN = function(i){
                                                        pmarker_i = eQTL_Albert2018$pmarker[i]
                                                        genotype_i = unlist(genotype_data[pmarker_i])

                                                        gene_i = eQTL_Albert2018$gene[i]
                                                        expression_i = unlist(expression_data[gene_i])

                                                        ct_expression = cor.test(expression_i, genotype_i)
                                                        r_expression = ct_expression$estimate
                                                        ci_expression = ct_expression$conf.int[2] - ct_expression$conf.int[1]

                                                        ct_i = cor.test(growth_i, genotype_i)
                                                        r_growth = ct_i$estimate
                                                        p_growth = ct_i$p.value
                                                        r_growth_ci_left = ct_i$conf.int[1]
                                                        r_growth_ci_right = ct_i$conf.int[2]
                                                        ci_growth = ct_i$conf.int[2] - ct_i$conf.int[1]
                                                        return(c("r_expression" = r_expression,
                                                                 "ci_expression" = ci_expression,
                                                                 "r_growth" = r_growth,
                                                                 "p_growth" = p_growth,
                                                                 "r_growth_ci_left" = r_growth_ci_left,
                                                                 "r_growth_ci_right" = r_growth_ci_right,
                                                                 "ci_growth" = ci_growth))
                                                      },
                                                      mc.cores = parallelly::availableCores())

                            # df_i = lapply(1:nrow(eQTL_Albert2018),
                            #                           FUN = function(i){
                            #                             pmarker_i = eQTL_Albert2018$pmarker[i]
                            #                             genotype_i = unlist(genotype_data[pmarker_i])
                            #                             
                            #                             gene_i = eQTL_Albert2018$gene[i]
                            #                             expression_i = unlist(expression_data[gene_i])
                            #                             
                            #                             ct_expression = cor.test(expression_i, genotype_i)
                            #                             r_expression = ct_expression$estimate
                            #                             ci_expression = ct_expression$conf.int[2] - ct_expression$conf.int[1]
                            #                             
                            #                             ct_i = cor.test(growth_i, genotype_i)
                            #                             r_growth = ct_i$estimate
                            #                             p_growth = ct_i$p.value
                            #                             r_growth_ci_left = ct_i$conf.int[1]
                            #                             r_growth_ci_right = ct_i$conf.int[2]
                            #                             ci_growth = ct_i$conf.int[2] - ct_i$conf.int[1]
                            #                             return(c("r_expression" = r_expression,
                            #                                      "ci_expression" = ci_expression,
                            #                                      "r_growth" = r_growth,
                            #                                      "p_growth" = p_growth,
                            #                                      "r_growth_ci_left" = r_growth_ci_left,
                            #                                      "r_growth_ci_right" = r_growth_ci_right,
                            #                                      "ci_growth" = ci_growth))
                            #                           })
                            df_i = as.data.frame(do.call(rbind, df_i))
                            eQTL_Albert2018 = cbind(eQTL_Albert2018, df_i)
                            
                            return(eQTL_Albert2018)
                          })
  names(computeRGrowth) = colnames(trait_data)
  
  # compute QEC
  computeQEC = function(eQTLTable, condition_i){
    qec_i = sapply(unique(eQTLTable$gene),
                   FUN = function(gene_i){
                     df_i = dplyr::filter(eQTLTable, gene == gene_i)
                     wts = (df_i$ci_expression*df_i$ci_growth)^-1
                     w = try({weights::wtd.cor(x = df_i$r_expression, y = df_i$r_growth.cor, 
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
  
  QEC = parallel::mclapply(colnames(trait_data),
               FUN = function(condition_i){
                 print(condition_i)
                 table_i = computeRGrowth[[condition_i]]
                 QEC_i = computeQEC(table_i, condition_i)
                 return(QEC_i)
               }, mc.cores = parallelly::availableCores())
  names(QEC) = colnames(trait_data)
  
  QEC = lapply(colnames(trait_data),
               FUN = function(condition_i){
                 table_i = QEC[[condition_i]]
                 q = p.adjust(table_i$p_QEC, method = "fdr")
                 table_i$q_QEC = q
                 return(table_i)
               })
                                       
  QEC_df = as.data.frame(do.call(rbind, QEC))
  return(QEC_df)
}

#compute GC ----
GC_i = computeGC(trait_GC, expression_GC)

#compute QEC ----
QEC_i = computeQEC(trait_QEC, expression_QEC, genotypes_QEC)

# collate into list and save
GCQEC = list("GC" = GC_i, "QEC" = QEC_i)

# save
save(GCQEC, file = paste0(results_dir, RObj_dir, "QECVsGC_5050/", "GC_QEC_", i, ".rda"))

