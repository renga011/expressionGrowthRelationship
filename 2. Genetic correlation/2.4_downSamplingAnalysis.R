#power analysis for phenotypic correlations----
library(dplyr)
library(ggpubr)

#global variables ---

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
plotting_dir = "plots_092522/"
RObj_dir = "RObjects_092522/"

#data -----

load(paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda")) #available in R objects folder on github
load(paste0(results_dir, RObj_dir, "colorPaletteForConditions.rda"))
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#plot histogram of the pearson Rs ----
phenotypicCorrelations = lapply(colnames(traitCommonSegregants_std), 
                         FUN = function(condition_i){
                           df_i = phenotypicCorrelations[[condition_i]]
                           df_i$condition = condition_i
                           return(df_i)
                         })
df_allPearsonRs = do.call(rbind, phenotypicCorrelations)
df_allPearsonRs$absR = abs(df_allPearsonRs$r)
df_allPearsonRs$significance = ifelse(df_allPearsonRs$p < 0.05 & df_allPearsonRs$q < 0.05, "TRUE", "FALSE")

pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations_distributionOfRs.pdf"), width = 15, height = 10)
gghistogram(df_allPearsonRs,
          x = "absR",
          facet.by = "condition") + theme_light() +
  theme_textProperties + xlim(c(0,0.3))

ggplot(df_allPearsonRs, aes(x = absR, color = condition, alpha = significance)) +
  stat_ecdf() + #facet_wrap(~condition) + 
  scale_color_manual(values = colorPalette) +
  ylab("Cumulative probability") +
  xlab("|r|") +
  #geom_hline(aes(yintercept = 0.5)) +
  theme_light() +
  theme_textProperties

dev.off()

#different sample sizes, get number of significant phenotypic correlations and the median abs r's for significant phenotypic correlations ------

sampleSizes = c(250, 500, 750, 979)
seeds = c(1:5)

r_expressionGrowth = function(gene, condition, expressionMatrix_i, traitMatrix_i){
  traitValues = unlist(traitMatrix_i[condition])
  expressionForThisGene = unlist(expressionMatrix_i[,gene])
  cor_expressionGrowth = cor.test(expressionForThisGene, traitValues)
  return(list(r = cor_expressionGrowth$estimate,
              p = cor_expressionGrowth$p.value))
}

computePhenotypicCorrelations = function(expressionMatrix, traitMatrix){
  phenotypicCorrelations_i = lapply(colnames(traitMatrix),
                                  FUN = function(x){
                                    print(x)
                                    df = parallel::mclapply(unique(eQTL_Albert2018$gene),
                                                            FUN = r_expressionGrowth,
                                                            condition = x,
                                                            expressionMatrix_i = expressionMatrix,
                                                            traitMatrix_i = traitMatrix,
                                                            mc.cores = parallelly::availableCores())
                                    df = as.data.frame(do.call(rbind, df))
                                    df = as.data.frame(apply(df, 2, FUN = unlist))
                                    df = as.data.frame(apply(df, 2, FUN = as.numeric))
                                    df$gene = unique(eQTL_Albert2018$gene)
                                    df$q = qvalue::qvalue(df$p)$qvalue
                                    return(df)
                                  })
  
  names(phenotypicCorrelations_i) = colnames(traitMatrix)
  return(phenotypicCorrelations_i)
}
pheno_diffSizes = lapply(sampleSizes, function(sampleSize){
  print(sampleSize)
  pheno_multipleSeeds = lapply(seeds, FUN = function(seed_i){
    print(seed_i)
    set.seed(seed_i)
    sample_i = sort(sample(1:nrow(expressionCommonSegregants_batchODCorrected_std), sampleSize))
    phenotypicCorrelations_i = computePhenotypicCorrelations(expressionMatrix = expressionCommonSegregants_batchODCorrected_std[sample_i,], traitMatrix = traitCommonSegregants_std[sample_i,])
    names(phenotypicCorrelations_i) = colnames(traitCommonSegregants_std)
    return(phenotypicCorrelations_i)
  })
  names(pheno_multipleSeeds) = seeds
 
  #save
  save(pheno_multipleSeeds, file = paste0(results_dir, RObj_dir,
                                             "phenotypicCorrelationTable_",
                                             sampleSize,"_pearson.rda"))
  
  #get summary statistics of the phenotypic correlations - nSignificantGenes, median(absR), max(absR), sd(absR)
  summary_i = lapply(seeds, FUN = function(seed_i){
    phenotypicCorrelations_i = pheno_multipleSeeds[[seed_i]]
    df_i = lapply(names(phenotypicCorrelations_i),
                  FUN = function(condition_i){
                    df_i = phenotypicCorrelations_i[[condition_i]]
                    df_i = dplyr::filter(df_i, p< 0.05 & q < 0.05) #filter significant genes
                    if(nrow(df_i) !=0){
                      median_i = median(abs(df_i$r))
                      min_i = min(abs(df_i$r))
                      max_i = max(abs(df_i$r))
                      sd_i = sd(abs(df_i$r))
                    } else {
                      median_i = NA
                      min_i = NA
                      max_i = NA
                      sd_i = NA
                    }
                    summary_df = data.frame(nSegregants = sampleSize,
                                            seed = seed_i,
                                            condition = condition_i,
                                            nSignificantGenes = nrow(df_i),
                                            absR_median = median_i,
                                            absR_min = min_i,
                                            absR_max = max_i,
                                            absR_sd = sd_i)
                    return(summary_df)
                  })
    df_i = as.data.frame(do.call(rbind, df_i))
    
    return(df_i)
  })
  summary_i = as.data.frame(do.call(rbind, summary_i))
  
  return(summary_i)
  
})

pheno_diffSizes = lapply(sampleSizes,
                         FUN = function(sampleSize){
                           load(paste0(results_dir, RObj_dir,
                                       "phenotypicCorrelationTable_",
                                       sampleSize,"_pearson.rda"))
                           summary_allSeeds = lapply(seeds, FUN = function(seed_i){
                             phenotypicCorrelations_i = pheno_multipleSeeds[[seed_i]]
                             
                             summary_i = lapply(names(phenotypicCorrelations_i),
                                                FUN = function(condition_i){
                                                  df_i = phenotypicCorrelations_i[[condition_i]]
                                                  df_i = dplyr::filter(df_i, p< 0.05 & q < 0.05) #filter significant genes
                                                  if(nrow(df_i) !=0){
                                                    median_i = median(abs(df_i$r))
                                                    min_i = min(abs(df_i$r))
                                                    max_i = max(abs(df_i$r))
                                                    sd_i = sd(abs(df_i$r))
                                                  } else {
                                                    median_i = NA
                                                    min_i = NA
                                                    max_i = NA
                                                    sd_i = NA
                                                  }
                                                  summary_df = data.frame(nSegregants = sampleSize,
                                                                          condition = condition_i,
                                                                          seed = seed_i,
                                                                          nSignificantGenes = nrow(df_i),
                                                                          absR_median = median_i,
                                                                          absR_min = min_i,
                                                                          absR_max = max_i,
                                                                          absR_sd = sd_i)
                                                  return(summary_df)
                                                })
                             summary_i = as.data.frame(do.call(rbind, summary_i))
                             return(summary_i)
                           })
                           summary_allSeeds = do.call(rbind, summary_allSeeds)
                           return(summary_allSeeds)
                         })

pheno_diffSizes = do.call(rbind,pheno_diffSizes)
pheno_diffSizes$seed = as.factor(pheno_diffSizes$seed)

plottingDf_diffSizes = reshape2::melt(dplyr::select(pheno_diffSizes, -absR_max, -absR_sd,
                                                    -nSignificantGenes),
                                      id.vars = c("nSegregants", "condition", "seed"))

pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations/", "phenotypicCorrelations_sampleSizeAnalyses_2.pdf"), width = 15, height = 10)

# lapply(colnames(pheno_diffSizes)[-c(1:3)],
#        FUN = function(measure){
#          plot = ggline(pheno_diffSizes,
#                 x = "nSegregants",
#                 y = measure,
#                 facet.by = "condition",
#                 color = "seed",
#                 title = measure,
#                 xlab = "Number of segregants") + theme_light() + theme_textProperties
#          
#        })

ggplot(data = pheno_diffSizes) +
  geom_point(mapping = aes(x = nSegregants, y = nSignificantGenes, shape = as.factor(seed))) +
  stat_summary(geom = "line", fun = median, aes(x = nSegregants, y = nSignificantGenes),
               color = "steel blue") +
  facet_wrap(~condition) + theme_bw() + theme_textProperties +
  scale_x_continuous(name = "No. of segregants",
                     breaks = c(250, 500,750, 979))

ggplot(data = plottingDf_diffSizes) +
  geom_point(mapping = aes(x = nSegregants, y = value, shape = as.factor(seed),
                           color = variable)) +
  stat_summary(geom = "line", fun = median, aes(x = nSegregants, y = value, 
                                                color = variable)) +
  scale_color_manual(values = c("absR_median" = "black", "absR_min" = "steel blue")) +
  facet_wrap(~condition) + theme_bw() + theme_textProperties +
  scale_x_continuous(name = "No. of segregants",
                     breaks = c(250, 500,750, 979))

dev.off()


a = pheno_diffSizes %>% group_by(condition, nSegregants) %>% summarise(median_associations = median(nSignificantGenes), median_medianAbsR = median(absR_median), median_min = median(absR_min))
#####

summary(pheno_diffSizes$absR_sd[which(pheno_diffSizes$nSegregants == 979)])
load(paste0(results_dir, "geneList_SGD.rda"))
geneList_SGD$essential = FALSE
geneList_SGD$essential[match(SGD_essentialGenes$X2, geneList_SGD$sysName)] = TRUE

essentialGeneList = geneList_SGD$sysName[which(geneList_SGD$essential == "TRUE")]
df_allPearsonRs$essential = FALSE
df_allPearsonRs$essential[match(essentialGeneList, df_allPearsonRs$gene)] = TRUE



