library(dplyr)
library(ggplot2)
library(ggpattern)

#global variables ---

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
plotting_dir = "plots_092522/"
RObj_dir = "RObjects_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda")) #made available on the github as R object

#make the table of trait and base-medium ----
trait_baseMediumTable = data.frame(trait = colnames(traitCommonSegregants_std))
trait_baseMediumTable$baseMedium = ifelse(grepl(trait_baseMediumTable$trait, pattern = "YPD"),
                                          "YPD", "YNB")

  ##save base medium table
ifelse(!(dir.exists(paste0(results_dir, otherFiles_dir, "supplementaryTables/"))),
  dir.create(paste0(results_dir, otherFiles_dir, "supplementaryTables/")), FALSE)

write.csv(trait_baseMediumTable, file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", "baseMediumTable.csv"), quote = FALSE, row.names = FALSE)

#functions ----
phenoR_baseMediumCorrection = function(trait_i, baseMedium){
  df = data.frame(trait = unlist(traitCommonSegregants_std[trait_i]),
                  baseMedium = unlist(traitCommonSegregants_std[baseMedium]))
  nonNA = which(complete.cases(df) == "TRUE")
  df = df[nonNA,]
  model = lm(trait~baseMedium, data = df)
  res = resid(model)
  
  #uncorrected
  phenoR_uncorrected = mclapply(unique(eQTL_Albert2018$gene),
                                FUN = function(gene_i){
                                  expression = expressionCommonSegregants_batchODCorrected_std[nonNA,gene_i]
                                  corTest = cor.test(expression, df$trait)
                                  return(data.frame(r = corTest$estimate,
                                                    p = corTest$p.value,
                                                    gene = gene_i))
                                },
                                mc.cores = parallelly::availableCores())
  phenoR_uncorrected = do.call(rbind, phenoR_uncorrected)
  phenoR_uncorrected$q = qvalue::qvalue(phenoR_uncorrected$p)$qvalues
  phenoR_uncorrected$significance = ifelse(phenoR_uncorrected$p < 0.05 & phenoR_uncorrected$q < 0.05, "TRUE", "FALSE")
  
  signifGenes_uncorrected = dplyr::filter(phenoR_uncorrected, significance == "TRUE")$gene
  
  #corrected 
  phenoR_corrected = mclapply(unique(eQTL_Albert2018$gene),
                      FUN = function(gene_i){
                      expression = expressionCommonSegregants_batchODCorrected_std[nonNA,gene_i]
                      corTest = cor.test(expression, res)
                      return(data.frame(r = corTest$estimate,
                                        p = corTest$p.value,
                                        gene = gene_i))
                  },
                      mc.cores = parallelly::availableCores())
  
  phenoR_corrected = do.call(rbind, phenoR_corrected)
  phenoR_corrected$q = qvalue::qvalue(phenoR_corrected$p)$qvalues
  phenoR_corrected$significance = ifelse(phenoR_corrected$p < 0.05 & phenoR_corrected$q < 0.05, "TRUE", "FALSE")
  
  signifGenes_corrected = dplyr::filter(phenoR_corrected, significance == "TRUE")$gene
  
  signifGenes_common = intersect(signifGenes_uncorrected, signifGenes_corrected)
  
  return(data.frame(uncorrected = length(signifGenes_uncorrected),
                    baseMedium = length(signifGenes_corrected),
                    commonGenes = length(signifGenes_common),
                    condition = trait_i))
}

#compute the number of genes before and after base medium correction ----

nGenes_beforeAfterBaseMediumCorrection = data.frame()
#nGenes_beforeAfterSelfCorrection = data.frame()

for(i in 1:nrow(trait_baseMediumTable)){
  print(i)
  trait_i = trait_baseMediumTable$trait[i]
  baseMedium_i = trait_baseMediumTable$baseMedium[i]
  
  nGenes_i = phenoR_baseMediumCorrection(trait_i, baseMedium_i)
  nGenes_beforeAfterBaseMediumCorrection = rbind(nGenes_beforeAfterBaseMediumCorrection,nGenes_i)
  
  #nGenes_self_i = phenoR_baseMediumCorrection(trait_i, trait_i)
  #nGenes_beforeAfterSelfCorrection = rbind(nGenes_beforeAfterSelfCorrection, nGenes_self_i)
}

nGenes_beforeAfterBaseMediumCorrection$proportionCommonGenes = nGenes_beforeAfterBaseMediumCorrection$commonGenes/nGenes_beforeAfterBaseMediumCorrection$uncorrected

#plot number of significant phenotypic correlations after correcting for base medium ----
melt_df = reshape2::melt(dplyr::select(nGenes_beforeAfterBaseMediumCorrection, uncorrected,
                                       baseMedium, commonGenes,condition), id.vars = c("condition", "commonGenes"))
melt_df$xend = rep(1:46, 2)


pdf(paste0(results_dir, plotting_dir, "phenotypicCorrelations/", "numberOfPhenoCorrelations_beforeAfterBaseMediumCorrection_3.pdf"), width = 15, height = 10)

ggplot() +
  geom_col(data = melt_df, aes(x = condition, y = value, fill = variable), width = 0.5, 
           position = position_dodge(), color = "black") +
  geom_col_pattern(data = melt_df, aes(x = condition, y = commonGenes, fill = variable), 
                   width = 0.5, position = position_dodge(), pattern_color = "white",
                   pattern_fill = "white", pattern = "stripe", pattern_spacing = 0.01,
                   pattern_size = 0.05, color = "black") +
  geom_segment(data = melt_df, aes(x = xend - 0.5, y = commonGenes, xend = xend +0.5, yend = commonGenes), color = "black")+
  geom_text(data = nGenes_beforeAfterBaseMediumCorrection,
            aes(x = condition, y = ifelse(uncorrected > baseMedium, uncorrected, baseMedium), 
                label = round(proportionCommonGenes, 2)),
            vjust = -1, size = 3) +
  theme_classic() +
  theme_textProperties +
  scale_fill_manual(values = c("uncorrected" = "#330000", "baseMedium" = "#CC0000")) + 
  ylab("Number of genes with \nsignificant phenotypic correlation") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()) +
  guides(fill = guide_legend(title = ""))

dev.off()


  
  
  
  
  



