# compare correlations across genetic, QTL, hotspot correlations
# directories are edited for consistency

### The output table from this code will have some extra columns corresponding to other analyses we performed (eg: FET on the signs of the correlations in the set of genes with sigificant genetic and QTL effect correlation; All analyses comparing genetic correlations and QTL effect correlations at a stricter FDR of 5%). These have been deleted from the published version of this table so that the latter is more reflective of the analyses discussed in the paper. However, the results displayed in the paper don't change in light of the extra analyses performed here and the overall conclusions remain consistent. We therefore publish the code with these extra analyses for more transparency.

library(tidyverse)
library(tidyr)
#library(readxl)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

allGeneResults <- read.table(file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", 
                                           "table_dataAcrossTests_updatedWithNewQECAndHECNumbers.csv"), sep = ",", header = TRUE)

# how many QTL effect correlations have ≥3 overlap points?
length(which(allGeneResults$nEQTL >= 3))
# 238,556

# how many unqiue genes in this?
length(unique(allGeneResults$gene[allGeneResults$nEQTL >= 3]))
# 5186

# distribution in this subset
summary(allGeneResults$nEQTL[allGeneResults$nEQTL >= 3])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.00    5.00    6.00    6.89    9.00   21.00 


# across everything, how many sig trait/gene pairs?
length(which(allGeneResults$nEQTL >= 3 & allGeneResults$q_QTLEffects <= 0.2))
# 50 at 5%FDR, 2038 at 20%FDR
# how many unique genes?
length(unique(allGeneResults[allGeneResults$nEQTL >= 3 & allGeneResults$q_QTLEffects <= 0.05,"gene"]))
# 42 genes at 5%; 1551 at 20%


# compare q-values at different minima of nOverlaps

# for each trait, how many genes have FDR=5% for the geneticQTL cors
QTLcorsFDR5_number <- sapply(unique(allGeneResults$condition), function(trait){
  length(unique(allGeneResults$gene[allGeneResults$q_QTLEffects <= 0.05 & allGeneResults$condition == trait & allGeneResults$nEQTL >= 3 & !(is.na(allGeneResults$q_QTLEffects))]))
})
summary(QTLcorsFDR5_number)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.000   0.000   1.087   1.000  17.000 

# and 20%
QTLcorsFDR20_number <- sapply(unique(allGeneResults$condition), function(trait){
  length(unique(allGeneResults$gene[allGeneResults$q_QTLEffects <= 0.2 & allGeneResults$condition == trait & allGeneResults$nEQTL >= 3 & !(is.na(allGeneResults$q_QTLEffects))]))
})
summary(QTLcorsFDR20_number)
# 20%
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0     0.0     5.0    44.3    15.0  1106.0 

length(which(QTLcorsFDR20_number > 0))
# 32 traits
summary(QTLcorsFDR20_number[QTLcorsFDR20_number > 0])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    4.50    9.50   63.69   27.25 1106.00 

# what are the correlation coefficients that go with this?
summary(sapply(unique(allGeneResults$condition), function(trait){
  median(abs(allGeneResults$r_QTLEffects[allGeneResults$q_QTLEffects <= 0.05 & allGeneResults$condition == trait & allGeneResults$nEQTL >= 3 & !(is.na(allGeneResults$q_QTLEffects))]))
}))
# 5% FDR:
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.8949  0.9421  0.9651  0.9577  0.9752  1.0000      30 
 #20% FDR:
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.8359  0.8950  0.9163  0.9211  0.9498  1.0000      14 
 #very high, as expected



# gather all the relevant comparisons between QTL effects and genetic cors
QTL_and_Genetic_Comparison <- sapply(unique(allGeneResults$condition), function(trait){
  subDat <- allGeneResults[allGeneResults$condition == trait & !(is.na(allGeneResults$q_QTLEffects)) & allGeneResults$nEQTL >= 3,]
  ret1 <- sapply(c("0.05", "0.2"), function(FDR){
    FDR <- as.numeric(FDR) # needed to preserve FDRs as column names
    selectTheseGenes <- subDat$gene[subDat$q_QTLEffects <= FDR]
    nGenesForCor <- length(selectTheseGenes)
    subDat4cor <- subDat[subDat$gene %in% selectTheseGenes,]
    QTL_genetic_Cor_result <- c(nrow(subDat), nGenesForCor, NA, NA)
    try({
      QTL_genetic_Cor <- cor.test(subDat4cor[, "r_genetic"], subDat4cor[, "r_QTLEffects"])
      QTL_genetic_Cor_result <- c(nrow(subDat), nGenesForCor, QTL_genetic_Cor$est, QTL_genetic_Cor$p.value)
      })
    names(QTL_genetic_Cor_result) <- c("nGenesWithData", "nGenesForCor", "QTLGeneticCor_r", "QTLGeneticCor_p")

    #FET for whether signs agree
    FET_matrix_sign <- cbind(
      c(length(which(subDat4cor$r_genetic > 0 & subDat4cor$r_QTLEffects > 0)), length(which(subDat4cor$r_genetic > 0 & subDat4cor$r_QTLEffects < 0))),
      c(length(which(subDat4cor$r_genetic < 0 & subDat4cor$r_QTLEffects > 0)), length(which(subDat4cor$r_genetic < 0 & subDat4cor$r_QTLEffects < 0)))
    )
    FET_sign <- fisher.test(FET_matrix_sign, alternative = "greater")
    FET_result_sign <- c(FET_matrix_sign, FET_sign$estimate, FET_sign$p.value)
    names(FET_result_sign) <- c("n_FET_sign_++", "n_FET_sign_+-", "n_FET_sign_-+", "n_FET_sign_--", "FET_sign_odds", "FET_sign_p")
    
    
    #FET for whether sig in both
    FET_matrix_sig <- cbind(
      c(length(which(subDat$q_QTLEffects <= FDR & subDat$q_genetic <= 0.05)), length(which(subDat$q_QTLEffects > FDR & subDat$q_genetic <= 0.05))),
      c(length(which(subDat$q_QTLEffects <= FDR & subDat$q_genetic > 0.05)), length(which(subDat$q_QTLEffects > FDR & subDat$q_genetic > 0.05)))
    )
    FET_sig <- fisher.test(FET_matrix_sig, alternative = "greater")
    FET_result_sig <- c(FET_matrix_sig, FET_sig$estimate, FET_sig$p.value)
    names(FET_result_sig) <- c("n_FET_sig_both", "n_FET_sig_G", "n_FET_sig_Q", "n_FET_sig_neither", "FET_sig_odds", "FET_sig_p")
    
    #ret <- (c("nGenes_QTL"=nGenes, QTL_genetic_Cor_result, FET_result))
    ret <- (c(QTL_genetic_Cor_result, FET_result_sign, FET_result_sig))
    return(ret)
  })
  # gymnastics below to get result into one matrix as output
  ret1 <- as.data.frame(ret1)
  ret1$resultType <- rownames(ret1)
  ret1 <- pivot_longer(ret1, cols = !resultType, values_to = "result")
  ret2 <- ret1$result
  names(ret2) <- paste(ret1$resultType, ret1$name, sep="_")
  ret2
})
#write.table(t(QTL_and_Genetic_Comparison), file="QTL_and_Genetic_Comparison_230731.txt", sep="\t", quote=FALSE)

# let's analyze this:
# how many sig tests (out of 46?)

sapply(c("QTLGeneticCor_p_0.05", "QTLGeneticCor_p_0.2", "FET_sign_p_0.05", "FET_sign_p_0.2", "FET_sig_p_0.05", "FET_sig_p_0.2"), function(thisTest){
  res <- table(QTL_and_Genetic_Comparison[thisTest,] < 0.05, useNA="always")
  # this is tortured because there are not always trues AND falses
  ret <- c("TRUE" = 0, "FALSE" = 0, "NA" = NA)
  if(!is.na(res["TRUE"])){ret["TRUE"] <- res["TRUE"]}
  if(!is.na(res["FALSE"])){ret["FALSE"] <- res["FALSE"]}
  ret["NA"] <- res["NA"]
  ret
})
#       QTLGeneticCor_p_0.05 QTLGeneticCor_p_0.2 FET_sign_p_0.05 FET_sign_p_0.2 FET_sig_p_0.05 FET_sig_p_0.2
#TRUE                     2                  24               1             18              3            21
#FALSE                    4                   4              45             28             43            25
#NA                      NA                  NA              NA             NA             NA            NA
# quite a bit of significance here

length(which(QTL_and_Genetic_Comparison["nGenesForCor_0.2",] > 2))
# 28 traits with at least 3 genes analyzed
length(which(QTL_and_Genetic_Comparison["nGenesForCor_0.2",] > 2 & QTL_and_Genetic_Comparison["FET_sig_p_0.2",] < 0.05 & !(is.na(QTL_and_Genetic_Comparison["FET_sig_p_0.2",]))))
# 21

# how many positive cors betwenn QEC and genetic?
length(which(QTL_and_Genetic_Comparison["nGenesForCor_0.2",] > 2 & QTL_and_Genetic_Comparison["QTLGeneticCor_p_0.2",] < 0.05 & QTL_and_Genetic_Comparison["QTLGeneticCor_r_0.2",] > 0))
# 22

# what if went across all traits?
##no filters ---
cor.test(allGeneResults[, "r_genetic"], allGeneResults[, "r_QTLEffects"], method = "spearman")
# rho = 0.69, p < 2e-16

##filter of >= 3 eQTLs and significance filter for QTL effects cor --- ---
subDatAcrossTraits <- allGeneResults[!(is.na(allGeneResults$q_QTLEffects)) & allGeneResults$nEQTL >= 3,]

#subDat4CorAcrossTraits <- subDatAcrossTraits[subDatAcrossTraits$q_QTLEffects <= 0.05,]
subDat4CorAcrossTraits <- subDatAcrossTraits[subDatAcrossTraits$q_QTLEffects <= 0.2,]
# 50 observations across traits for 5%, 2038 for 20%
cor.test(subDat4CorAcrossTraits[, "r_genetic"], subDat4CorAcrossTraits[, "r_QTLEffects"])
# 5%: r=0.91, p = 2e-11
#20%: r=0.88, p < 2e-16

###########################
# hotspot vs genetic

# how many sig genes?
length(which(allGeneResults$q_hotspotEffects <= 0.05))
# 64,023
# how many unique genes?
length(unique(allGeneResults[allGeneResults$q_hotspotEffects <= 0.05,"gene"]))
# 5284

# per trait
hotspotCorsFDR5_number <- sapply(unique(allGeneResults$condition), function(trait){
  length(unique(allGeneResults$gene[allGeneResults$q_hotspotEffects <= 0.2 & allGeneResults$condition == trait & !(is.na(allGeneResults$q_hotspotEffects))]))
})
summary(hotspotCorsFDR5_number)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.0   287.5  1650.5  1391.8  2149.2  3650.0


# comparison to genetic cors
# loop is copied from above; just ignore the 20% FDR
hotspot_and_Genetic_Comparison <- sapply(unique(allGeneResults$condition), function(trait){
  subDat <- allGeneResults[allGeneResults$condition == trait & !(is.na(allGeneResults$q_hotspotEffects)),]
  #ret1 <- sapply(c("0.05", "0.2"), function(FDR){
  ret1 <- sapply(c("0.05"), function(FDR){
    FDR <- as.numeric(FDR) # needed to preserve FDRs as column names
    selectTheseGenes <- subDat$gene[subDat$q_hotspotEffects <= FDR]
    nGenesForCor <- length(selectTheseGenes)
    subDat4cor <- subDat[subDat$gene %in% selectTheseGenes,]
    QTL_genetic_Cor_result <- c(nrow(subDat), nGenesForCor, NA, NA)
    try({
      QTL_genetic_Cor <- cor.test(subDat4cor[, "r_genetic"], subDat4cor[, "r_hotspotEffects"])
      QTL_genetic_Cor_result <- c(nrow(subDat), nGenesForCor, QTL_genetic_Cor$est, QTL_genetic_Cor$p.value)
    })
    names(QTL_genetic_Cor_result) <- c("nGenesWithData", "nGenesForCor", "hotspotGeneticCor_r", "hotspotGeneticCor_p")
    
    # ALL genes, not the subset with q < 0.05 for hotspots
    QTL_genetic_Cor_result_ALL <- c(NA, NA)
    try({
      QTL_genetic_Cor_ALL <- cor.test(subDat[, "r_genetic"], subDat[, "r_hotspotEffects"])
      QTL_genetic_Cor_result_ALL <- c(QTL_genetic_Cor_ALL$est, QTL_genetic_Cor_ALL$p.value)
    })
    names(QTL_genetic_Cor_result_ALL) <- c("hotspotGeneticCor_r_ALL", "hotspotGeneticCor_p_ALL")
    
        
    #FET for whether sig in both
    FET_matrix_sig <- cbind(
      c(length(which(subDat$q_hotspotEffects <= FDR & subDat$q_genetic <= 0.05)), length(which(subDat$q_hotspotEffects > FDR & subDat$q_genetic <= 0.05))),
      c(length(which(subDat$q_hotspotEffects <= FDR & subDat$q_genetic > 0.05)), length(which(subDat$q_hotspotEffects > FDR & subDat$q_genetic > 0.05)))
    )
    FET_sig <- fisher.test(FET_matrix_sig, alternative = "greater")
    FET_result_sig <- c(FET_matrix_sig, FET_sig$estimate, FET_sig$p.value)
    names(FET_result_sig) <- c("n_FET_sig_both", "n_FET_sig_G", "n_FET_sig_Q", "n_FET_sig_neither", "FET_sig_odds", "FET_sig_p")
    
    #ret <- (c(QTL_genetic_Cor_result, FET_result_sig))
    ret <- (c(QTL_genetic_Cor_result, QTL_genetic_Cor_result_ALL, FET_result_sig))
    return(ret)
  })
  # gymnastics below to get result into one matrix as output
  ret1 <- as.data.frame(ret1)
  ret1$resultType <- rownames(ret1)
  ret1 <- pivot_longer(ret1, cols = !resultType, values_to = "result")
  ret2 <- ret1$result
  names(ret2) <- paste(ret1$resultType, ret1$name, sep="_")
  ret2
})
#write.table(t(hotspot_and_Genetic_Comparison), file="hotspot_and_Genetic_Comparison_230728.txt", sep="\t", quote=FALSE)

sapply(c("hotspotGeneticCor_p_0.05", "hotspotGeneticCor_p_ALL_0.05", "FET_sig_p_0.05"), function(thisTest){
  res <- table(hotspot_and_Genetic_Comparison[thisTest,] < 0.05, useNA="always")
  # this is tortured because there are not always trues AND falses
  ret <- c("TRUE" = 0, "FALSE" = 0, "NA" = NA)
  if(!is.na(res["TRUE"])){ret["TRUE"] <- res["TRUE"]}
  if(!is.na(res["FALSE"])){ret["FALSE"] <- res["FALSE"]}
  ret["NA"] <- res["NA"]
  ret
})
#      hotspotGeneticCor_p_0.05 hotspotGeneticCor_p_ALL_0.05 FET_sig_p_0.05
#TRUE                        45                           46             45
#FALSE                        1                            0              1
#NA                          NA                           NA             NA

# very strong agreement!


# summary of p-values for the cors
summary(hotspot_and_Genetic_Comparison["hotspotGeneticCor_p_ALL_0.05",])
# all reported as 0

summary(hotspot_and_Genetic_Comparison["hotspotGeneticCor_r_ALL_0.05",])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5636  0.8074  0.8597  0.8370  0.8927  0.9372

# what if went across all traits?
cor.test(allGeneResults[, "r_genetic"], allGeneResults[, "r_hotspotEffects"], method = "spearman")
# rho = 0.86, p < 2e-16

##
QTL_and_Genetic_Comparison = as.data.frame(t(QTL_and_Genetic_Comparison))
  QTL_and_Genetic_Comparison$condition = row.names(QTL_and_Genetic_Comparison)
hotspot_and_Genetic_Comparison = as.data.frame(t(hotspot_and_Genetic_Comparison))
  hotspot_and_Genetic_Comparison$condition = row.names(hotspot_and_Genetic_Comparison)
sheet3 = list(select(QTL_and_Genetic_Comparison, condition, everything()),
              select(hotspot_and_Genetic_Comparison, condition, everything()))

names(sheet3) = c("QTL_and_Genetic_Comparison", "hotspot_and_Genetic_Comparison")

library(openxlsx)
write.xlsx(sheet3, 
           file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", "Supplementary Table 4 - QTLEffectCorVsGeneticCor_HotspotEffectCorVsGeneticCor.xlsx"))




