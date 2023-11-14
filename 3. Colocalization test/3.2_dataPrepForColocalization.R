#install.packages("qtl2pleio")
#install.packages("qtl2")
#install.packages("qtl")
#install.packages("rlist")

#load packages
library(qtl)
library(qtl2)
library(qtl2pleio)
library(ggplot2)
library(rlist)
library(stringr)
library(dplyr)
library(tidyr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "genotypesCommonSegregants_LD95_forColoc.rda"))

#load the dummy control conditions data 
load(paste0(results_dir, RObj_dir, "dummyConditionsForColocalizationBenchmarking.rda"))
traitCommonSegregants_std = cbind(traitCommonSegregants_std, dummyConditions)


##try making the genotype matrix with our data split chromosome wise ----
#prepare data to save as a CSVs file format to prepare a cross class object
#GENOTYPE:(A)Column 1 - IDs; (B)ROW 1 - Id and marker names; (C) ROW2: blank and chrID
#PHENOTYPE: (A) Column 1 - IDs; (B) ROW 1 - Id and phenotypes; NOTE: No blanks anywhere

getGenotypesMatrix = function(genotypesCommonSegregants_LD95){
  chromosomes = paste0("chr", as.roman(c(1:16)))
  genotype_allChr = list()
  markers = colnames(genotypesCommonSegregants_LD95)
  for(chr in chromosomes){
    markers_thisChr = markers[grep(pattern = paste0(chr, ":"),
                                   markers)]
    genotypes_thisChr = dplyr::select(genotypesCommonSegregants_LD95, markers_thisChr)
    genotype_allChr = list.append(genotype_allChr, as.matrix(genotypes_thisChr))
  }
  names(genotype_allChr) = chromosomes
  
markers = str_split(colnames(genotypesCommonSegregants_LD95), pattern = ":", simplify = TRUE)
markers = as.data.frame(markers)
colnames(markers) = c("chr", "marker")
markers$bp = extract_numeric(markers$marker)
markers$cm_init = 0.4*markers$bp/1000

markers = markers %>% group_by(chr) %>% mutate(minCM = min(cm_init))
markers$cm_final = markers$cm_init - markers$minCM

# ---
genotypes = genotypesCommonSegregants_LD95
genotypes = apply(genotypes, 2,as.factor)
genotypes = as.data.frame(genotypes)
genotypes$ID = paste0("S",c(1:nrow(genotypes)))
genotypes = dplyr::select(genotypes, ID, everything())

row1 = c("", as.integer(as.factor(str_remove(markers$chr, pattern = "chr"))))
row2 = c("", markers$cm_final)
genotypes = rbind(row1, row2, genotypes)
return(genotypes)
}

#genotypes matrix for LD = 0.95 -- for coarse scan
genotypes_LD95 = getGenotypesMatrix(genotypesCommonSegregants_LD95)

#phenotypes matrix
phenotypes = cbind(paste0("S", c(1:nrow(genotypesCommonSegregants_LD95))),
                   traitCommonSegregants_std,
                   expressionCommonSegregants_batchODCorrected_std)
colnames(phenotypes)[1] = "ID"


#save ----
  #create any reqd directories
ifelse(!dir.exists(file.path(paste0(results_dir, otherFiles_dir))), 
       dir.create(file.path(paste0(results_dir, otherFiles_dir))),
       FALSE)
  #genotypes - LD95
  write.csv(genotypes_LD95, file = paste0(results_dir, otherFiles_dir, "genotypes_LD95.csvs"), 
          quote = FALSE, row.names = FALSE)
  #genotypes - LD98
#  write.csv(genotypes_LD98, file = paste0(results_dir, otherFiles_dir, "genotypes_LD98.csvs"), 
#         quote = FALSE, row.names = FALSE)

  #phenotypes 
  write.csv(phenotypes, file = paste0(results_dir, otherFiles_dir, "phenotype.csvs"), 
          quote = FALSE, row.names = FALSE)


