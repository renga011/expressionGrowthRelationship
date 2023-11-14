#code to prep the data. Source data from Albert et al. (2018) & Bloom et al. (2013)
library(dplyr) 
library(readr)
library(readxl)
library(tidyr)
library(rlang)
library(stringr)
library(gtools)
library(bindata)
library(caret)
library(gridExtra)

#global variables ---

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#check/create directories ---
ifelse(!dir.exists(file.path(results_dir)), 
       dir.create(file.path(results_dir)), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir))), 
       dir.create(file.path(results_dir, RObj_dir)), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir))), 
       dir.create(file.path(results_dir, plotting_dir)), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, otherFiles_dir))), 
       dir.create(file.path(results_dir, otherFiles_dir)), FALSE)

#load datasets ----
genotypes_all <- read_delim(paste0(results_dir, otherFiles_dir, "SI_Data_03_genotypes.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)  #AB18 - source data 3
  colnames(genotypes_all)[1] = "wellID"

traitData_all <- read_delim(paste0(results_dir, otherFiles_dir,"BYxRM_PhenoData_2013.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)  #from Bloom et al. 2013 - supplementary data
  colnames(traitData_all)[1] = "wellID"

expression1000segregants <- read_delim(paste0(results_dir, otherFiles_dir, "SI_Data_01_expressionValues.txt"), 
                                       "\t", escape_double = FALSE, trim_ws = TRUE) # AB18 - source data 1
 colnames(expression1000segregants)[1] = "wellID"

covariates_2018 <- read_csv(paste0(results_dir, otherFiles_dir,"covariates_2018.csv")) #AB18 - source data 2
  colnames(covariates_2018)[1] = "wellID"

BloomQTL = read_csv(paste0(results_dir, otherFiles_dir, "BloomQTL.csv")) #Bloom 2013 - supplement 5
eQTL_Albert2018 = read_csv(paste0(results_dir, otherFiles_dir,"1000SegregantsExpressionQTLs.csv")) #AB18 - source data 4

load(paste0(results_dir, RObj_dir, "scanoneLODS_OD_stranded.rdata")) #contains the reduced set of 11530 markers

###process datasets -----

#retrieve segregants common to both studies and subset ---

expression1000segregants$wellID = stringr::str_split_fixed(expression1000segregants$wellID, pattern = "-", n = 2)[,1]

covariates_2018$wellID = stringr::str_split_fixed(covariates_2018$wellID, pattern = "-", n = 2)[,1]

commonWellIDs = intersect(expression1000segregants$wellID, traitData_all$wellID)

expressionCommonSegregants = dplyr::filter(expression1000segregants, wellID %in% commonWellIDs)
  expressionCommonSegregants = dplyr::select(expressionCommonSegregants, -wellID)

traitCommonSegregants = dplyr::filter(traitData_all, wellID %in% commonWellIDs)
  traitCommonSegregants = dplyr::select(traitCommonSegregants, -wellID)

genotypesCommonSegregants = dplyr::filter(genotypes_all, wellID %in% commonWellIDs)
  genotypesCommonSegregants = dplyr::select(genotypesCommonSegregants, -wellID)
  
covariatesCommonSegregants = dplyr::filter(covariates_2018, wellID %in% commonWellIDs)
  covariatesCommonSegregants = dplyr::select(covariatesCommonSegregants, -wellID)

# genotypes - reduce number of markers from 42052 to 11530 markers which aren't in perfect LD 
genotypesCommonSegregants = genotypesCommonSegregants[which(colnames(genotypesCommonSegregants) %in% colnames(scanoneLODS.OD$r))]

#batch and OD correct expression ---
batch = covariatesCommonSegregants$batch
OD = covariatesCommonSegregants$OD_covariate

expressionCommonSegregants_batchODCorrected = sapply(expressionCommonSegregants,
                                        FUN = function(x){
                                          model = lm(x ~ batch + OD)
                                          x_batchODCorrected = residuals(model)
                                          return(x_batchODCorrected)
                                        })
expressionCommonSegregants_batchODCorrected = as.data.frame(expressionCommonSegregants_batchODCorrected)

# standardize expression and trait---
preProcFunction_expression = caret::preProcess(expressionCommonSegregants_batchODCorrected, method = c("center", "scale"))
expressionCommonSegregants_batchODCorrected_std = predict(preProcFunction_expression, expressionCommonSegregants_batchODCorrected)

preProcFunction_trait= caret::preProcess(traitCommonSegregants, method = c("center", "scale"))
traitCommonSegregants_std = predict(preProcFunction_trait, traitCommonSegregants)

#change Magnesium sulfate to manganese sulfate on the traitCommonSegregants table
colnames(traitCommonSegregants_std)[which(colnames(traitCommonSegregants_std) == "Magnesium_Sulfate")] = "Manganese_Sulfate"

#add the 'a' beside trait names that begin with a number 
x = c("4-Hydroxybenzaldehyde", "5-Fluorocytosine", "4NQO", "6-Azauracil", "5-Fluorouracil")
BloomQTL$Trait[which(BloomQTL$Trait %in% x)] = paste0("a", BloomQTL$Trait[which(BloomQTL$Trait %in% x)])

colnames(traitCommonSegregants_std)[which(colnames(traitCommonSegregants_std) %in% paste0("x", x))] = paste0("a", x)

traitCommonSegregants_std = traitCommonSegregants_std[, order(colnames(traitCommonSegregants_std))]


#get the closest peak marker for the different growth QTLs and put in bloomQTL table
pMarkerPositions = extract_numeric(colnames(genotypesCommonSegregants))
pMarkerChromosomes = unlist(str_split(colnames(genotypesCommonSegregants), pattern = ":"))[c(TRUE, FALSE)]
pMarkerChromosomes =  roman2int(str_replace(pMarkerChromosomes, pattern = "chr", replacement = ""))

pMarkers = data.frame(chromosome = pMarkerChromosomes, position = pMarkerPositions)

for(i in 1:nrow(BloomQTL)){
  chromosome = BloomQTL$Chromosome[i]
  peakPosition = BloomQTL$`Peak Position (bp)`[i]
  sameChromosomeMarkers = which(pMarkerChromosomes == chromosome)
  sameChromosomeMarkerPositions = pMarkerPositions[sameChromosomeMarkers]
  closestMarker = sameChromosomeMarkerPositions[which(abs(sameChromosomeMarkerPositions - peakPosition) == min(abs(sameChromosomeMarkerPositions-peakPosition)))]
  BloomQTL$closestPeakMarkerPosition[i] = closestMarker
  BloomQTL$growthQTL[i] = paste0("chr", as.roman(chromosome), ":", closestMarker)
}

#change BloomQTL column names to not have the spaces because it is inconvenient to refer to them otherwise.

colnames(BloomQTL) = c("Trait", "FractionPhenotypicVariance", "Chromosome", "PeakPosition",
                       "LOD", "LOD_leftCI", "LOD_rightCI", "genesUnderCI_commonName",
                       "genesUnderCI_standardName", "closestPeakMarkerPosition", "growthQTL")


#get the coordinate positions of the left and right of all eQTL CIs -------
eQTL_Albert2018$CILeft = extract_numeric(eQTL_Albert2018$CI.l)
eQTL_Albert2018$CIRight = extract_numeric(eQTL_Albert2018$CI.r)
eQTL_Albert2018$peak = extract_numeric(eQTL_Albert2018$pmarker)
eQTL_Albert2018$chromosome = roman2int(str_remove(eQTL_Albert2018$chr, pattern = "chr"))


#compute the plotting positions for the eQTLs so that you can plot all chromosomes on a common axis ----
chromosomeSizes = data.frame(chromosome = c(1:16, 0),
                             size = c(230218,
                                      813184,
                                      316620,
                                      1531933,
                                      576874,
                                      270161,
                                      1090940,
                                      562643,
                                      439888,
                                      745751,
                                      666816,
                                      1078177,
                                      924431,
                                      784333,
                                      1091291,
                                      948066,
                                      85779))

spacingBetweenChromosomes = 100000
chromosomeSizes$chromosomeStart_plotting = 0
chromosomeSizes$chromosomeEnd_plotting = 230218 #this is size of chromosome 1

for(i in 2:nrow(chromosomeSizes)){
  previousChromosomeSize = chromosomeSizes$size[i-1]
  currentChromosomeSize = chromosomeSizes$size[i]
  currentChromosome = chromosomeSizes$chromosome[i]
  previousAdditiveFactor = chromosomeSizes$chromosomeStart_plotting[i-1]
  
  chromosomeSizes$chromosomeStart_plotting[i] = previousAdditiveFactor + 
    previousChromosomeSize + spacingBetweenChromosomes
  
  chromosomeSizes$chromosomeEnd_plotting[i] = chromosomeSizes$chromosomeStart_plotting[i] + currentChromosomeSize
}

#save the chromosome sizes as a separate file

save(chromosomeSizes, file = paste0(results_dir, RObj_dir, "chromosomeSizes.rda"))

for(i in 1:nrow(eQTL_Albert2018)){
  chromosome = eQTL_Albert2018$chromosome[i]
  additiveFactorForThisChromosome = chromosomeSizes$chromosomeStart_plotting[which(chromosomeSizes$chromosome == chromosome)]
  eQTL_Albert2018$position_plotting[i] = eQTL_Albert2018$peak[i] + additiveFactorForThisChromosome
  eQTL_Albert2018$start_plotting[i] = eQTL_Albert2018$CILeft[i] + additiveFactorForThisChromosome
  eQTL_Albert2018$end_plotting[i] = eQTL_Albert2018$CIRight[i] + additiveFactorForThisChromosome
}

for(i in 1:nrow(BloomQTL)){
  chromosome = BloomQTL$Chromosome[i]
  additiveFactorForThisChromosome = chromosomeSizes$chromosomeStart_plotting[which(chromosomeSizes$chromosome == chromosome)]
  BloomQTL$position_plotting[i] = BloomQTL$PeakPosition[i] + additiveFactorForThisChromosome
  BloomQTL$start_plotting[i] = BloomQTL$LOD_leftCI[i] + additiveFactorForThisChromosome
  BloomQTL$end_plotting[i] = BloomQTL$LOD_rightCI[i] + additiveFactorForThisChromosome
}


#save----

save(eQTL_Albert2018, BloomQTL, 
     genotypesCommonSegregants, 
     expressionCommonSegregants_batchODCorrected_std, traitCommonSegregants_std,
     file = paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))
