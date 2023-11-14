library(rlist)
library(stringr)
library(dplyr)
library(tidyr)

#global variables ---
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

#reduce the number of markers so that the colocalization becomes computationally less intensive.
## We need more markers for high-res QTL mapping, not for doing the colocalization. So decrease stingency of LD Threshold and recompute which markers to use ----

#function
computeRForMarkers = function(marker_i, starterMarker){
  genotype_i = unlist(genotypesCommonSegregants[grep(colnames(genotypesCommonSegregants), pattern = paste0(marker_i, "_"))])
  genotype_start = unlist(genotypesCommonSegregants[grep(colnames(genotypesCommonSegregants), pattern = paste0(starterMarker, "_"))])
  return(cor(genotype_i, genotype_start))
}

markersList = str_split(colnames(genotypesCommonSegregants), pattern = ":", simplify = TRUE)
markersList = as.data.frame(markersList)
colnames(markersList) = c("chr", "position")
markersList$position = extract_numeric(markersList$position)
markersList$markerName = paste0(markersList$chr, ":", markersList$position)

LDThreshold = 0.95 #change this per requirement

markers_new = c()
for(chr_i in unique(markersList$chr)){
  print(chr_i)
  markersList_thisChr = dplyr::filter(markersList, chr == chr_i)
  i = 1
  while(i < nrow(markersList_thisChr)){
    print(markersList_thisChr$markerName[i])
    startChr = markersList_thisChr$markerName[i]
    Rs_allOtherMarkers = lapply(markersList_thisChr$markerName[(i+1): nrow(markersList_thisChr)],
                                FUN = computeRForMarkers, starterMarker = startChr)
    Rs_allOtherMarkers = unlist(Rs_allOtherMarkers) 
    names(Rs_allOtherMarkers) = markersList_thisChr$markerName[(i+1): nrow(markersList_thisChr)]
    rows_moreThanLDThreshold = which(Rs_allOtherMarkers >= LDThreshold) 
    if(!rlang::is_empty(rows_moreThanLDThreshold)){
      markersList_thisChr = markersList_thisChr[-(rows_moreThanLDThreshold+i),]
    }
    i = i+1
  }
  print(markersList_thisChr$markerName)
  markers_new = append(markers_new, markersList_thisChr$markerName)
  print(length(markers_new))
}

markersList_new = dplyr::filter(markersList, markerName %in% markers_new)

nMarkers = markersList %>% group_by(chr) %>% summarise(LD100 = n())
nMarkers_new = markersList_new %>% group_by(chr) %>% summarise(LD95 = n())
nMarkers$LD95 = nMarkers_new$LD95

#reduce the number of markers in genotypesCommonSegregants----

genotypesCommonSegregants_LD95 = genotypesCommonSegregants
colnames(genotypesCommonSegregants_LD95) = str_split(colnames(genotypesCommonSegregants),
                                                     pattern = "_", simplify = TRUE)[,1]
genotypesCommonSegregants_LD95 = dplyr::select(genotypesCommonSegregants_LD95,
                                               markers_new)

save(genotypesCommonSegregants_LD95, file = paste0(results_dir, RObj_dir,"genotypesCommonSegregants_LD95_forColoc.rda"))
