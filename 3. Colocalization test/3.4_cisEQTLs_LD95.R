#!/usr/bin/env Rscript --vanilla
#commands required for MSI 
setwd("~/ExpressionPhenotypeProject/")

#load packages
library(qtl)
library(qtl2)
library(qtl2pleio)
library(ggplot2)
library(rlist)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(GenomicRanges)
library(readr)

#global variables ----
home_dir = '~/ExpressionPhenotypeProject/' 
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
otherFiles_dir = "otherFiles_101522/"
plotting_dir = "plots_092522/"
unfinished_dir = "preLD95_ctrl/"
destination_dir = "LD95_ctrl/"
paddingInKb = 5

#load data ----

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

load(paste0(results_dir, RObj_dir, "genotypesCommonSegregants_LD95_forColoc.rda"))
  
  #load filenames
load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "fileList_QTLLODsOver10_final.rda"))

#variables from command line -------
filenameIndex = commandArgs(trailingOnly = TRUE)

#get the condition, gene, filename for this index ----
filename_this = filename_list[as.numeric(filenameIndex)]
filename_this = stringr::str_remove(filename_this, pattern = ".rda")
filename_brokenDown = stringr::str_split(filename_this, pattern = "__", n = 3, simplify = TRUE)
colnames(filename_brokenDown) = c("condition", "gene", "eQTL")

gene_this = filename_brokenDown[1, "gene"]
condition_this = filename_brokenDown[1, "condition"]

rm(filename_list, filename_brokenDown) #free up memory

##get markers list ----

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

#read tables in cross format ----
DOEx_myData = read.cross(format = "csvs",
                         genfile = paste0(results_dir, otherFiles_dir, "genotypes_LD95.csvs"),
                         phefile = paste0(results_dir, otherFiles_dir, "phenotype.csvs"),
                         genotypes = c("1", "-1"),
                         crosstype = "haploid")

DOEx_myData_2 = convert2cross2(DOEx_myData)

rm(DOEx_myData) #free up memory

probs = calc_genoprob(DOEx_myData_2)
pr = genoprob_to_alleleprob(probs)
names(pr) = sort(unique(markers$chr))

kinship = calc_kinship(probs = pr, type = "loco") #calculate kinship matrix


#compute the pValue for the gene -----

#take one gene, its eQTLs overlapping with the growthQTLs, and do a colocalization scan with the smallest interval that contains both the growth and eQTLs (including the CIs)

findingTheSmallestRegionContainingEQTLAndGrowthQTL = function(condition_i, eQTLList){
  #eQTLsOverlappingGrowthQTLs_i = eQTLsOverlappingGrowthQTLs[[condition_i]]
  #eQTLRow = which(eQTLList$gene == gene_i & eQTLList$pmarker == eQTL_i)
  eQTL_interval = c(eQTLList$eQTL_CILeft[1],
                    eQTLList$eQTL_CIRight[1])
  gQTL_interval = c(eQTLList$gQTL_CILeft[1],
                    eQTLList$gQTL_CIRight[1])
  
  #growthQTLsOverlappingEQTLs = unlist(str_split(eQTLList$overlappingGrowthQTLs[1],
  #                                              pattern = " ; "))
  
  #growthQTLs_df = filter(BloomQTL, Trait == condition_i & growthQTL %in% growthQTLsOverlappingEQTLs)
  #growthQTL_interval = c(growthQTLs_df$LOD_leftCI, growthQTLs_df$LOD_rightCI)
  shortestInterval = c("chromosome" = eQTLList$chromosome[1],
                       "pmarker" = eQTLList$eQTL[1],
                       "left" = min(c(eQTL_interval, gQTL_interval)),
                       "right" = max(c(eQTL_interval, gQTL_interval)))
  
  return(shortestInterval)
}

findNearestMarkerPosition = function(position, chromosome){
  markersList = str_split(colnames(genotypesCommonSegregants_LD95), pattern = ":", simplify = TRUE)
  markersList = as.data.frame(markersList)
  colnames(markersList) = c("chr", "position")
  markersList$position = extract_numeric(markersList$position)
  #markersList$markerName = paste0(markersList$chr, ":", markersList$position)
  
  markers_thisChr = markersList$position[which(markersList$chr == paste0("chr", as.roman(chromosome)))]
  diff = markers_thisChr - position
  minAbsDiff = min(abs(diff))
  
  diff_i =which(abs(diff)== minAbsDiff)
  
  return(markers_thisChr[min(diff_i)]) #we are returning the min in case where two positions have same diff on account of being on equidistant on either side of this position
}

calculatePValueForGene = function(inputVector, gene_i, condition_i, paddingBases){
  print(gene_i)
  gene_i = str_replace(gene_i, pattern = "-", ".") #this is done because the gene name sometimes contains a hyphen in case of paralogs
  condition_i = str_replace(condition_i, pattern = ":", ".")
  condition_i = str_replace(condition_i, pattern = "-", ".")
  
  chr_i = as.numeric(inputVector[grep(names(inputVector), pattern = "chromosome")])
  nearestMarker_leftCI = findNearestMarkerPosition(position = (as.numeric(inputVector["left"]) - paddingBases), chromosome = chr_i)
  nearestMarker_rightCI = findNearestMarkerPosition(position = (as.numeric(inputVector["right"]) + paddingBases), chromosome = chr_i)
  pmarker = inputVector[grep(names(inputVector), pattern = "pmarker")]
  pmarker = str_replace_all(pmarker, pattern = "[[:punct:]]", replacement = "_")
  
  #continue the colocalization test
  pp = pr[[paste0("chr",as.roman(chr_i))]]
  traitMatrix = matrix(0, nrow= nrow(pp), ncol = 2)
  rownames(traitMatrix) = rownames(pp)
  colnames(traitMatrix) = c(condition_i, gene_i)
  traitMatrix[,condition_i] = DOEx_myData_2$pheno[,condition_i]
  traitMatrix[,gene_i] = DOEx_myData_2$pheno[,gene_i]
  
  x = dimnames(pp)[[3]] #markers for pp
  
  startSNP_i = grep(x, pattern = nearestMarker_leftCI)
  endSNP_i = grep(x, pattern = nearestMarker_rightCI)
  nSNP_i = endSNP_i - startSNP_i + 1
  
  print(nSNP_i)
  
  out = scan_pvl(probs = pp,
                 pheno = traitMatrix,
                 kinship = kinship[[as.roman(chr_i)]],
                 start_snp = startSNP_i,
                 n_snp = nSNP_i,
                 cores = parallelly::availableCores())
  
  out_lods = calc_profile_lods(out)
  
  out_lods$marker_position = extract_numeric(out_lods$marker)
  out_lods$cM = DOEx_myData_2$gmap[[paste0("chr", as.roman(chr_i))]][out_lods$marker] #add centimorgan measures
  out_lods$trait[which(out_lods$trait == "tr1")] = condition_i
  out_lods$trait[which(out_lods$trait == "tr2")] = gene_i
  
  lrt = calc_lrt_tib(out)
  
  y = grep(x, pattern = min(out_lods$marker_position))
  pleio_index = find_pleio_peak_tib(out,
                                    start_snp = y)
  
  
  set.seed(456) # set for reproducibility purposes.
  b_out = suppressMessages(boot_pvl(probs = pp,
                                    pheno = traitMatrix,
                                    pleio_peak_index = pleio_index,
                                    kinship = kinship[[as.roman(chr_i)]],
                                    nboot = 1000,
                                    start_snp = startSNP_i,
                                    n_snp = nSNP_i,
                                    cores = parallelly::availableCores()
  ))
  
  pvalue = mean(b_out >= lrt)
  
  #plot the profile
  ##title
  title = paste0("chr", as.roman(chr_i), ":",
                 nearestMarker_leftCI, "-", nearestMarker_rightCI,
                 "  (nMarkers = ", nSNP_i, ")")
  subtitle = paste0("coloc_test p = ", round(pvalue, digits = 3))
  
  plot_trace =  ggplot(out_lods) + 
    geom_line(aes(x = marker_position, y = profile_lod, colour = trait)) + 
    geom_vline(xintercept = extract_numeric(inputVector["pmarker"]), color = "black") +
    ggtitle(title, subtitle = subtitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces"))), 
         dir.create(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces"))), FALSE)
  
  ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces/cisEQTLs"))), 
         dir.create(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces/cisEQTLs"))), FALSE)
  
  ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces/cisEQTLs/", destination_dir))), 
         dir.create(file.path(paste0(results_dir, plotting_dir, "colocalizationLODTraces/cisEQTLs/", destination_dir))),
         FALSE)
  
  pdf(paste0(results_dir, plotting_dir, 
             "colocalizationLODTraces/cisEQTLs/", destination_dir,
             filename_this, ".pdf"))
  print(plot_trace)
  dev.off()
  
  return(pvalue)
}

#load the results table
load(file = paste0(results_dir, RObj_dir,
                   "colocalizationTest/cisEQTLs/",
                   unfinished_dir, "/",
                   filename_this, ".rda"))

inputVector_i = findingTheSmallestRegionContainingEQTLAndGrowthQTL(condition_i = condition_this,
                                                                   eQTLList = eQTLTable_QTLThis)

coloc_pValue = calculatePValueForGene(inputVector = inputVector_i,
                                      gene_i = gene_this,
                                      condition_i = condition_this,
                                      paddingBases = paddingInKb*1000)

eQTLTable_QTLThis$coloc_pValue = coloc_pValue
eQTLTable_QTLThis$coloc_significance = ifelse(eQTLTable_QTLThis$coloc_pValue < 0.05, TRUE, FALSE)


ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))), FALSE)

ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", destination_dir))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", destination_dir))), FALSE)

save(eQTLTable_QTLThis, file = paste0(results_dir, RObj_dir,
                                      "colocalizationTest/cisEQTLs/",
                                      destination_dir, "/",
                                      filename_this, ".rda"))

