#plot specific colocalization examples

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
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "colocDone_QTLLODGreaterThan10.rda")) ## eQTL LOD >= 10, gQTL LOD >= 10
QTLs_colocDone$filename_ID = row.names(QTLs_colocDone)
QTLs_colocDone$gQTL_ID = paste0(QTLs_colocDone$condition, "__", QTLs_colocDone$growthQTL)

load(paste0(results_dir, RObj_dir, "genotypesCommonSegregants_LD95_forColoc.rda"))
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))

#get objects for running the coloc functions --- very similar to the cisEQTLs_LD95.r code

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

#do the colocalization test - plot the different traces for the given example

findingTheSmallestRegionContainingEQTLAndGrowthQTL = function(condition_i, eQTLList){
  #eQTLsOverlappingGrowthQTLs_i = eQTLsOverlappingGrowthQTLs[[condition_i]]
  #eQTLRow = which(eQTLList$gene == gene_i & eQTLList$pmarker == eQTL_i)
  eQTL_interval = c(min(eQTLList$eQTL_CILeft),
                    max(eQTLList$eQTL_CIRight))
  gQTL_interval = c(min(eQTLList$gQTL_CILeft),
                   max(eQTLList$gQTL_CIRight))
  
  #growthQTLsOverlappingEQTLs = unlist(str_split(eQTLList$overlappingGrowthQTLs[1],
  #                                              pattern = " ; "))
  
  #growthQTLs_df = filter(BloomQTL, Trait == condition_i & growthQTL %in% growthQTLsOverlappingEQTLs)
  #growthQTL_interval = c(growthQTLs_df$LOD_leftCI, growthQTLs_df$LOD_rightCI)
  shortestInterval = c("chromosome" = eQTLList$chromosome[1],
                       #"pmarker" = eQTLList$eQTL[1],
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

#QTLs list for a given gQTL
gQTL = "Tunicamycin__chrX:250151"
  
getPlottingCoordinates = function(i, QTLs_list){
  condition_i = QTLs_list$condition[i]
  print(condition_i)
  inputVector = findingTheSmallestRegionContainingEQTLAndGrowthQTL(condition_i, QTLs_list)
  
  gene_i = str_replace(QTLs_list$geneNames[i], pattern = "-", ".") #this is done because the gene name sometimes contains a hyphen in case of paralogs
  condition_i = str_replace(QTLs_list$condition[i], pattern = ":", ".")
  condition_i = str_replace(condition_i, pattern = "-", ".")
  print(condition_i)
  
  chr_i = as.numeric(inputVector[grep(names(inputVector), pattern = "chromosome")])
  nearestMarker_leftCI = findNearestMarkerPosition(position = (as.numeric(inputVector["left"]) - 5000), chromosome = chr_i)
  nearestMarker_rightCI = findNearestMarkerPosition(position = (as.numeric(inputVector["right"]) + 5000), chromosome = chr_i)
  #pmarker = inputVector[grep(names(inputVector), pattern = "pmarker")]
  #pmarker = str_replace_all(pmarker, pattern = "[[:punct:]]", replacement = "_")
  
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
  
  return(out_lods)
  
}
plot_LODTrace = function(gQTL){
QTLs_list = filter(QTLs_colocDone, gQTL_ID == gQTL & nMarkers ==61)

plottingCoords = lapply(1:nrow(QTLs_list), FUN = getPlottingCoordinates, QTLs_list)

#get the coordinates for plotting the overlap of QTL cartoon 
gQTL_start = min(QTLs_list$gQTL_CILeft)
gQTL_end = max(QTLs_list$gQTL_CIRight)
QTLs_list$gene_left = geneList_SGD$GenomePosition_start[match(QTLs_list$geneNames, geneList_SGD$sysName)]
QTLs_list$gene_right = geneList_SGD$GenomePosition_end[match(QTLs_list$geneNames, geneList_SGD$sysName)]

pdf(paste0(results_dir, plotting_dir, "colocalizationLODTraceExample__", gQTL, ".pdf"))

lapply(1:length(plottingCoords), function(i){
  plottingCoords_i = plottingCoords[[i]]
  plottingCoords_i = plottingCoords_i[!duplicated(plottingCoords_i),]
  plottingCoords_i = dplyr::filter(plottingCoords_i, trait != "pleiotropy")
  
  QTL_list_i = QTLs_list[i,]
  
  plot_trace =  ggplot() + 
    geom_line(plottingCoords_i, 
              mapping = aes(x = marker_position, y = profile_lod, colour = trait)) + 
    geom_rect(aes(xmin = gQTL_start - 5000, xmax = gQTL_end + 5000, ymin = -0.05, ymax = 0.05),
              color = "darkblue", fill = "darkblue", alpha =0.4) +
    geom_rect(QTL_list_i, 
              mapping = aes(xmin = gene_left, xmax = gene_right, ymin = 0.2 - 0.05, ymax = 0.2 + 0.05),color = "black", fill = "black", alpha = 0.1) +
    geom_rect(QTL_list_i, 
              mapping = aes(xmin = eQTL_CILeft - 5000, xmax = eQTL_CIRight + 5000, ymin = 0.2 - 0.05, ymax = 0.2 + 0.05), color = "black", alpha = 0) +
    ggtitle(gQTL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  #draw plot
  print(plot_trace)
  
})

dev.off()
}

  ##example 1 - "Tunicamycin__chrX:250151"
plot_LODTrace("Tunicamycin__chrX:250151")

## get QTL overlap layout for candidate examples

QTLs_list = filter(QTLs_colocDone, gQTL_ID == gQTL)


pdf(paste0(results_dir, plotting_dir, "colocalizationQTLOverlapExample__", gQTL, ".pdf"))
lapply(1:nrow(QTLs_list), 
       FUN = function(i){
         inputVector = findingTheSmallestRegionContainingEQTLAndGrowthQTL(QTLs_list$condition[i], QTLs_list)
         QTL_list_i = QTLs_list[i,]
          ggplot() +
  geom_rect(aes(xmin = gQTL_start - 5000, xmax = gQTL_end + 5000, ymin = -0.05, ymax = 0.05),
            color = "darkblue", fill = "darkblue", alpha =0.4) +
  geom_rect(QTL_list_i, 
            mapping = aes(xmin = gene_left, xmax = gene_right, ymin = 0.2 - 0.05, ymax = 0.2 + 0.05, color = geneNames, fill = geneNames), alpha = 0.1) +
  geom_rect(QTL_list_i, 
            mapping = aes(xmin = eQTL_CILeft - 5000, xmax = eQTL_CIRight + 5000, ymin = 0.2 - 0.05, ymax = 0.2 + 0.05, color = geneNames), alpha = 0) +
  ylim(-0.2,0.5) + 
  xlim(inputVector["left"] - 5000, inputVector["right"] + 5000) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "grey") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
})
dev.off()
