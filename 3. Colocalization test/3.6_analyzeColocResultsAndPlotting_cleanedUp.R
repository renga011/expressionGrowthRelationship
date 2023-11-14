#code for analysing the cis-eQTL colocalization results

library(dplyr)
library(ggpubr)
library(stringr)
library(rlang)
library(tidyr)
library(readr)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "colocDone_QTLLODGreaterThan10_final.rda")) ## eQTL LOD >= 10, gQTL LOD >= 10
  QTLs_colocDone$filename_ID = row.names(QTLs_colocDone)
  QTLs_colocDone$gQTL_ID = paste0(QTLs_colocDone$condition, "__", QTLs_colocDone$growthQTL)
  QTLs_colocDone$eQTL_ID = paste0(QTLs_colocDone$geneNames, "__", QTLs_colocDone$eQTL)
  
load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/",
            "fileList_QTLLODsOver10_final.rda"))
  
load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))
load(paste0(results_dir, RObj_dir, "colorPaletteForConditions.rda"))
load(paste0(results_dir, RObj_dir, "chromosomeSizes.rda"))
#remove mitochondrial chromosome row
chromosomeSizes= filter(chromosomeSizes, chromosome!=0)
#create directories ---
ifelse(!dir.exists(file.path(paste0(results_dir, plotting_dir, "cisEQTLs_coloc/"))), 
       dir.create(file.path(paste0(results_dir, plotting_dir, "cisEQTLs_coloc/"))), FALSE)

#get no of LD95 markers in overlapping intervals of gQTL and eQTL ----
  load(paste0(results_dir, RObj_dir, "genotypesCommonSegregants_LD95_forColoc.rda"))
  #nMarkers in region
  findingTheSmallestRegionContainingEQTLAndGrowthQTL = function(eQTLList){
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
  
  paddingBases = 5000
  
  for(i in 1:nrow(QTLs_colocDone)){
    condition_i = QTLs_colocDone$condition[i]
    eQTLList_i = QTLs_colocDone[i,]
    gene_i = QTLs_colocDone$geneNames[i]
    inputVector_i = findingTheSmallestRegionContainingEQTLAndGrowthQTL(eQTLList_i)
    
    #
    gene_i = stringr::str_replace(gene_i, pattern = "-", ".") #this is done because the gene name sometimes contains a hyphen in case of paralogs
    condition_i = stringr::str_replace(condition_i, pattern = ":", ".")
    condition_i = stringr::str_replace(condition_i, pattern = "-", ".")
    
    chr_i = as.numeric(inputVector_i[grep(names(inputVector_i), pattern = "chromosome")])
    nearestMarker_leftCI = findNearestMarkerPosition(position = (as.numeric(inputVector_i["left"]) - paddingBases), chromosome = chr_i)
    nearestMarker_rightCI = findNearestMarkerPosition(position = (as.numeric(inputVector_i["right"]) + paddingBases), chromosome = chr_i)
    
    pmarker = inputVector_i[grep(names(inputVector_i), pattern = "pmarker")]
    pmarker = stringr::str_replace_all(pmarker, pattern = "[[:punct:]]", replacement = "_")
    
    startSNP_i = grep(colnames(genotypesCommonSegregants_LD95), pattern = paste0("chr", as.roman(chr_i), ":",nearestMarker_leftCI))
    endSNP_i = grep(colnames(genotypesCommonSegregants_LD95), pattern = paste0("chr", as.roman(chr_i), ":",nearestMarker_rightCI))
    nSNP_i = endSNP_i - startSNP_i + 1
    
    QTLs_colocDone$nMarkers[i] = nSNP_i
  }

#set the coloc significance of those with <= 3 markers as TRUE even if the colocalization test returned coloc significance as FALSE
  QTLs_colocDone$pleiotropy = FALSE
  QTLs_colocDone$pleiotropy[which(QTLs_colocDone$coloc_significance == "FALSE" & QTLs_colocDone$nMarkers > 3)] = TRUE
  
    ##save for supplementary table
save(QTLs_colocDone, 
     file = paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "colocalizationResults_forSupplement.rda"))
  
#get no. of pleiotropic loci ----
  n = which(QTLs_colocDone$filename_ID %in% filename_list)
  QTLs_colocDone = QTLs_colocDone[n,]
  
  nPleiotropic = length(which(QTLs_colocDone$coloc_significance == FALSE)) #before three marker filter
  nPleiotropic_moreThan3Markers = length(which(QTLs_colocDone$coloc_significance == FALSE & QTLs_colocDone$nMarkers > 3))
  
  nGQTLs = length(unique(QTLs_colocDone$gQTL_ID))
  nGenes = length(unique(QTLs_colocDone$geneNames))
  nConditions = length(unique(QTLs_colocDone$condition))
  nLocalEQTLs = length(unique(QTLs_colocDone$eQTL_ID))

#get the number of colocalizing eQTLs per gQTL ---
  nColocalizingLocalEQTLsPerGQTL = sapply(unique(QTLs_colocDone$gQTL_ID), FUN = function(gQTL_i){
    df = dplyr::filter(QTLs_colocDone, gQTL_ID == gQTL_i)
    nOverlapping = length(unique(df$eQTL_ID))
    nPleiotropic_overThreeMarkers = length(unique(df$eQTL_ID[which(df$coloc_significance == FALSE & df$nMarkers > 3)]))
    
    return(c("nOverlapping" = nOverlapping, "nPleiotropic_moreThan3Markers" = nPleiotropic_overThreeMarkers))
  })
  
  nColocalizingLocalEQTLsPerGQTL = as.data.frame(t(nColocalizingLocalEQTLsPerGQTL))
  nColocalizingLocalEQTLsPerGQTL$gQTL_ID = row.names(nColocalizingLocalEQTLsPerGQTL)
  
##plot number of genes and number of eQTLs with 0, 1, >1 colocalizing eQTLs
nColocalizingLocalEQTLsPerGQTL$QTL_label = "0"
nColocalizingLocalEQTLsPerGQTL$QTL_label[which(nColocalizingLocalEQTLsPerGQTL$nPleiotropic_moreThan3Markers == 1)] = "1"
nColocalizingLocalEQTLsPerGQTL$QTL_label[which(nColocalizingLocalEQTLsPerGQTL$nPleiotropic_moreThan3Markers > 1)] = "> 1"

nQTLs = nColocalizingLocalEQTLsPerGQTL %>% group_by(QTL_label) %>% summarise(nGQTLs = n())

#plot stacked bar plot of nQTLs
melt_nQTLs = reshape2::melt(nQTLs, id.vars= "QTL_label")
melt_nQTLs$QTL_label = factor(melt_nQTLs$QTL_label, levels = c("> 1", "1", "0"))

  pdf(paste0(results_dir, plotting_dir, "cisEQTLs_coloc/", "nQTLs_nGenes.pdf"))
ggbarplot(data = melt_nQTLs,
          x = "variable",
          y = "value",
          fill = "QTL_label",
          width = 0.25) + scale_fill_manual(values = c("0" = "dimgrey", "1" = "darkblue", "> 1" = "gold")) 
  dev.off()

## save number of colocalized local eQTLs for different gQTLs 
  nColocalizingLocalEQTLsPerGQTL_df = cbind(nColocalizingLocalEQTLsPerGQTL, as.data.frame(str_split(nColocalizingLocalEQTLsPerGQTL$gQTL_ID, pattern = "__", simplify = TRUE)))
  nColocalizingLocalEQTLsPerGQTL_df = select(nColocalizingLocalEQTLsPerGQTL_df, -gQTL_ID, -QTL_label)
  colnames(nColocalizingLocalEQTLsPerGQTL_df) = c("nOverlappingLocalEQTLs", "nPleiotropicLocalEQTLs", "condition", "gQTL_pMarker")
  
  write.csv(nColocalizingLocalEQTLsPerGQTL_df, file = paste0(results_dir, otherFiles_dir, "nColocalizingLocalEQTLsPerGQTL_tableS3.csv"), quote = FALSE, row.names = FALSE)
  

#histogram of number of colocalizing gQTLs with N localizing eQTLs
  pdf(paste0(results_dir, plotting_dir, "cisEQTLs_coloc/", "histogramOfNColocalizingEQTLsAtGQTLs.pdf"))
  gghistogram(nColocalizingLocalEQTLsPerGQTL,
              x = "nPleiotropic_moreThan3Markers",
              bins = 22, fill = "grey") + ylab("Number of gQTLs") + xlab("Number of colocalizing local-eQTLs")
  dev.off()
    
#Plot genome browser view of colocalization examples -------

plotLayout = function(gQTL_i, condition_i, filename, causalGene){

#gQTL intervals 
gQTL_CILeft = BloomQTL$LOD_leftCI[which(BloomQTL$Trait == condition_i & BloomQTL$growthQTL == gQTL_i)]
gQTL_CIRight = BloomQTL$LOD_rightCI[which(BloomQTL$Trait == condition_i & BloomQTL$growthQTL == gQTL_i)]
gQTL_chr = BloomQTL$Chromosome[which(BloomQTL$Trait == condition_i & BloomQTL$growthQTL == gQTL_i)]

#get the genes around the interval
genes_withinGQTL = geneList_SGD[which(geneList_SGD$GenomePosition_start >= gQTL_CILeft - 20000 & geneList_SGD$GenomePosition_end <= gQTL_CIRight + 20000 & geneList_SGD$chromosome == gQTL_chr),]

genes_withinGQTL$hasLocalOverLOD10 = ifelse(genes_withinGQTL$sysName %in% QTLs_colocDone$geneNames[which(QTLs_colocDone$condition == condition_i)], TRUE, FALSE)

genes_withinGQTL$hasColocalizedEQTL = ifelse(genes_withinGQTL$sysName %in% QTLs_colocDone$geneNames[which(QTLs_colocDone$condition == condition_i & QTLs_colocDone$coloc_significance == "FALSE")], TRUE, FALSE)

genes_withinGQTL$causalGene = ifelse(grepl(genes_withinGQTL$stdName, pattern = causalGene), TRUE, FALSE)

#random y-positions
set.seed(13)
genes_withinGQTL$y = sample(1:nrow(genes_withinGQTL), size = nrow(genes_withinGQTL))

#plot layout

pdf(paste0(results_dir, plotting_dir, "cisEQTLs_coloc/", "exampleOfColocalizedCis_linkageArtifacts_", filename, ".pdf"), width = 12, height = 10)

plot = ggplot(genes_withinGQTL) +
  geom_rect(aes(xmin = GenomePosition_start, xmax = GenomePosition_end,
                ymin = y*0.5-0.2, ymax = y*0.5+0.2, size = hasLocalOverLOD10, 
                color = hasColocalizedEQTL, fill = causalGene)) +
  geom_vline(xintercept =gQTL_CILeft, color = "black", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = gQTL_CIRight, color = "black", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "grey") +
  theme_classic() +
  geom_text(aes(x = 0.5*(GenomePosition_start + GenomePosition_end),label = stdName, y = y*0.5+ 0.4), angle = 0, size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "colocalized local- eQTL?") +
  scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.5), name = "local-eQTL?") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "light grey")) +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = sort(genes_withinGQTL$y*0.5)) +
  theme_textProperties +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.ticks.y = element_blank())

print(plot)
dev.off()
}

##Example 1 - ENA1 - chromosome 4 gQTL in LiCl
plotLayout(gQTL_i = "chrIV:523964", condition_i = "Lithium_Chloride",
           filename = "ENA1_LithiumChloride", causalGene = "ENA1")

##Example 2 - CUP1 - chromosome 8 gQTL in Copper
plotLayout(gQTL_i = "chrVIII:208566", condition_i = "Copper",
           filename = "CUP1_Copper", causalGene = "CUP1")

##Example 3 - PMR1 - chromosome 7 gQTL in MnSO4
plotLayout(gQTL_i = "chrVII:187372", condition_i = "Manganese_Sulfate",
           filename = "PMR1_Manganese_Sulfate", causalGene = "PMR1")

#make supplementary table 4 - number of overlapping and pleiotropic eQTLs for each gQTL ----

gQTL_ID = stringr::str_split(nColocalizingLocalEQTLsPerGQTL$gQTL_ID, pattern = "__", simplify = TRUE)
colnames(gQTL_ID) = c("trait", "gQTL_peak_position")

s4 = as.data.frame(cbind(gQTL_ID, dplyr::select(nColocalizingLocalEQTLsPerGQTL, nOverlapping, nPleiotropic_moreThan3Markers)))

colnames(s4) = c(colnames(gQTL_ID), "n_overlapping_eQTLs", "n_pleiotropic_eQTLs")

  ##save
openxlsx::write.xlsx(s4, file = paste0(results_dir, otherFiles_dir, "supplementaryTables/", "Supplementary Table 3 - number of overlapping and pleiotropic eQTLs for gQTLs used in the colocalization test.xlsx"))


