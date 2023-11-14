library(dplyr)
library(GenomicRanges)
library(magrittr)

#global variables ---
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
paddingInKb = 5

#load data ----

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))

#overlap growthQTLs and eQTLs ---

growth_gRanges = GRanges(seqnames = BloomQTL$Chromosome,
                         IRanges(start = BloomQTL$LOD_leftCI - paddingInKb*1000, 
                                 end = BloomQTL$LOD_rightCI + paddingInKb*1000),
                         pMarker = BloomQTL$closestPeakMarkerPosition,
                         LOD = BloomQTL$LOD,
                         growthQTL = BloomQTL$growthQTL,
                         condition = BloomQTL$Trait,
                         r_growth = BloomQTL$r_growth,
                         se_growth = BloomQTL$se_growth,
                         gQTL_CILeft = BloomQTL$LOD_leftCI,
                         gQTL_CIRight = BloomQTL$LOD_rightCI)

eQTL_gRanges = GRanges(seqnames = eQTL_Albert2018$chromosome,
                       IRanges(start = eQTL_Albert2018$CILeft - paddingInKb*1000, 
                               end = eQTL_Albert2018$CIRight + paddingInKb*1000),
                       geneNames = eQTL_Albert2018$gene,
                       eQTL = eQTL_Albert2018$pmarker,
                       eQTL_peak = eQTL_Albert2018$peak,
                       cis = eQTL_Albert2018$cis,
                       LOD = eQTL_Albert2018$LOD,
                       r_expression = eQTL_Albert2018$r_expression,
                       se_expression = eQTL_Albert2018$se_expression,
                       eQTL_CILeft = eQTL_Albert2018$CILeft,
                       eQTL_CIRight = eQTL_Albert2018$CIRight)

hits = findOverlaps(growth_gRanges, eQTL_gRanges)

eQTLs_overlapping = as.data.frame(eQTL_gRanges[subjectHits(hits)])
growthQTLs_overlapping = as.data.frame(growth_gRanges[queryHits(hits)])

eQTLs_overlapping$growthQTL = growthQTLs_overlapping$growthQTL
eQTLs_overlapping$condition = growthQTLs_overlapping$condition
eQTLs_overlapping$gQTL_LOD = growthQTLs_overlapping$LOD
eQTLs_overlapping$r_growth = growthQTLs_overlapping$r_growth
eQTLs_overlapping$se_growth = growthQTLs_overlapping$se_growth
eQTLs_overlapping$gQTL_peak = growthQTLs_overlapping$pMarker
eQTLs_overlapping$gQTL_CILeft = growthQTLs_overlapping$gQTL_CILeft
eQTLs_overlapping$gQTL_CIRight = growthQTLs_overlapping$gQTL_CIRight
eQTLs_overlapping$interPeakDistance = abs(eQTLs_overlapping$eQTL_peak- eQTLs_overlapping$gQTL_peak)

#get the LD between gQTL and eQTLs: You maybe don't need to do coloc for markers that are already super close-by.
cisEQTLs_overlapping = dplyr::filter(eQTLs_overlapping, cis == "TRUE")

cisEQTLs_overlapping$gQTL_LOD = BloomQTL$LOD[match(paste0(cisEQTLs_overlapping$growthQTL,cisEQTLs_overlapping$condition), paste0(BloomQTL$growthQTL, BloomQTL$Trait))]

cisEQTLs_overlapping = dplyr::select(cisEQTLs_overlapping, -start, -end, -width, -strand)
colnames(cisEQTLs_overlapping)[1] = "chromosome"

for(i in 1:nrow(cisEQTLs_overlapping)){
  eQTL_marker = cisEQTLs_overlapping$eQTL[i]
  gQTL_marker = cisEQTLs_overlapping$growthQTL[i]
  
  gQTL_col = grep(colnames(genotypesCommonSegregants), 
                  pattern = paste0(gQTL_marker, "_"))
  
  LD = cor(unlist(genotypesCommonSegregants[eQTL_marker]),
           unlist(genotypesCommonSegregants[gQTL_col]))
  
  cisEQTLs_overlapping$LD[i] = round(as.numeric(LD),2)
}

#remove cases where the cis-eQTL overlaps with 2 gQTLs - keep only the overlap with strongest gQTL

cisEQTLs_overlapping$ID = paste0(cisEQTLs_overlapping$condition,
                                 "__", cisEQTLs_overlapping$geneNames,
                                 "__", stringr::str_replace_all(cisEQTLs_overlapping$eQTL, "[[:punct:]]", "_"))

cisEQTLs_overlapping = cisEQTLs_overlapping[order(cisEQTLs_overlapping$ID, -cisEQTLs_overlapping$LOD), ] 
cisEQTLs_overlapping = cisEQTLs_overlapping[!duplicated(cisEQTLs_overlapping$ID),]
cisEQTLs_overlapping = dplyr::select(cisEQTLs_overlapping, -ID)

#break into QTLs and save -----
  #create directories for saving ---
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest"))),
       FALSE)
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))),
       FALSE)
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "preLD95"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "preLD95"))),
       FALSE)

for(i in 1:nrow(cisEQTLs_overlapping)){
  gene_i = cisEQTLs_overlapping$geneNames[i]
  condition_i = cisEQTLs_overlapping$condition[i]
  eQTL_i = stringr::str_replace_all(cisEQTLs_overlapping$eQTL[i], "[[:punct:]]", "_")
  
  eQTLTable_QTLThis = cisEQTLs_overlapping[i,]
  
  filename_i = paste0(condition_i, "__", gene_i, "__",
                      eQTL_i, ".rda")
  
  save(eQTLTable_QTLThis, 
       file = paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/preLD95/", filename_i))
  
}

#make the files list to determine the order of running jobs on MSI
  #CHANGE AS PER RUN ORDER
  #df = filter(cisEQTLs_overlapping, LOD >= 10 & gQTL_LOD >=10 & LD > 0.98)

  df = filter(cisEQTLs_overlapping, LOD >= 10 & gQTL_LOD >=10)
  df = arrange(df, desc(LOD))
  
filename_list = paste0(df$condition, "__", 
                  df$geneNames, "__",
                  stringr::str_replace_all(df$eQTL, "[[:punct:]]", "_"),
                  ".rda")

#duplicates_cisEQTLs = paste0(duplicates_cisEQTLs, ".rda")
#filename_list = filename_list[which(filename_list %in% duplicates_cisEQTLs)]

  #filename_list - not done
#filename_list_1 = filename_list[-which(filename_list %in% colocResults_df$filename[!is.na(colocResults_df$coloc_significance)])]

  #append filename_list for eQTL LOD >=10 and gQTL LOD <10, all LDs to filename_list_1

#save the file-list
save(filename_list, file = paste0(results_dir, RObj_dir, 
                             "colocalizationTest/cisEQTLs/",
                             "fileList_QTLLODsOver10_final.rda"))
     
##################


#controls - benchmarking expectations ----
#pure negative control ---
nonOverlapping = as.data.frame(eQTL_gRanges[eQTL_gRanges %outside% growth_gRanges,]) #these get eQTLs that don't overlap with any of the gQTLs across conditions.

nonOverlapping$growthQTL = nonOverlapping$eQTL
nonOverlapping$gQTL_LOD = nonOverlapping$LOD
nonOverlapping$r_growth = nonOverlapping$r_expression
nonOverlapping$se_growth = nonOverlapping$se_expression
nonOverlapping$gQTL_peak = nonOverlapping$eQTL_peak
nonOverlapping$gQTL_CILeft = nonOverlapping$eQTL_CILeft
nonOverlapping$gQTL_CIRight = nonOverlapping$eQTL_CIRight
nonOverlapping$interPeakDistance = abs(nonOverlapping$eQTL_peak- nonOverlapping$gQTL_peak)

nonOverlapping$condition = sample(BloomQTL$Trait, size = nrow(nonOverlapping), replace = TRUE)
#randomly assign conditions to the eQTLs

#shortlist nonOverlapping which has eQTL LOD >=10 and eQTl CI Left != eQTL CI Right
nonOverlapping = dplyr::filter(nonOverlapping, LOD >=10 & eQTL_CILeft != eQTL_CIRight)

#negative control - get targets where the confidence intervals don't overlap (are right next to each other)
eQTL_gr = eQTL_gRanges[which(eQTL_gRanges$cis == "TRUE" & eQTL_gRanges$LOD >= 10)]

nonOverlapping2 = as.data.frame(eQTL_gr)

for(i in 1:length(eQTL_gr)){
  dist_i = GenomicRanges::distance(eQTL_gr[i], growth_gRanges, select = "all")
  minDist = min(dist_i[which(dist_i!=0)], na.rm= TRUE)
  
  gQTL_closest = growth_gRanges[which(dist_i == minDist)]
  
  nonOverlapping2$distanceToGQTL[i] = minDist
  nonOverlapping2$gQTL_LOD[i] = gQTL_closest$LOD
  nonOverlapping2$growthQTL[i] = gQTL_closest$growthQTL
  nonOverlapping2$condition[i] = gQTL_closest$condition
  nonOverlapping2$r_growth[i] = gQTL_closest$r_growth
  nonOverlapping2$se_growth[i] = gQTL_closest$se_growth
  nonOverlapping2$gQTL_peak[i] = gQTL_closest$pMarker
  nonOverlapping2$gQTL_CILeft[i] = gQTL_closest$gQTL_CILeft
  nonOverlapping2$gQTL_CIRight[i] = gQTL_closest$gQTL_CIRight
  nonOverlapping2$interPeakDistance[i] = abs(nonOverlapping2$eQTL_peak[i] - nonOverlapping2$gQTL_peak[i])
  
}

nonOverlapping2 = dplyr::filter(nonOverlapping2, LOD >=10 & gQTL_LOD >= 10 )

#third set of nonOverlappings - this time - pick from QTLs that don't overlap with any gQTLs but in the overlapping use a dummy growthQTL condition where growth = K*expressionOfGene (This is when I expect the eQTL to overlap the gQTL because the growth and expression are pretty much multiples)

##first select 100 eQTls
set.seed(123)
rows_sampled= sample(1:nrow(nonOverlapping), 100, replace = FALSE)
nonOverlapping3 = nonOverlapping[rows_sampled,]
nonOverlapping3$condition = paste0("dummy_", 1:100)

dummyConditions = as.data.frame(matrix(NA, nrow = nrow(traitCommonSegregants_std),
                                       ncol = 100))
colnames(dummyConditions) = nonOverlapping3$condition

set.seed(234)
multiplicativeFactor = rnorm(100)

for(i in 1:nrow(nonOverlapping3)){
  condition_i = nonOverlapping3$condition[i]
  gene_i = nonOverlapping3$geneNames[i]
  growth_i = multiplicativeFactor[i]*unlist(expressionCommonSegregants_batchODCorrected_std[gene_i])
  
  #put some noise on the data: growth = multiplicativeFactor*expression + noise
  set.seed(i)
  noise = rnorm(length(growth_i), mean =0, sd = 1)
  
  #growth_i = growth_i + noise
  dummyConditions[condition_i] = growth_i
}

  #save the dummy conditions values - shall be appended to the traitCommonSegregants_std table
save(dummyConditions, file = paste0(results_dir, RObj_dir,
                                  "dummyConditionsForColocalizationBenchmarking.rda"))


#make final set of QTL-benchmarking set:
  ## 100 from the first set
  ## QTLs with distance between 0-1500 in second for gQTls with LOD >=10 - this is less than 100
  ## 100 from the third set
b1_rows = setdiff(1:nrow(nonOverlapping), rows_sampled)
set.seed(345)
b1 = nonOverlapping[b1_rows,]

b2 = dplyr::filter(nonOverlapping2, distanceToGQTL <= 1500 & gQTL_LOD >=10) %>% dplyr::select(-distanceToGQTL)

benchmarkingSet = rbind(b1, b2, nonOverlapping3)
colnames(benchmarkingSet)[1] = "chromosome"

#break into QTLs and save ---
#create directories for saving ---
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest"))),
       FALSE)
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/"))),
       FALSE)
ifelse(!dir.exists(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "preLD95_ctrl"))), 
       dir.create(file.path(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "preLD95_ctrl"))),
       FALSE)

for(i in 1:nrow(benchmarkingSet)){
  gene_i = benchmarkingSet$geneNames[i]
  condition_i = benchmarkingSet$condition[i]
  eQTL_i = stringr::str_replace_all(benchmarkingSet$eQTL[i], "[[:punct:]]", "_")
  
  eQTLTable_QTLThis = benchmarkingSet[i,]
  
  filename_i = paste0(condition_i, "__", gene_i, "__",
                      eQTL_i, ".rda")
  
  save(eQTLTable_QTLThis, 
       file = paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/preLD95_ctrl/", filename_i))
}

# make and save file list ---
#make the files list to determine the order of running jobs on MSI
#CHANGE AS PER RUN ORDER

filename_list = paste0(benchmarkingSet$condition, "__", 
                       benchmarkingSet$geneNames, "__",
                       stringr::str_replace_all(benchmarkingSet$eQTL, "[[:punct:]]", "_"),
                       ".rda")

#save
save(filename_list, file = paste0(results_dir, RObj_dir, 
                                  "colocalizationTest/cisEQTLs/",
                                  "fileList_benchmarkingControls.rda"))
load(paste0(results_dir, RObj_dir, 
            "colocalizationTest/cisEQTLs/",
            "fileList_benchmarkingControls.rda"))
