#get number of markers in LD with hotspot markers

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(file = paste0(results_dir, RObj_dir, "randomMarkersForRandomHeritabilityEstimate.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

#functions to compute LD -----
computeRForMarkers = function(marker_i, starterMarker){
  genotype_i = unlist(genotypesCommonSegregants[grep(colnames(genotypesCommonSegregants), pattern = paste0(marker_i))])
  genotype_start = unlist(genotypesCommonSegregants[grep(colnames(genotypesCommonSegregants), pattern = paste0(starterMarker))])
  return(cor(genotype_i, genotype_start))
}

# add the actual set to random marker set so that you can process all of them at the same time ----
randomMarkers[[1001]] = hotspotData_Albert2018$hotspotMarker
names(randomMarkers) = c(1:1000, "actual")

# compute the LD with all markers
allMarkers = colnames(genotypesCommonSegregants)
getLDForMarkerSet = lapply(randomMarkers,
                           FUN = function(markerSet){
  LDforMarkers = parallel::mclapply(markerSet,
                        FUN = function(hotspot_i){
                          print(hotspot_i)
                          #get the chromosome and all markers in the same chromosome
                          chr_i = stringr::str_split_i(hotspot_i, pattern = ":", i =1)
                          markersInThisChr = allMarkers[grep(allMarkers, pattern = paste0(chr_i, ":"))]
                          LD_i = lapply(markersInThisChr,
                                        FUN = function(mkr){
                                          print(mkr)
                                          ld = computeRForMarkers(marker_i = mkr, starterMarker = hotspot_i)
                                          return(c("hotspot_mkr" = hotspot_i,
                                                   "marker" = mkr,
                                                   "LD" = ld))
                                        })
                          LD_df = as.data.frame(do.call(rbind, LD_i))
                          LD_df$LD = as.numeric(LD_df$LD)
                          
                          LD_df = LD_df[which(LD_df$LD != 1),]
                          
                          # Categorize values into bins
                          bins = cut(LD_df$LD, breaks = c(-Inf, 0.6, 0.7, 0.8, 0.9, Inf),
                                     labels = c("<0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", ">0.9"))
                          # Count values in each bin
                          counts = table(bins)
                          return(counts)
                          # proportions = counts/length(markersInThisChr)
                          # 
                          # return(proportions)
                        },
                        mc.cores = parallelly::availableCores())
  
  LD_proportions = as.data.frame(do.call(rbind, LDforMarkers))
  LD_proportions$hotspotMarker = markerSet
  return(LD_proportions)
})

names(getLDForMarkerSet) = names(randomMarkers)

  #save
save(getLDForMarkerSet, file = paste0(results_dir, RObj_dir, "reviewerAnalysis_nMarkersAtDifferentLD.rda"))

# analyze data  -----
  load(paste0(results_dir, RObj_dir, "reviewerAnalysis_nMarkersAtDifferentLD.rda"))

actualLDData = getLDForMarkerSet[["actual"]]
actualLDData = dplyr::select(actualLDData, -hotspotMarker, -`<0.6`)
actualLDData = as.data.frame(t(apply(actualLDData, 1, function(x) rev(cumsum(rev(x))))))
colnames(actualLDData) = c(">0.6", ">0.7", ">0.8", ">0.9")
actualLDData = as.data.frame(actualLDData)
actualLDData$set = "actualSet"

#
getLDForMarkerSet = getLDForMarkerSet[1:1000]
avgLDForEachBin = lapply(getLDForMarkerSet,
                         function(x){
                           x = dplyr::select(x, - `<0.6`, -hotspotMarker)
                           x = as.data.frame(t(apply(x, 1, function(i) rev(cumsum(rev(i))))))
                           colnames(x) = c(">0.6", ">0.7", ">0.8", ">0.9")
                           x_medians = apply(x, 2, FUN = median)
                           return(x_medians)
                         })

avgLDForEachBin_randomMarkers = as.data.frame(do.call(rbind,avgLDForEachBin))
avgLDForEachBin_randomMarkers$set = "randomMarkers"

#
allLDData = rbind(actualLDData, avgLDForEachBin_randomMarkers)

allLDData_melt = reshape2::melt(allLDData, id.vars = "set")
colnames(allLDData_melt) = c("set", "LDRange", "nMarkers")

pdf(paste0(results_dir, plotting_dir, "revision_figures_CellGenomics/", "nMarkersInLDWithHotspots_actualVsRandomSets.pdf"))

ggplot(allLDData_melt) +
  geom_boxplot(aes(x = LDRange, y = nMarkers, color = set)) +
  scale_color_manual(values = c("randomMarkers" = "darkgray", "actualSet" = "black")) +
  stat_compare_means(aes(x = LDRange, y = nMarkers, group = set), method = "wilcox.test", label = "p.format") +
  ylab("Number of markers in LD with hotspot marker") +
  xlab("LD range") +
  theme_bw() +
  theme(legend.position = "top")

dev.off()

# number of markers in random sets for table in review (R2,C8) ------
getLDForMarkerSet = getLDForMarkerSet[1:1000]
LDMarkers = lapply(getLDForMarkerSet,
                   function(df){
                     x = dplyr::select(df, - `<0.6`, -hotspotMarker)
                     x = as.data.frame(t(apply(x, 1, function(i) rev(cumsum(rev(i))))))
                     colnames(x) = c(">0.6", ">0.7", ">0.8", ">0.9")
                     sumMkr = apply(x, 2, sum)
                     return(sumMkr)
                   })

LDMarkers_all = as.data.frame(do.call(rbind, LDMarkers))
nMarkers = apply(LDMarkers_all, 2, median) 
fractionMarkers = nMarkers/ncol(genotypesCommonSegregants)

getLDForMarkerSet_df = as.data.frame(do.call(rbind, getLDForMarkerSet))
allMarkers = colnames(genotypesCommonSegregants)
getLDForMarkerSet_df$nMarkersInChromosome = sapply(getLDForMarkerSet_df$hotspotMarker,
       FUN = function(mkr){
         chr_i = stringr::str_split_i(mkr, pattern = ":", i =1)
         markersInThisChr = allMarkers[grep(allMarkers, pattern = paste0(chr_i, ":"))]
         return(length(markersInThisChr) )
       })

getLDForMarkerSet_df[,1:5] = getLDForMarkerSet_df[,1:5]/getLDForMarkerSet_df$nMarkersInChromosome
x = dplyr::select(getLDForMarkerSet_df, - `<0.6`, -hotspotMarker, -nMarkersInChromosome)
x = as.data.frame(t(apply(x, 1, function(i) rev(cumsum(rev(i))))))
colnames(x) = c(">0.6", ">0.7", ">0.8", ">0.9")
propMarkersInLD = apply(x, 2, median, simplify = TRUE)