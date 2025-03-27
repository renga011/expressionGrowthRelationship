library(ggpubr)
library(dplyr)
library(stringr)
library(readr)
library(GenomicRanges)

#GLOBAL VARIABLES -----

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

 #global functions ----
addPlottingCoords=  function(QTLs){
  for(i in 1:nrow(QTLs)){
    chr_i = QTLs$chromosome[i]
    additiveFactor = chromosomeSizes$chromosomeStart_plotting[which(chromosomeSizes$chromosome == chr_i)]
    CILeft_i = QTLs$CILeft[i]
    CIRight_i = QTLs$CIRight[i]
    QTLs$plotting_start[i] = CILeft_i + additiveFactor
    QTLs$plotting_end[i] = CIRight_i + additiveFactor
  }
  return(QTLs)
}

#load data -----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "geneList_SGD.rda"))
  geneList_SGD$plotting_position = 0.5*(geneList_SGD$plottingPosition_start + geneList_SGD$plottingPosition_end)

load(paste0(results_dir, RObj_dir, "chromosomeSizes.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))
  hotspotData_Albert2018$chromosome = stringr::str_remove(hotspotData_Albert2018$chromosome, pattern = "chr")
  hotspotData_Albert2018$chromosome = as.integer(as.roman(hotspotData_Albert2018$chromosome))
  
hotspotEffects_growth = readr::read_csv(paste0(results_dir, otherFiles_dir, "hotspotEffectsOnGrowth.csv"))

load(paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))

#prep eQTL table----
eQTLs = dplyr::select(eQTL_Albert2018, gene, r_expression, position_plotting, cis, chromosome)
eQTLs$gene_chr = geneList_SGD$chromosome[match(eQTLs$gene, geneList_SGD$sysName)]
eQTLs$gene_plottingPosition = geneList_SGD$plotting_position[match(eQTLs$gene, geneList_SGD$sysName)]
eQTLs$gene_stdName = geneList_SGD$stdName[match(eQTLs$gene, geneList_SGD$sysName)]

#prep hotspot table ----
hotspotTable = dplyr::select(hotspotData_Albert2018, hotspotMarker, chromosome, 
                             bootstrapIntervalLeft, bootstrapIntervalRight,
                             numberEQTLInHotspot, numberNonzeroEffects)
colnames(hotspotTable) = c("hotspotMarker", "chromosome", "CILeft", "CIRight", 
                           "nEQTLs", "nNonzeroEffects")

hotspotTable = addPlottingCoords(hotspotTable)

  #get nNonZeroEffects. for growth
nNonZeroEffects_growth = sapply(hotspotEffects_growth[-1], FUN = function(hotspot){
  nNonZeroEffects = sum(hotspot!=0)
  return(nNonZeroEffects)
})
hotspotTable$nNonzeroEffects_growth = nNonZeroEffects_growth[match(hotspotTable$hotspotMarker,
                                                                   names(nNonZeroEffects_growth))]

#prep gQTL table ----
gQTLs = dplyr::select(BloomQTL, Trait, Chromosome, position_plotting, r_growth)
gQTLs$rank = dplyr::dense_rank(gQTLs$Trait)

#plot----
chromosomeSizes$plotting_locus= 0.5* (chromosomeSizes$chromosomeStart_plotting + chromosomeSizes$chromosomeEnd_plotting)
chromosomeSizes$label= as.character(as.roman(chromosomeSizes$chromosome))
chromosomeSizes = dplyr::filter(chromosomeSizes, chromosome != 0)

ifelse(!(dir.exists(paste0(results_dir, plotting_dir, "QTLOverlap/"))),
       dir.create(paste0(results_dir, plotting_dir, "QTLOverlap/")),
       FALSE)
pdf(paste0(results_dir, plotting_dir, "QTLOverlap/", "eQTLGQTLChromosomeLayout.pdf"), width = 10, height = 15)

plot1 = ggplot() +
  geom_col(data = hotspotTable,
           mapping = aes(x = 0.5*(plotting_start + plotting_end),
               y = sqrt(nNonzeroEffects)),
           color = "black", fill = "dark gray") +
  geom_point(data = gQTLs, 
             aes(x = position_plotting, 
                 y = -4*rank,
                 size = abs(r_growth),
                 color = ifelse(r_growth >= 0, "Higher in RM", "Higher in BY")), 
             alpha = 0.5) +
  scale_size_binned(breaks = seq(0,1,0.2),
                    range = c(0, 4)) +
  scale_color_manual(values = c("Higher in BY" = "dark blue", "Higher in RM" = "maroon")) +
  geom_hline(yintercept = 0) +
  geom_vline(data = chromosomeSizes,
             aes(xintercept = chromosomeStart_plotting-50000),
             color = "light gray",
             size = 0.5) +
  geom_vline(data = chromosomeSizes,
             aes(xintercept = chromosomeEnd_plotting + 50000),
             color = "light gray",
             size = 0.5) +
  theme_classic() +
  scale_x_continuous(breaks = chromosomeSizes$plotting_locus, 
                     labels = chromosomeSizes$label,
                     name = "QTL position (chromosome)") + 
  scale_y_continuous(breaks = c(-4*gQTLs$rank, 20, 40, 60, sqrt(5643)), 
                     labels = c(gQTLs$Trait, 400, 1600, 3600, 5643),
                     name = "Number of genes \naffected by hotspot") +
  guides(size = guide_legend(title = "gQTL effect\n magnitude",
                             ncol = 4),
         color = guide_legend(title = "gQTL effect\n direction")) +
  theme_textProperties +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
  
print(plot1) 

dev.off()

##########

#overlap gQTL and eQTL hotspots ----
hotspotGRanges = GRanges(seqnames= hotspotData_Albert2018$chromosome,
                         ranges = IRanges(start = hotspotData_Albert2018$bootstrapIntervalLeft,
                                          end = hotspotData_Albert2018$bootstrapIntervalRight),
                         hotspotMarker = hotspotData_Albert2018$hotspotMarker,
                         nNonZeroEffects = hotspotData_Albert2018$numberNonzeroEffects)

growth_gRanges = GRanges(seqnames = BloomQTL$Chromosome,
                         IRanges(start = BloomQTL$LOD_leftCI, 
                                 end = BloomQTL$LOD_rightCI),
                         pMarker = BloomQTL$closestPeakMarkerPosition,
                         LOD = BloomQTL$LOD,
                         growthQTL = BloomQTL$growthQTL,
                         condition = BloomQTL$Trait,
                         r_growth = BloomQTL$r_growth,
                         se_growth = BloomQTL$se_growth,
                         gQTL_CILeft = BloomQTL$LOD_leftCI,
                         gQTL_CIRight = BloomQTL$LOD_rightCI)

gQTLs$overlapsWithHotspots = countOverlaps(growth_gRanges, hotspotGRanges)
gQTLs$growthQTL = BloomQTL$growthQTL
  #number of gQTLs overlapping hotspot
length(which(gQTLs$overlapsWithHotspots != 0))

## get expectation of number of gQTL/hotspot overlaps by chance - raindown random gQTLs and see how many overlap hotspots
load(paste0(results_dir, RObj_dir, "randomGQTLs.rda"))

expectedNoOfOverlappingGQTLs = lapply(1:1000,
                          FUN = function(n){
                            gQTL = randomGQTLs_all[[n]]
                            growth_gRanges = GRanges(seqnames = gQTL$chromosome,
                                                     IRanges(start = gQTL$CILeft, 
                                                             end = gQTL$CIRight),
                                                     condition = gQTL$trait)
                            overlaps_i = countOverlaps(growth_gRanges, hotspotGRanges)
                            return(length(which(overlaps_i != 0)))
                          })
expectedNoOfOverlappingGQTLs = unlist(expectedNoOfOverlappingGQTLs)

overlappingGQTLs_Actual = countOverlaps(growth_gRanges, hotspotGRanges)
numberOverlappingGQTLs_actual = length(which(overlappingGQTLs_Actual != 0))

#get p-value
n0 = length(which(expectedNoOfOverlappingGQTLs >= numberOverlappingGQTLs_actual))
p = n0/1000

#plot histogram
pdf(paste0(results_dir, plotting_dir, "QTLOverlap/", "expectedGQTLOverlapWithHotspots_histogram.pdf"))
gghistogram(data = data.frame(n = expectedNoOfOverlappingGQTLs),
            x = "n",
            xlab = "number of gQTLs overlapping trans eQTL hotspots") + 
  geom_vline(xintercept = numberOverlappingGQTLs_actual, color = "red")
dev.off()

#######
## Reviewer comment - what proportion of 102 trans eQTL hotspots is overlapping a gQTL
expectedNoOfOverlappingHotspots = lapply(1:1000,
                                      FUN = function(n){
                                        gQTL = randomGQTLs_all[[n]]
                                        growth_gRanges = GRanges(seqnames = gQTL$chromosome,
                                                                 IRanges(start = gQTL$CILeft, 
                                                                         end = gQTL$CIRight),
                                                                 condition = gQTL$trait)
                                        overlaps_i = countOverlaps(hotspotGRanges, growth_gRanges)
                                        return(length(which(overlaps_i != 0)))
                                      })
expectedNoOfOverlappingHotspots = unlist(expectedNoOfOverlappingHotspots)

overlappingHotspots_Actual = countOverlaps(hotspotGRanges, growth_gRanges)
numberOverlappingHotspots_actual = length(which(overlappingHotspots_Actual != 0))

  # plot histogram

pdf(paste0(results_dir, plotting_dir, "QTLOverlap/", "expectedHotspotOverlapWithGQTL_histogram.pdf"))
gghistogram(data = data.frame(n = expectedNoOfOverlappingHotspots),
            x = "n",
            xlab = "number of trans eQTL hotspots overlapping gQTLs") + 
  geom_vline(xintercept = numberOverlappingHotspots_actual, color = "red")+ xlim(80,102)
dev.off()

  #p value
n0_hotspots = length(which(expectedNoOfOverlappingHotspots <= numberOverlappingHotspots_actual))
p = n0_hotspots/1000

#######
## Reviewer comment: difference between the gQTLs that overlap a hotspot vs not

gQTLs$hotspotOverlapStatus = ifelse(gQTLs$overlapsWithHotspots !=0, TRUE, FALSE)
diffBetweenGQTLs = gQTLs %>% group_by(hotspotOverlapStatus) %>% 
  summarise(avg_effectSizeMagnitude = median(abs(r_growth)), 
            n_plus = sum(r_growth > 0), n_minus = sum(r_growth < 0), 
            n_plus_normalized = n_plus/n(), n_minus_normalized = n_minus/n())

  #do some statistical tests to establish significant difference or not.
  ## difference in magnitude of effect size
  
  wilcox.test(x = abs(gQTLs$r_growth[which(gQTLs$hotspotOverlapStatus == "TRUE")]), 
              y = abs(gQTLs$r_growth[which(gQTLs$hotspotOverlapStatus == "FALSE")]))
  
  # W = 33750, p = 0.017
  
  ## difference in direction of effect size
  
  fisher.test(as.matrix(dplyr::select(diffBetweenGQTLs, n_plus, n_minus)))
  
  # odds ratio = 0.97, p = 0.93
  
  ## count number per condition
  diffBetweenGQTLs_conditionWise = gQTLs %>% 
    group_by(Trait, hotspotOverlapStatus) %>%
    summarize(n = n())
  
  pdf(paste0(results_dir, plotting_dir, "revision_figures_cellGenomics/", "gQTLs_overlappingHotspotsVsNot_conditionWiseNumber.pdf"), width = 10, height = 10)
  ggbarplot(diffBetweenGQTLs_conditionWise,
            x = "Trait",
            y = "n",
            fill = "hotspotOverlapStatus") + theme_textProperties + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Number of gQTLs")

    dev.off()


