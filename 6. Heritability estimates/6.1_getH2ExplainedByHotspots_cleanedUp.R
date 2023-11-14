#get the h2 explained by hotspots for each of the 46 conditions 

devtools::install_github("variani/lme4qtl")

library(lme4qtl)
library(parallel)
library(parallelly)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load data

load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
hotspotData_Albert2018 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"hotspotData_Albert2018.xlsx"))

#MAKE GRMs ------
#grm for hotspot markers
grm_h = cov(t(dplyr::select(genotypesCommonSegregants, hotspotData_Albert2018$hotspotMarker)))
colnames(grm_h) = 1:979
row.names(grm_h) = 1:979

#grm all markers
grm = cov(t(genotypesCommonSegregants))
colnames(grm)= 1:979
row.names(grm)= 1:979

#MAKE MIXED MODEL AND GET PROPORTION OF VAR EXPLAINED BY THE GRMs -----
getPropVarExplained = lapply(colnames(traitCommonSegregants_std), 
                             function(condition_i){
                               print(condition_i)
                               
                               data1 = data.frame(trait = unlist(traitCommonSegregants_std[condition_i]),
                                                  A = 1:979,
                                                  H = 1:979,
                                                  G = 1:979,
                                                  T3 = 1:979)
                               
                               #all genotype markers ---
                               m_all = relmatLmer(trait ~ (1|A), data1, relmat = list(A = grm))
                               
                                 #get proportion of variance explained
                                 A = VarProp(m_all)
                               
                               #all hotspots ---
                               m_hotspots = relmatLmer(trait ~ (1|H), data1, relmat = list(H = grm_h))
                               
                                 #get proportion of variance explained
                                 H = VarProp(m_hotspots)
                               
                               #all gQTLs ---
                                ##make a grm
                                 gQTLs = dplyr::filter(BloomQTL, Trait == condition_i)
                                 print(nrow(gQTLs))
                                 genotypeMarkers = stringr::str_split(colnames(genotypesCommonSegregants), pattern = "_", simplify = TRUE)[,1]
                                 col_i = which(genotypeMarkers %in% gQTLs$growthQTL)
                                 genotypeMatrix = genotypesCommonSegregants[col_i]
                                 
                                 grm_g = cov(t(genotypeMatrix))
                                 colnames(grm_g) = 1:979
                                 row.names(grm_g) = 1:979
                                 
                                ##mixed model
                                 m_gQTLs = relmatLmer(trait ~ (1|G), data1, relmat = list(G = grm_g))
                                  ##proportion variance explained
                                  G = VarProp(m_gQTLs)
                                  
                                ##top 3 hotspots ---
                                  ##make grm
                                grm_top3Hotspots = cov(t(dplyr::select(genotypesCommonSegregants,
                                                                       `chrXIV:466588_T/G`, #MKT1
                                                                       `chrXII:657022_T/C`, #HAP1
                                                                       `chrXV:171150_T/C`, #IRA2
                                                                       )))
                                  colnames(grm_top3Hotspots) = 1:979
                                  row.names(grm_top3Hotspots) = 1:979
                                  
                                  ##mixed model
                                  m_T3 = relmatLmer(trait ~ (1|T3), data1, relmat = list(T3 = grm_top3Hotspots))
                                 
                                    ##proportion variance explained
                                  T3 = VarProp(m_T3)
                               return(c("H" = H$prop[which(H$grp == "H")],
                                        "A" = A$prop[which(A$grp == "A")],
                                        "G" = G$prop[which(G$grp == "G")],
                                        "T3" = T3$prop[which(T3$grp == "T3")]))    

})

propVarExplained_allConditions = as.data.frame(do.call(rbind, getPropVarExplained))
propVarExplained_allConditions$condition  = colnames(traitCommonSegregants_std)

#save and plot ---
save(propVarExplained_allConditions, file = paste0(results_dir, RObj_dir, "h2Estimates_allVsHotspot_allConditions.rda"))
