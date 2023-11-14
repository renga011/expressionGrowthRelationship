library(dplyr)

#GLOBAL variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#data ----
load(paste0(results_dir, RObj_dir, "phenotypicCorrelationTable_pearson.rda"))

#get the average phenotypic R for genes that belong to the UP, DOWN & UNRESPONSIVE groups in Brauer et al. 2008 -----
Brauer_S1 = readxl::read_excel(paste0(results_dir, otherFiles_dir,"Brauer_S1.xls"),
                               col_types = c("text", "text", "text", "text", "numeric", 
                                             "numeric","numeric", "text", "text", "text",
                                             "numeric", "numeric", "numeric", "numeric", 
                                             "numeric", "text")) #supplementary table 1 from Brauer et al. 2008

up = Brauer_S1$ORF[which(Brauer_S1$`Growth Rate Response @1.5SD`=="up")]
down = Brauer_S1$ORF[which(Brauer_S1$`Growth Rate Response @1.5SD`=="down")]
unresponsive = Brauer_S1$ORF[which(Brauer_S1$`Unresponsive?` == "yes")]

#get the average phenotypic R for genes that are up and down regulated in the ESR response - Gasch et al. 2000 ------

up_ESR = Brauer_S1$ORF[which(Brauer_S1$ESR =="up")]
down_ESR = Brauer_S1$ORF[which(Brauer_S1$ESR =="down")]

#
aveGeneticR_brauerGroups = sapply(colnames(traitCommonSegregants_std),
                                  FUN = function(condition_i){
                                    r_i = phenotypicCorrelations[[condition_i]]
                                    #r_i = filter(r_i, p < 0.05)
                                    up_r = dplyr::filter(r_i, gene %in% up)$r
                                    down_r = dplyr::filter(r_i, gene %in% down)$r
                                    unresponsive_r = dplyr::filter(r_i, gene %in% unresponsive)$r
                                    return(c("up in fast growth" = mean(up_r, na.rm = TRUE),
                                             "up in slow growth" = mean(down_r, na.rm = TRUE),
                                             "unresponsive" = mean(unresponsive_r, na.rm = TRUE)))
                                  })
names(aveGeneticR_brauerGroups) = colnames(traitCommonSegregants_std)

aveGeneticR_ESRGroups = sapply(colnames(traitCommonSegregants_std),
                               FUN = function(condition_i){
                                 r_i = phenotypicCorrelations[[condition_i]]
                                 up_r = dplyr::filter(r_i, gene %in% up_ESR)$r
                                 down_r = dplyr::filter(r_i, gene %in% down_ESR)$r
                                 return(c("up in ESR" = mean(up_r),
                                          "down in ESR" = mean(down_r)
                                 ))
                               })

names(aveGeneticR_ESRGroups) = colnames(traitCommonSegregants_std)

# save
save(aveGeneticR_brauerGroups, aveGeneticR_ESRGroups, file = paste0(results_dir, RObj_dir, "aveGeneticR.rda"))


