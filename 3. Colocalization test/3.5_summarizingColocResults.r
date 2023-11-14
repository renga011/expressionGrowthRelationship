setwd("~/ExpressionPhenotypeProject/") #for working on MSI

#global variables
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
coloc_dir = "LD95/"
fileList = "fileList_QTLLODsOver10_final.rda"

#load data ----
load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", fileList))
#make list of all coloc done ----
list_colocDone = list()
colocResults_df = data.frame(filename = filename_list,
                             coloc_p = NA,
                             coloc_significance = NA)
for(filename in filename_list){
  if(file.exists(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", coloc_dir, filename))){
  load(paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", coloc_dir, filename))
  list_colocDone[[filename]]= eQTLTable_QTLThis
  colocResults_df$coloc_p[which(colocResults_df$filename == filename)] = eQTLTable_QTLThis$coloc_pValue
  colocResults_df$coloc_significance[which(colocResults_df$filename == filename)] = eQTLTable_QTLThis$coloc_significance
  }
}

QTLs_colocDone = do.call(rbind, list_colocDone)

#save ----
save(QTLs_colocDone, colocResults_df, 
     file = paste0(results_dir, RObj_dir, "colocalizationTest/cisEQTLs/", "colocDone_QTLLODGreaterThan10_final.rda"))


