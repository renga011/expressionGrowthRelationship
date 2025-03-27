#compile the bootstrapping results and get standard error of h2 estimates

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"
plotting_dir = "plots_092522/"
otherFiles_dir = "otherFiles_101522/"

#load datasets----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "h2Estimates_allVsHotspot_allConditions.rda"))

#compile the bootstrapped datasets first
bootH2Estimates = data.frame()
for(i in 0:49){
  load(paste0(results_dir, RObj_dir, "h2_bootEstimates/", "h2Estimates_boot_", i, ".rda"))
  df = as.data.frame(do.call(rbind, propVarExplained_allConditions_boot))
  df = df %>%
    mutate(across(-all_of("condition"), as.numeric))
  bootH2Estimates = rbind(bootH2Estimates, df)
}

bootH2_molten = reshape2::melt(bootH2Estimates, id.vars = "condition")
se_h2 = bootH2_molten %>% group_by(variable, condition) %>% summarise(se_h2 = sd(value, na.rm = TRUE))

#save
save(se_h2, file = paste0(results_dir, RObj_dir, "h2_se.rda"))


