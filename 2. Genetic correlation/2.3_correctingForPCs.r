library(dplyr)
library(ggpubr)
library(psych)


results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
plotting_dir = "plots_092522/"
RObj_dir = "RObjects_092522/"

#load datasets ----
load(paste0(results_dir, RObj_dir, "dataAfterPrep_120821.rda"))

mydata = traitCommonSegregants_std

#remove columns with > 20% missing data

missings = which(colSums(is.na(mydata))/nrow(mydata) > 0.20)
mydata = mydata[,-missings]

#correlation and remove overcorrelated (rho > 0.8)
  #mydata = dplyr::select(mydata, -Trehalose, -Lactate)

completeRows = complete.cases(mydata)

mydata = mydata[completeRows,]

#principal components analysis
fit <- prcomp(mydata, center = TRUE, scale = TRUE)
Y = summary(fit)

x = data.frame(PC1_loadings = fit$rotation[1],
               conditions = colnames(mydata))



#PLOT SCREE PLOT
pdf(paste0(results_dir, plotting_dir, "PC-ScreePlot.pdf"))
scree(mydata, factors = FALSE)
plot(cumsum(fit$sdev^2 / sum(fit$sdev^2)), type="b", ylab = "cumulative variance explained",
     xlab = "principal components")
dev.off()
