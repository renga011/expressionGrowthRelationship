library(ggplot2)

#global variables ----
results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"

#custom color palette for different conditions ----


colorPalette = c("peachpuff","antiquewhite3","rosybrown",
                 "aquamarine", "aquamarine4",
                 "red", "burlywood",
                 "black", "blue", "yellowgreen",
                 "blue4", "blueviolet", "thistle",
                 "brown", "darksalmon", "dodgerblue",
                 "darkcyan","chartreuse", "gold", "magenta",
                 "powderblue", "slateblue", "violetred",
                 "coral1", "cornsilk4", "plum",
                 "cyan", "cyan3",
                 "forestgreen","darkgoldenrod", "chocolate",
                 "darkgoldenrod1", "darkgray", "deepskyblue4",
                 "darkgreen", "darkkhaki", "deepskyblue3",
                 "darkmagenta", "darkolivegreen1", "deepskyblue",
                 "darkorange", "yellow2", "deeppink",
                 "darkseagreen", "darkseagreen1", "olivedrab"
)

names(colorPalette) = colnames(traitCommonSegregants_std)

save(colorPalette, file = paste0(results_dir, RObj_dir, "colorPaletteForConditions.rda"))

#theme for axes and legend ----

theme_textProperties = theme(legend.position = "top", legend.text = element_text(size = 10),
                             legend.title = element_text(size = 12, face = "bold"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             axis.title.x = element_text(size = 14, face = "bold"),
                             axis.title.y = element_text(size = 14, face = "bold"))

save(theme_textProperties, file = paste0(results_dir, RObj_dir, "theme_legendAndAxes.rda"))
