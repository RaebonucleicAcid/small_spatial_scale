##TO-DO:
### 1) figure out if these should be generated for chloro + cyano
#--------------------------------------------------------------------------
###ALPHA DIVERSITY PLOTS
#Create alpha diversity plot outputs directory
dir.alpha <- "./OUTPUTS/Alpha_Diversity"
dir.create(paste0(dir.alpha),
           showWarnings = FALSE, recursive = TRUE)

for (group in groups){
  
  #get object
  physeq_name <- paste0("physeq_", group) 
  physeq_temp <- get(physeq_name)
  
  # Create PDF file of each boxplot
  pdf(file = paste0(
    dir.alpha, "/", group, "_alpha_diversity_plots.pdf"),
    height = 7, width = 8)
  
  ##Alpha Diversity:
  alpha_plot_temp <- plot_richness(physeq_temp, measures = 
                                    c("Observed", "Shannon", "Simpson")) + 
    geom_boxplot()
  
  alpha_plot_gulf_grouping_temp <- plot_richness(physeq_temp, x="gulf_group", measures = 
                                                  c("Observed", "Shannon", "Simpson")) + 
    geom_boxplot()
  
  print(alpha_plot_temp)
  print(alpha_plot_gulf_grouping_temp)
  dev.off()
}
  
