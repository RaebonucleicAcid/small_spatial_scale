##To-do: make plots look nicer across all levels 
##and maybe change from a minimum abundance threshold 
##to a specific number of top features
#--------------------------------------------------------------------------
###RELATIVE ABUNDANCE PLOTS
#Create stacked bar plot outputs directory
dir.stacked <- "./OUTPUTS/Stacked_Bar_Plots/"
dir.create(dir.stacked, showWarnings = FALSE, recursive = TRUE)

for (group in groups){
  
  #get object
  rel_physeq_name <- paste0("rel_physeq_", group) 
  rel_physeq <- get(rel_physeq_name)
  
  # Create multi-page PDF file of boxplots
  pdf(file = paste0(
    dir.stacked, group, "_Stacked_Bar_Plots.pdf"),
    height = 7, width = 8)
  
  for (taxa_level in taxa_levels){
    #get plot df
    taxa_plot_temp_name <- paste0(group,"_",taxa_level,"_plot")
    taxa_plot_temp <- get(taxa_plot_temp_name)

    ##create sub-directory for specific tax-level
    dir.create(paste0(dir.stacked, "/", taxa_level),
               showWarnings = FALSE, recursive = TRUE)
    
    ##PLOTTING
    # Make abundance plot:
    fill_var <- rlang::sym(taxa_level)
    
    plot_temp <- ggplot(data=taxa_plot_temp, 
                        aes_string(x = "Sample",y = "Abundance", fill = taxa_level)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(values = stacked_colors) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab(paste0("Relative Abundance ( >", filter*100 ,"%) \n")) +
      theme(legend.key.size = unit(0.5, "cm"),
            legend.title = element_text(size = 12, face = "bold"), 
            legend.text = element_text(size = 10)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(shape = guide_legend(override.aes = list(size = 0.2))) +
      guides(color = guide_legend(override.aes = list(size = 0.2))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size =10)) +
      ggtitle(paste0(group, " ", taxa_level,"-level stacked bar plot"))
    
    print(plot_temp) 
  }
  dev.off()
}
