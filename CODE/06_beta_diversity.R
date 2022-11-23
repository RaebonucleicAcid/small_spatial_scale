### TO-DO: add loop for all physeq objects and maybe 2 or three taxonomic levels EACH
#-------------------------------------------------
###BETA DIVERSITY PLOTS
#Create beta diversity plot outputs directory
dir.beta <- "./OUTPUTS/Beta_Diversity"
dir.create(paste0(dir.beta),
           showWarnings = FALSE, recursive = TRUE)

for (group in groups){
  
  #get physeq object
  physeq_name <- paste0("physeq_", group) 
  physeq_temp <- get(physeq_name)
  
  #get otumat
  otumat_name <- paste0("otumat_", group) 
  otumat_temp <- get(otumat_name)
  
  # Create PDF file of each boxplot
  pdf(file = paste0(
    dir.beta, "/", group, "_beta_diversity_plots.pdf"),
    height = 7, width = 8)
  
  ###NMDS
  OTU_temp <- otu_table(physeq_temp)
  otumat_temp_t <- t(otumat_temp)
  bray <- vegdist(otumat_temp_t, method = "bray")
  distMat = as.dist(bray)

  set.seed(23984)
  NMDS1 = metaMDS(otumat_temp_t, k = 2, trymax=999, autotransform = FALSE, 
                  distance = "bray")
  goodness(NMDS1)
  print(stressplot(NMDS1))
  summary(NMDS1)
  str(NMDS1)
  
  ##PLOTTING NMDS
  p1 <- plot_ordination(physeq_temp, NMDS1, color = "gulf_group",
                        label = "Sample_ID") + geom_point(size = 2) + 
    theme_bw() + 
    scale_colour_brewer(type = "qual", palette = "Set1")

  print(p1)
  
  ###conducting homogeneity of dispersal test + PERMANOVA
  ##see https://rpubs.com/DKCH2020/587758:
  veganotu = function(physeq_temp) {
    require("vegan")
    OTU = otu_table(physeq_temp)
    if (taxa_are_rows(OTU)) {
      OTU = t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  # export data from phyloseq to vegan-compatible object
  physeq_vegan <- veganotu(physeq_temp)
  
  # make a data frame that can be used in vegan from the sample_data in the
  # phyloseq object
  sampledf <- data.frame(sample_data(physeq_temp))
  physeq_BRAY <- vegdist(wisconsin(sqrt(physeq_vegan)), method = "bray")
  
  betadisp_physeq <- betadisper(physeq_BRAY, sampledf$gulf_group)
  betadisper_plot <- boxplot(betadisp_physeq, xlab = "", las = 2, cex.axis = 0.8)
  print(betadisper_plot)
  
  #are the variances the same?
  anova(betadisp_physeq)
  
  ###are the centroids of the clusters seen on the NMDS distinct?
  sampledf <- data.frame(sample_data(physeq_temp))
  adonis(distMat ~ sampledf$gulf_group, data = sampledf)
  
  dev.off()
  
}
