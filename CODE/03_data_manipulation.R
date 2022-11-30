
###Making ALL taxa-specific physeq objects:
groups <- c("prok", "euk", "chloro")
SAMPLE_DATA <- sample_data(metadata)
taxa_levels <- c("Phylum", "Family", "Genus", "Species")
filter <- 0.15

for (group in groups) {
  ##get otu and tax tables
  otumat_temp_name <- paste0("otumat_", group)
  taxmat_temp_name <- paste0("taxmat_", group)
  otumat_temp <- get(otumat_temp_name)
  taxmat_temp <- get(taxmat_temp_name)
  
  ##make phyloseq object
  OTU_temp <- otu_table(otumat_temp, taxa_are_rows = TRUE);
  TAX_temp <- tax_table(taxmat_temp);
  physeq_temp <- phyloseq(OTU_temp, TAX_temp, SAMPLE_DATA)
  physeq_name <- paste0("physeq_", group)
  assign(physeq_name, physeq_temp)
  
  #transform to relative abundance + make separate Phyloseq object:
  rel_physeq_temp <- transform_sample_counts(physeq_temp, function(x) x / sum(x))
  rel_physeq_name <- paste0("rel_physeq_", group)
  assign(rel_physeq_name, rel_physeq_temp)
  
  #remove singletons + make separate phyloseq object with pruned OTU table
  #so that both are available for subsequent analyses
  pruned_otu_temp <- prune_taxa(taxa_sums(OTU_temp) > 1, OTU_temp)
  pruned_physeq_name <- paste0("pruned_physeq_", group)
  pruned_physeq_temp <- phyloseq(pruned_otu_temp, TAX_temp, SAMPLE_DATA)
  assign(pruned_physeq_name,pruned_physeq_temp)  
  
  ###MAKE ALL SUBSETTED TAXONOMIC-LEVEL PHYSEQ OBJECTS FOR EACH GROUP
  for (taxa_level in taxa_levels){
    
    #filter out low abundance taxa for visualization purposes
    physeq_taxa <- filter_taxa(rel_physeq_temp, function(x) sum(x) > filter, TRUE) 
    
    #agglomerate at taxa level and melt df for plotting
    taxa_plot_temp <- physeq_taxa %>%
      tax_glom(taxrank = taxa_level) %>%         # agglomerate at phylum level
      psmelt() %>%                             # Melt to long format
      arrange(Sample)                   # Sort data frame alphabetically by phylum
    
    #assign df to perm object 
    taxa_plot_temp_name <- paste0(group,"_",taxa_level,"_plot")
    assign(taxa_plot_temp_name, taxa_plot_temp)
  }
}
### SUBSETTING OTHER TAXA into objects
#alphaproteobacteria
physeq_Alphaproteobacteria <- subset_taxa(physeq_prok, Class==" Alphaproteobacteria")
rel_physeq_Alphaproteobacteria <- transform_sample_counts(physeq_Alphaproteobacteria, function(x) x / sum(x))

#flavo
physeq_Flavobacteriales <- subset_taxa(physeq_prok, Order==" Flavobacteriales")
rel_physeq_Flavobacteriales <- transform_sample_counts(physeq_Flavobacteriales, function(x) x / sum(x))

#flavo
physeq_cyano <- subset_taxa(physeq_prok, Phylum==" Cyanobacteria")
rel_physeq_cyano <- transform_sample_counts(physeq_cyano, function(x) x / sum(x))

