rm(list=ls())

##setwd
setwd("~/Desktop/SPATIAL/small_spatial_scale")

##source needed functions
source("./FUNCTIONS/distance_decay.R")
  
##load packages
library("ggplot2");
library("vegan");
library("dplyr");
library("phyloseq");
library("data.table");
library("tidyr")

# Set colors for plotting
stacked_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "steelblue2","#618f5e",
  "#6e5e8f", "#eb102e", "#673770","#D14285", "#652926", "#C84248"
)
