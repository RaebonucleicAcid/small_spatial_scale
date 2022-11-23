### Show most interesting plot in the paper and put the rest in supplementary figures
#-------------------------------------------------
###DISTANCE DECAY PLOTS
#Create distance decay outputs directory
dir.dd <- "./OUTPUTS/Distance_Decay"
dir.create(paste0(dir.dd),
           showWarnings = FALSE, recursive = TRUE)

groups <- c("prok", "euk", "cyano", "chloro", "Alphaproteobacteria", "Flavobacteriales")

for (group in groups){
  
  #get physeq object
  physeq_name <- paste0("rel_physeq_", group) 
  physeq_temp <- get(physeq_name)
  
  # Create PDF file of each boxplot
  pdf(file = paste0(
    dir.dd, "/", group, "_distance_decay_plots.pdf"),
    height = 7, width = 8)
  
  ##making plot for distance decay
  metadata_temp <- data.frame(sample_data(physeq_temp))
  comm_physeq_temp <- t(otu_table(physeq_temp))
  dist <- data.frame(row.names(comm_physeq_temp))
  geo_physeq <- data.frame("lat" = metadata_temp$Latitude, 
                           "long" = metadata_temp$Lontitude)
  
  dist$lat<-geo_physeq$lat
  dist$lon<-geo_physeq$long
  colnames(dist)[colnames(dist)=="row.names.comm_physeq_temp."] <- "name"
  
  dist_m <-GeoDistanceInMetresMatrix(dist) 
  #this does not round and keeps in meters
  dist_m[upper.tri(dist_m)] <- NA
  
  OTU1 = as(otu_table(physeq_temp), "matrix")
  OTU1<-t(OTU1)
  commdist <- vegdist(OTU1, method="bray",na.rm=TRUE,binary=FALSE)
  #binary = TRUE makes it sorenson
  
  ###getting average and stdev of Bray Curtis values 
  commdist_onerow <- as.matrix(commdist)
  commdist_onerow <- commdist_onerow[1,]
  commdist_onerow <- as.matrix(commdist_onerow)
  commdist_onerow$dissimilarity <- 1-commdist_onerow
  mean <- mean(commdist_onerow$dissimilarity)
  sd <- sd(commdist_onerow$dissimilarity)
  
  #this makes this a 3 column dataframe instead of matrix, 
  # melting the data
  comm.dist.ls <- liste(commdist, entry="comm")
  coord.dist.ls <- liste(dist_m, entry="dist")
  coord.dist.ls <- na.omit(coord.dist.ls)
  coord.dist.ls <- coord.dist.ls[!(coord.dist.ls$NBX==coord.dist.ls$NBY),]
  coord.dist.ls <- as.data.frame(coord.dist.ls)
  dim(as.data.frame(coord.dist.ls, rownames = FALSE))
  comm.dist.ls$comm <-comm.dist.ls$comm-0.0000001
  

  ###log model
  model_log<-lm(log(1-comm.dist.ls$comm)~log(coord.dist.ls$dist))
  summary(model_log)
  
  ### non-log model
  model <- lm((1-comm.dist.ls$comm)~(coord.dist.ls$dist))
  p_value <- summary(model)$coefficients[8]
  slope <- summary(model)$coefficients[2]
  summary(model)
  
  ##plotting
  df <- data.frame(sim = (1-comm.dist.ls$comm), dist = (coord.dist.ls$dist))
  
  p <- ggplot(df, aes(x = dist, y =sim)) + geom_point()+
    geom_smooth(method = "lm", se=FALSE, color="red", aes(group=1))+
    labs(y = "Bray-Curtis Similarity")  +
    labs(x = "Spatial distance (m)")+ 
    theme_bw() + theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"), 
                       axis.title.x = element_text(face="bold", size=16),
                       axis.title.y = element_text(face="bold", size=16)) + 
    ylim(0, 0.9) + 
    labs(subtitle = paste0("slope: ", slope, "\n p-value: ", p_value),
         title = paste0("Distance decay: ",group)) +
    theme(plot.subtitle=element_text(size=10, hjust=0)) +
    theme(plot.title=element_text(size=16, hjust=0))
  
  
  print(p)
  dev.off()
}
