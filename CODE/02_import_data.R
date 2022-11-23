###DATA MANIPULATION

#setwd
setwd("~/Desktop/SPATIAL/")

##prok: READ IN OTU + TAX + metadata & re-format:
#OTU table
otumat_prok <- read.csv("./DATA/prok/OTU_table.tsv", header = TRUE, sep="\t");
row.names(otumat_prok) <- otumat_prok$OTU
otumat_prok <- subset(otumat_prok, select = 
                   -c(OTU, F028, F062, F103, F010, F012, F040, F051, F085));
otumat_prok <- as.matrix(otumat_prok)

#Taxonomy table
taxmat_prok <- read.csv("./DATA/prok/TAX_table.tsv", header = TRUE, sep = "\t");
row.names(taxmat_prok) <- taxmat_prok$OTU ;
taxmat_prok <- subset(taxmat_prok, select = -c(OTU));
taxmat_prok <- as.matrix(taxmat_prok)

#Metadata, also adding gulf stream grouping 
metadata <- read.csv("./DATA/prok/Florida_subset_metadata.tsv", header = TRUE, sep="\t");
row.names(metadata) <- metadata$Sample_ID;
gulf_group <- c((rep(x="A", times=38)),(rep(x="B", times=11)))
gulf_group <- as.character(gulf_group)
metadata <- cbind(metadata, gulf_group)

############################################################
##euk: READ IN OTU + TAX & re-format:
#OTU table
otumat_euk <- read.csv("./DATA/euk/OTU_table.tsv", header = TRUE, sep="\t");
row.names(otumat_euk) <- otumat_euk$OUT_ID
otumat_euk <- subset(otumat_euk, select = -c(OUT_ID, F108));
otumat_euk <- as.matrix(otumat_euk)

taxmat_euk <- read.csv("./DATA/euk/TAX_table.tsv", header = TRUE, sep = "\t");
row.names(taxmat_euk) <- taxmat_euk$OTU ;
taxmat_euk <- subset(taxmat_euk, select = -c(OTU));
taxmat_euk <- as.matrix(taxmat_euk)

############################################################
##CYANOBACTERIA: READ IN OTU + TAX & re-format:
#OTU table
otumat_cyano <- read.csv("./DATA/cyano/cyano_OTU_table.tsv", header = TRUE, sep="\t");
row.names(otumat_cyano) <- otumat_cyano$OTU
otumat_cyano <- subset(otumat_cyano, select = -c(OTU, F028, F062, F103, F010, F012, F040, F041, F051, F085, F099, F102));
otumat_cyano <- as.matrix(otumat_cyano)

#Taxonomy table
taxmat_cyano <- read.csv("./DATA/cyano/cyano_TAX_table.tsv", header = TRUE, sep = "\t");
row.names(taxmat_cyano) <- taxmat_cyano$OTU ;
taxmat_cyano <- subset(taxmat_cyano, select = -c(OTU));
taxmat_cyano <- as.matrix(taxmat_cyano)

############################################################
##CHLOROPLASTS: READ IN OTU + TAX & re-format:
#OTU table
otumat_chloro <- read.csv("./DATA/chloro/OTU_table.tsv", header = TRUE, sep="\t");
row.names(otumat_chloro) <- otumat_chloro$OTU
#removing some low sequencing depth samples
otumat_chloro <- subset(otumat_chloro, select = -c(OTU, F028, F040, F085, F103)); 
otumat_chloro <- as.matrix(otumat_chloro)

taxmat_chloro <- read.csv("./DATA/chloro/TAX_table.tsv", header = TRUE, sep = "\t");
row.names(taxmat_chloro) <- taxmat_chloro$OTU ;
taxmat_chloro <- subset(taxmat_chloro, select = -c(OTU));
taxmat_chloro <- as.matrix(taxmat_chloro)

