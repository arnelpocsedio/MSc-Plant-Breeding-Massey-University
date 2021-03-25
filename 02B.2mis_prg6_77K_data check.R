#'Data quality check of the filtered vcf file
#'violin plots of read depth
#'heat map of missingness and allele 
#'PCoA looking at (non)clustering of pop/gen, batch and checks 
#'

#'Missing data
#'
#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

prg6.vcf <- read.vcfR("01_prg6_snptaxafiltered_pol.vcf.gz") 
show(prg6.vcf)#77,261 snps; 180 samples

#'Prepare data 
library(poppr)
prg6.gl<- vcfR2genlight(prg6.vcf) #7789 non-biallelic snps removed
ploidy(prg6.gl)<-2



#'subset batch and pop
#'
#'population info
pop(prg6.gl)<-substr(indNames(prg6.gl),1,3)

#'batch information
batch1<-c("P42_A1",	"P42_A2",	"P42_A3", "P42_A4",	"P42_A5",
          "P42_A6", "P46_A7",	"P46_A8",	"P46_A9",	"P46_A10", 
          "P46_A11", "P42_E3_02", "P42_B1",	"P42_B2", "P42_B4",	
          "P42_B5",	"P42_B6",	"P46_B7",	"P46_B8",	"P46_B9",	
          "P46_B10", "P46_B11", "P46_B12", "P42_C1", "P42_C2",	
          "P42_C3",	"P42_C4",	"P42_C5",	"P42_C6",	"P46_C7",
          "P46_C8",	"P46_C9",	"P46_C10",	"P46_C11", "P46_C12", 
          "P42_D1","P42_D2",	"P42_D3",	"P42_D4",	"P42_D5",	
          "P42_D6",	"P46_D7", "P46_D8",	"P46_D9",	"P46_A6_02",	
          "P46_D11",	"GA66_SQ2718", "P46_A6_01", "P42_E2", 
          "P42_E3_01",	"P42_E4",	"P42_E5",	"P42_E6", "P46_E7",	
          "P46_E8",	"P46_E9",	"P46_E10",	"P46_E11",	"P46_E12",
          "P42_F1",	"P42_F2",	"P42_F3",	"P42_F4",	"P42_F5",	"P42_F6",
          "P46_F7",	"P46_F8",	"P46_F9",	"P46_F10", "P46_F11",
          "P46_F12", "P42_G1", "P42_G2",	"P42_G3",	"P42_G4",	
          "P42_G5",	"P42_G6", "P46_G7",	"P46_G8",	"P46_G9",	
          "P46_G10",	"P46_G11",	"P46_G12", "P42_H1", "P42_H2",	
          "P42_H3",	"P42_H4",	"P42_H5",	"P42_H6", "P46_H7",
          "P46_H8",	"P46_H9",	"P46_H10",	"P46_H11",	"P46_H12")

#'treat batch info as population

prg6v2.gl<-prg6.gl
pop(prg6v2.gl)<-ifelse(indNames(prg6v2.gl) %in% batch1,"B1","B2")

#'separate data into population and batches
#'
prg_pop<-seppop(prg6.gl)
names(prg_pop)

prg_bat<-seppop(prg6v2.gl)
names(prg_bat)

#missing data heatmap
glPlot(prg6.gl, col = c("black", "black", "black"), legend = FALSE)#whole data
glPlot(prg_pop$P42, col = c("black", "black", "black"), legend = FALSE)
glPlot(prg_pop$P46, col = c("black", "black", "black"), legend = FALSE)
glPlot(prg_bat$B1, col = c("black", "black", "black"), legend = FALSE)
glPlot(prg_bat$B2, col = c("black", "black", "black"), legend = FALSE)

#alleles heatmap
library(RColorBrewer)
mycol<-brewer.pal(9,"Set1")
glPlot(prg_pop$P46, col = c("#77AADD", "#4477AA", "#114477"))


#clean mem
rm(gt)
rm(prg_b1.gl)
rm(prg_b2.gl)
rm(prg_b1.VCF)
rm(prg_b2.VCF)
rm(prg_p42.gl)
rm(prg_p46.gl)
rm(prg_p42.VCF)
rm(prg_p46.VCF)
rm(prg.VCF)
gc()