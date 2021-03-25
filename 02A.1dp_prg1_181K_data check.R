#'Data quality check of the initial vcf file
  #'violin plots of read depth
  #'heat map of missingness and allele 
  #'PCoA looking at (non)clustering of pop/gen, batch and checks 

#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

prg.vcf <- read.vcfR("arnel_out.recode.vcf") #189K snps
show(prg.vcf)#189K


#'Prepare data for plotting
#'Read depths
#'
#'extract read depth
dp1 <- extract.gt(prg.vcf,  element = "DP", as.numeric = TRUE)

rm(prg.vcf) #resolve low memory issues
gc()

#'population information
pop<-substr(colnames(dp1),1,3)

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



batch<-vector()
batch<-ifelse(colnames(dp1) %in% batch1,"B1","B2")
rm(batch1)


#'violin plots
library(ggplot2)
library(reshape2)

#'per batch
dp.bat<-dp1
colnames(dp.bat)<-batch

#'per population
dp.pop<-dp1
colnames(dp.pop)<-pop

dp.pop<-dp.pop[,which(colnames(dp.pop)!="GA6")]

#'whole dataset
dpw<-dp1
colnames(dpw)<-NULL

#'per batch/per pop/whole dataset
#'rerun script below, modifying object to melt (dpf) everytime
  #'dp.bat
  #'dp.pop
  #'dpw

dpf <- melt(dpw, varnames = c("Index", "Sample"),
            value.name = "Depth", na.rm = TRUE)
dpf <- dpf[ dpf$Depth > 0, ]

p <- ggplot(dpf, aes(x = Sample, y = Depth))
p <- p + geom_violin(fill = "#C0C0C0", adjust = 1.0,
                     scale = "count", trim = TRUE)
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + scale_y_continuous(trans = scales::log2_trans(),
                            breaks = c(1, 10, 100, 800),
                            minor_breaks = c(1:10, 2:10 * 10, 2:8 * 100))
p <- p + theme(panel.grid.major.y = element_line(color = "#A9A9A9", size = 0.6))
p <- p + theme(panel.grid.minor.y = element_line(color = "#C0C0C0", size = 0.2))
p <- p + ylab("Depth (DP)")
p

#for whole data set
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks = element_blank())
p <- p + scale_x_continuous(breaks = NULL)

#' Save plots everytime
#' Clean memory for a new run
rm(p)
rm(dpf)
rm(dpw)
gc()
