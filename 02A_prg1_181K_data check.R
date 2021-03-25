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

rm(prg.vcf) #low memory
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

dpf <- melt(dp.pop, varnames = c("Index", "Sample"),
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

#' Save plots everytime
#' Clean memory for a new run
rm(p)
rm(dpf)
rm(dp.bat)
gc()


#'Missing data
#'
#'Extract genotype matrix (gt)
gt <- extract.gt(prg.VCF,  element = "GT")#extract gt

#subset batch and pop
batch<-ifelse(colnames(gt) %in% batch1,"B1","B2")
pop<-substr(colnames(gt),1,3)
rm(batch1)

prg_b1.VCF<-prg.VCF
prg_b1.VCF@gt <- prg_b1.VCF@gt[, c(TRUE, batch=="B1")]

prg_b2.VCF<-prg.VCF
prg_b2.VCF@gt <- prg_b2.VCF@gt[, c(TRUE, batch=="B2")]

prg_p42.VCF<-prg.VCF
prg_p42.VCF@gt <- prg_p42.VCF@gt[, c(TRUE, pop=="P42")]

prg_p46.VCF<-prg.VCF
prg_p46.VCF@gt <- prg_p46.VCF@gt[, c(TRUE, pop=="P46")]

#convert to genlight obj
#run glplot for allele&missing heatmap
library(poppr)

#per batch
prg_b1.gl<- vcfR2genlight(prg_b1.VCF) #7789 non-biallelic snps removed
glPlot(prg_b1.gl, col = c("green", "blue" ,"red"))
glPlot(prg_b1.gl, col = c("black", "black" ,"black"), legend = FALSE)

prg_b2.gl<- vcfR2genlight(prg_b2.VCF) #7789 non-biallelic snps removed
glPlot(prg_b2.gl, col = c("green", "blue" ,"red"))
glPlot(prg_b2.gl, col = c("black", "black" ,"black"), legend = FALSE)


#per pop
prg_p42.gl<- vcfR2genlight(prg_p42.VCF) #7789 non-biallelic snps removed
glPlot(prg_p42.gl, col = c("green", "blue" ,"red"))
glPlot(prg_p42.gl, col = c("black", "black" ,"black"), legend = FALSE)

prg_p46.gl<- vcfR2genlight(prg_p46.VCF) #7789 non-biallelic snps removed
glPlot(prg_p46.gl, col = c("green", "blue" ,"red"))
glPlot(prg_p46.gl, col = c("black", "black" ,"black"), legend = FALSE)

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

#pca
#using fast PCA implementation in SNPrelate

library(SNPRelate)
vcf.fn <- "C:/Users/PB-ARNEL/Documents/Thesis2/GBS SNP files/arnel_out.recode.vcf"
snpgdsVCF2GDS_R(vcf.fn, "prg1.gds", method="biallelic.only")
snpgdsSummary("prg1.gds") #non-biallelic snps not included

genofile <- snpgdsOpen("prg1.gds", readonly = FALSE)

#maybe not needed, subset in the final output
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code<-substr(sample.id,1,3)
samp.annot <- data.frame(pop.group = pop_code )
add.gdsn(genofile, "sample.annot", samp.annot)
pop <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))


#snpgdsClose(genofile)


set.seed(999)
#pca
pca1 <- snpgdsPCA(genofile, algorithm = "randomized")

#pop and batch
pop<-substr(pca1$sample.id,1,3)
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
batch<-ifelse(pca1$sample.id %in% batch1,"B1","B2")

#combine to df
tab <- data.frame(sample.id = pca1$sample.id,
                  pop, batch,
                  EV1 = pca1$eigenvect[,1],    # the first eigenvector
                  EV2 = pca1$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

tab<-tab[order(tab$pop,tab$sample.id),]

biorep<- c("P42_E3_01","P42_E3_02", "P42_E3_03",
           "P42_F5","P42_F5_02", "P42_F5_03", 
           "P46_A6_01","P46_A6_02", 
           "P46_A10","P46_A10_02", "P46_A10_03",  
           "P46_D11","P46_D11_02", 
           "GA66_SQ2761", "GA66_SQ2718")

#P42_F5 NOT P42_F5_01 
#P46_A10 NOT P46_A10_01 
#P46_D11 NOT P46_D11_01

tab2<-tab[-which(tab$sample.id %in% biorep),]
tabrep<-tab[which(tab$sample.id %in% biorep),]

plot(x=tab2$EV1, y=tab2$EV2,
     ylim = c(-0.20,0.20), xlim = c(-0.20, 0.20),
     xlab="eigenvector 1", ylab="eigenvector 2",
     col=as.factor(tab2$pop), pch=as.numeric(substr(tab2$batch,2,2))+15)

sym<-c("o","+",
       "o","+","x",
       "o","+","x",
       "o","+","x",
       "o","+",
       "o","+")


points(x=tabrep$EV1, y=tabrep$EV2, 
       col=c(rep("green",2),rep("azure4",6), rep("brown4",7)),
       pch=sym)

repname<-c("GA66","P42_E3", "P42_F5", "P46_A10",
           "P46_A6","P46_D11")

text(x=tabrep[c(1,4,7,11,12,14),]$EV1, 
     y=tabrep[c(1,4,7,11,12,14),]$EV2,
     labels =repname,
     pos = c(1,1,1,2,1,4),
     col = c("green", rep("azure4",2), rep("brown4",3))
     )

#optional legend
leg<-c("P42 batch1", "P46 batch1",
       "P42 batch2", "P46 batch2",
       "Reps")

legend("bottomright", legend=leg, pch="o", col=1:nlevels(tab$pop))

#pca run after removing inconsistent snps

library(SNPRelate)
vcf.fn <- "C:/Users/PB-ARNEL/Documents/Thesis2/GBS SNP files/prg_bi_repmat.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "prg2.gds", method="biallelic.only")
#it removed some snps, why? from 123004 to 122973
snpgdsSummary("prg2.gds") 

genofile <- snpgdsOpen("prg2.gds", readonly = FALSE)

#snpgdsClose(genofile)


set.seed(101010)
#pca
pca2 <- snpgdsPCA(genofile, algorithm = "randomized")

#pop and batch
pop<-substr(pca2$sample.id,1,3)
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
batch<-ifelse(pca2$sample.id %in% batch1,"B1","B2")

#combine to df
tab <- data.frame(sample.id = pca2$sample.id,
                  pop, batch,
                  EV1 = pca2$eigenvect[,1],    # the first eigenvector
                  EV2 = pca2$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

tab<-tab[order(tab$pop,tab$sample.id),]

biorep<- c("P42_E3_01","P42_E3_02", "P42_E3_03",
           "P42_F5","P42_F5_02", "P42_F5_03", 
           "P46_A6_01","P46_A6_02", 
           "P46_A10","P46_A10_02", "P46_A10_03",  
           "P46_D11","P46_D11_02", 
           "GA66_SQ2761", "GA66_SQ2718")

#P42_F5 NOT P42_F5_01 
#P46_A10 NOT P46_A10_01 
#P46_D11 NOT P46_D11_01

tab2<-tab[-which(tab$sample.id %in% biorep),]
tabrep<-tab[which(tab$sample.id %in% biorep),]

plot(x=tab2$EV1, y=tab2$EV2,
     ylim = c(-0.20,0.20), xlim = c(-0.20, 0.20),
     xlab="eigenvector 1", ylab="eigenvector 2",
     col=as.factor(tab2$pop), pch=as.numeric(substr(tab2$batch,2,2))+15)

sym<-c("o","+",
       "o","+","x",
       "o","+","x",
       "o","+","x",
       "o","+",
       "o","+")


points(x=tabrep$EV1, y=tabrep$EV2, 
       col=c(rep("green",2),rep("azure4",6), rep("brown4",7)),
       pch=sym)

repname<-c("GA66","P42_E3", "P42_F5", "P46_A10",
           "P46_A6","P46_D11")

text(x=tabrep[c(1,4,7,11,12,14),]$EV1, 
     y=tabrep[c(1,4,7,11,12,14),]$EV2,
     labels =repname,
     pos = c(1,1,1,2,1,4),
     col = c("green", rep("azure4",2), rep("brown4",3))
)

#pca run after removing 20% missing snps

library(SNPRelate)
vcf.fn <- "C:/Users/PB-ARNEL/Documents/Thesis2/GBS SNP files/prg_bi_repmat_mis.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "prg3.gds", method="biallelic.only")
#it removed some snps, why? from 80326 to 80298
snpgdsSummary("prg3.gds") 

genofile <- snpgdsOpen("prg3.gds", readonly = FALSE)

#snpgdsClose(genofile)


set.seed(888)
#pca
pca3 <- snpgdsPCA(genofile, algorithm = "randomized")

#pop and batch
pop<-substr(pca3$sample.id,1,3)
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
batch<-ifelse(pca3$sample.id %in% batch1,"B1","B2")

#combine to df
tab <- data.frame(sample.id = pca3$sample.id,
                  pop, batch,
                  EV1 = pca3$eigenvect[,1],    # the first eigenvector
                  EV2 = pca3$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

tab<-tab[order(tab$pop,tab$sample.id),]

biorep<- c("P42_E3_01","P42_E3_02", "P42_E3_03",
           "P42_F5","P42_F5_02", "P42_F5_03", 
           "P46_A6_01","P46_A6_02", 
           "P46_A10","P46_A10_02", "P46_A10_03",  
           "P46_D11","P46_D11_02", 
           "GA66_SQ2761", "GA66_SQ2718")

#P42_F5 NOT P42_F5_01 
#P46_A10 NOT P46_A10_01 
#P46_D11 NOT P46_D11_01

tab2<-tab[-which(tab$sample.id %in% biorep),]
tabrep<-tab[which(tab$sample.id %in% biorep),]

plot(x=tab2$EV1, y=tab2$EV2,
     ylim = c(-0.20,0.20), xlim = c(-0.20, 0.20),
     xlab="eigenvector 1", ylab="eigenvector 2",
     col=as.factor(tab2$pop), pch=as.numeric(substr(tab2$batch,2,2))+15)

sym<-c("o","+",
       "o","+","x",
       "o","+","x",
       "o","+","x",
       "o","+",
       "o","+")


points(x=tabrep$EV1, y=tabrep$EV2, 
       col=c(rep("green",2),rep("azure4",6), rep("brown4",7)),
       pch=sym)

repname<-c("GA66","P42_E3", "P42_F5", "P46_A10",
           "P46_A6","P46_D11")

text(x=tabrep[c(1,4,7,11,12,14),]$EV1, 
     y=tabrep[c(1,4,7,11,12,14),]$EV2,
     labels =repname,
     pos = c(3,3,3,3,4,4),
     col = c("green", rep("azure4",2), rep("brown4",3))
)
