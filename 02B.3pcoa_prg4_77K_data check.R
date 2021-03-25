#'Data quality check of the initial vcf file
#'violin plots of read depth
#'heat map of missingness and allele 
#'PCoA looking at (non)clustering of pop/gen, batch and checks 

#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

#'load data with taxa not filtered (i.e. reps 2 and 3 are retained)
prg4.vcf <- read.vcfR("01_prg4_snpfiltered.vcf.gz") 
show(prg4.vcf) #77263 snps; 190 samples

#'monomorphic snps filtered out
prg4v2.vcf<-prg4.vcf[is.polymorphic(prg4.vcf, na.omit = TRUE),]
show(prg4v2.vcf)#removed 2 snps




#'Prepare data for plotting
#'
library(poppr)
prg4.gl<- vcfR2genlight(prg4v2.vcf)
ploidy(prg4.gl)<-2

#'PCoA by distance first
#'
#'bitwise distance
prg4.bwdist<-bitwise.dist(prg4.gl) #defaults
prg4.bweuc<-bitwise.dist(prg4.gl, euclidean = TRUE) #euclidean distance

mds<-cmdscale(prg4.bweuc)

#'subset batch and pop
#'
#'population info
pop<-substr(row.names(mds),1,3)

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

batch<-ifelse(row.names(mds) %in% batch1,"B1","B2")

#combine to df
tab <- data.frame(pop, batch,
                  C1 = mds[,1],
                  C2 = mds[,2],
                  stringsAsFactors = FALSE)
head(tab)

tab<-tab[order(tab$pop,row.names(tab)),]

biorep<- c("P42_E3_01","P42_E3_02", "P42_E3_03",
           "P42_F5","P42_F5_02", "P42_F5_03", 
           "P46_A6_01","P46_A6_02", 
           "P46_A10","P46_A10_02", "P46_A10_03",  
           "P46_D11","P46_D11_02", 
           "GA66_SQ2761", "GA66_SQ2718")

#P42_F5 NOT P42_F5_01 
#P46_A10 NOT P46_A10_01 
#P46_D11 NOT P46_D11_01

tab2<-tab[-which(row.names(tab) %in% biorep),]
tabrep<-tab[which(row.names(tab) %in% biorep),]

plot(x=tab2$C1, y=tab2$C2,cex=2,
     #xlim = c(-0.10, 0.10), ylim = c(-0.10, 0.10),
     xlim = c(-80,80), ylim = c(-80, 80),
     xlab="PCo 1", ylab="PCo 2",
     col=as.factor(tab2$pop), pch=as.numeric(substr(tab2$batch,2,2))+15)

sym<-c("o","+",
       "o","+","x",
       "o","+","x",
       "o","+","x",
       "o","+",
       "o","+")


points(x=tabrep$C1, y=tabrep$C2, cex=2,
       col=c(rep("green",2),rep("azure4",6), rep("brown4",7)),
       pch=sym)

repname<-c("GA66","P42_E3", "P42_F5", "P46_A10",
           "P46_A6","P46_D11")

text(x=tabrep[c(1,4,7,11,12,14),]$C1, 
     y=tabrep[c(1,4,7,11,12,14),]$C2,
     labels =repname,
     pos = c(1,4,2,3,3,1),
     col = c("green", rep("azure4",2), rep("brown4",3))
)

#optional legend
leg<-c("P42 batch1", "P46 batch1",
       "P42 batch2", "P46 batch2",
       "Reps")

legend("bottomright", legend=leg, pch="o", col=1:nlevels(tab$pop))

#'PCoA based on proportion of shared allele distance
#'still need to figure this out

library(pegas)
prg1.loci<-vcfR2loci(prg.vcf) 
prg1.asd<-dist.asd(prg1.loci, na.) #NA's
save(prg1.loci, file = "~/Thesis2/Data analysis/Data/02A.3pegas_prg1_loci.Rdata")

#'vcf to genind
prg1.gi<-vcfR2genind(prg.vcf)
save(prg1.gi, file = "~/Thesis2/Data analysis/Data/02A.3adgnt_prg1_genind.Rdata")
