#'VCF filtering
#'scheme:
#'
  #'snps filter first
  #'VCF0 =>[>5% maf; <50% mis]=> VCF1 (this is "arnel_recode.vcf")
  #'VCF1 =>[mean DP>=5; biallelic]=> VCF2
  #'VCF2 =>[<20% mis; No "N" base calls]=> VCF3
  #'VCF3 =>[consitent across bioreps]=> VCF4
  #'VCF4 =>[taxa filter: remove reps 2 and 3; <20% mis]=> VCF5
  #'VCF5 =>[remove monomorphic snps]=> VCF6
  #'VCF6 is for the pop genetic analysis
  #'
#'
#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

prg.vcf <- read.vcfR("arnel_out.recode.vcf") #189K snps
show(prg.vcf)#189K

#'  #'VCF1 =>[mean DP>=5; biallelic]=> VCF2
#'  
#'  extract depth
dp1 <- extract.gt(prg.vcf,  element = "DP", as.numeric = TRUE)
#'select dp of at least 5
prg_dp5.vcf<- prg.vcf[rowMeans(dp1)>=5, ]
show(prg_dp5.vcf)#142K

#'retain biallelic snps only
sum(is.biallelic(prg_dp5.vcf)) #how many biallelic
prg_bi.vcf<-prg_dp5.vcf[is.biallelic(prg_dp5.vcf),]
show(prg_bi.vcf)#136K
#same as csv file from mingshu

#'VCF2 =>[<20% mis; No "N" base calls]=> VCF3
#'
#'Missing data first
#'extract genotype matrix
GT1 <- extract.gt(prg_bi.vcf,  element = "GT") #extract GT matrix 
sum(is.na(GT1)) #how many missing

#'Replace NAs with NA, since it is masked in vcfR format
#'[,-1] do not change the first col, "FORMAT"
prg_bi2.vcf<-prg_bi.vcf
prg_bi2.vcf@gt[,-1][is.na(GT1) == TRUE] <- NA  
show(prg_bi2.vcf)#stil 136K but 10.53% missing data

#'create rowwise count/percent of missing/NAs
myMiss <- apply(GT1, MARGIN = 1, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / ncol(GT1)

#'include SNPs with less than 20% missing 
prg_mis.vcf<- prg_bi2.vcf[myMiss < 0.2, ]
show(prg_mis.vcf)#107K; 4.33% missing

#'N base calls
#'check for "N" bases in the "fix" part of vcf file
prg_mis.vcf@fix[prg_mis.vcf@fix[,5]=="N",]

N_snp<-prg_mis.vcf@fix[prg_mis.vcf@fix[,5]=="N",]#subset
N_snp[,3] #snp names

GT2 <- extract.gt(prg_mis.vcf,  element = "GT") #update GT matrix
prg_N.vcf<-prg_mis.vcf[-which(row.names(GT2) %in% N_snp[,3]),]
show(prg_N.vcf)#107K; 4.33% missing; 30 snps removed

#'VCF3 =>[consitent across bioreps]=> VCF4
#'
#'filter inconsistent snps among reps
#'
#'extract new/updated GT
GT3 <- extract.gt(prg_N.vcf,  element = "GT")

#' convert geno into 0,1,2
GT3[GT3=="0/0"]<-0
GT3[GT3=="0/1"]<-1
GT3[GT3=="1/0"]<-1
GT3[GT3=="1/1"]<-2

#'Subset bioreps 
#'check for inconsitencies
#' 
#' P42_E3
P42_E3<-GT3[,grepl("P42_E3", colnames(GT3))]
P42_E3.1<-P42_E3[complete.cases(P42_E3),]
P42_E3mis12<-P42_E3.1[rownames(P42_E3.1[P42_E3.1[,1] != P42_E3.1[,2],]),]
P42_E3mis13<-P42_E3.1[rownames(P42_E3.1[P42_E3.1[,1] != P42_E3.1[,3],]),]
P42_E3mis23<-P42_E3.1[rownames(P42_E3.1[P42_E3.1[,2] != P42_E3.1[,3],]),]

#' P42_F5
P42_F5<-GT3[,grepl("P42_F5", colnames(GT3))]
P42_F5.1<-P42_F5[complete.cases(P42_F5),]
P42_F5mis12<-P42_F5.1[rownames(P42_F5.1[P42_F5.1[,1] != P42_F5.1[,2],]),]
P42_F5mis13<-P42_F5.1[rownames(P42_F5.1[P42_F5.1[,1] != P42_F5.1[,3],]),]
P42_F5mis23<-P42_F5.1[rownames(P42_F5.1[P42_F5.1[,2] != P42_F5.1[,3],]),]

#' P46_A6
P46_A6<-GT3[,grepl("P46_A6", colnames(GT3))]
P46_A6.1<-P46_A6[complete.cases(P46_A6),]
P46_A6mis12<-P46_A6.1[rownames(P46_A6.1[P46_A6.1[,1] != P46_A6.1[,2],]),]

#' P46_A10
P46_A10<-GT3[,grepl("P46_A10", colnames(GT3))]
P46_A10.1<-P46_A10[complete.cases(P46_A10),]
P46_A10mis12<-P46_A10.1[rownames(P46_A10.1[P46_A10.1[,1] != P46_A10.1[,2],]),]
P46_A10mis13<-P46_A10.1[rownames(P46_A10.1[P46_A10.1[,1] != P46_A10.1[,3],]),]
P46_A10mis23<-P46_A10.1[rownames(P46_A10.1[P46_A10.1[,2] != P46_A10.1[,3],]),]

#' P46_D11
P46_D11<-GT3[,grepl("P46_D11", colnames(GT3))]
P46_D11.1<-P46_D11[complete.cases(P46_D11),]
P46_D11mis12<-P46_D11.1[rownames(P46_D11.1[P46_D11.1[,1] != P46_D11.1[,2],]),]

#' GA66
GA66<-GT3[,grepl("GA66", colnames(GT3))]
GA66.1<-GA66[complete.cases(GA66),]
GA66mis12<-GA66.1[rownames(GA66.1[GA66.1[,1] != GA66.1[,2],]),]

#'check subsets
head(P42_E3)
head(P46_A10)
head(GA66)

#'inconsistent snps across reps
mis<- unique(c(rownames(P42_E3mis12), rownames(P42_E3mis13), rownames(P42_E3mis23),
               rownames(P42_F5mis12), rownames(P42_F5mis13), rownames(P42_F5mis23),
               rownames(P46_A6mis12),
               rownames(P46_A10mis12), rownames(P46_A10mis13), rownames(P46_A10mis23),
               rownames(P46_D11mis12),
               rownames(GA66mis12)
))

#remove inconsistent snps
prg_repmis.vcf<-prg_N.vcf[which(row.names(GT3) %in% mis),]
prg_repmat.vcf<-prg_N.vcf[-which(row.names(GT3) %in% mis),]

show(prg_repmis.vcf) #29K snps removed
show(prg_repmat.vcf) #77k snps REMAINING!

write.vcf(prg_repmat.vcf, file="01_prg4_snpfiltered.vcf.gz", mask=TRUE)
save(prg_bi.vcf,file = "~/Thesis2/Data analysis/Results/01_VCF2.RData")
save(prg_N.vcf,file = "~/Thesis2/Data analysis/Results/01_VCF3.RData")

#' 
#' Taxa filter
#' VCF4 =>[remove reps 2 and 3; <20% mis]=> VCF5
#' 
#' Remove biological reps 2 and 3, as well as GA66


biorep2<- c("P42_E3_02", "P42_E3_03",
            "P42_F5_02", "P42_F5_03", 
            "P46_A6_02", 
            "P46_A10_02", "P46_A10_03",  
            "P46_D11_02", 
            "GA66_SQ2761", "GA66_SQ2718")


prg_repmat2.vcf<-prg_repmat.vcf #new vcfr file
GT4 <- extract.gt(prg_repmat2.vcf,  element = "GT") #update GT matrix

#'Select col/individuals not in biorep; include 1st col="FORMAT"
prg_repmat2.vcf@gt<-prg_repmat2.vcf@gt[,c("FORMAT",colnames(GT4[,-which(colnames(GT4) %in% biorep2)]))]
show(prg_repmat2.vcf)#77K
colnames(prg_repmat2.vcf@gt)
head(prg_repmat2.vcf@gt[1:5,1:5])

#'Remove taxa with 20% missing
#'
#'Update mymiss
GT5 <- extract.gt(prg_repmat2.vcf,  element = "GT") #update GT matrix

myMiss <- apply(GT5, MARGIN = 2, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / nrow(GT5)

#'Include individual with less than 20% missing; TRUE=1st col="FORMAT"
prg_mistaxa.vcf<-prg_repmat2.vcf
prg_mistaxa.vcf@gt<- prg_mistaxa.vcf@gt[, c(TRUE, myMiss < 0.2 )] 
show(prg_mistaxa.vcf)#no change

write.vcf(prg_mistaxa.vcf, file="01_prg5_snptaxafiltered.vcf.gz", mask=TRUE)


#'VCF5 =>[remove monomorphic snps]=> VCF6
#'remove monomorphic snps
prg_pol.vcf<-prg_mistaxa.vcf[is.polymorphic(prg_mistaxa.vcf, na.omit = TRUE),]
show(prg_pol.vcf)#removed 2 snps

write.vcf(prg_pol.vcf, file="01_prg6_snptaxafiltered_pol.vcf.gz", mask=TRUE)

#'Data ready for genetic diversity and population structure analysis
#'
#'END