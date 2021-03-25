#'VCF additional filtering
#'scheme:
#'
#'VCF6 =>[ DP = c(5,150) ]=> VCF7
#'VCF7 is for the selection signature analysis
#'
#'
#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

prg6.vcf <- read.vcfR("01_prg6_snptaxafiltered_pol.vcf.gz") 
show(prg6.vcf)



#'select dp of at least 5
dp6 <- extract.gt(prg6.vcf,  element = "DP", as.numeric = TRUE)
mindp <- apply(dp6, MARGIN = 1, function(x){ min(x, na.rm = T) } )
prg_minDP5.vcf<- prg6.vcf[mindp>=5, ]
show(prg_minDP5.vcf) #30K

#'select dp of at most 150
dp6.2 <- extract.gt(prg_minDP5.vcf,  element = "DP", as.numeric = TRUE)
maxdp <- apply(dp6.2, MARGIN = 1, function(x){ max(x, na.rm = T) } )
prg_DP_5to150.vcf<- prg_minDP5.vcf[maxdp>=150, ]
show(prg_DP_5to150.vcf) #17K

prg7.vcf<-prg_DP_5to150.vcf

write.vcf(prg7.vcf, file="04_prg7_snptaxafiltered_pol_5dp150.vcf.gz", mask=TRUE)
