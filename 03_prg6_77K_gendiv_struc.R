#'Set up and load the data

#'Set up and load the data
setwd("~/Thesis2/Data analysis/Data")
library(vcfR)

prg6.vcf <- read.vcfR("01_prg6_snptaxafiltered_pol.vcf.gz") 
show(prg6.vcf)


#'convert to genlight
library(poppr)
prg6.gl<-vcfR2genlight(prg6.vcf)
ploidy(prg6.gl) <- 2 #specify ploidy

#'population info
pop(prg6.gl)<-substr(indNames(prg6.gl),1,3)

#'subset per pop
prg6_pop.gl<-seppop(prg6.gl)


#'Allele frequency spectrum
#'
#'whole data set
myFreq <- glMean(prg6.gl, alleleAsUnit = TRUE)
hist(myFreq, proba=TRUE, col="blue", xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
lines(density(myFreq), lwd=3)

#'per pop
myFreq_p42 <- glMean(prg6_pop.gl$P42, alleleAsUnit = TRUE)
hist(myFreq_p42, proba=TRUE, col="blue", xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
lines(density(myFreq_p42), lwd=3)

myFreq_p46 <- glMean(prg6_pop.gl$P46, alleleAsUnit = TRUE)
hist(myFreq_p46, proba=TRUE, col="blue", xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
lines(density(myFreq_p46), lwd=3)

#'overlapping
#'grey
hist(myFreq_p42, proba=TRUE, col=rgb(0.1,0.1,0.1,0.5), xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
hist(myFreq_p46, proba=TRUE, col=rgb(0.8,0.8,0.8,0.5), add=T)
legend("topright", legend=c("P42", "P46"), fill = c(rgb(0.1,0.1,0.1,0.5), rgb(0.8,0.8,0.8,0.5)))

#'black and red
hist(myFreq_p42, proba=TRUE, col=adjustcolor("black", alpha=0.7), xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
hist(myFreq_p46, proba=TRUE, col=adjustcolor("red", alpha=0.9), add=T)
legend("topright", legend=c("P42", "P46"), 
       fill = c(adjustcolor("black", alpha=0.7), adjustcolor("red", alpha=0.9)))


#box()

#'both alleles
myFreq2 <- c(myFreq, 1-myFreq)
hist(myFreq2, proba=TRUE, col="blue", xlab="Allele frequencies",
     main="Frequency distribution of the alternative allele")
lines(density(myFreq2), lwd=3)


#'Summary statistics
#'convert to genind
#'
#install.packages('dartR')
library(dartR)
prg6.gi<-gl2gi(prg6.gl)
save(prg6.gi, file = "prg6.gi.rdata")
genstat_pop<- summary(prg6.gi) #Hobs and Hexp
save(genstat_pop, file = "prg6.gisum_het.rdata")

#'per pop
prg6_perpop.gi<-seppop(prg6.gi)
genstat_popp42<- summary(prg6_perpop.gi$P42) #Hobs and Hexp
genstat_popp46<- summary(prg6_perpop.gi$P46) #Hobs and Hexp

boxplot(genstat_pop$Hexp, genstat_popp42$Hexp, genstat_popp46$Hexp)
boxplot(genstat_pop$Hobs, genstat_popp42$Hobs, genstat_popp46$Hobs)
#'looks same for the 3


plot(genstat_popp42$Hexp,genstat_popp42$Hobs, 
     col=adjustcolor("black", alpha=0.5))
points(genstat_popp46$Hexp,genstat_popp46$Hobs, 
     col=adjustcolor("red", alpha=0.7))

plot(genstat_popp42$Hobs,genstat_popp42$Hexp, 
     col=adjustcolor("black", alpha=0.1))
points(genstat_popp46$Hobs,genstat_popp46$Hexp, 
       col=adjustcolor("red", alpha=0.1))
#'I need to figure out how to present Hexp and Hobs

library(poppr)
#prg_diversity<-poppr(prg6.gi)
#'other summary stats
#'
#'Population Structure
#'
#'tree matrix
#'
tree_bwdist <- aboot(prg6.gl, tree = "upgma", distance = bitwise.dist, sample = 500, showtree = F, cutoff = 90, quiet = F)
tree_euc <- aboot(prg6.gl, tree = "upgma", distance = function(x){bitwise.dist(x, euclidean = T)}, sample = 500, showtree = F, cutoff = 90, quiet = F)

library(RColorBrewer)
cols <- brewer.pal(n = nPop(prg6.gl), name = "Accent")

library(ape)
plot.phylo(tree_bwdist, cex = 0.6, font = 1.5, adj = 0, 
           tip.color =  cols[pop(prg6.gl)],
           use.edge.length = F)
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
legend('topleft', legend = c("P42","P46"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

plot.phylo(tree_euc, cex = 0.6, font = 1.5, adj = 0, 
           tip.color =  cols[pop(prg6.gl)],
           use.edge.length = F)
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
legend('topleft', legend = c("P42","P46"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (Euclidean)")

#' PCA
prg6.pca <- glPca(prg6.gl, nf = 10)
save(prg6.pca, file="prg6.pca.Rdata")
barplot(100*prg6.pca$eig/sum(prg6.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

prg6.pca.scores <- as.data.frame(prg6.pca$scores)
prg6.pca.scores$pop <- pop(prg6.gl)

library(ggplot2)
set.seed(999)
p <- ggplot(prg6.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p


#'amova

prg6_strata.gl<-prg6.gl
strata(prg6_strata.gl)<-data.frame(pop = pop(prg6_strata.gl))
prg6_strata.gl

poppr.amova(prg6.gl, ~pop)

