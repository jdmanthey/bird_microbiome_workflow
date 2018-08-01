# load package for processing data
library(RColorBrewer)
cols <- brewer.pal(12,"Paired")
cols <- c(cols[1:4], cols[7:12])

i_phylum <- read.table("~/Desktop/microbe/figures_new/i_phylum_output.txt", sep="\t", header=T)
i_phylum2 <- read.table("~/Desktop/microbe/figures_new/i_phylum_output2.txt", sep="\t", header=T)
l_phylum <- read.table("~/Desktop/microbe/figures_new/l_phylum_output.txt", sep="\t", header=T)
l_phylum2 <- read.table("~/Desktop/microbe/figures_new/l_phylum_output2.txt", sep="\t", header=T)

common_bacteria <- as.character(i_phylum[i_phylum[,3] > 0, 1])
test <- paste("p__", common_bacteria, sep="")
common_abund <- i_phylum2[ ,colnames(i_phylum2) %in% test]
less_abund <- i_phylum2[ ,!(colnames(i_phylum2) %in% test)]
less_abund <- less_abund[,2:ncol(less_abund)]
other <- apply(less_abund, 1, sum)
total <- cbind(common_abund, other)

l_phylum2 <- rbind(l_phylum2[11:19,], l_phylum2[1:10,])
common_abund <- l_phylum2[ ,colnames(l_phylum2) %in% test]
less_abund <- l_phylum2[ ,!(colnames(l_phylum2) %in% test)]
less_abund <- less_abund[,2:ncol(less_abund)]
other <- apply(less_abund, 1, sum)
total2 <- cbind(common_abund, other)

par(mar=c(0.25,0,0.4,0))
layout(matrix(c(1,1,1,2, 1,1,1,2), nrow=2, ncol=4, byrow=T))
barplot(t(as.matrix(total)), col=cols, yaxt="n")
barplot(t(as.matrix(total2)), col=cols, yaxt="n")