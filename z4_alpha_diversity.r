# three datasets each for individuals and localities
ps1 <- ps
ps2 <- tip_glom(ps, h=0.03)
ps3 <- tip_glom(ps, h=0.06)

psl1 <- ps.locales
psl2 <- tip_glom(ps.locales, h=0.03)
psl3 <- tip_glom(ps.locales, h=0.06)

# calculate five measures of alpha diversity and write to tables
# use rarefaction (100 reps) to estimate relative diversity
# do this for each of the three agglomerated datasets
otus1 <- c()
chao11 <- c()
shannon1 <- c()
simpson1 <- c()
phy.div1 <- c()
otus2 <- c()
chao12 <- c()
shannon2 <- c()
simpson2 <- c()
phy.div2 <- c()
otus3 <- c()
chao13 <- c()
shannon3 <- c()
simpson3 <- c()
phy.div3 <- c()
for(a in 1:100) {
	print(a)
	a.rep1 <- rarefy_even_depth(ps1, 15000)
	otus1 <- cbind(otus1, estimate_richness(a.rep1)[,1])
	chao11 <- cbind(chao11, estimate_richness(a.rep1)[,2])
	shannon1 <- cbind(shannon1, estimate_richness(a.rep1)[,6])
	simpson1 <- cbind(simpson1, estimate_richness(a.rep1)[,7])
	phy.div1 <- cbind(phy.div1, pd((otu_table(a.rep1)@.Data), ps1@phy_tree, include.root=F)$PD)
	a.rep2 <- rarefy_even_depth(ps2, 15000)
	otus2 <- cbind(otus2, estimate_richness(a.rep2)[,1])
	chao12 <- cbind(chao12, estimate_richness(a.rep2)[,2])
	shannon2 <- cbind(shannon2, estimate_richness(a.rep2)[,6])
	simpson2 <- cbind(simpson2, estimate_richness(a.rep2)[,7])
	phy.div2 <- cbind(phy.div2, pd((otu_table(a.rep2)@.Data), ps1@phy_tree, include.root=F)$PD)
	a.rep3 <- rarefy_even_depth(ps3, 15000)
	otus3 <- cbind(otus3, estimate_richness(a.rep3)[,1])
	chao13 <- cbind(chao13, estimate_richness(a.rep3)[,2])
	shannon3 <- cbind(shannon3, estimate_richness(a.rep3)[,6])
	simpson3 <- cbind(simpson3, estimate_richness(a.rep3)[,7])
	phy.div3 <- cbind(phy.div3, pd((otu_table(a.rep3)@.Data), ps1@phy_tree, include.root=F)$PD)
}
output <- data.frame(Sample.ID=ps@sam_data$Sample.ID,
			OTU_Median_Original=apply(otus1, 1, function(x) median(x)),
			OTU_SD_Original=apply(otus1, 1, function(x) sd(x)),
			CHAO1_Median_Original=apply(chao11, 1, function(x) median(x)),
			CHAO1_SD_Original=apply(chao11, 1, function(x) sd(x)),
			Shannon_Median_Original=apply(shannon1, 1, function(x) median(x)),
			Shannon_SD_Original=apply(shannon1, 1, function(x) sd(x)),
			Simpson_Median_Original=apply(simpson1, 1, function(x) median(x)),
			Simpson_SD_Original=apply(simpson1, 1, function(x) sd(x)),
			PD_Median_Original=apply(phy.div1, 1, function(x) median(x)),
			PD_SD_Original=apply(phy.div1, 1, function(x) sd(x)),
			OTU_Median_H03=apply(otus2, 1, function(x) median(x)),
			OTU_SD_H03=apply(otus2, 1, function(x) sd(x)),
			CHAO1_Median_H03=apply(chao12, 1, function(x) median(x)),
			CHAO1_SD_H03=apply(chao12, 1, function(x) sd(x)),
			Shannon_Median_H03=apply(shannon2, 1, function(x) median(x)),
			Shannon_SD_H03=apply(shannon2, 1, function(x) sd(x)),
			Simpson_Median_H03=apply(simpson2, 1, function(x) median(x)),
			Simpson_SD_H03=apply(simpson2, 1, function(x) sd(x)),
			PD_Median_H03=apply(phy.div2, 1, function(x) median(x)),
			PD_SD_H03=apply(phy.div2, 1, function(x) sd(x)),
			OTU_Median_H06=apply(otus3, 1, function(x) median(x)),
			OTU_SD_H06=apply(otus3, 1, function(x) sd(x)),
			CHAO1_Median_H06=apply(chao13, 1, function(x) median(x)),
			CHAO1_SD_H06=apply(chao13, 1, function(x) sd(x)),
			Shannon_Median_H06=apply(shannon3, 1, function(x) median(x)),
			Shannon_SD_H06=apply(shannon3, 1, function(x) sd(x)),
			Simpson_Median_H06=apply(simpson3, 1, function(x) median(x)),
			Simpson_SD_H06=apply(simpson3, 1, function(x) sd(x)),
			PD_Median_H06=apply(phy.div3, 1, function(x) median(x)),
			PD_SD_H06=apply(phy.div3, 1, function(x) sd(x)))

write.table(output, file="alpha_diversity_individuals.txt", sep="\t",quote=F, row.names=F)

otus1 <- c()
chao11 <- c()
shannon1 <- c()
simpson1 <- c()
phy.div1 <- c()
otus2 <- c()
chao12 <- c()
shannon2 <- c()
simpson2 <- c()
phy.div2 <- c()
otus3 <- c()
chao13 <- c()
shannon3 <- c()
simpson3 <- c()
phy.div3 <- c()
for(a in 1:100) {
	print(a)
	a.rep1 <- rarefy_even_depth(psl1, 65000)
	otus1 <- cbind(otus1, estimate_richness(a.rep1)[,1])
	chao11 <- cbind(chao11, estimate_richness(a.rep1)[,2])
	shannon1 <- cbind(shannon1, estimate_richness(a.rep1)[,6])
	simpson1 <- cbind(simpson1, estimate_richness(a.rep1)[,7])
	phy.div1 <- cbind(phy.div1, pd((otu_table(a.rep1)@.Data), ps1@phy_tree, include.root=F)$PD)
	a.rep2 <- rarefy_even_depth(psl2, 65000)
	otus2 <- cbind(otus2, estimate_richness(a.rep2)[,1])
	chao12 <- cbind(chao12, estimate_richness(a.rep2)[,2])
	shannon2 <- cbind(shannon2, estimate_richness(a.rep2)[,6])
	simpson2 <- cbind(simpson2, estimate_richness(a.rep2)[,7])
	phy.div2 <- cbind(phy.div2, pd((otu_table(a.rep2)@.Data), ps1@phy_tree, include.root=F)$PD)
	a.rep3 <- rarefy_even_depth(psl3, 65000)
	otus3 <- cbind(otus3, estimate_richness(a.rep3)[,1])
	chao13 <- cbind(chao13, estimate_richness(a.rep3)[,2])
	shannon3 <- cbind(shannon3, estimate_richness(a.rep3)[,6])
	simpson3 <- cbind(simpson3, estimate_richness(a.rep3)[,7])
	phy.div3 <- cbind(phy.div3, pd((otu_table(a.rep3)@.Data), ps1@phy_tree, include.root=F)$PD)
}
output <- data.frame(Sample.ID=psl1@sam_data$Species.Location,
			OTU_Median_Original=apply(otus1, 1, function(x) median(x)),
			OTU_SD_Original=apply(otus1, 1, function(x) sd(x)),
			CHAO1_Median_Original=apply(chao11, 1, function(x) median(x)),
			CHAO1_SD_Original=apply(chao11, 1, function(x) sd(x)),
			Shannon_Median_Original=apply(shannon1, 1, function(x) median(x)),
			Shannon_SD_Original=apply(shannon1, 1, function(x) sd(x)),
			Simpson_Median_Original=apply(simpson1, 1, function(x) median(x)),
			Simpson_SD_Original=apply(simpson1, 1, function(x) sd(x)),
			PD_Median_Original=apply(phy.div1, 1, function(x) median(x)),
			PD_SD_Original=apply(phy.div1, 1, function(x) sd(x)),
			OTU_Median_H03=apply(otus2, 1, function(x) median(x)),
			OTU_SD_H03=apply(otus2, 1, function(x) sd(x)),
			CHAO1_Median_H03=apply(chao12, 1, function(x) median(x)),
			CHAO1_SD_H03=apply(chao12, 1, function(x) sd(x)),
			Shannon_Median_H03=apply(shannon2, 1, function(x) median(x)),
			Shannon_SD_H03=apply(shannon2, 1, function(x) sd(x)),
			Simpson_Median_H03=apply(simpson2, 1, function(x) median(x)),
			Simpson_SD_H03=apply(simpson2, 1, function(x) sd(x)),
			PD_Median_H03=apply(phy.div2, 1, function(x) median(x)),
			PD_SD_H03=apply(phy.div2, 1, function(x) sd(x)),
			OTU_Median_H06=apply(otus3, 1, function(x) median(x)),
			OTU_SD_H06=apply(otus3, 1, function(x) sd(x)),
			CHAO1_Median_H06=apply(chao13, 1, function(x) median(x)),
			CHAO1_SD_H06=apply(chao13, 1, function(x) sd(x)),
			Shannon_Median_H06=apply(shannon3, 1, function(x) median(x)),
			Shannon_SD_H06=apply(shannon3, 1, function(x) sd(x)),
			Simpson_Median_H06=apply(simpson3, 1, function(x) median(x)),
			Simpson_SD_H06=apply(simpson3, 1, function(x) sd(x)),
			PD_Median_H06=apply(phy.div3, 1, function(x) median(x)),
			PD_SD_H06=apply(phy.div3, 1, function(x) sd(x)))

write.table(output, file="alpha_diversity_localities.txt", sep="\t",quote=F, row.names=F)
