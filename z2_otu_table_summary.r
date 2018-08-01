# load package for processing data
library(dada2); packageVersion("dada2")
library(phangorn)
library(DECIPHER)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(reshape)
library(RColorBrewer)
library(plyr)
library(picante)
cols <- brewer.pal(8,"Set2")
cols <- c(cols[2], cols[1], cols[3:8])

load("microbe_workflow6.Rdata")

# define output tables
i_kingdom <- c()
reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Kingdom"]))
reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Kingdom"]
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
	i_kingdom <- cbind(i_kingdom, a.rep)
}
colnames(i_kingdom) <- reps
inds_names <- rownames(i_kingdom)

i_phylum <- c()
reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Phylum"]))
reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Phylum"]
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
	i_phylum <- cbind(i_phylum, a.rep)
}
colnames(i_phylum) <- reps

i_class <- c()
reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Class"]))
reps <- na.omit(reps)
reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Class"]
reps_matrix <- na.omit(reps_matrix)
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
	i_class <- cbind(i_class, a.rep)
}
colnames(i_class) <- reps

i_genus <- c()
reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Genus"]))
reps <- na.omit(reps)
reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Genus"]
reps_matrix <- na.omit(reps_matrix)
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
	i_genus <- cbind(i_genus, a.rep)
}
colnames(i_genus) <- reps

# kingdom output
i_tmp <- c()
for(a in 1:nrow(i_kingdom)) {
	a.rep <- i_kingdom[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	i_tmp <- rbind(i_tmp, a.rep)
}
i_kingdom <- i_tmp
i_kingdom_output <- cbind(colnames(i_kingdom), apply(i_kingdom, 2, mean), apply(i_kingdom, 2, min), apply(i_kingdom, 2, max))
i_kingdom_output <- data.frame(Kingdom=as.vector(sapply(strsplit(i_kingdom_output[,1], "__"), "[[", 2)), Mean=(as.numeric(i_kingdom_output[,2])), Min=(as.numeric(i_kingdom_output[,3])), Max=(as.numeric(i_kingdom_output[,4])))
					 
# phylum output					 
i_tmp <- c()
for(a in 1:nrow(i_phylum)) {
	a.rep <- i_phylum[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	i_tmp <- rbind(i_tmp, a.rep)
}
i_phylum <- i_tmp
i_phylum_output <- cbind(colnames(i_phylum), apply(i_phylum, 2, mean), apply(i_phylum, 2, min), apply(i_phylum, 2, max))
i_phylum_output <- data.frame(Phylum=as.vector(sapply(strsplit(i_phylum_output[,1], "__"), "[[", 2)), Mean=(as.numeric(i_phylum_output[,2])),Min=(as.numeric(i_phylum_output[,3])), Max=(as.numeric(i_phylum_output[,4])))

# class output
i_tmp <- c()
for(a in 1:nrow(i_class)) {
	a.rep <- i_class[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	i_tmp <- rbind(i_tmp, a.rep)
}
i_class <- i_tmp
i_class_output <- cbind(colnames(i_class), apply(i_class, 2, mean), apply(i_class, 2, min), apply(i_class, 2, max))
i_class_output <- data.frame(Class=substr(i_class_output[,1], 4, nchar(i_class_output[,1])), Mean=(as.numeric(i_class_output[,2])),Min=(as.numeric(i_class_output[,3])), Max=(as.numeric(i_class_output[,4])))

# genus output
i_tmp <- c()
for(a in 1:nrow(i_genus)) {
	a.rep <- i_genus[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	i_tmp <- rbind(i_tmp, a.rep)
}
i_genus <- i_tmp
i_genus_output <- cbind(colnames(i_genus), apply(i_genus, 2, mean), apply(i_genus, 2, min), apply(i_genus, 2, max))
i_genus_output <- data.frame(Genus=substr(i_genus_output[,1], 4, nchar(i_genus_output[,1])), Mean=(as.numeric(i_genus_output[,2])),Min=(as.numeric(i_genus_output[,3])), Max=(as.numeric(i_genus_output[,4])))

writeLines(paste("Summary across individuals"))
i_kingdom_output
i_phylum_output
i_class_output
i_genus_output
write.table(i_kingdom_output, file="i_kingdom_output.txt", sep="\t", quote=F, row.names=F)
write.table(i_phylum_output, file="i_phylum_output.txt", sep="\t", quote=F, row.names=F)
write.table(i_class_output, file="i_class_output.txt", sep="\t", quote=F, row.names=F)
write.table(i_genus_output, file="i_genus_output.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(inds_names, i_phylum), file="i_phylum_output2.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(inds_names, i_class), file="i_class_output2.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(inds_names, i_genus), file="i_genus_output2.txt", sep="\t", quote=F, row.names=F)


l_phylum <- c()
reps <- as.vector(unique(ps.locales@tax_table[ ,colnames(ps.locales@tax_table) == "Phylum"]))
reps_matrix <- ps.locales@tax_table[,colnames(ps.locales@tax_table) == "Phylum"]
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps.locales@otu_table[ ,match(a.rep, colnames(ps.locales@otu_table))], 1, sum)
	l_phylum <- cbind(l_phylum, a.rep)
}
colnames(l_phylum) <- reps
locs_names <- rownames(l_phylum)

l_class <- c()
reps <- as.vector(unique(ps.locales@tax_table[ ,colnames(ps.locales@tax_table) == "Class"]))
reps <- na.omit(reps)
reps_matrix <- ps.locales@tax_table[,colnames(ps.locales@tax_table) == "Class"]
reps_matrix <- na.omit(reps_matrix)
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps.locales@otu_table[ ,match(a.rep, colnames(ps.locales@otu_table))], 1, sum)
	l_class <- cbind(l_class, a.rep)
}
colnames(l_class) <- reps

l_genus <- c()
reps <- as.vector(unique(ps.locales@tax_table[ ,colnames(ps.locales@tax_table) == "Genus"]))
reps <- na.omit(reps)
reps_matrix <- ps.locales@tax_table[,colnames(ps.locales@tax_table) == "Genus"]
reps_matrix <- na.omit(reps_matrix)
for(a in 1:length(reps)) {
	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
	a.rep <- apply(ps.locales@otu_table[ ,match(a.rep, colnames(ps.locales@otu_table))], 1, sum)
	l_genus <- cbind(l_genus, a.rep)
}
colnames(l_genus) <- reps

# phylum output					 
l_tmp <- c()
for(a in 1:nrow(l_phylum)) {
	a.rep <- l_phylum[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	l_tmp <- rbind(l_tmp, a.rep)
}
l_phylum <- l_tmp
l_phylum_output <- cbind(colnames(l_phylum), apply(l_phylum, 2, mean), apply(l_phylum, 2, min), apply(l_phylum, 2, max))
l_phylum_output <- data.frame(Phylum=as.vector(sapply(strsplit(l_phylum_output[,1], "__"), "[[", 2)), Mean=(as.numeric(l_phylum_output[,2])),Min=(as.numeric(l_phylum_output[,3])), Max=(as.numeric(l_phylum_output[,4])))

# class output
l_tmp <- c()
for(a in 1:nrow(l_class)) {
	a.rep <- l_class[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	l_tmp <- rbind(l_tmp, a.rep)
}
l_class <- l_tmp
l_class_output <- cbind(colnames(l_class), apply(l_class, 2, mean), apply(l_class, 2, min), apply(l_class, 2, max))
l_class_output <- data.frame(Class=substr(l_class_output[,1], 4, nchar(l_class_output[,1])), Mean=(as.numeric(l_class_output[,2])),Min=(as.numeric(l_class_output[,3])), Max=(as.numeric(l_class_output[,4])))

# genus output
l_tmp <- c()
for(a in 1:nrow(l_genus)) {
	a.rep <- l_genus[a,]
	a.rep <- a.rep / sum(a.rep) * 100
	l_tmp <- rbind(l_tmp, a.rep)
}
l_genus <- l_tmp
l_genus_output <- cbind(colnames(l_genus), apply(l_genus, 2, mean), apply(l_genus, 2, min), apply(l_genus, 2, max))
l_genus_output <- data.frame(Genus=substr(l_genus_output[,1], 4, nchar(l_genus_output[,1])), Mean=(as.numeric(l_genus_output[,2])),Min=(as.numeric(l_genus_output[,3])), Max=(as.numeric(l_genus_output[,4])))

writeLines(paste("Summary across localities"))
l_phylum_output
l_class_output
l_genus_output
write.table(l_phylum_output, file="l_phylum_output.txt", sep="\t", quote=F, row.names=F)
write.table(l_class_output, file="l_class_output.txt", sep="\t", quote=F, row.names=F)
write.table(l_genus_output, file="l_genus_output.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(locs_names, l_phylum), file="l_phylum_output2.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(locs_names, l_class), file="l_class_output2.txt", sep="\t", quote=F, row.names=F)
write.table(cbind(locs_names, l_genus), file="l_genus_output2.txt", sep="\t", quote=F, row.names=F)