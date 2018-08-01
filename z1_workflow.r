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

# set the path to where the fastq files are
path <- "~/Desktop/microbe/samples"
list.files(path)

# sort the order of the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

# extract sample names from files
sample.names <- paste(sapply(strsplit(fnFs, "_"), `[`, 1), sapply(strsplit(fnFs, "_"), `[`, 2), sapply(strsplit(fnFs, "_"), `[`, 3), sapply(strsplit(fnFs, "_"), `[`, 4), sapply(strsplit(fnFs, "_"), `[`, 5), sep="_")
sample.names

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# look at quality of sequences (test some different numbers below)
plotQualityProfile(fnFs[1:20])
plotQualityProfile(fnRs[1:20])
#save plots as examples


# file paths for putting the filtered reads
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


# filter and trim the samples
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)


# learn the error rates for the sequencing
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


# look at plots of errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# dereplicate all reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# rename the dereplicated reads files
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# run the main file to call all of the unique sequences
dadaFs <- dada(derepFs, err=errF, pool=T,multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=T,multithread=TRUE)

save.image(file="microbe_workflow1.Rdata")
# look at the output for a few of the samples
dadaFs
dadaRs


# merge the forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# keep all mergers with length near the mode (253 + or - 4)
seqtab2 <- seqtab
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(249,257)]


# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # percent of sequences not chimeric


# summarize the filtering
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
write.table(track,"tracked_filtering.txt",sep="\t", quote=F)

# assign taxa
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/microbe/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
unname(head(taxa))

save.image(file="microbe_workflow2.Rdata")

# create phylogenetic tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs))
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=T, optGamma=T, rearrangement="stochastic", control=pml.control(trace=0))

save.image(file="microbe_workflow3.Rdata")





# Make a data.frame holding the sample data
samples.out <- rownames(seqtab.nochim)
samp.number <- sapply(strsplit(samples.out, "_"), `[`, 1)
sample.id <- sapply(strsplit(samples.out, "_"), `[`, 2)
genus.id <- sapply(strsplit(samples.out, "_"), `[`, 3)
species.id <- sapply(strsplit(samples.out, "_"), `[`, 4)
location <- sapply(strsplit(samples.out, "_"), `[`, 5)
species.location <- paste(species.id, location, sep=".")
micro.df <- data.frame(Sample.number=samp.number, Sample.ID=sample.id, Genus.ID=genus.id,
						Species.ID=species.id, Location=location, Species.Location=species.location)
rownames(micro.df) <- samples.out


# construct a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),
		sample_data(micro.df),
		tax_table(taxa),
		phy_tree(fitGTR$tree))
ps
ps <- subset_taxa(ps, Phylum != "p__Cyanobacteria")
ps.locales <- merge_samples(ps, ps@sam_data$Species.Location)
ps.locales@sam_data$Genus.ID <- as.factor(c(rep("Certhia", 10), rep("Sitta", 9)))
ps.locales@sam_data$Species.ID <- as.factor(c(rep("Certhia.albescens", 3), rep("Certhia.americana", 7), rep("Sitta.carolinensis", 9)))
ps.locales@sam_data$Location <- as.factor(c("Chiricahua", "Huachuca", "SantaRita", "MogRimE", "MogRimW1", "MogRimW2", "Pinal", "Pinaleno", "Prescott", "SantaCatalina", "Chiricahua", "Huachuca", "MogRimE", "MogRimW", "Pinal", "Pinaleno", "Prescott", "SantaCatalina","SantaRita"))
ps.locales@sam_data$Species.Location <- as.factor(paste(ps.locales@sam_data$Species.ID, ps.locales@sam_data$Location, sep="."))
ps.locales <- subset_taxa(ps.locales, Phylum != "p__Cyanobacteria")
save.image(file="microbe_workflow4.Rdata")


# summarize counts for each sample
median(sample_sums(ps))
sd(sample_sums(ps))
min(sample_sums(ps))
max(sample_sums(ps))
median(sample_sums(ps.locales))
sd(sample_sums(ps.locales))
min(sample_sums(ps.locales))
max(sample_sums(ps.locales))

# summarize observed taxa
sort(get_taxa_unique(ps, "Phylum"))
sort(get_taxa_unique(ps, "Class"))
sort(get_taxa_unique(ps, "Order"))

# run z_alpha_diversity.r



#plot beta diversity individuals
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
wu.df <- data.frame(
	x = as.vector(ord.nmds.bray$points[,1]),
	y = as.vector(ord.nmds.bray$points[,2]),
	Species.ID = species.id,
	Species.Location = species.location)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps, ord.nmds.bray, color="Species.ID", shape="Species.ID", title="Bray NMDS") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)

#plot beta diversity locales
ord.nmds.bray <- ordinate(ps.locales, method="NMDS", distance="bray")
wu.df <- data.frame(
	x = as.vector(ord.nmds.bray$points[,1]),
	y = as.vector(ord.nmds.bray$points[,2]),
	Species.ID = ps.locales@sam_data$Species.ID)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps.locales, ord.nmds.bray, color="Species.ID", shape="Species.ID", title="Bray NMDS") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)


# weighted unifrac plot individuals
ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)
wu.df <- data.frame(
	x = ord.weighted.unifrac$vectors[,1],
	y = ord.weighted.unifrac$vectors[,2],
	Species.ID = species.id,
	Species.Location = species.location)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps, ord.weighted.unifrac, color="Species.ID", shape="Species.ID",title="Weighted Unifrac PCoA") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)

# weighted unifrac plot localities
ord.weighted.unifrac <- ordinate(rarefy_even_depth(ps.locales, 65000), "PCoA", "unifrac", weighted=T)
wu.df <- data.frame(
	x = ord.weighted.unifrac$vectors[,1],
	y = ord.weighted.unifrac$vectors[,2],
	Species.ID = ps.locales@sam_data$Species.ID)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps.locales, ord.weighted.unifrac, color="Species.ID", shape="Species.ID",title="Weighted Unifrac PCoA") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)


# unweighted unifrac plot individuals
ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)
wu.df <- data.frame(
	x = ord.unweighted.unifrac$vectors[,1],
	y = ord.unweighted.unifrac$vectors[,2],
	Species.ID = species.id,
	Species.Location = species.location)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps, ord.unweighted.unifrac, color="Species.ID", shape="Species.ID",title="Unweighted Unifrac PCoA") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)

# unweighted unifrac plot localities
ord.unweighted.unifrac <- ordinate(rarefy_even_depth(ps.locales, 65000), "PCoA", "unifrac", weighted=F)
wu.df <- data.frame(
	x = ord.unweighted.unifrac$vectors[,1],
	y = ord.unweighted.unifrac$vectors[,2],
	Species.ID = ps.locales@sam_data$Species.ID)
chulls <- ddply(wu.df, .(Species.ID), function(df) df[chull(df$x, df$y), ])
plot_ordination(ps.locales, ord.unweighted.unifrac, color="Species.ID", shape="Species.ID",title="Unweighted Unifrac PCoA") + geom_polygon(data=chulls, aes(x=x, y=y, fill=Species.ID), alpha=0.2) + geom_point(size=3)




top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:500]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample.number", fill="Phylum") #+ facet_wrap(~Species.ID, scales="free_x")

ps2 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
plot_bar(ps2, x="Sample.number", fill="Phylum") 

gpt <- subset_taxa(ps, Kingdom=="k__Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:300]), gpt)
plot_heatmap(gpt, sample.label="Species.Location", low="#000033", high="#FF3300")



