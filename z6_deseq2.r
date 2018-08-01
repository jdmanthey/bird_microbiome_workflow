load("microbe_workflow5.Rdata")
library(phyloseq); packageVersion("phyloseq")
library(DESeq2)
library(ggplot2)
ps@sam_data$Genus.ID
# genus
bird_genus <- phyloseq_to_deseq2(ps, ~ Genus.ID)
bird_genus <- DESeq(bird_genus, test="Wald", fitType="parametric")

# certhia
certhia_ps <- subset_samples(ps, Genus.ID != "Sitta")
certhia <- phyloseq_to_deseq2(certhia_ps, ~ Species.ID)
certhia <- DESeq(certhia, test = "Wald", fitType="parametric")

#format outputs
res_genus <- results(bird_genus, cooksCutoff=F)
alpha <- 0.01
sigtab_genus <- res_genus[which(res_genus$padj < alpha), ]
sigtab_genus <- cbind(as(sigtab_genus, "data.frame"), as(tax_table(ps)[rownames(sigtab_genus), ], "matrix"))
head(sigtab_genus)
# high log2FoldChange = more in Sitta

res_certhia <- results(certhia, cooksCutoff=F)
alpha <- 0.01
sigtab_certhia <- res_certhia[which(res_certhia$padj < alpha), ]
sigtab_certhia <- cbind(as(sigtab_certhia, "data.frame"), as(tax_table(certhia_ps)[rownames(sigtab_certhia), ], "matrix"))
head(sigtab_certhia)
# high log2FoldChange = more in americana

# plots

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}


# Phylum order
x = tapply(sigtab_genus$log2FoldChange, sigtab_genus$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_genus$Phylum = factor(as.character(sigtab_genus$Phylum), levels=names(x))
# Order order
x = tapply(sigtab_genus$log2FoldChange, sigtab_genus$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_genus$Genus = factor(as.character(sigtab_genus$Genus), levels=names(x))
ggplot(sigtab_genus, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# Phylum order
x = tapply(sigtab_certhia$log2FoldChange, sigtab_certhia $Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_certhia$Phylum = factor(as.character(sigtab_certhia$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_certhia$log2FoldChange, sigtab_certhia$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_certhia$Order = factor(as.character(sigtab_certhia$Genus), levels=names(x))
ggplot(sigtab_certhia, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

write.table(sigtab_genus, file="deseq2_certhia_sitta.txt",quote=F)
write.table(sigtab_certhia, file="deseq2_certhia_albescens_americana.txt",quote=F)




# extract OTU tables of microbes different between groups
genus_diff <- otu_table(certhia_ps)[ ,rownames(sigtab_certhia)]
genus_diff <- cbind(certhia_ps@sam_data$Species.ID, genus_diff)
genus_diff <- as.data.frame(genus_diff)
colnames(genus_diff)


