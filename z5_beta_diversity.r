# three datasets each for individuals and localities
ps1 <- ps
ps2 <- tip_glom(ps, h=0.03)
ps3 <- tip_glom(ps, h=0.06)

psl1 <- ps.locales
psl2 <- tip_glom(ps.locales, h=0.03)
psl3 <- tip_glom(ps.locales, h=0.06)

bray1 <- c()
wuni1 <- c()
uuni1 <- c()
bray2 <- c()
wuni2 <- c()
uuni2 <- c()
bray3 <- c()
wuni3 <- c()
uuni3 <- c()
for(a in 1:100) {
	print(a)
	a.combs <- combn(nrow(ps1@sam_data),2) # all combinations
	# obtain coordinates for each of the three methods
	a.rep1 <- rarefy_even_depth(ps1, 15000)
	bray.tmp1 <- ordinate(a.rep1, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp1 <- ordinate(a.rep1, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp1 <- ordinate(a.rep1, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	a.rep2 <- rarefy_even_depth(ps2, 15000)
	bray.tmp2 <- ordinate(a.rep2, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp2 <- ordinate(a.rep2, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp2 <- ordinate(a.rep2, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	a.rep3 <- rarefy_even_depth(ps3, 15000)
	bray.tmp3 <- ordinate(a.rep3, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp3 <- ordinate(a.rep3, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp3 <- ordinate(a.rep3, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	bt1 <- c()
	wt1 <- c()
	ut1 <- c()
	bt2 <- c()
	wt2 <- c()
	ut2 <- c()
	bt3 <- c()
	wt3 <- c()
	ut3 <- c()
	for(b in 1:ncol(a.combs)) { # obtain distances between coordinates
		bt1 <- c(bt1, sqrt((bray.tmp1[a.combs[1,b], 1] - bray.tmp1[a.combs[2,b], 1])^2 + 
							(bray.tmp1[a.combs[1,b], 2] - bray.tmp1[a.combs[2,b], 2])^2))
		wt1 <- c(wt1, sqrt((wuni.tmp1[a.combs[1,b], 1] - wuni.tmp1[a.combs[2,b], 1])^2 + 
							(wuni.tmp1[a.combs[1,b], 2] - wuni.tmp1[a.combs[2,b], 2])^2))		
		ut1 <- c(ut1, sqrt((uuni.tmp1[a.combs[1,b], 1] - uuni.tmp1[a.combs[2,b], 1])^2 + 
							(uuni.tmp1[a.combs[1,b], 2] - uuni.tmp1[a.combs[2,b], 2])^2))		
		bt2 <- c(bt2, sqrt((bray.tmp2[a.combs[1,b], 1] - bray.tmp2[a.combs[2,b], 1])^2 + 
							(bray.tmp2[a.combs[1,b], 2] - bray.tmp2[a.combs[2,b], 2])^2))
		wt2 <- c(wt2, sqrt((wuni.tmp2[a.combs[1,b], 1] - wuni.tmp2[a.combs[2,b], 1])^2 + 
							(wuni.tmp2[a.combs[1,b], 2] - wuni.tmp2[a.combs[2,b], 2])^2))		
		ut2 <- c(ut2, sqrt((uuni.tmp2[a.combs[1,b], 1] - uuni.tmp2[a.combs[2,b], 1])^2 + 
							(uuni.tmp2[a.combs[1,b], 2] - uuni.tmp2[a.combs[2,b], 2])^2))	
		bt3 <- c(bt3, sqrt((bray.tmp3[a.combs[1,b], 1] - bray.tmp3[a.combs[2,b], 1])^2 + 
							(bray.tmp3[a.combs[1,b], 2] - bray.tmp3[a.combs[2,b], 2])^2))
		wt3 <- c(wt3, sqrt((wuni.tmp3[a.combs[1,b], 1] - wuni.tmp3[a.combs[2,b], 1])^2 + 
							(wuni.tmp3[a.combs[1,b], 2] - wuni.tmp3[a.combs[2,b], 2])^2))		
		ut3 <- c(ut3, sqrt((uuni.tmp3[a.combs[1,b], 1] - uuni.tmp3[a.combs[2,b], 1])^2 + 
							(uuni.tmp3[a.combs[1,b], 2] - uuni.tmp3[a.combs[2,b], 2])^2))						
	}
	bray1 <- cbind(bray1, bt1)
	wuni1 <- cbind(wuni1, wt1)
	uuni1 <- cbind(uuni1, ut1)
	bray2 <- cbind(bray2, bt2)
	wuni2 <- cbind(wuni2, wt2)
	uuni2 <- cbind(uuni2, ut2)
	bray3 <- cbind(bray3, bt3)
	wuni3 <- cbind(wuni3, wt3)
	uuni3 <- cbind(uuni3, ut3)

}
sample1 <- c()
sample2 <- c()
for(a in 1:ncol(a.combs)) {
	sample1 <- c(sample1, as.character(ps1@sam_data$Sample.ID[a.combs[1,a]]))
	sample2 <- c(sample2, as.character(ps1@sam_data$Sample.ID[a.combs[2,a]]))
}
output <- data.frame(Sample1=sample1, Sample2=sample2,
			Bray_Original=apply(bray1, 1, function(x) mean(x)),
			Weighted_Unifrac_Original=apply(wuni1, 1, function(x) mean(x)),
			Unweighted_Unifrac_Original=apply(uuni1, 1, function(x) mean(x)),
			Bray_H03=apply(bray2, 1, function(x) mean(x)),
			Weighted_Unifrac_H03=apply(wuni2, 1, function(x) mean(x)),
			Unweighted_Unifrac_H03=apply(uuni2, 1, function(x) mean(x)),
			Bray_H06=apply(bray3, 1, function(x) mean(x)),
			Weighted_Unifrac_H06=apply(wuni3, 1, function(x) mean(x)),
			Unweighted_Unifrac_H06=apply(uuni3, 1, function(x) mean(x)))
write.table(output, file="beta_diversity_individuals.txt", sep="\t",quote=F, row.names=F)






bray1 <- c()
wuni1 <- c()
uuni1 <- c()
bray2 <- c()
wuni2 <- c()
uuni2 <- c()
bray3 <- c()
wuni3 <- c()
uuni3 <- c()
for(a in 1:100) {
	print(a)
	a.combs <- combn(nrow(psl1@sam_data),2) # all combinations
	# obtain coordinates for each of the three methods
	a.rep1 <- rarefy_even_depth(psl1, 65000)
	bray.tmp1 <- ordinate(a.rep1, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp1 <- ordinate(a.rep1, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp1 <- ordinate(a.rep1, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	a.rep2 <- rarefy_even_depth(psl2, 65000)
	bray.tmp2 <- ordinate(a.rep2, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp2 <- ordinate(a.rep2, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp2 <- ordinate(a.rep2, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	a.rep3 <- rarefy_even_depth(psl3, 65000)
	bray.tmp3 <- ordinate(a.rep3, method="NMDS", distance="bray")$points[,1:2]
	wuni.tmp3 <- ordinate(a.rep3, "PCoA", "unifrac", weighted=T)$vectors[,1:2]
	uuni.tmp3 <- ordinate(a.rep3, "PCoA", "unifrac", weighted=F)$vectors[,1:2]
	bray.tmp1 <- bray.tmp1[c(11:19, 1:3, 5:6, 4, 7:10),]
	wuni.tmp1 <- wuni.tmp1[c(11:19, 1:3, 5:6, 4, 7:10),]
	uuni.tmp1 <- uuni.tmp1[c(11:19, 1:3, 5:6, 4, 7:10),]
	bray.tmp2 <- bray.tmp2[c(11:19, 1:3, 5:6, 4, 7:10),]
	wuni.tmp2 <- wuni.tmp2[c(11:19, 1:3, 5:6, 4, 7:10),]
	uuni.tmp2 <- uuni.tmp2[c(11:19, 1:3, 5:6, 4, 7:10),]
	bray.tmp3 <- bray.tmp3[c(11:19, 1:3, 5:6, 4, 7:10),]
	wuni.tmp3 <- wuni.tmp3[c(11:19, 1:3, 5:6, 4, 7:10),]
	uuni.tmp3 <- uuni.tmp3[c(11:19, 1:3, 5:6, 4, 7:10),]
	bt1 <- c()
	wt1 <- c()
	ut1 <- c()
	bt2 <- c()
	wt2 <- c()
	ut2 <- c()
	bt3 <- c()
	wt3 <- c()
	ut3 <- c()
	for(b in 1:ncol(a.combs)) { # obtain distances between coordinates
		bt1 <- c(bt1, sqrt((bray.tmp1[a.combs[1,b], 1] - bray.tmp1[a.combs[2,b], 1])^2 + 
							(bray.tmp1[a.combs[1,b], 2] - bray.tmp1[a.combs[2,b], 2])^2))
		wt1 <- c(wt1, sqrt((wuni.tmp1[a.combs[1,b], 1] - wuni.tmp1[a.combs[2,b], 1])^2 + 
							(wuni.tmp1[a.combs[1,b], 2] - wuni.tmp1[a.combs[2,b], 2])^2))		
		ut1 <- c(ut1, sqrt((uuni.tmp1[a.combs[1,b], 1] - uuni.tmp1[a.combs[2,b], 1])^2 + 
							(uuni.tmp1[a.combs[1,b], 2] - uuni.tmp1[a.combs[2,b], 2])^2))		
		bt2 <- c(bt2, sqrt((bray.tmp2[a.combs[1,b], 1] - bray.tmp2[a.combs[2,b], 1])^2 + 
							(bray.tmp2[a.combs[1,b], 2] - bray.tmp2[a.combs[2,b], 2])^2))
		wt2 <- c(wt2, sqrt((wuni.tmp2[a.combs[1,b], 1] - wuni.tmp2[a.combs[2,b], 1])^2 + 
							(wuni.tmp2[a.combs[1,b], 2] - wuni.tmp2[a.combs[2,b], 2])^2))		
		ut2 <- c(ut2, sqrt((uuni.tmp2[a.combs[1,b], 1] - uuni.tmp2[a.combs[2,b], 1])^2 + 
							(uuni.tmp2[a.combs[1,b], 2] - uuni.tmp2[a.combs[2,b], 2])^2))	
		bt3 <- c(bt3, sqrt((bray.tmp3[a.combs[1,b], 1] - bray.tmp3[a.combs[2,b], 1])^2 + 
							(bray.tmp3[a.combs[1,b], 2] - bray.tmp3[a.combs[2,b], 2])^2))
		wt3 <- c(wt3, sqrt((wuni.tmp3[a.combs[1,b], 1] - wuni.tmp3[a.combs[2,b], 1])^2 + 
							(wuni.tmp3[a.combs[1,b], 2] - wuni.tmp3[a.combs[2,b], 2])^2))		
		ut3 <- c(ut3, sqrt((uuni.tmp3[a.combs[1,b], 1] - uuni.tmp3[a.combs[2,b], 1])^2 + 
							(uuni.tmp3[a.combs[1,b], 2] - uuni.tmp3[a.combs[2,b], 2])^2))						
	}
	bray1 <- cbind(bray1, bt1)
	wuni1 <- cbind(wuni1, wt1)
	uuni1 <- cbind(uuni1, ut1)
	bray2 <- cbind(bray2, bt2)
	wuni2 <- cbind(wuni2, wt2)
	uuni2 <- cbind(uuni2, ut2)
	bray3 <- cbind(bray3, bt3)
	wuni3 <- cbind(wuni3, wt3)
	uuni3 <- cbind(uuni3, ut3)

}
sample1 <- c()
sample2 <- c()
for(a in 1:ncol(a.combs)) {
	sample1 <- c(sample1, as.character(rownames(bray.tmp1)[a.combs[1,a]]))
	sample2 <- c(sample2, as.character(rownames(bray.tmp1)[a.combs[2,a]]))
}
output <- data.frame(Sample1=sample1, Sample2=sample2,
			Bray_Original=apply(bray1, 1, function(x) mean(x)),
			Weighted_Unifrac_Original=apply(wuni1, 1, function(x) mean(x)),
			Unweighted_Unifrac_Original=apply(uuni1, 1, function(x) mean(x)),
			Bray_H03=apply(bray2, 1, function(x) mean(x)),
			Weighted_Unifrac_H03=apply(wuni2, 1, function(x) mean(x)),
			Unweighted_Unifrac_H03=apply(uuni2, 1, function(x) mean(x)),
			Bray_H06=apply(bray3, 1, function(x) mean(x)),
			Weighted_Unifrac_H06=apply(wuni3, 1, function(x) mean(x)),
			Unweighted_Unifrac_H06=apply(uuni3, 1, function(x) mean(x)))
write.table(output, file="beta_diversity_localities.txt", sep="\t",quote=F, row.names=F)




















