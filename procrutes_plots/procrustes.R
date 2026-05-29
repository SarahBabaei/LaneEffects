library(vegan)
library(adegenet)
library(corrplot)
setwd("C:/Users/Court/Documents/LaneEffectMS/procrustes")
color <- c("#a50026", "purple", "#a6d96a", "blue")

##########read in lane effect free dataset################
#once this is read in and everything, leave it cuz all other datasets are being compared to this one.
pl <- "A_L4_LD_no2005_popmap5_NonNuc_LG_consensusoutliers_rem.gen"
read <- read.genepop(pl, ncode=3)
read

#population names: rename
print(levels(pop(read)))
levels(pop(read)) <- c("3Ps", "Northern", "GSL", "Southern")

#get genetic distance matrix
tab <- dist(tab(read, NA.method = "mean"))

#get PCA of genetic data
dist <- tab(read, NA.method="mean")
pca1 <- dudi.pca(dist, scannf = FALSE, scale = FALSE, nf=4) #nf=#of PC axes to retain

#sanitycheckplot
s.class(pca1$li,pop(read),col=transp(color),xax=1,yax=2,axesel=FALSE, 
        clabel = FALSE, cstar=0, cpoint = 2, grid=FALSE, cellipse = 0)
legend("right", legend = levels(pop(read)), col=color, pch = 19, cex = 0.8)

##########read in mitigated dataset##############
pl2 <- "B_NEUTRAL_5pop_testmissingindv_0.2_library1rem_NonNucrem.gen"
read2 <- read.genepop(pl2, ncode=3)
read2

#popchange
print(levels(pop(read2)))
levels(pop(read2))  <- c("3Ps", "Southern", "Northern", "GSL")

#get genetic distance matrix
tab2 <- dist(tab(read2, NA.method = "mean"))

#get PCA of genetic data
dist2 <- tab(read2, NA.method="mean")
pca2 <- dudi.pca(dist2, scannf = FALSE, scale = FALSE, nf=4) #nf=#of PC axes to retain

#sanitycheckplot
s.class(pca2$li,pop(read2),col=transp(color),xax=1,yax=2,axesel=FALSE, 
        clabel = FALSE, cstar=0, cpoint = 2, grid=FALSE, cellipse = 0)
legend("right", legend = levels(pop(read2)), col=color, pch = 19, cex = 0.8)


###############Procrustes Rotation: Population LEvel##############
df_pca1 <- data.frame(pca1$li[, 1:2], Pop = pop(read))
df_pca2 <- data.frame(pca2$li[, 1:2], Pop = pop(read2))

centroids1 <- aggregate(cbind(Axis1, Axis2) ~ Pop, data = df_pca1, FUN = mean)
centroids2 <- aggregate(cbind(Axis1, Axis2) ~ Pop, data = df_pca2, FUN = mean)

rownames(centroids1) <- centroids1$Pop
rownames(centroids2) <- centroids2$Pop

shared_pops <- intersect(rownames(centroids1), rownames(centroids2))
shared_pops

mat_A_mitigated <- centroids1[shared_pops, c("Axis1", "Axis2")]
mat_B_free      <- centroids2[shared_pops, c("Axis1", "Axis2")]

#run procrustes
proc <- procrustes(X = mat_A_mitigated, Y = mat_B_free, scale = TRUE, symmetric = TRUE)
summary(proc)
#significnace test
prot <- protest(X = mat_A_mitigated, Y = mat_B_free, permutations = 999, symmetric = TRUE)
print(prot)
#see pop movement
plot(proc, kind = 1, main = "Procrustes: Mitigated A vs. Lane-Effect Free Populations")
text(proc, display = "target", labels = rownames(mat_A_mitigated), cex = 0.8, pos = 3, col = "blue") #add pop names
