###########setup#############
library(vegan) #RDA function
library(adegenet) #read in genepop file

setwd("C:/Users/Court/Documents/LaneEffectMS") #set working directory to where all input files are

##########Read in and impute genepop file############
pl <- "NEUTRAL_testmissingindv_0.2_library1rem_NonNucrem.gen"
read <- read.genepop(pl, ncode=3) #reads in genepop file

lane <- "Ind_Lane.txt" #2 columns: one with individual ID and one with lane #(1or2)
                       #Make sure that this file has individuals in the same order as the genepop file
l <- read.table(lane)
lane <- data.frame(l)
lane_data <- factor(lane$V2)
lane_data

#convert genind object to a genotype matrix
genotype_matrix <- tab(read, freq=TRUE)

#function for imputation
impute_mode <- function(x) {
  mode_value <- as.numeric(names(sort(table(x), decreasing = TRUE)[1])) #mode: most frequent value
  x[is.na(x)] <- mode_value
  return(x)
}

#impute using the global model imputation function built earlier
genotype_matrix_imputed <- impute_mode(genotype_matrix)


##################RDA######################
#run RDA with imputed genotype matrix and lane data from earlier
rda_lane <- rda(genotype_matrix_imputed~lane_data, na.action=na.exclude)

#plot RDA
plot(rda_lane, scaling = 2, display = c("sites", "bp"))

# add colors by lane
points(rda_lane, display = "sites", col = as.numeric(lane_data), pch = 19, scaling = 2)

# add legend
legend("topright", legend = levels(lane_data),
       col = 1:length(levels(lane_data)), pch = 19)

#plot <- ordiplot(rda_lane, type = "text", cex = 0.7)

#test if lane is still a significant predictor in the RDA
anova_terms <- anova.cca(rda_lane, by = "term", permutations = 999)
anova_terms


###########################Remove lane effect using RDA##################
#first, run the RDA using the non mitigated dataset (NEUTRAL_testmissingindv_0.2_library1rem_NonNucrem.gen)
#do this above

# Get species (SNP) scores
snp_loadings <- scores(rda_lane, display = "species", scaling = 2)

# For example, loadings on RDA1
axis1 <- snp_loadings[, 1]   # numeric vector of loadings per SNP

mean_axis1 <- mean(axis1, na.rm = TRUE)
sd_axis1 <- sd(axis1, na.rm = TRUE)

# 4 standard deviation threshold
upper <- mean_axis1 + 4 * sd_axis1
lower <- mean_axis1 - 4 * sd_axis1

#find outliers based on the 4 STDEV thresholds
outliers <- which(axis1 > upper | axis1 < lower)
length(outliers)   # how many SNPs
outlier_snps <- rownames(snp_loadings)[outliers]

#remove outlier snps from the genotype matrix
genotype_matrix_filtered <- genotype_matrix_imputed[, -outliers]

#run another RDA based on lane to see if it was removed
rda_stdev <- rda(genotype_matrix_filtered~lane_data, na.action=na.exclude)

#test if lane is still a significant predictor in the RDA
anova_terms <- anova.cca(rda_stdev, by = "term", permutations = 999)
anova_terms

#plot RDA
plot(rda_stdev, scaling = 2, display = c("sites", "bp"))

# add colors by lane
points(rda_stdev, display = "sites", col = as.numeric(lane_data), pch = 19, scaling = 2)

# add legend
legend("topright", legend = levels(lane_data),
       col = 1:length(levels(lane_data)), pch = 19)

#PCA using this dataset to see if lane effect is gone
library(adegenet)

# Filtered genind (keep everything consistent)
genind_filtered <- read[, -outliers]

#get PCA
tab <- tab(genind_filtered, NA.method="mean")
pca1 <- dudi.pca(tab, scannf = FALSE, scale = FALSE, nf=4) #nf=#of PC axes to retain
temp <- as.integer(pop(read))

#plot
goodcolors <- c("black","#b2182b" ,"#fddbc7", "#92c5de", "#2166ac")
s.class(pca1$li,pop(read),col=transp(goodcolors),xax=1,yax=2,axesel=FALSE, clabel = FALSE, cstar=0, cpoint = 2, grid=FALSE)
legend("right", legend = levels(pop(read)), col=goodcolors, pch = 19, cex = 0.8)
add.scatter.eig(pca1$eig[1:20],3,1,2, ratio=.2)

# Basic PCA scatterplot without points
plot(pca1$li[,1], pca1$li[,2], type = "n", xlab = "PC1", ylab = "PC2")

# Add sample names instead of points
text(pca1$li[,1], pca1$li[,2], labels = rownames(pca1$li), col = transp(goodcolors)[pop(read)])

#variation and more PC info
pca1$eig[1]
eig.perc <- 100*pca1$eig/sum(pca1$eig)
head(eig.perc)
