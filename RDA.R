#########RDA##########
#setup: load libraries we need and prep the work space so it has everything we need
library("ape")
library("ggplot2") #Plotting functions! we love ggplot
library ("adegenet") #Does most of the PCA, DAPC, and plotting functions
library("gplots")
library("devtools")
library(rgl)
library(plyr)
library(vegan)
library(radiator)

setwd("C:/Users/Court/Documents/LaneEffectMS") #set working directory to where all input files are

#########load data in, genepop file (found in GitHub repo)#############
View(read)
table <- tidy_genepop(pl, strata = lane, tidy = FALSE, filename = "table.txt")
View (table)

G <- read.table("table.txt")
lane <- "Ind_Lane.txt"

#######################################################################
pl <- "BYLANE_NEUTRAL_testmissingindv_0.2_library1rem_NonNucrem.gen"
read <- read.genepop(pl, ncode=3)
lane <- "Ind_Lane.txt"
l <- read.table(lane)
lane2


a.comm <- data.frame(as.matrix(read$tab))
str(a.comm)

a.env <- data.frame(as.matrix(table$POP_ID))
str(a.env)
a.env

a.rda <- vegan::rda(a.comm ~ read$pop, data = a.comm, na.action = na.exclude)

help(package = "vegan")
#############################################
pl <- "BYLANE_NEUTRAL_testmissingindv_0.2_library1rem_NonNucrem.gen"
read <- read.genepop(pl, ncode=3)
lane2 <- "Ind_Lane.txt"
l <- read.table(lane2)
popmap <- read.table("popmap.txt")
popmap <- factor(popmap$V2)

genotype_matrix <- tab(read, freq=TRUE)
sum(is.na(genotype_matrix))
genotype_matrix_clean <- na.omit(genotype_matrix)

lane <- data.frame(l)
lane_data <- factor(lane$V2)
lane_data
sum(is.na(lane_data)) #no missing data

impute_mode <- function(x) {
  mode_value <- as.numeric(names(sort(table(x), decreasing = TRUE)[1])) # Most frequent value
  x[is.na(x)] <- mode_value
  return(x)
}

# Apply this function to each column (locus) of the genotype matrix
genotype_matrix_clean <- apply(genotype_matrix, 2, impute_mode)
nrow(genotype_matrix_clean)
sum(is.na(genotype_matrix_clean))

rda_lane <- rda(genotype_matrix_clean~lane_data, na.action=na.exclude)

plot <- ordiplot(rda_lane, type = "text", cex = 0.7)
plot <- ordiplot(rda_lane, type = "points", cex = 0.7)
pop <- ordiellipse(plot, groups = lane)
rda_lane

colors <- rainbow(length(levels(popmap)))  # Create a color palette for each day
group_colors <- colors[as.numeric(popmap)]
plot(rda_lane, type = "none")
points(rda_lane, display = "sites", col = group_colors, pch = 16)
# Optionally, add ellipses for each lane group (or another grouping factor)
pop <- ordiellipse(rda_lane, groups = lane_data, col = "black", lwd = 2)
# Display RDA results
rda_lane
legend("topright", 
       legend = levels(popmap),  # Use the levels of the day_data factor
       col = colors,               # Use the colors for the legend
       pch = 16,                   # Same point type used in the plot
       title = "pop",    # Title for the legend
       bty = "n")

rda_lanenrow(genotype_matrix)
nrow(lane_data)
length(lane_data)
####################IMPUTE USING LEA##################
pl <- "NEUTRAL_BAYESCANOUTLIERS_testmissingindv0.2_library1rem_NonNucrem.gen"
genepop <- read.genepop(pl, ncode=3)

library(dartR)
genlight <- gi2gl(genepop, parallel=FALSE, verbose=NULL)

#genlight to geno
path <- "C:/Users/Court/Documents/LaneEffectMS/impute"
geno <- gl2geno(genlight, outfile="bayescan.geno", outpath=path, verbose=NULL)

#change working directory to impute folder
setwd("C:/Users/Court/Documents/LaneEffectMS/impute") #set working directory to where all input files are

#run snmf 
project.snmf = snmf("bayescan.geno.lfmm", K = 4, alpha = 10,
                    entropy = TRUE, repetitions = 25, iterations = 200,
                    project = "new")

plot(project.snmf, lwd = 5, col = "red", pch=1)

ce = cross.entropy(project.snmf, K = 4)
best = which.min(cross.entropy(project.snmf, K = 4))

impute(project.snmf, "bayescan.geno.lfmm", method = 'mode', K = 4, run = best)

###########################
library(LEA)
library(vegan)
geno <- read.lfmm("bayescan.geno.lfmm")

rda <- rda(geno ~ lane_data)

plot <- ordiplot(rda, type = "text", cex = 0.7)
plot <- ordiplot(rda, type = "points", cex = 0.7)
pop <- ordiellipse(plot, groups = lane)
