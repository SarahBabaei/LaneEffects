#########PCA##########
#setup: load libraries we need and prep the work space so it has everything we need
library ("adegenet") #Does most of the PCA, DAPC, and plotting functions

setwd("C:/Users/Court/Documents/LaneEffectMS")

###load data in, genepop file (found in GitHub repo)
pl <- "DAPCtop10%rem_NEUTRAL_testmissingindv_0.2_library1rem_NonNucrem.gen"
read <- read.genepop(pl, ncode=3)
read

#get PCA
tab <- tab(read, NA.method="mean")
pca1 <- dudi.pca(tab, scannf = FALSE, scale = FALSE, nf=10) #nf=#of PC axes to retain
temp <- as.integer(pop(read))

#plot
#goodcolors<- c("black", "#a1d99b", "#41ab5d", "#00441b", "#6a51a3", "#9e9ac8", "#dadaeb", "#6baed6", "#08306b")
goodcolors<- c("black", "#a1d99b", "#00441b","#a1d99b", "#00441b", "#dadaeb", "#6a51a3", "#92c5de", "#08306b")
#goodcolors<- c("grey", "#a1d99b", "#00441b","#a1d99b", "#00441b", "#dadaeb", "#542788", "#92c5de", "#08306b")
s.class(pca1$li,pop(read),col=transp(goodcolors),xax=1,yax=2,axesel=FALSE, clabel = FALSE, cstar=0, cpoint = 2, grid=FALSE, cellipse = 0)
legend("right", legend = levels(pop(read)), col=goodcolors, pch = 19, cex = 0.8)
add.scatter.eig(pca1$eig[1:20],3,1,2, ratio=.1)

#variation per PC axis and more PC info
pca1$eig[1]
eig.perc <- 100*pca1$eig/sum(pca1$eig)
head(eig.perc)
