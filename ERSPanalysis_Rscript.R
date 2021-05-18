library(ape)
library(adegenet)
library(ade4)
library(vegan)
library(vcfR)
library(poppr)
library(strataG)
library(pegas)
library(ggplot2)
library(tidyverse)
library(mmod)
library(genepop)
library(doBy)
library(dplyr)
##################################################################################################################
#####
ERSPvcf <- read.vcfR(file = "ERSP_pop4_haps_011921.vcf", verbose = FALSE)
ERSPgenind <- vcfR2genind(ERSPvcf)
ERSPgenind
dim(ERSPgenind@tab)
ERSPgenind@tab[1:194,1:5]
populations <- read.csv(file = 'pop_names.csv', header = TRUE)
pop(ERSPgenind) <- populations$pops
######Add spatial coordinates in the @other$xy slot of the genind file
xy <-cbind(c("BE1","CO1","CR1","DC1","DF8","FL2","FL7","HC1","LC1","MC1","MF1","MF5","MF8","PA1","PC1","PR1","PR2","SC1","SP1","SQ1","UW4","WC1","WH1"),
           c("S","S","S","S","S","S","S","N","S","N","S","S","S","N","N","N","N","N","S","N","S","S","S"),
           c(42.324851,43.503071,40.65012,38.2531,38.01667,38.498851,38.48345,43.69868,41.821355,45.479338,39.13333,39.58354,39.86674,44.94229,45.1857,45.171797,45.716715,44.80359,38.4996,44.381692,40.5411,41.228727,42.979552),
           c(-113.6026558,-110.9673183,-109.500053,-106.6725,-112.20011,-112.008531,-111.450168,-114.188931,-111.4731907,-112.0787368,-111.60007,-111.21674,-111.26687,-115.94232,-114.1484,-108.4314097,-113.3419806,-111.40581,-106.4539,-114.497194,-111.0322167,-111.8539976,-109.7567062))
colnames(xy) <- c("pop","kpop","x","y")
other(ERSPgenind) <- xy
xy <- as.data.frame(xy)
ERSPgenind$other$xy <-cbind(c(42.324851,43.503071,40.65012,38.2531,38.01667,38.498851,38.48345,43.69868,41.821355,45.479338,39.13333,39.58354,39.86674,44.94229,45.1857,45.171797,45.716715,44.80359,38.4996,44.381692,40.5411,41.228727,42.979552),
c(-113.6026558,-110.9673183,-109.500053,-106.6725,-112.20011,-112.008531,-111.450168,-114.188931,-111.4731907,-112.0787368,-111.60007,-111.21674,-111.26687,-115.94232,-114.1484,-108.4314097,-113.3419806,-111.40581,-106.4539,-114.497194,-111.0322167,-111.8539976,-109.7567062))


#######PCA
pca1 <- tab(ERSPgenind, freq=TRUE, NA.method="mean")
pca2 <- dudi.pca(pca1, center=TRUE, scale=FALSE)
s.label(pca2$li)
###PCA with ellipses
s.class(pca2$li, fac=pop(ERSPgenind), col=funky(23),clabel = 0.8, sub = "PC1/PC2", possub = "topleft")
###Eigenvector percent variation
eig.perc <- 100*pca2$eig/sum(pca2$eig)
head(eig.perc)
###Relationship between Eigenvectors and Latitude
Eig_vect <- as.data.frame((pca2$li))
Eigen_lat <- read.csv(file = 'Eigenvectors.csv', header = TRUE)
Eigen_lat_sum <- summaryBy(Axis1+Axis2~Pop, data = Eigen_lat, FUN = c(mean))
xy$x <- as.numeric(as.character(xy$x))
xy$y <- as.numeric(as.character(xy$y))
Eigen_lat_sum <- cbind(Eigen_lat_sum,xy)
South_lat <- subset(Eigen_lat_sum, kpop=="S")
North_lat <- subset(Eigen_lat_sum, kpop=="N")
#hh <- ggplot(South_lat, aes(Axis2.mean,x))
#hh + geom_point() + geom_smooth(method = lm, se = FALSE) + theme_bw() +
south_reg <- lm(South_lat$x~South_lat$Axis2.mean)
summary(south_reg)
coeff = coefficients(south_reg)
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))

###Range-wide IBD
ERSPgenpop <- genind2genpop(ERSPgenind)
ERSPdist <- dist.genpop(ERSPgenpop,method = 2)
ERSPspatial <- dist(ERSPgenind$other$xy)
ibd <- mantel.randtest(ERSPdist,ERSPspatial)
ibd
plot(ibd)
plot(ERSPdist,ERSPspatial)
#gg <- ggplot(data = NULL, aes(ERSPdist,ERSPspatial))
#gg + geom_point() + geom_smooth(method = "lm", se = FALSE)
km_fst <- read.csv(file = 'km_fst.csv', header = TRUE)
cor.test(km_fst$KM,km_fst$FST)
ff <- ggplot(data = km_fst, aes(KM,FST))
ff + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_bw()
             
spat_lm <- lm(ERSPdist~ERSPspatial)
summary(spat_lm)

##Diversity and geographic distance
diversity <- read.csv(file = "diversity_sum_stat.csv", header = TRUE)
diversity <- cbind(diversity, Eigen_lat_sum$Axis1.mean,Eigen_lat_sum$Axis2.mean)
div_south <- dplyr::filter(diversity, Pop == "South")
colnames(div_south)[30] <- "Vect1"
colnames(div_south)[31] <- "Vect2"
##Plot South Pop lat with Pi
gg <- ggplot(data = div_south, aes(Pi,latitude, label = Site_ID))
gg + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_bw() + geom_text(hjust = 0,nudge_x = 0.0005)
cor.test(div_south$latitude,div_south$Pi)

#Plot South Pop lat with Eigenvectors
hh <- ggplot(data = div_south, aes(Vect2,latitude, label = Site_ID))
hh + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_bw() + geom_text(hjust = 0,nudge_x = 0.05)
cor.test(div_south$latitude,div_south$Vect2)

#Plot south population for Fig3
plot(South_lat$Axis2.mean,South_lat$x.1)
abline(south_reg,col="blue")
text(South_lat$Axis2.mean,South_lat$x.1,labels = South_lat$Pop,cex=0.8, pos=2)
cor.test(South_lat$Axis2.mean,South_lat$x.1)
#Plot north population
plot(North_lat$Axis2.mean,North_lat$x.1)
cor.test(North_lat$Axis1.mean,North_lat$x.1)


###UPGMA Tree w/ outgroups
ERSPvcf_full <- read.vcfR(file = "populations.haps.vcf", verbose = FALSE)
ERSPgenind_full <- vcfR2genind(ERSPvcf_full)
ERSPgenind_full
dim(ERSPgenind_full@tab)
ERSPgenind_full@tab[1:214,1:5]
#Subset out DF2 samples
ERSPgenind_full@pop
ERSP_mod <- popsub(ERSPgenind_full, sublist = c("BE1", "CO1", "CR1", "DC1","DF8", "FL2", "FL7", "HC1", "LC1", "MC1", "MF1", "MF5", "MF8", "ML1", "PA1","PC1","PR1","PR2","SC1","SP1","SQ1","UW4","WC1","WH1"))

populations_full <- read.csv(file = 'pop_names_full.csv', header = TRUE)
pop(ERSPgenind_full) <- populations_full$pops

gpop1 <- genind2genpop(ERSP_mod)
poptree <- aboot(gpop1)

