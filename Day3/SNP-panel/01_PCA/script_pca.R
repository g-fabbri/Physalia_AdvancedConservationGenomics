library(ade4)
library(adegenet)
library(ggplot2)
library(mclust)

# Set your working directory (or navigate to it in R Studio)
setwd("D:/Side projects/Physalia course/Day3/tutorial_snpchip")



## Prepare the dataset
# Read in the data
testudo <- read.structure("testudo_dataset.stru",
                          n.ind = 70,
                          n.loc = 3182,
                          onerowperind = FALSE,
                          col.lab = 1,
                          col.pop = 2,
                          row.marknames = 1,
                          NA.char = "0")

# Check the structure of the input file
testudo

# Add metadata file
pops <- read.table("Info_dataset.txt", header = T)

# Make the Pop_ID column a factor
pops$Sampling_location <- as.factor(pops$Sampling_location)
str(pops)



## Check missingness patterns
# By sample
prop_ind <- propTyped(testudo, by = "ind")
hist(prop_ind)

# Who has more than 20% missing data?
prop_ind[which(prop_ind < 0.8)]

# Filter dataset for sample missingness
gen_filt <- testudo[prop_ind >= 0.8]
str(gen_filt)

# Extract population info for kept samples
pops_filt <- pops[pops$Sample %in% rownames(gen_filt$tab),]

# By marker
prop_loc <- propTyped(gen_filt, by = "loc")
hist(prop_loc)

# Filter dataset for locus missingness
x.tab <- tab(gen_filt, freq=TRUE, NA.method="asis")
dim(x.tab)
x.tab <- x.tab[, seq(1, ncol(x.tab), by = 2)]
x.tab.filt <- x.tab[, colMeans(is.na(x.tab)) < 0.20]
dim(x.tab.filt)

# Transform the missing data
for (j in seq_len(ncol(x.tab.filt))) {
  x.tab.filt[is.na(x.tab.filt[, j]), j] <- mean(x.tab.filt[, j], na.rm = TRUE)
}

# Check if there are still missing data
sum(is.na(x.tab.filt))



## Perform PCA
pca.testudo <- dudi.pca(x.tab.filt, center=TRUE, scale=FALSE)

# Extract eigenvalues
eig.perc <- 100*pca.testudo$eig/sum(pca.testudo$eig)
head(eig.perc)
plot(eig.perc[0:6], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot by the first 6 PCs")



## Plot
# Basic PCA representation with sample labels, useful to check where each sample falls
s.label(pca.testudo$li,clabel = 0.35)

# PC1-PC2
p12 <- ggplot(
  data = pca.testudo$li,
  aes(x = Axis1, y = Axis2, color = pops_filt$Sampling_location)) +
  geom_point(shape = 16, size = 3, alpha = 0.7) +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  labs(
    y = paste("PC2", "(", round(eig.perc[2], digits = 2), "%", ")"),
    x = paste("PC1", "(", round(eig.perc[1], digits = 2), "%", ")"),
    color = "Populations"
  ) +
  scale_color_brewer(palette = "Paired")

p12

# PC3-PC4
p34 <- ggplot(
  data = pca.testudo$li,
  aes(x = Axis3, y = Axis4, color = pops_filt$Sampling_location)) +
  geom_point(shape = 16, size = 3, alpha = 0.7) +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  labs(
    y = paste("PC4", "(", round(eig.perc[4], digits = 2), "%", ")"),
    x = paste("PC3", "(", round(eig.perc[3], digits = 2), "%", ")"),
    color = "Populations"
  ) +
  scale_color_brewer(palette = "Paired")

p34



## Decide clustering
# Extract PC coordinates
pca_coord <- as.data.frame(pca.testudo$li)

# Classification based on the first 4 PCs
nclust <- Mclust(pca_coord[,1:4])
nclust
summary(nclust)

# Re-plot with grouping-based colors
# PC1-PC2
ggplot(data = pca_coord, aes(x = Axis1, y = Axis2, color = as.factor(nclust$classification)))+
  geom_point(shape = 16, size = 3, alpha = 0.7) +
  scale_color_manual(values = c("1" = "deepskyblue", "2" = "purple", "3" = "black", "4" = "brown", "5" = "orange")) +
  theme_classic() +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  xlab(paste0("PC1 (", round(eig.perc[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(eig.perc[2], digits = 2), "%)")) +
  labs(color = "Clusters")

# PC3-PC4
ggplot(data = pca_coord, aes(x = Axis3, y = Axis4, color = as.factor(nclust$classification)))+
  geom_point(shape = 16, size = 3, alpha = 0.7) +
  scale_color_manual(values = c("1" = "deepskyblue", "2" = "purple", "3" = "black", "4" = "brown", "5" = "orange")) +
  theme_classic() +
  geom_hline(aes(yintercept=0)) +
  xlab(paste0("PC3 (", round(eig.perc[3], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(eig.perc[4], digits = 2), "%)")) +
  labs(color = "Clusters")



## Save metadata with new information on classification
meta <- cbind(pops_filt, nclust$classification)
write.table(meta, file = "Info_dataset_classification_PCA.txt", sep = "\t", quote = F, row.names = F)



#######################################################
# Decide threshold for markers (depends on your dataset)
table(pops$Classification) # this is to check the smallest groups so you set a threshold to avoi having missing data just there
prop_loc_keep <- (1 - length(which(pops$Classification == 4 | pops$Classification == 5))/nrow(pops))
