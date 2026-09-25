library(ade4)
library(adegenet)
library(ggplot2)

# Set your working directory (or navigate to it in R Studio)
setwd("D:/Side projects/Physalia course/Day3/tutorial_snpchip/01_pca")

# Read in the data
testudo <- read.structure("../batch_1.stru",
                          n.ind = 70,
                          n.loc = 3182,
                          onerowperind = FALSE,
                          col.lab = 1,
                          col.pop = 2,
                          row.marknames = 1,
                          NA.char = "0")

# Add metadata file
pops <- read.table("Info_dataset_classification_PCA.txt", header = T)

# Make the Pop_ID column a factor
pops$Pop_ID <- as.factor(pops$Pop_ID)
str(pops)



## Check missingness patterns
# By sample
prop_ind <- propTyped(testudo, by = "ind")
hist(prop_ind)

# Who has more than 20% missing data?
prop_ind[which(prop_ind < 0.8)]

# By marker
prop_loc <- propTyped(testudo, by = "loc")
hist(prop_loc)

# Decide threshold for markers (depends on your dataset)
table(pops$Classification)
prop_loc_keep <- (1 - length(which(pops$Classification == 4 | pops$Classification == 5))/nrow(pops))

# Filter
gen_filt <- testudo[
  prop_ind >= 0.8 & indNames(testudo) != "CAL_RO1",
  loc = prop_loc >= prop_loc_keep
]



## Re-run PCA
# Transform the missing data
x.test <- tab(gen_filt, freq=TRUE, NA.method="mean")
x.test <- x.test[, seq(1, ncol(x.test), by = 2)]

# Perform PCA
pca.testudo <- dudi.pca(x.test, center=TRUE, scale=FALSE)

# Basic PCA representation with sample labels, useful to check where each sample falls
s.label(pca.testudo$li,clabel = 0.35)

# Extract eigenvalues
eig.perc <- 100*pca.testudo$eig/sum(pca.testudo$eig)
head(eig.perc)
plot(eig.perc[0:6], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot by the first 6 PCs")

# Extract population info for kept samples
pops_filt <- pops[pops$Sample %in% rownames(gen_filt$tab),]
pops_filt$Classification <- as.factor(pops_filt$Classification)

# Plot pc1 - pc2
p12 <- ggplot(
  data = pca.testudo$li,
  aes(x = Axis1, y = Axis2, color = pops_filt$Classification)
) +
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

p34 <- ggplot(
  data = pca.testudo$li,
  aes(x = Axis3, y = Axis4, color = pops_filt$Classification)
) +
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



## Loadings
# Loading plot for PC1, PC2, PC3, and PC4.
# threshold at top 99%*
# * only for Axis1 we set it as 98.9% as all the highest values correspond to the threshold itself and this gives a bug that doesn't allow to save the markers
load1 <- loadingplot(pca.testudo$c1[,1]^2, threshold=quantile(pca.testudo$c1[,1]^2,0.989),
                     axis=1, cex.lab=0.3, cex.fac=1, lab.jitter=0,
                     main="Loading plot", xlab="Variables", ylab="Loadings",
                     srt = 0, adj = NULL, lab = rownames(pca.testudo$c1))
topload1 <- data.frame(load1$var.names, load1$var.values)
topload1 <- topload1[order(topload1[,2], decreasing = TRUE), ]

load2 <- loadingplot(pca.testudo$c1[,2]^2, threshold=quantile(pca.testudo$c1[,2]^2,0.99),
                     axis=1, cex.lab=0.3, cex.fac=1, lab.jitter=0,
                     main="Loading plot", xlab="Variables", ylab="Loadings",
                     srt = 0, adj = NULL, lab = rownames(pca.testudo$c1))
topload2 <- data.frame(load2$var.names, load2$var.values)
topload2 <- topload2[order(topload2[,2], decreasing = TRUE), ]

load3 <- loadingplot(pca.testudo$c1[,3]^2, threshold=quantile(pca.testudo$c1[,3]^2,0.99),
                     axis=1, cex.lab=0.3, cex.fac=1, lab.jitter=0,
                     main="Loading plot", xlab="Variables", ylab="Loadings",
                     srt = 0, adj = NULL, lab = rownames(pca.testudo$c1))
topload3 <- data.frame(load3$var.names, load3$var.values)
topload3 <- topload3[order(topload3[,2], decreasing = TRUE), ]

load4 <- loadingplot(pca.testudo$c1[,4]^2, threshold=quantile(pca.testudo$c1[,4]^2,0.99),
                     axis=1, cex.lab=0.3, cex.fac=1, lab.jitter=0,
                     main="Loading plot", xlab="Variables", ylab="Loadings",
                     srt = 0, adj = NULL, lab = rownames(pca.testudo$c1))
topload4 <- data.frame(load4$var.names, load4$var.values)
topload4 <- topload4[order(topload4[,2], decreasing = TRUE), ]



## Save results
write.table(topload1, file = "loadings_PC1.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload2, file = "loadings_PC2.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload3, file = "loadings_PC3.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload4, file = "loadings_PC4.txt", sep = "\t", quote = F, row.names = T, col.names = F)
