## Attention!!!
## The script below works if you have run script_pca.R first

library(adegenet)
library(hierfstat)
library(pheatmap)



pops <- read.table("Info_dataset_classification_PCA.txt", header = T)

# Assign population info to the relevant genind slot
testudo@pop <- as.factor(pops$Classification)

# Compute Fst
fst <- pairwise.WCfst(testudo)
fst

# FST matrix
fst_mat <- as.matrix(fst)
rownames(fst_mat) <- colnames(fst_mat) <- c("ITP+NCA", "SCA+SIC+SAR", "CCA", "MEC", "GRE")

# Clustering for plot
hc <- hclust(as.dist(fst_mat), method = "average")

# Heatmap
pheatmap(
  fst_mat,
  cluster_rows = hc,
  cluster_cols = hc,
  display_numbers = TRUE,
  number_format = "%.3f",
  border_color = "grey80",
  main = "Pairwise FST among genetic clusters",
  angle_col = 45
)
