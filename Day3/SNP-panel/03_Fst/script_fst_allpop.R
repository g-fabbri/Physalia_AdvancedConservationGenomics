library(adegenet)
library(hierfstat)
library(pheatmap)

# Set your working directory (or navigate to it in R Studio)
setwd("D:/Side projects/Physalia course/Day3/tutorial_snpchip/03_fst")

# Read in the data
testudo <- read.structure("../testudo_dataset.stru",
                                n.ind = 70,
                                n.loc = 3182,
                                onerowperind = FALSE,
                                col.lab = 1,
                                col.pop = 2,
                                row.marknames = 1,
                                NA.char = "0")

pops <- read.table("../01_pca/Info_dataset_classification_PCA.txt", header = T)

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
