library(adegenet)
library(hierfstat)

# Set your working directory (or navigate to it in R Studio)
setwd("D:/Side projects/Physalia course/Day3/tutorial_snpchip/03_fst")

# Read in the data
testudo_node1 <- read.structure("../testudo_dataset.stru",
                                n.ind = 70,
                                n.loc = 3182,
                                onerowperind = FALSE,
                                col.lab = 1,
                                col.pop = 2,
                                row.marknames = 1,
                                NA.char = "0")



## Check missingness patterns
# By sample
prop_ind <- propTyped(testudo_node1, by = "ind")
hist(prop_ind)

# Who has more than 20% missing data?
prop_ind[which(prop_ind < 0.8)]

# By marker
prop_loc <- propTyped(testudo_node1, by = "loc")
hist(prop_loc)

# Decide threshold for markers (depends on your dataset)
table(pops$Classification)
prop_loc_keep <- (1 - length(which(pops$Classification == 4 | pops$Classification == 5))/nrow(pops))

# Filter
gen_filt <- testudo_node1[
  prop_ind >= 0.8,
  loc = prop_loc >= prop_loc_keep
]

## Node 1
# Assign populations to node 1
pops <- read.table("../01_pca/Info_dataset_classification_PCA.txt", header = T)
pops_filt <- pops[pops$Sample %in% rownames(gen_filt$tab),]
pops_filt$node1 <- ifelse(
  pops_filt$Pop_ID %in% c("GRE", "MEC"),
  1,
  2
)

# Prepare dataset for Fst between the two groups at node 1
gen_node1 <- gen_filt
pop(gen_node1) <- pops_filt$node1

# Fst by marker for node 1
fst_node1 <- wc(gen_node1)

fst_node1_df <- data.frame(
  locus = names(gen_node1@all.names),
  FST = fst_node1$per.loc$FST
)

fst_node1_df <- fst_node1_df[order(fst_node1_df$FST, decreasing = TRUE), ]

# Save results
write.table(fst_node1_df, file = "fst_node1.txt", sep = "\t", quote = F, row.names = F, col.names = F)


## Node 2
# Extract samples 
gen_node2 <- gen_filt[pop(gen_node1) == 1, ]

# Assign populations to node 1
pops_filt_node2 <- pops_filt[pops_filt$Sample %in% rownames(gen_node2$tab),]
pops_filt_node2$node2 <- ifelse(
  pops_filt_node2$Pop_ID %in% c("GRE"),
  1,
  2
)

# Prepare dataset for Fst between the two groups at node 2
pop(gen_node2) <- pops_filt_node2$node2

# Fst by marker for node 2
fst_node2 <- wc(gen_node2)

fst_node2_df <- data.frame(
  locus = names(gen_node2@all.names),
  FST = fst_node2$per.loc$FST
)

fst_node2_df <- fst_node2_df[order(fst_node2_df$FST, decreasing = TRUE), ]



## Save results
write.table(fst_node2_df, file = "fst_node2.txt", sep = "\t", quote = F, row.names = F, col.names = F)

