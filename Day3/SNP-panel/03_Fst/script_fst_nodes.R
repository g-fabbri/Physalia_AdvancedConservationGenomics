## Attention!!!
## The script below works if you have run script_pca.R first

library(adegenet)
library(hierfstat)


## Node 1
# Assign populations to node 1
pops <- read.table("Info_dataset_classification_PCA.txt", header = T)
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
write.table(fst_node1_df, file = "03_fst/fst_node1.txt", sep = "\t", quote = F, row.names = F, col.names = F)


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
write.table(fst_node2_df, file = "03_fst/fst_node2.txt", sep = "\t", quote = F, row.names = F, col.names = F)


