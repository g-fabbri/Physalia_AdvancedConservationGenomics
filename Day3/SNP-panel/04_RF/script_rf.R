library(MASS)
library(randomForest)
library(ade4)
library(adegenet)
library(plyr)

# Set your working directory (or navigate to it in R Studio)
setwd("D:/Side projects/Physalia course/Day3/tutorial_snpchip/04_rf")



## Read in the data
testudo_node1 <- read.structure("../testudo_dataset.stru",
                                n.ind = 70,
                                n.loc = 3182,
                                onerowperind = FALSE,
                                col.lab = 1,
                                col.pop = 2,
                                row.marknames = 1,
                                NA.char = "0")

pops <- read.table("../01_pca/Info_dataset_classification_PCA.txt", header = T)


## Check missingness patterns
# By sample
prop_ind <- propTyped(testudo_node1, by = "ind")

# By marker
prop_loc <- propTyped(testudo_node1, by = "loc")


# Decide threshold for markers (depends on your dataset)
prop_loc_keep <- (1 - length(which(pops$Classification == 4 | pops$Classification == 5))/nrow(pops))

# Filter
gen_filt <- testudo_node1[
  prop_ind >= 0.8,
  loc = prop_loc >= prop_loc_keep
]

pops_filt <- pops[pops$Sample %in% rownames(gen_filt$tab),]



## Format dataset for RF
x.test <- tab(gen_filt, freq = TRUE, NA.method = "asis")
x.test <- x.test[, seq(1, ncol(x.test), by = 2)]
x.test.fixed <- na.roughfix(x.test)

testudo_factor <- as.data.frame(x.test.fixed, row.names = pops_filt$Sample)
testudo_factor$pop <- as.factor(pops_filt$Classification)
testudo_factor <- testudo_factor[, c("pop", setdiff(names(testudo_factor), "pop"))]



## Find the best value for parameter mtry
set.seed(3549)
bestmtry <- tuneRF(testudo_factor[,-1], testudo_factor$pop, ntreeTry = 8000, stepFactor = 1.2, improve = 0.01, trace = T, plot = T)
bestmtry
best_mtry <- bestmtry[which.min(bestmtry[, 2]), "mtry"]



## Run RF
# Store importance results over the 10 runs
n_runs <- 10

importance_list <- vector("list", n_runs)

for (i in 1:n_runs) {
  
  cat("\nRunning RF", i, "of", n_runs, "\n")
  
  # Different seed for each run
  set.seed(3549 + i)
  
  test_rf <- randomForest(
    testudo_factor[,-1],
    y = testudo_factor$pop,
    mtry = best_mtry,
    ntree = 8000,
    importance = TRUE
  )
  
  # Extract importance
  impo <- importance(test_rf, class = NULL)
  
  # Keep only Mean Decrease Accuracy
  importance_list[[i]] <- impo[, "MeanDecreaseAccuracy"]
  
  cat(
    "SNPs with MDA >= 1:",
    sum(impo[, "MeanDecreaseAccuracy"] >= 1),
    "\n"
  )
}


# Combine the 10 runs into one matrix
mda_matrix <- do.call(cbind, importance_list)
colnames(mda_matrix) <- paste0("Run", 1:n_runs)

# Summarize SNP importance across runs
mda_summary <- data.frame(
  SNP = rownames(mda_matrix),
  mda_matrix,
  Mean_MDA = rowMeans(mda_matrix),
  Median_MDA = apply(mda_matrix, 1, median),
  SD_MDA = apply(mda_matrix, 1, sd),
  Runs_MDA_ge1 = rowSums(mda_matrix >= 1),
  All_10_MDA_ge1 = rowSums(mda_matrix >= 1) == n_runs
)

# Rank SNPs according to their mean MDA
mda_summary <- mda_summary[order(-mda_summary$Mean_MDA),]

# Save results of the 10 runs together
write.table(mda_summary, file = "RF_MDA_10runs_summary.txt", sep = "\t", quote = F, row.names = FALSE)
