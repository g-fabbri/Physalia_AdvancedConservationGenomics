## Attention!!!
## The script below works if you have run script_pca.R first

library(MASS)
library(randomForest)
library(ade4)
library(adegenet)
library(plyr)

# Format dataset
testudo_factor <- as.data.frame(x.tab.filt, row.names = pops_filt$Sample)
testudo_factor$pop <- as.factor(pops_filt$nclust.classification)
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
write.table(mda_summary, file = "04_rf/RF_MDA_10runs_summary.txt", sep = "\t", quote = F, row.names = FALSE)
