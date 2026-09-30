library(adegenet)
library(recipes)
library(assignPOP)

## Load datasets
panel1 <- read.Genepop("05_test/panel1_pca_dataset_5pop.gen", pop.names = c("GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))
panel2 <- read.Genepop("05_test/panel2_fst_dataset_5pop.gen", pop.names = c("GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))
panel3 <- read.Genepop("05_test/panel3_rf_dataset_5pop.gen", pop.names = c("GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))


# Put together
panels <- list(
  panel1 = panel1,
  panel2 = panel2,
  panel3 = panel3
)

# Self-assignment with Monte Carlo cross-validation on all three panels
for (i in seq_along(panels)) {
  
  panel_name <- names(panels)[i]
  
  # Cross-validation
  assign.MC(
    panels[[i]],
    train.inds = 0.7,
    train.loci = 1,
    loci.sample = "random",
    iterations = 50,
    dir = paste0("05_test/mc_train_svm_", panel_name, "/"),
    model = "svm",
    svm.kernel = "linear",
    multiprocess = F
  )
  
  # Summarize reps
  accuMC <- accuracy.MC(dir = paste0("05_test/mc_train_svm_", panel_name, "/"))
  assign(paste0("accuMC_", panel_name), accuMC)

}

# Plot assignment accuracy across reps
accuracy.plot(accuMC_panel1, pop=c("all", "GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))
accuracy.plot(accuMC_panel2, pop=c("all", "GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))
accuracy.plot(accuMC_panel3, pop=c("all", "GRE", "MEC", "ITP", "CCA", "SCA+SIC+SAR"))

# To visualize mean and standard deviaton of assignments
assign.matrix( dir="05_test/mc_train_svm_panel1/", train.inds=c(0.7), train.loci=c(1))



#############################################
# Self-assignment with KFOLD cross-validation
#assign.kfold(genin, k.fold=c(2, 3), train.loci=1, loci.sample="random", dir="kfold_train_lda/", model="lda", multiprocess = F)
#accuKF <- accuracy.kfold(dir = "kfold_train_lda/")
#membership.plot(dir = "kfold_train_lda/")
