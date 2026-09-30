## Attention!!!
## The script below works if you have run script_pca.R first

library(ade4)
library(adegenet)
library(ggplot2)


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
write.table(topload1, file = "01_pca/loadings_PC1.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload2, file = "01_pca/loadings_PC2.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload3, file = "01_pca/loadings_PC3.txt", sep = "\t", quote = F, row.names = T, col.names = F)
write.table(topload4, file = "01_pca/loadings_PC4.txt", sep = "\t", quote = F, row.names = T, col.names = F)

