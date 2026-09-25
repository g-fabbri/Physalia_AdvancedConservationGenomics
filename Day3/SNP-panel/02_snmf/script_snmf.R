library(LEA)

# Convert from STRUCTURE to Geno file
testudo_node1_geno <- struct2geno(stru_file, ploidy = 2, FORMAT = 2, extra.row = 1, extra.column = 2)

# Load Geno file
geno_file <- "D:/Side projects/Physalia course/Day3/tutorial_snpchip/testudo_dataset_mod.stru.geno"

# Run snmf
testudo.snmf <- snmf(
geno_file,
K = 1:10,
repetitions = 20,
entropy = TRUE,
alpha = 100,
iterations = 200,
project = "new",
CPU = 1,
seed = 2384672
)

# Look at the cross entropy
plot(
testudo.snmf,
cex = 1.2,
pch = 19
)

# Plot K = 3-6
Q3 <- Q(testudo.snmf, K = 3, run = which.min(cross.entropy(testudo.snmf, K = 3)))
Q4 <- Q(testudo.snmf, K = 4, run = which.min(cross.entropy(testudo.snmf, K = 4)))
Q5 <- Q(testudo.snmf, K = 5, run = which.min(cross.entropy(testudo.snmf, K = 5)))

rownames(Q3) <- names(mclust_class)
rownames(Q4) <- names(mclust_class)
rownames(Q5) <- names(mclust_class)

make_Q_df <- function(Qmat, K, mclust_class, sample_order) {
	df <- as.data.frame(Qmat)
	colnames(df) <- paste0("Cluster", 1:K)
	df$sample <- rownames(Qmat)
	df$Mclust <- mclust_class[df$sample]
	df <- df %>%
	mutate(
	sample = factor(sample, levels = sample_order)
	) %>%
	arrange(sample)
	df_long <- df %>%
	tidyr::pivot_longer(
	cols = starts_with("Cluster"),
	names_to = "Cluster",
	values_to = "Ancestry"
	)
	df_long$K <- K
	return(df_long)
}

df3 <- make_Q_df(Q3, 3, mclust_class, sample_order)
df4 <- make_Q_df(Q4, 4, mclust_class, sample_order)
df5 <- make_Q_df(Q5, 5, mclust_class, sample_order)
df <- bind_rows(df3, df4, df5)

# Plot
ggplot(df, aes(x = sample, y = Ancestry, fill = Cluster)) +
geom_bar(stat = "identity", width = 1) +
facet_wrap(~K, ncol = 1, scales = "free_x") +
theme_classic() +
theme(
axis.text.x = element_text(
angle = 90,
hjust = 1,
vjust = 0.5,
size = 6
),
axis.title.x = element_blank(),
panel.spacing = unit(0.8, "lines")
) +
ylab("Ancestry proportion")
