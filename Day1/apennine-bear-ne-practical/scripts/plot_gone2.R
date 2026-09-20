#!/usr/bin/env Rscript

outdir <- "results/gone2"
args <- commandArgs(trailingOnly = TRUE)

# Defaults can be overridden with:
# Rscript scripts/plot_gone2.R MIN_GENERATION MAX_GENERATION
min_generation <- if (length(args) >= 1) as.numeric(args[1]) else 1
max_generation <- if (length(args) >= 2) as.numeric(args[2]) else 100

if (!is.finite(min_generation) || !is.finite(max_generation) ||
    min_generation < 0 || max_generation <= min_generation) {
  stop("Use a non-negative minimum and a larger maximum generation")
}

read_gone2 <- function(path) {
  if (!file.exists(path)) stop("Missing GONE2 result: ", path)
  x <- read.table(path, header = TRUE)

  if (!"Generation" %in% names(x)) {
    stop("Generation column not found in: ", path)
  }

  ne_columns <- grep("^Ne", names(x), value = TRUE)
  if (!length(ne_columns)) stop("Ne column not found in: ", path)

  preferred <- if ("Ne" %in% ne_columns) "Ne" else ne_columns[1]
  y <- data.frame(Generation = x$Generation, Ne = x[[preferred]])
  y <- y[is.finite(y$Generation) &
           y$Generation >= min_generation &
           y$Generation <= max_generation &
           is.finite(y$Ne) & y$Ne > 0, , drop = FALSE]
  y[order(y$Generation), , drop = FALSE]
}

abb <- read_gone2(file.path(outdir, "ABB_GONE2_Ne"))
sbb <- read_gone2(file.path(outdir, "SBB_GONE2_Ne"))

if (!nrow(abb) || !nrow(sbb)) {
  stop("No positive, finite estimates in generations ",
       min_generation, " to ", max_generation)
}

range_label <- paste0(
  format(min_generation, scientific = FALSE, trim = TRUE), "_",
  format(max_generation, scientific = FALSE, trim = TRUE)
)
output_pdf <- file.path(
  outdir,
  paste0("GONE2_ABB_SBB_generations_", range_label, ".pdf")
)

pdf(output_pdf, width = 7, height = 5)
plot(abb$Generation, abb$Ne,
     type = "l", log = "y", lwd = 2, col = "firebrick",
     xlim = c(min_generation, max_generation),
     ylim = range(c(abb$Ne, sbb$Ne), finite = TRUE),
     xlab = "Generations before present",
     ylab = "Effective population size")
lines(sbb$Generation, sbb$Ne, lwd = 2, col = "steelblue")
legend("topright", legend = c("ABB", "SBB"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")
dev.off()

cat("Plotted generations:", min_generation, "to", max_generation, "\n")
cat("Plot saved to:", output_pdf, "\n")
