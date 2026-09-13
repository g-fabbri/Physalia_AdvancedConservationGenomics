#!/usr/bin/env Rscript

outdir <- "results/msmc2"
mu <- 1.82e-8
generation_time <- 11

scale_msmc <- function(path) {
  if (!file.exists(path)) stop("Missing MSMC2 result: ", path)
  x <- read.table(path, header = TRUE)
  rate_column <- intersect(c("lambda_00", "lambda"), names(x))
  if (!length(rate_column) ||
      !all(c("left_time_boundary", "right_time_boundary") %in% names(x))) {
    stop("Unexpected MSMC2 columns in: ", path)
  }
  midpoint <- sqrt(x$left_time_boundary * x$right_time_boundary)
  y <- data.frame(
    years = midpoint / mu * generation_time,
    Ne = 1 / (2 * mu * x[[rate_column[1]]])
  )
  y[is.finite(y$years) & y$years > 0 &
      is.finite(y$Ne) & y$Ne > 0, , drop = FALSE]
}

abb <- scale_msmc(file.path(outdir, "ABB_4573.final.txt"))
sbb <- scale_msmc(file.path(outdir, "SBB_U1916.final.txt"))
if (!nrow(abb) || !nrow(sbb)) stop("No positive, finite values to plot")

output_pdf <- file.path(outdir, "MSMC2_ABB_4573_SBB_U1916.pdf")
pdf(output_pdf, width = 8, height = 6)
plot(abb$years, abb$Ne,
     type = "s", log = "xy", lwd = 2, col = "firebrick",
     xlim = range(c(abb$years, sbb$years)),
     ylim = range(c(abb$Ne, sbb$Ne)),
     xlab = "Years before present", ylab = "Effective population size")
lines(sbb$years, sbb$Ne, type = "s", lwd = 2, col = "steelblue")
legend("topleft", legend = c("ABB: 4573", "SBB: U1916"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")
dev.off()

cat("Plot saved to:", output_pdf, "\n")
