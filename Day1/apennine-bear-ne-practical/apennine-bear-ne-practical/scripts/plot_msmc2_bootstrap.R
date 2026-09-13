#!/usr/bin/env Rscript

outdir <- "results/msmc2"
bootdir <- file.path(outdir, "bootstrap")
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
abb_files <- Sys.glob(file.path(bootdir, "ABB_*", "*.final.txt"))
sbb_files <- Sys.glob(file.path(bootdir, "SBB_*", "*.final.txt"))
if (!length(abb_files) || !length(sbb_files)) {
  stop("Bootstrap results not found under: ", bootdir)
}
abb_boot <- lapply(abb_files, scale_msmc)
sbb_boot <- lapply(sbb_files, scale_msmc)
all_curves <- c(list(abb, sbb), abb_boot, sbb_boot)
if (any(!vapply(all_curves, nrow, integer(1)))) {
  stop("At least one result has no positive, finite values")
}

output_pdf <- file.path(outdir, "MSMC2_ABB_SBB_bootstrap.pdf")
pdf(output_pdf, width = 8, height = 6)
plot(abb$years, abb$Ne, type = "n", log = "xy",
     xlim = range(unlist(lapply(all_curves, `[[`, "years"))),
     ylim = range(unlist(lapply(all_curves, `[[`, "Ne"))),
     xlab = "Years before present", ylab = "Effective population size")
for (x in abb_boot) {
  lines(x$years, x$Ne, type = "s",
        col = adjustcolor("firebrick", alpha.f = 0.20))
}
for (x in sbb_boot) {
  lines(x$years, x$Ne, type = "s",
        col = adjustcolor("steelblue", alpha.f = 0.20))
}
lines(abb$years, abb$Ne, type = "s", lwd = 3, col = "firebrick")
lines(sbb$years, sbb$Ne, type = "s", lwd = 3, col = "steelblue")
legend("topleft", legend = c("ABB: 4573", "SBB: U1916"),
       col = c("firebrick", "steelblue"), lwd = 3, bty = "n")
dev.off()

cat("ABB bootstrap runs:", length(abb_files), "\n")
cat("SBB bootstrap runs:", length(sbb_files), "\n")
cat("Plot saved to:", output_pdf, "\n")
