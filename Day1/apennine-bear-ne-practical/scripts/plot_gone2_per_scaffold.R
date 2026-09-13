# Instructor sensitivity analysis: one PDF per scaffold and -u, ABB/SBB overlaid.
# Run: Rscript scripts/plot_gone2_per_scaffold.R
root <- "results/gone2_per_scaffold"
input <- file.path(root, "output")
plots <- file.path(root, "plots")
dir.create(plots, recursive = TRUE, showWarnings = FALSE)

read_ne <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  x <- tryCatch(read.table(path, header = TRUE), error = function(e) NULL)
  if (is.null(x) || !"Generation" %in% names(x)) return(NULL)
  ne_col <- grep("^Ne", names(x), value = TRUE)[1]
  if (is.na(ne_col)) return(NULL)
  data.frame(generation = x$Generation, ne = x[[ne_col]])
}

made <- 0L
for (scaffold in 1:36) for (u in c("002", "005")) {
  prefix <- function(pop) file.path(input, sprintf("%s_s%d_u%s_GONE2_Ne", pop, scaffold, u))
  abb <- read_ne(prefix("ABB"))
  sbb <- read_ne(prefix("SBB"))
  if (is.null(abb) || is.null(sbb)) {
    message("No comparison plot for Scaffold_", scaffold, " u=0.", substring(u, 2), ": missing/invalid output")
    next
  }
  all_ne <- c(abb$ne, sbb$ne)
  positive <- all_ne[is.finite(all_ne) & all_ne > 0]
  if (!length(positive)) next
  filename <- file.path(plots, sprintf("GONE2_ABB_SBB_Scaffold_%d_u%s.pdf", scaffold, u))
  pdf(filename, width = 7, height = 5)
  plot(abb$generation, abb$ne, type = "l", col = "firebrick", lwd = 2,
       log = "y", xlim = rev(range(c(abb$generation, sbb$generation))),
       ylim = range(positive), xlab = "Generations before present",
       ylab = "Effective population size (Ne)",
       main = sprintf("Scaffold_%d; -u 0.%s", scaffold, substring(u, 2)))
  lines(sbb$generation, sbb$ne, col = "steelblue", lwd = 2)
  legend("topright", c("ABB", "SBB"), col = c("firebrick", "steelblue"), lwd = 2)
  dev.off()
  made <- made + 1L
}
message("Created ", made, " comparison PDFs in ", plots, " (maximum 72).")
