dat <- read.table("comparison_assignment_panels.txt",
                  header = TRUE)

# Convert to long format
library(tidyr)
library(dplyr)

dat_long <- dat %>%
  pivot_longer(
    cols = starts_with("panel"),
    names_to = "panel",
    values_to = "accuracy"
  )

# Make panel a factor
dat_long$panel <- factor(dat_long$panel)

# Repeated-measures ANOVA
anova_model <- aov(accuracy ~ panel + Error(rep/panel),
                   data = dat_long)

summary(anova_model)

##
pairwise.wilcox.test(
  dat_long$accuracy,
  dat_long$panel,
  paired = TRUE,
  p.adjust.method = "holm"
)
