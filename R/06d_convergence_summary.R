# 06d_convergence_summary.R: one convergence table for every model, for the slides

# Turn on the packages for tables and MCMC output
library(tidyverse)
library(coda)

# List every model's saved draws
samples_files <- c(
  "Gibbs sampler (M1)"           = "output/m1_gibbs_samples.rds",
  "M1"                           = "output/m1_samples.rds",
  "M2"                           = "output/m2_samples.rds",
  "M3"                           = "output/m3_samples.rds",
  "M4 (main model)"              = "output/m4_samples.rds",
  "Second model"                 = "output/nav_samples.rds",
  "Second model + VIX change"    = "output/nav_robust_samples.rds",
  "Second model, Student-t"      = "output/nav_t_samples.rds",
  "Second model + VIX change, t" = "output/nav_robust_t_samples.rds")

# For each model, find its worst Gelman-Rubin value and smallest effective sample size
conv_summary <- imap_dfr(samples_files, function(file, model) {
  s <- readRDS(file)                                                       # load the saved draws
  gd <- gelman.diag(s, multivariate = FALSE, autoburnin = FALSE)$psrf      # Gelman-Rubin value for each parameter
  ess <- effectiveSize(s)                                                  # effective draws for each parameter
  tibble(model = model,
         parameters = nrow(gd),                                            # how many parameters were checked
         worst_rhat = round(max(gd[, 1]), 3),                              # the largest Gelman-Rubin value
         smallest_ess = round(min(ess)),                                   # the fewest effective draws
         all_converged = max(gd[, 1]) <= 1.02 & min(ess) > 400)            # TRUE if every parameter passes
})

# Show the table and save it for the slides
conv_summary
dir.create("output/tables", showWarnings = FALSE)
write_csv(conv_summary, "output/tables/convergence_summary.csv")

# Save trace plots of the main model's key parameters as one image for the slides
dir.create("output/figures", showWarnings = FALSE)
s4 <- readRDS("output/m4_samples.rds")
png("output/figures/fig_convergence_m4.png", width = 2400, height = 2400, res = 300)
plot(s4[, c("mu_b[1]", "mu_b[2]", "mu_b[3]", "b2[2]", "b2[3]", "nu")], density = FALSE)  # trace plots only
dev.off()