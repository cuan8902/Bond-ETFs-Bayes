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

# Normal QQ plot of M3's residuals, to show why we switched to Student-t errors (M4)
train <- readRDS("data/processed/etf_weekly.rds") |> filter(set == "train")
type_levels <- c("Treasury", "Investment grade", "High yield")
etf_levels  <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_id  <- as.integer(factor(train$etf, levels = etf_levels))
type_id <- as.integer(factor(train$type, levels = type_levels))

pm3 <- colMeans(as.matrix(readRDS("output/m3_samples.rds")))   # M3's posterior means (loaded, not refitted)
get_par <- function(names) pm3[[intersect(names, names(pm3))[1]]]
b0_m3 <- sapply(1:7, function(j) get_par(c(paste0("b0[", j, "]"), paste0("beta0[", j, "]"))))
b1_m3 <- sapply(1:7, function(j) get_par(c(paste0("b1[", j, "]"), paste0("beta1[", j, "]"))))
b2_m3 <- sapply(1:3, function(k) get_par(c(paste0("b2[", k, "]"), paste0("beta2[", k, "]"))))
sigma_m3 <- if ("sigma" %in% names(pm3)) pm3[["sigma"]] else 1 / sqrt(pm3[["tau"]])

fitted_m3 <- b0_m3[etf_id] + b1_m3[etf_id] * train$vix_c + b2_m3[type_id] * train$lag_premium
resid_m3 <- tibble(std_resid = (train$premium - fitted_m3) / sigma_m3) |> filter(!is.na(std_resid))

qq_m3 <- ggplot(resid_m3, aes(sample = std_resid)) +
  stat_qq() + stat_qq_line() +
  labs(title = "Normal QQ plot of M3 residuals", x = "Normal values", y = "Standardised residual") +
  theme_minimal()
ggsave("output/figures/fig_m3_resid_qq.png", qq_m3, width = 7, height = 5, dpi = 300, bg = "white")
