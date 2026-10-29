# 04_conjugate.R: calm vs stressed weeks with Normal-Gamma priors (Lectures 7-8)

# Turn on the tidyverse tools for tables and plots
library(tidyverse)

# Load the weekly table and keep only the training years, 2019-2021
model_data <- readRDS("data/processed/etf_weekly.rds") |>
  filter(set == "train")

# Update a Normal-Gamma prior with data, using the formulas from Lecture 8
ng_update <- function(x, mu0, kappa0, alpha0, beta0) {
  n <- length(x)                     # count the observations
  xbar <- mean(x)                    # work out their average
  kappa_n <- kappa0 + n              # add the data's weight to the prior's weight
  c(mu = (kappa0 * mu0 + n * xbar) / kappa_n,  # posterior mean: weighted average of prior mean and data average
    kappa = kappa_n,                           # posterior weight on the mean
    alpha = alpha0 + n / 2,                    # posterior shape for the precision
    beta = beta0 + sum((x - xbar)^2) / 2 +
      kappa0 * n * (xbar - mu0)^2 / (2 * kappa_n))  # posterior rate for the precision
}

# Draw random values of the mean and precision from a Normal-Gamma posterior
ng_sim <- function(p, S) {
  phi <- rgamma(S, shape = p["alpha"], rate = p["beta"])              # draw precisions from the Gamma part
  theta <- rnorm(S, mean = p["mu"], sd = 1 / sqrt(p["kappa"] * phi))  # draw a mean for each precision
  cbind(theta = theta, phi = phi)                                     # put the draws side by side
}

# Split the corporate bond funds' weeks into calm (VIX below 30) and stressed (VIX 30 or above)
corp <- model_data |> filter(type != "Treasury")  # keep only investment-grade and high-yield funds
calm <- corp$premium[corp$vix < 30]               # premiums in calm weeks
stress <- corp$premium[corp$vix >= 30]            # premiums in stressed weeks

# Update the prior with each group's data (prior mean 0, prior guess for the spread of about 0.5 points)
post_calm <- ng_update(calm, mu0 = 0, kappa0 = 0.25, alpha0 = 1, beta0 = 0.25)
post_stress <- ng_update(stress, mu0 = 0, kappa0 = 0.25, alpha0 = 1, beta0 = 0.25)

# Draw 20,000 values from each posterior
set.seed(417)  # make the random draws the same every time
S <- 20000
sims_calm <- ng_sim(post_calm, S)
sims_stress <- ng_sim(post_stress, S)
diff <- sims_stress[, "theta"] - sims_calm[, "theta"]  # stressed-week mean minus calm-week mean

# Results
length(calm); length(stress)       # number of fund-weeks in each group
mean(diff < 0)                     # probability that premiums are lower in stressed weeks
quantile(diff, c(0.025, 0.975))    # 95% credible interval for the difference
quantile(1 / sqrt(sims_calm[, "phi"]), c(0.025, 0.5, 0.975))    # spread (SD) of premiums in calm weeks
quantile(1 / sqrt(sims_stress[, "phi"]), c(0.025, 0.5, 0.975))  # spread (SD) of premiums in stressed weeks

# Plot the posterior of the difference in mean premium
fig_diff <- ggplot(tibble(diff = diff), aes(x = diff)) +
  geom_density() +                                    # draw the smoothed shape of the simulated differences
  geom_vline(xintercept = 0, linetype = "dashed") +   # draw a dashed line at no difference
  labs(x = "Stressed minus calm mean premium (%)", y = "Density",
       title = "Posterior of the difference in mean premium") +  # add axis labels and a title
  theme_minimal()                                     # use a clean, plain style

# Show the figure and save it for the slides
fig_diff
ggsave("output/figures/fig5_calm_vs_stress.png", fig_diff, width = 6, height = 4, dpi = 300)

