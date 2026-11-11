# 07b_prior_sensitivity.R: refit M4 with different priors and save the draws

# Turn on the packages for tables, JAGS and MCMC output
library(tidyverse)
library(rjags)
library(coda)

# Load the weekly table and keep only the training years, 2019-2021
train <- readRDS("data/processed/etf_weekly.rds") |>
  filter(set == "train")

# Number the bond types and funds the same way as M4
type_levels <- c("Treasury", "Investment grade", "High yield")
type_id     <- as.integer(factor(train$type, levels = type_levels))   # 1 = Treasury, 2 = Investment grade, 3 = High yield
etf_levels  <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_id      <- as.integer(factor(train$etf, levels = etf_levels))     # 1 = GOVT up to 7 = USHY

# ---- No-pooling version of M4: each fund's VIX slope estimated on its own ----

# Write out M4 with one change: the fund slopes no longer share a bond-type distribution
m4_nopool_code <- "
model {
  for (i in 1:N) {
    premium[i] ~ dt(mu[i], tau, nu)                  # each premium follows a Student-t around its predicted value
    mu[i] <- b0[etf[i]] + b1[etf[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]   # same prediction as M4
  }
  for (j in 1:J) {
    b0[j] ~ dnorm(0, 4)                              # prior for each fund's intercept: centred on 0, SD 0.5
    b1[j] ~ dnorm(0, 100)                            # each fund's slope gets its own prior, centred on 0, SD 0.1: no sharing
  }
  for (k in 1:K) {
    b2[k] ~ dnorm(0.5, 16)                           # prior for each bond type's persistence: same as M4
  }
  tau ~ dgamma(1, 0.25)                              # prior for the precision of the errors: same as M4
  sigma <- 1 / sqrt(tau)                             # turn the precision into a scale so it is easier to read
  nu ~ dgamma(2, 0.1)                                # prior for the degrees of freedom: same as M4
}"
writeLines(m4_nopool_code, "models/m4_nopool.jags")   # save the model as a file in the models folder

# Put the data JAGS needs into a list
nopool_data <- list(premium = train$premium, vix_c = train$vix_c,
                    lag_premium = train$lag_premium, etf = etf_id, type = type_id,
                    N = nrow(train), J = 7, K = 3)

# Give each of the 3 chains a different starting point and a fixed random seed
nopool_inits <- list(
  list(b0 = rep(-1, 7), b1 = rep(-0.05, 7), b2 = rep(0, 3),   tau = 0.5, nu = 3,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(b0 = rep(0, 7),  b1 = rep(0, 7),     b2 = rep(0.5, 3), tau = 4,   nu = 10,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(b0 = rep(1, 7),  b1 = rep(0.05, 7),  b2 = rep(0.9, 3), tau = 20,  nu = 30,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))

# Set up the model, let the chains settle, then keep 5,000 draws per chain
m4_nopool <- jags.model("models/m4_nopool.jags", data = nopool_data,
                        inits = nopool_inits, n.chains = 3)
update(m4_nopool, 2000)                                                  # run 2,000 draws and throw them away
s4_nopool <- coda.samples(m4_nopool, c("b1", "b2", "nu"), n.iter = 5000)  # run 5,000 more draws per chain and keep them

# Save the draws so the shrinkage plot can load them
saveRDS(s4_nopool, "output/m4_nopool_samples.rds")

# Compare the 3 chains for each parameter: values close to 1 mean they agree
gelman.diag(s4_nopool, multivariate = FALSE)

# Count how many effectively independent draws each parameter has
effectiveSize(s4_nopool)

# Keep only the LQD and HYG weeks where the VIX was at least 20 points above its usual level
stress_weeks <- train |> filter(etf %in% c("LQD", "HYG"), vix_c >= 20)
# For each fund, count those weeks, how many had a discount bigger than 2%, and what share that is
stress_weeks |> group_by(etf) |> summarise(weeks = n(), big_discounts = sum(premium < -2), share = mean(premium < -2))
# List those weeks with their date, centred VIX and premium, oldest first
stress_weeks |> select(etf, date, vix_c, premium) |> arrange(etf, date)

# Prior sensitivity refits of M4 ----
# List the 7 funds in the model's order
etf_levels <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
# Give each fund its bond type: 1 = Treasury, 2 = Investment grade, 3 = High yield
etf_type <- c(1, 1, 1, 2, 2, 3, 3)
# Put the training data into the list JAGS needs: number of weeks, premiums, centred VIX, last week's premium, fund number, bond type number, and each fund's type
sens_data <- list(N = nrow(train),
                  premium = train$premium,
                  vix_c = train$vix_c,
                  lag_premium = train$lag_premium,
                  etf = match(train$etf, etf_levels),
                  type = etf_type[match(train$etf, etf_levels)],
                  etf_type = etf_type)

# Write the vague-prior version of M4 to a file in the models folder
writeLines("
model {
  # Go through every training week
  for (i in 1:N) {
    # Say each premium is Student-t around its expected value
    premium[i] ~ dt(mu[i], tau, nu)
    # Work out the expected premium: fund level + fund VIX slope x VIX + type persistence x last week's premium
    mu[i] <- b0[etf[i]] + b1[etf[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]
  }
  # Go through the 7 funds
  for (j in 1:7) {
    # Fund level prior: Normal, mean 0, SD 5 (M4 has SD 0.5)
    b0[j] ~ dnorm(0, 0.04)
    # Fund VIX slope: drawn around its bond type's average slope
    b1[j] ~ dnorm(mu_b[etf_type[j]], tau_b)
  }
  # Go through the 3 bond types
  for (k in 1:3) {
    # Bond-type average VIX slope prior: Normal, mean 0, SD 1 (M4 has SD 0.1)
    mu_b[k] ~ dnorm(0, 1)
    # Persistence prior: Normal, mean 0.5, SD 2.5 (M4 has SD 0.25)
    b2[k] ~ dnorm(0.5, 0.16)
  }
  # Spread of fund slopes within a type: anywhere from 0 to 1 (M4 has 0 to 0.1)
  sigma_b ~ dunif(0, 1)
  # Turn that spread into the precision JAGS uses
  tau_b <- 1 / pow(sigma_b, 2)
  # Noise precision prior: Gamma(0.01, 0.0025), same average as M4's Gamma(1, 0.25) but far more spread
  tau ~ dgamma(0.01, 0.0025)
  # Tail-heaviness prior: same as M4
  nu ~ dgamma(2, 0.1)
}
", "models/m4_sens_vague.jags")

# Write the informative-prior version of M4 to a file in the models folder
writeLines("
model {
  # Go through every training week
  for (i in 1:N) {
    # Say each premium is Student-t around its expected value
    premium[i] ~ dt(mu[i], tau, nu)
    # Work out the expected premium: fund level + fund VIX slope x VIX + type persistence x last week's premium
    mu[i] <- b0[etf[i]] + b1[etf[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]
  }
  # Go through the 7 funds
  for (j in 1:7) {
    # Fund level prior: same as M4
    b0[j] ~ dnorm(0, 4)
    # Fund VIX slope: drawn around its bond type's average slope
    b1[j] ~ dnorm(mu_b[etf_type[j]], tau_b)
  }
  # Go through the 3 bond types
  for (k in 1:3) {
    # Bond-type average VIX slope prior: Normal, mean -0.05, SD 0.05 (the only change from M4)
    mu_b[k] ~ dnorm(-0.05, 400)
    # Persistence prior: same as M4
    b2[k] ~ dnorm(0.5, 16)
  }
  # Spread of fund slopes within a type: same as M4
  sigma_b ~ dunif(0, 0.1)
  # Turn that spread into the precision JAGS uses
  tau_b <- 1 / pow(sigma_b, 2)
  # Noise precision prior: same as M4
  tau ~ dgamma(1, 0.25)
  # Tail-heaviness prior: same as M4
  nu ~ dgamma(2, 0.1)
}
", "models/m4_sens_inform.jags")

# Make starting values for one chain, given a random seed and a starting tail value
make_inits <- function(seed, nu_start) {
  list(b0 = rep(0, 7), b1 = rep(0, 7), mu_b = rep(0, 3), sigma_b = 0.05,
       b2 = rep(0.3, 3), tau = 100, nu = nu_start,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
}
# Make starting values for the 3 chains, each with a different seed and starting tail value
sens_inits <- list(make_inits(1, 3), make_inits(2, 10), make_inits(3, 30))

# Set up the vague-prior model in JAGS with 3 chains
fit_vague <- jags.model("models/m4_sens_vague.jags", data = sens_data, inits = sens_inits, n.chains = 3)
# Run 2,000 warm-up draws and throw them away
update(fit_vague, 2000)
# Run 10,000 more draws per chain and keep them for the listed parameters
s_vague <- coda.samples(fit_vague, c("b1", "b2", "mu_b", "sigma_b", "tau", "nu"), n.iter = 10000)
# Save the draws to a file so we never have to rerun this
saveRDS(s_vague, "output/m4_sens_vague_samples.rds")

# Set up the informative-prior model in JAGS with 3 chains
fit_inform <- jags.model("models/m4_sens_inform.jags", data = sens_data, inits = sens_inits, n.chains = 3)
# Run 2,000 warm-up draws and throw them away
update(fit_inform, 2000)
# Run 10,000 more draws per chain and keep them for the listed parameters
s_inform <- coda.samples(fit_inform, c("b1", "b2", "mu_b", "sigma_b", "tau", "nu"), n.iter = 10000)
# Save the draws to a file so we never have to rerun this
saveRDS(s_inform, "output/m4_sens_inform_samples.rds")

# Prior sensitivity results ----
# Load the data-handling tools (tables, pivoting, text matching)
library(tidyverse)
# Load the saved draws for the vague fit
s_vague <- readRDS("output/m4_sens_vague_samples.rds")
# Load the saved draws for the informative fit
s_inform <- readRDS("output/m4_sens_inform_samples.rds")

# Show the worst (largest) Gelman-Rubin value across all parameters in the vague fit
max(gelman.diag(s_vague, multivariate = FALSE)$psrf[, "Point est."])
# Show the 3 parameters with the smallest effective sample size in the vague fit
sort(effectiveSize(s_vague))[1:3]
# Show the worst (largest) Gelman-Rubin value across all parameters in the informative fit
max(gelman.diag(s_inform, multivariate = FALSE)$psrf[, "Point est."])
# Show the 3 parameters with the smallest effective sample size in the informative fit
sort(effectiveSize(s_inform))[1:3]

# Turn each fit's draws into one big table: one row per draw, one column per parameter
d_m4 <- as.matrix(readRDS("output/m4_samples.rds"))
d_vague <- as.matrix(s_vague)
d_inform <- as.matrix(s_inform)

# Make a summary of the key numbers for one fit: median, 95% interval and chance of being below zero
summarise_fit <- function(d, fit_name) {
  # Use the saved noise size if the fit has one; otherwise work it out from tau
  sigma <- if ("sigma" %in% colnames(d)) d[, "sigma"] else 1 / sqrt(d[, "tau"])
  # Pick out the draws we compare, with readable names; VIX slopes are multiplied by 40 to show the effect of a 40-point VIX rise
  draws <- list("VIX effect x40: Treasury"     = d[, "mu_b[1]"] * 40,
                "VIX effect x40: Inv. grade"   = d[, "mu_b[2]"] * 40,
                "VIX effect x40: High yield"   = d[, "mu_b[3]"] * 40,
                "Persistence: Treasury"        = d[, "b2[1]"],
                "Persistence: Inv. grade"      = d[, "b2[2]"],
                "Persistence: High yield"      = d[, "b2[3]"],
                "Fund spread sigma_b x40"      = d[, "sigma_b"] * 40,
                "Noise size sigma"             = sigma,
                "Tails nu"                     = d[, "nu"])
  # Work out the median, the 95% interval and the chance of being below zero for each one, and put them in a table
  tibble(fit = fit_name, parameter = names(draws),
         median = sapply(draws, median),
         lo95 = sapply(draws, quantile, 0.025),
         hi95 = sapply(draws, quantile, 0.975),
         p_below_0 = sapply(draws, function(x) mean(x < 0)))
}

# Summarise all three fits and stack them into one long table
sens_long <- bind_rows(summarise_fit(d_m4, "M4 (main)"),
                       summarise_fit(d_vague, "Vague"),
                       summarise_fit(d_inform, "Informative"))
# Turn each result into one readable cell like "-0.190 (-0.350 to -0.030)", with one column per fit
sens_table <- sens_long |>
  mutate(cell = sprintf("%.3f (%.3f to %.3f)", median, lo95, hi95)) |>
  select(parameter, fit, cell) |>
  pivot_wider(names_from = fit, values_from = cell)
# Show the comparison table in full
print(as.data.frame(sens_table))
# Show the chance each bond type's VIX effect is below zero, under each prior
sens_long |> filter(str_detect(parameter, "VIX")) |> select(fit, parameter, p_below_0) |>
  pivot_wider(names_from = fit, values_from = p_below_0)
# Save all the numbers as a table file for the slides
write_csv(sens_long, "output/tables/prior_sensitivity.csv")

# Work out the Gelman-Rubin value for every parameter in the vague fit
psrf_vague <- gelman.diag(s_vague, multivariate = FALSE)$psrf[, "Point est."]
# Show the 5 parameters with the largest values
sort(psrf_vague, decreasing = TRUE)[1:5]
# Save the names of the 3 worst parameters
worst <- names(sort(psrf_vague, decreasing = TRUE))[1:3]
# Show the median of those 3 parameters in each chain separately (one column per chain)
sapply(s_vague, function(chain) apply(as.matrix(chain)[, worst, drop = FALSE], 2, median))
# Draw trace plots of those 3 parameters, with the 3 chains overlaid
traceplot(s_vague[, worst])
