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
