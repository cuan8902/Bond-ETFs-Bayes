# 06_fit_jags.R: fit the Bayesian regression models in JAGS

# Turn on the packages for tables, JAGS and MCMC output
library(tidyverse)
library(rjags)
library(coda)

# Load the weekly table and keep only the training years, 2019-2021
train <- readRDS("data/processed/etf_weekly.rds") |>
  filter(set == "train")

# Write out model M1: one intercept and one VIX slope shared by every fund
m1_code <- "
model {
  for (i in 1:N) {
    premium[i] ~ dnorm(mu[i], tau)    # each premium is Normal around its predicted value
    mu[i] <- b0 + b1 * vix_c[i]       # predicted premium: intercept plus slope times centred VIX
  }
  b0  ~ dnorm(0, 4)                   # prior for the intercept: centred on 0, SD 0.5
  b1  ~ dnorm(0, 100)                 # prior for the VIX slope: centred on 0, SD 0.1
  tau ~ dgamma(1, 0.25)               # prior for the precision: guess of an SD around 0.5
  sigma <- 1 / sqrt(tau)              # turn the precision into an SD so it is easier to read
}"
writeLines(m1_code, "models/m1.jags")   # save the model as a file in the models folder

# Put the data JAGS needs into a list
m1_data <- list(premium = train$premium, vix_c = train$vix_c, N = nrow(train))

# Give each of the 3 chains a different starting point and a fixed random seed
m1_inits <- list(
  list(b0 = -1, b1 = -0.1, tau = 0.5, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(b0 =  0, b1 =  0,   tau = 4,   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(b0 =  1, b1 =  0.1, tau = 20,  .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))

# Set up M1 in JAGS with the data and starting points
m1 <- jags.model("models/m1.jags", data = m1_data, inits = m1_inits, n.chains = 3)

# Run 2,000 draws to let the chains settle, then discard them
update(m1, 2000)

# Run 10,000 more draws per chain and keep them
s1 <- coda.samples(m1, c("b0", "b1", "sigma"), n.iter = 10000)

# Save the draws so later scripts can use them
saveRDS(s1, "output/m1_samples.rds")

# Show the posterior mean, SD and 95% interval for each parameter
summary(s1)

# Compare the 3 chains: values close to 1 mean they agree
gelman.diag(s1)

# Count how many effectively independent draws each parameter has
effectiveSize(s1)

# Draw the trace plots and posterior densities
plot(s1)

# Load the hand-written Gibbs sampler's draws
g1 <- readRDS("output/m1_gibbs_samples.rds")

# Put the Gibbs and JAGS results side by side
pars <- c("b0", "b1", "sigma")
tibble(parameter  = pars,
       gibbs_mean = colMeans(as.matrix(g1))[pars],        # average of the Gibbs draws
       jags_mean  = colMeans(as.matrix(s1))[pars],        # average of the JAGS draws
       gibbs_sd   = apply(as.matrix(g1)[, pars], 2, sd),  # spread of the Gibbs draws
       jags_sd    = apply(as.matrix(s1)[, pars], 2, sd))  # spread of the JAGS draws

# Work out the DIC for M1 so the models can be compared later
dic1 <- dic.samples(m1, n.iter = 5000)
saveRDS(dic1, "output/m1_dic.rds")
dic1

# Write out model M2: an intercept, VIX slope and persistence for each bond type
m2_code <- "
model {
  for (i in 1:N) {
    premium[i] ~ dnorm(mu[i], tau)    # each premium is Normal around its predicted value
    mu[i] <- b0[type[i]] + b1[type[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]   # prediction uses the fund's bond type
  }
  for (k in 1:K) {
    b0[k] ~ dnorm(0, 4)               # prior for each bond type's intercept: centred on 0, SD 0.5
    b1[k] ~ dnorm(0, 100)             # prior for each bond type's VIX slope: centred on 0, SD 0.1
    b2[k] ~ dnorm(0.5, 16)            # prior for each bond type's persistence: centred on 0.5, SD 0.25
  }
  tau ~ dgamma(1, 0.25)               # prior for the precision: guess of an SD around 0.5
  sigma <- 1 / sqrt(tau)              # turn the precision into an SD so it is easier to read
}"
writeLines(m2_code, "models/m2.jags")   # save the model as a file in the models folder

# Number the bond types: 1 = Treasury, 2 = Investment grade, 3 = High yield
type_levels <- c("Treasury", "Investment grade", "High yield")
type_id <- as.integer(factor(train$type, levels = type_levels))

# Put the data JAGS needs into a list
m2_data <- list(premium = train$premium, vix_c = train$vix_c,
                lag_premium = train$lag_premium, type = type_id,
                N = nrow(train), K = 3)

# Give each of the 3 chains a different starting point, with one value per bond type
m2_inits <- list(
  list(b0 = rep(-1, 3), b1 = rep(-0.1, 3), b2 = rep(0, 3),   tau = 0.5, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(b0 = rep(0, 3),  b1 = rep(0, 3),    b2 = rep(0.5, 3), tau = 4,   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(b0 = rep(1, 3),  b1 = rep(0.1, 3),  b2 = rep(0.9, 3), tau = 20,  .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))

# Set up M2 in JAGS, let the chains settle, then keep 10,000 draws per chain
m2 <- jags.model("models/m2.jags", data = m2_data, inits = m2_inits, n.chains = 3)
update(m2, 2000)
s2 <- coda.samples(m2, c("b0", "b1", "b2", "sigma"), n.iter = 10000)

# Save the draws so later scripts can use them
saveRDS(s2, "output/m2_samples.rds")

# Show the posterior mean, SD and 95% interval for each parameter
summary(s2)

# Compare the 3 chains: values close to 1 mean they agree
gelman.diag(s2)

# Count how many effectively independent draws each parameter has
effectiveSize(s2)

# Draw the trace plots and posterior densities
plot(s2)

# Work out the DIC for M2 and compare it with M1
dic2 <- dic.samples(m2, n.iter = 5000)
saveRDS(dic2, "output/m2_dic.rds")
dic2
diffdic(dic1, dic2)   # M1's DIC minus M2's: a positive number means M2 fits better

# Write out model M3: each fund gets its own intercept and VIX slope, with slopes grouped by bond type
m3_code <- "
model {
  for (i in 1:N) {
    premium[i] ~ dnorm(mu[i], tau)    # each premium is Normal around its predicted value
    mu[i] <- b0[etf[i]] + b1[etf[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]   # prediction uses the fund and its bond type
  }
  for (j in 1:J) {
    b0[j] ~ dnorm(0, 4)                       # prior for each fund's intercept: centred on 0, SD 0.5
    b1[j] ~ dnorm(mu_b[etf_type[j]], tau_b)   # each fund's VIX slope is drawn from its bond type's distribution
  }
  for (k in 1:K) {
    mu_b[k] ~ dnorm(0, 100)           # prior for each bond type's average VIX slope: centred on 0, SD 0.1
    b2[k]   ~ dnorm(0.5, 16)          # prior for each bond type's persistence: centred on 0.5, SD 0.25
  }
  sigma_b ~ dunif(0, 0.1)             # prior for how much fund slopes differ within a bond type
  tau_b <- 1 / pow(sigma_b, 2)        # turn that SD into a precision for JAGS
  tau ~ dgamma(1, 0.25)               # prior for the precision: guess of an SD around 0.5
  sigma <- 1 / sqrt(tau)              # turn the precision into an SD so it is easier to read
}"
writeLines(m3_code, "models/m3.jags")   # save the model as a file in the models folder

# Number the funds 1 to 7, and record each fund's bond type (1 = Treasury, 2 = Investment grade, 3 = High yield)
etf_levels <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_id     <- as.integer(factor(train$etf, levels = etf_levels))
etf_type   <- c(1, 1, 1, 2, 2, 3, 3)

# Put the data JAGS needs into a list
m3_data <- list(premium = train$premium, vix_c = train$vix_c,
                lag_premium = train$lag_premium, etf = etf_id, type = type_id,
                etf_type = etf_type, N = nrow(train), J = 7, K = 3)

# Give each of the 3 chains a different starting point
m3_inits <- list(
  list(b0 = rep(-1, 7), b1 = rep(-0.05, 7), mu_b = rep(-0.05, 3), b2 = rep(0, 3),   sigma_b = 0.01, tau = 0.5,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(b0 = rep(0, 7),  b1 = rep(0, 7),     mu_b = rep(0, 3),     b2 = rep(0.5, 3), sigma_b = 0.05, tau = 4,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(b0 = rep(1, 7),  b1 = rep(0.05, 7),  mu_b = rep(0.05, 3),  b2 = rep(0.9, 3), sigma_b = 0.09, tau = 20,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))

# Set up M3 in JAGS, let the chains settle, then keep 20,000 draws per chain
m3 <- jags.model("models/m3.jags", data = m3_data, inits = m3_inits, n.chains = 3)
update(m3, 5000)
# Run M3 for longer so the within-type spread (sigma_b) settles properly
s3 <- coda.samples(m3, c("mu_b", "b1", "b2", "b0", "sigma_b", "sigma"),
                   n.iter = 60000, thin = 3)          # run 60,000 more draws per chain and keep every 3rd one

# Save the longer run over the earlier one
saveRDS(s3, "output/m3_samples.rds")

# Show the posterior mean, SD and 95% interval for each parameter
summary(s3)

# Compare the 3 chains for each parameter: values close to 1 mean they agree
gelman.diag(s3, multivariate = FALSE)

# Count how many effectively independent draws each parameter has
effectiveSize(s3)
