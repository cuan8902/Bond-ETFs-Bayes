# 06c_fit_m4.R: model M4, which is M3 with Student-t errors

# Turn on the packages for tables, JAGS and MCMC output
library(tidyverse)
library(rjags)
library(coda)

# Load the weekly table and keep only the training years, 2019-2021
train <- readRDS("data/processed/etf_weekly.rds") |>
  filter(set == "train")

# Number the bond types: 1 = Treasury, 2 = Investment grade, 3 = High yield
type_levels <- c("Treasury", "Investment grade", "High yield")
type_id <- as.integer(factor(train$type, levels = type_levels))

# Number the funds 1 to 7, and record each fund's bond type
etf_levels <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_id <- as.integer(factor(train$etf, levels = etf_levels))
etf_type <- c(1, 1, 1, 2, 2, 3, 3)

# Write out model M4: the same as M3, but with Student-t errors instead of Normal ones
m4_code <- "
model {
  for (i in 1:N) {
    premium[i] ~ dt(mu[i], tau, nu)                                              # each premium follows a Student-t around its predicted value
    mu[i] <- b0[etf[i]] + b1[etf[i]] * vix_c[i] + b2[type[i]] * lag_premium[i]  # prediction uses the fund and its bond type
  }
  for (j in 1:J) {
    b0[j] ~ dnorm(0, 4)                        # prior for each fund's intercept: centred on 0, SD 0.5
    b1[j] ~ dnorm(mu_b[etf_type[j]], tau_b)    # each fund's VIX slope is drawn from its bond type's distribution
  }
  for (k in 1:K) {
    mu_b[k] ~ dnorm(0, 100)                    # prior for each bond type's average VIX slope: centred on 0, SD 0.1
    b2[k] ~ dnorm(0.5, 16)                     # prior for each bond type's persistence: centred on 0.5, SD 0.25
  }
  sigma_b ~ dunif(0, 0.1)                      # prior for how much fund slopes differ within a bond type
  tau_b <- 1 / pow(sigma_b, 2)                 # turn that SD into a precision for JAGS
  tau ~ dgamma(1, 0.25)                        # prior for the precision of the errors
  sigma <- 1 / sqrt(tau)                       # turn the precision into a scale so it is easier to read
  nu ~ dgamma(2, 0.1)                          # prior for the degrees of freedom: small values mean fat tails
}"
dir.create("models", showWarnings = FALSE)     # make the models folder if it isn't there yet
dir.create("output", showWarnings = FALSE)     # make the output folder if it isn't there yet
writeLines(m4_code, "models/m4.jags")          # save the model as a file in the models folder

# Put the data JAGS needs into a list
m4_data <- list(premium = train$premium, vix_c = train$vix_c,
                lag_premium = train$lag_premium, etf = etf_id, type = type_id,
                etf_type = etf_type, N = nrow(train), J = 7, K = 3)

# Give each of the 3 chains a different starting point
m4_inits <- list(
  list(b0 = rep(-1, 7), b1 = rep(-0.05, 7), mu_b = rep(-0.05, 3), b2 = rep(0, 3), sigma_b = 0.01, tau = 0.5, nu = 3,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(b0 = rep(0, 7), b1 = rep(0, 7), mu_b = rep(0, 3), b2 = rep(0.5, 3), sigma_b = 0.05, tau = 4, nu = 10,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(b0 = rep(1, 7), b1 = rep(0.05, 7), mu_b = rep(0.05, 3), b2 = rep(0.9, 3), sigma_b = 0.09, tau = 20, nu = 30,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))

# Set up M4 in JAGS, let the chains settle, then run 60,000 draws per chain and keep every 3rd
m4 <- jags.model("models/m4.jags", data = m4_data, inits = m4_inits, n.chains = 3)
update(m4, 5000)
s4 <- coda.samples(m4, c("mu_b", "b1", "b2", "b0", "sigma_b", "sigma", "nu"),
                   n.iter = 60000, thin = 3)

# Save the draws so later scripts can use them
saveRDS(s4, "output/m4_samples.rds")

# Show the posterior mean, SD and 95% interval for each parameter
summary(s4)

# Compare the 3 chains for each parameter: values close to 1 mean they agree
gelman.diag(s4, multivariate = FALSE)

# Count how many effectively independent draws each parameter has
effectiveSize(s4)

# Draw trace plots and densities for the key parameters only
plot(s4[, c("b2[1]", "b2[2]", "b2[3]", "mu_b[3]", "sigma_b", "nu")])

# Work out the DIC for M4 and compare it with M3's saved DIC
dic4 <- dic.samples(m4, n.iter = 5000)
saveRDS(dic4, "output/m4_dic.rds")
dic4
dic3 <- readRDS("output/m3_dic.rds")           # load M3's DIC from file
diffdic(dic3, dic4)                            # M3's DIC minus M4's: a positive number means M4 fits better

