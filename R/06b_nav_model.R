# 06b_nav_model.R: second model (which side closes the gap), robustness check and convergence checks

# Turn on the tools for JAGS, MCMC output and tables
library(rjags)
library(coda)
library(tidyverse)

set.seed(417)  # make the random draws the same every time

# Load the weekly table, add next week's VIX change and keep only the training years, 2019-2021
model_data <- readRDS("data/processed/etf_weekly.rds") |>
  arrange(etf, date) |>                                  # sort by fund, then by date
  group_by(etf) |>                                       # work on each fund separately
  mutate(vix_change_next = lead(vix) - vix) |>           # change in the VIX over the same week as nav_change_next
  ungroup() |>
  filter(set == "train")

types <- c("Treasury", "Investment grade", "High yield")  # order the bond types appear in
B <- length(types)

# Compile a JAGS model, run the burn-in and draw from the posterior
fit_jags <- function(model_string, data, inits, params, burn = 5000, iter = 10000) {
  m <- jags.model(textConnection(model_string), data = data, inits = inits,
                  n.chains = length(inits), quiet = TRUE)
  update(m, burn, progress.bar = "none")                 # throw away the burn-in draws
  s <- coda.samples(m, params, n.iter = iter, progress.bar = "none")
  list(model = m, samples = s)
}

# Fixed random-number seed for each chain so results can be reproduced
chain_seed <- function(k) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 100 + k)

# Posterior mean, 95% interval and P(gamma > 0) for each bond type
gamma_table <- function(samples) {
  g <- as.matrix(samples)[, paste0("gamma[", 1:B, "]"), drop = FALSE]
  tibble(type     = types,
         mean     = round(colMeans(g), 3),
         lower    = round(apply(g, 2, quantile, 0.025), 3),
         upper    = round(apply(g, 2, quantile, 0.975), 3),
         prob_pos = round(colMeans(g > 0), 3))
}


# ---- Second model: next week's NAV change on this week's premium ----
# gamma near 1: the NAV catches up to the price (NAV was stale)
# gamma near 0: the price moves back to the NAV (price was wrong)

nav_data <- model_data |> filter(!is.na(nav_change_next), !is.na(premium))

data_nav <- list(
  dnav    = nav_data$nav_change_next,
  premium = nav_data$premium,
  type    = as.integer(factor(nav_data$type, levels = types)),
  N = nrow(nav_data), B = B
)

nav_model <- "
model {
  for (n in 1:N) {
    dnav[n] ~ dnorm(mu[n], tau_N)
    mu[n] <- alpha[type[n]] + gamma[type[n]] * premium[n]
  }
  for (b in 1:B) {
    alpha[b] ~ dnorm(0, 4)
    gamma[b] ~ dnorm(0, 4)
  }
  tau_N ~ dgamma(1, 1)
  sigma_N <- 1 / sqrt(tau_N)
}"

# Three chains started from spread-out values
inits_nav <- lapply(1:3, function(k) c(list(
  alpha = rnorm(B, 0, 1), gamma = rnorm(B, 0, 1), tau_N = runif(1, 0.2, 5)
), chain_seed(k)))

fit_nav <- fit_jags(nav_model, data_nav, inits_nav, c("alpha", "gamma", "sigma_N"))
print(summary(fit_nav$samples))
saveRDS(fit_nav$samples, "output/nav_samples.rds")          # save the draws for the convergence summary

gamma_nav <- gamma_table(fit_nav$samples)
print(gamma_nav)


# ---- Robustness check: add the change in VIX to the second model ----
rob_data <- nav_data |> filter(!is.na(vix_change_next))

data_rob <- list(
  dnav       = rob_data$nav_change_next,
  premium    = rob_data$premium,
  vix_change = rob_data$vix_change_next,
  type       = as.integer(factor(rob_data$type, levels = types)),
  N = nrow(rob_data), B = B
)

rob_model <- "
model {
  for (n in 1:N) {
    dnav[n] ~ dnorm(mu[n], tau_N)
    mu[n] <- alpha[type[n]] + gamma[type[n]] * premium[n] + delta[type[n]] * vix_change[n]
  }
  for (b in 1:B) {
    alpha[b] ~ dnorm(0, 4)
    gamma[b] ~ dnorm(0, 4)
    delta[b] ~ dnorm(0, 4)
  }
  tau_N ~ dgamma(1, 1)
  sigma_N <- 1 / sqrt(tau_N)
}"

inits_rob <- lapply(1:3, function(k) c(list(
  alpha = rnorm(B, 0, 1), gamma = rnorm(B, 0, 1), delta = rnorm(B, 0, 1),
  tau_N = runif(1, 0.2, 5)
), chain_seed(k)))

fit_rob <- fit_jags(rob_model, data_rob, inits_rob, c("alpha", "gamma", "delta", "sigma_N"))
print(summary(fit_rob$samples))
saveRDS(fit_rob$samples, "output/nav_robust_samples.rds")   # save the draws for the convergence summary

gamma_rob <- gamma_table(fit_rob$samples)

# Compare gamma with and without the VIX change
gamma_compare <- tibble(
  type       = types,
  base       = gamma_nav$mean,
  base_95    = paste0("(", gamma_nav$lower, ", ", gamma_nav$upper, ")"),
  with_vix   = gamma_rob$mean,
  with_vix_95 = paste0("(", gamma_rob$lower, ", ", gamma_rob$upper, ")")
)
print(gamma_compare)


# ---- Convergence checks: trace plots, Gelman-Rubin, effective sample size ----
check_convergence <- function(samples, label) {
  pdf(paste0("output/figures/trace_", label, ".pdf"))    # trace and density plots
  plot(samples)
  dev.off()
  
  pdf(paste0("output/figures/gelman_", label, ".pdf"))   # Gelman-Rubin shrink plots
  gelman.plot(samples, autoburnin = FALSE)
  dev.off()
  
  gd  <- gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE)$psrf
  ess <- effectiveSize(samples)[rownames(gd)]
  
  tibble(parameter = rownames(gd),
         rhat      = round(gd[, 1], 3),
         rhat_upr  = round(gd[, 2], 3),
         ess       = round(ess),
         converged = gd[, 1] <= 1.02 & ess > 400)         # Gelman-Rubin 1.02 or below and at least 400 effective draws
}

conv_nav <- check_convergence(fit_nav$samples, "nav_model")
conv_rob <- check_convergence(fit_rob$samples, "robustness")
print(conv_nav)
print(conv_rob)

