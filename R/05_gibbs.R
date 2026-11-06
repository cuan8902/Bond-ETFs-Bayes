# Turn on the packages for tables and MCMC output
library(tidyverse)
library(coda)

# Load the weekly table and keep only the training years, 2019-2021
train <- readRDS("data/processed/etf_weekly.rds") |>
  filter(set == "train")

# Set up the data: the premiums, and a table with a column of 1s and the centred VIX
y <- train$premium
X <- cbind(intercept = 1, vix_c = train$vix_c)

# Gibbs sampler for M1, using the same priors as the proposal
gibbs_lm <- function(y, X, n_iter, tau_start,
                     m0 = c(0, 0),              # prior means for the intercept and slope
                     P0 = diag(c(4, 100)),      # prior precisions: SD 0.5 for the intercept, SD 0.1 for the slope
                     a = 1, b = 0.25) {         # Gamma prior for the precision
  n   <- length(y)                                      # count the observations
  XtX <- t(X) %*% X                                     # work out X'X once, since it never changes
  Xty <- t(X) %*% y                                     # work out X'y once, since it never changes
  draws <- matrix(NA, n_iter, 4,
                  dimnames = list(NULL, c("b0", "b1", "tau", "sigma")))   # empty table to store every draw
  tau <- tau_start                                      # starting value for the precision
  for (s in 1:n_iter) {
    V    <- solve(tau * XtX + P0)                       # spread of the coefficients, given the current precision
    m    <- V %*% (tau * Xty + P0 %*% m0)               # centre of the coefficients, given the current precision
    beta <- as.vector(m + t(chol(V)) %*% rnorm(2))      # draw both coefficients from that Normal distribution
    ssr  <- sum((y - X %*% beta)^2)                     # add up the squared errors with the new coefficients
    tau  <- rgamma(1, shape = a + n / 2, rate = b + ssr / 2)   # draw the precision from its Gamma distribution
    draws[s, ] <- c(beta, tau, 1 / sqrt(tau))           # store this round's draws, plus the SD
  }
  draws                                                 # hand back the full table of draws
}

# Run 3 chains from different starting precisions, 12,000 draws each
set.seed(417)                                           # make the random draws the same every time
chains <- lapply(c(0.5, 4, 20), function(t0) gibbs_lm(y, X, n_iter = 12000, tau_start = t0))

# Drop the first 2,000 draws of each chain and store the rest as MCMC output
g1 <- mcmc.list(lapply(chains, function(ch) mcmc(ch[-(1:2000), ])))

# Save the draws so they can be compared with JAGS later
saveRDS(g1, "output/m1_gibbs_samples.rds")

# Show the posterior mean, SD and 95% interval for each parameter
summary(g1)

# Compare the 3 chains: values close to 1 mean they agree
gelman.diag(g1)

# Count how many effectively independent draws each parameter has
effectiveSize(g1)

# Draw the trace plots and posterior densities
plot(g1)

