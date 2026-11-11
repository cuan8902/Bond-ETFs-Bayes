# 07_compare_check.R: compare the models and check the chosen one

# Turn on the packages for tables, JAGS and MCMC output
library(tidyverse)
library(rjags)
library(coda)

# Load Ryan's saved M4 draws and both models' DIC
s4   <- readRDS("output/m4_samples.rds")
dic3 <- readRDS("output/m3_dic.rds")
dic4 <- readRDS("output/m4_dic.rds")

# Show the 95% intervals for persistence, the bond-type VIX slopes and the degrees of freedom
summary(s4)$quantiles[c("b2[1]", "b2[2]", "b2[3]", "mu_b[1]", "mu_b[2]", "mu_b[3]", "nu"), ]

# Show M4's DIC, then M3's DIC minus M4's: a positive number means M4 fits better
dic4
diffdic(dic3, dic4)

# Compare the 3 chains for each M4 parameter: values close to 1 mean they agree
gelman.diag(s4, multivariate = FALSE)

# Count how many effectively independent draws each M4 parameter has
effectiveSize(s4)

# Make a folder for results tables
dir.create("output/tables", showWarnings = FALSE)

# ---- DIC comparison table ----

# Load all four models' DIC
dics <- list(M1 = readRDS("output/m1_dic.rds"), M2 = readRDS("output/m2_dic.rds"),
             M3 = dic3, M4 = dic4)

# Build a table of fit, complexity and DIC for each model
dic_table <- tibble(
  model         = names(dics),
  description   = c("Pooled: one VIX slope for all funds",
                    "Bond-type slopes and persistence",
                    "Hierarchical: fund slopes grouped by bond type",
                    "M3 with Student-t errors"),
  mean_deviance = map_dbl(dics, ~ sum(.x$deviance)),   # how badly the model fits the data (lower is better)
  penalty       = map_dbl(dics, ~ sum(.x$penalty)),    # the effective number of parameters
  DIC           = mean_deviance + penalty) |>          # fit plus a charge for complexity (lower is better)
  mutate(across(where(is.numeric), ~ round(.x, 1)))    # round the numbers to 1 decimal place

# Show the table and save it for the slides
dic_table
write_csv(dic_table, "output/tables/dic_table.csv")

# ---- Posterior of the bond-type VIX slopes (mu_b) ----

# Stack M4's 3 chains into one table of draws
draws4 <- as.matrix(s4)
type_names <- c("Treasury", "Investment grade", "High yield")
mu_cols    <- c("mu_b[1]", "mu_b[2]", "mu_b[3]")

# Summarise each bond type's VIX slope
mu_b_table <- tibble(
  bond_type     = type_names,
  mean          = colMeans(draws4[, mu_cols]),                       # average slope
  lower_95      = apply(draws4[, mu_cols], 2, quantile, 0.025),      # bottom of the 95% interval
  upper_95      = apply(draws4[, mu_cols], 2, quantile, 0.975),      # top of the 95% interval
  prob_negative = colMeans(draws4[, mu_cols] < 0),                   # share of draws below zero
  effect_vix_40 = mean * 40)                                         # predicted premium change for a 40-point VIX rise

# Show the table and save it for the slides
mu_b_table
write_csv(mu_b_table, "output/tables/mu_b_table.csv")

# Work out high yield's slope minus investment grade's slope for every draw
hy_minus_ig <- draws4[, "mu_b[3]"] - draws4[, "mu_b[2]"]
c(mean = mean(hy_minus_ig),
  quantile(hy_minus_ig, c(0.025, 0.975)),
  prob_hy_lower = mean(hy_minus_ig < 0))                             # share of draws where high yield's slope is lower

# ---- Half-life of a mispricing ----

# Turn each draw of persistence into a half-life in weeks, for the two corporate bond types
half_life <- tibble(
  `Investment grade` = log(0.5) / log(draws4[, "b2[2]"]),
  `High yield`       = log(0.5) / log(draws4[, "b2[3]"]))

# Summarise the half-lives in weeks and in trading days
half_life |>
  pivot_longer(everything(), names_to = "bond_type", values_to = "weeks") |>
  group_by(bond_type) |>
  summarise(median_weeks = median(weeks),
            lower_95     = quantile(weeks, 0.025),
            upper_95     = quantile(weeks, 0.975),
            median_days  = median(weeks) * 5)                        # 5 trading days per week

# ---- Figure: posterior of each bond type's VIX slope ----

# Put the slope draws into one long table, one row per draw per bond type
mu_b_long <- tibble(Treasury           = draws4[, "mu_b[1]"],
                    `Investment grade` = draws4[, "mu_b[2]"],
                    `High yield`       = draws4[, "mu_b[3]"]) |>
  pivot_longer(everything(), names_to = "bond_type", values_to = "slope") |>
  mutate(bond_type = factor(bond_type, levels = type_names))         # set the order the bond types appear in

# Draw the smoothed shape of each bond type's slope draws
fig_mu_b <- ggplot(mu_b_long, aes(x = slope, fill = bond_type)) +
  geom_density(alpha = 0.4) +                                        # one partly see-through curve per bond type
  geom_vline(xintercept = 0, linetype = "dashed") +                  # draw a dashed line at no effect
  labs(x = "Change in premium (%) per 1-point rise in the VIX", y = "Density", fill = NULL,
       title = "Posterior of the stress effect by bond type (M4)") + # add axis labels and a title
  theme_minimal() +                                                  # use a clean, plain style
  theme(legend.position = "bottom")                                  # move the colour key below the chart

# Show the figure and save it for the slides
fig_mu_b
ggsave("output/figures/fig6_mu_b_posteriors.png", fig_mu_b, width = 8, height = 4, dpi = 300)

# ---- Shrinkage plot: each fund's VIX slope next to its bond type's average ----

# Step 1: summarise each fund's VIX slope and each bond type's average slope from M4

# Load M4's draws and stack the 3 chains into one table
draws4 <- as.matrix(readRDS("output/m4_samples.rds"))

# List the funds in the order JAGS numbered them, and each fund's bond type
etf_levels <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_type   <- c(1, 1, 1, 2, 2, 3, 3)                       # 1 = Treasury, 2 = Investment grade, 3 = High yield
type_names <- c("Treasury", "Investment grade", "High yield")

# Name the columns that hold the fund slopes and the bond-type averages
b1_cols <- paste0("b1[", 1:7, "]")                          # b1[1] to b1[7]: one slope per fund
mu_cols <- paste0("mu_b[", 1:3, "]")                        # mu_b[1] to mu_b[3]: one average per bond type

# For each fund, work out the middle value and 95% interval of its slope, for a 40-point VIX rise
fund_slopes <- tibble(
  etf       = etf_levels,
  bond_type = type_names[etf_type],                                  # the fund's bond type, in words
  median    = apply(draws4[, b1_cols], 2, median) * 40,              # middle value of the slope draws
  lo95      = apply(draws4[, b1_cols], 2, quantile, 0.025) * 40,     # bottom of the 95% interval
  hi95      = apply(draws4[, b1_cols], 2, quantile, 0.975) * 40)     # top of the 95% interval

# Do the same for each bond type's average slope
type_slopes <- tibble(
  bond_type = type_names,
  median    = apply(draws4[, mu_cols], 2, median) * 40,              # middle value of the average slope
  lo95      = apply(draws4[, mu_cols], 2, quantile, 0.025) * 40,     # bottom of the 95% interval
  hi95      = apply(draws4[, mu_cols], 2, quantile, 0.975) * 40)     # top of the 95% interval

# Show both tables
fund_slopes
type_slopes

# Show how much fund slopes differ within a bond type, for a 40-point VIX rise
quantile(draws4[, "sigma_b"] * 40, c(0.025, 0.5, 0.975))

# Step 3: compare each fund's slope on its own with its slope in M4, then draw the shrinkage plot

# Load the no-pooling draws and stack the 3 chains into one table
draws_nopool <- as.matrix(readRDS("output/m4_nopool_samples.rds"))

# Put each fund's slope on its own, its slope in M4 and its bond type's average side by side
shrink <- fund_slopes |>
  select(etf, bond_type, hierarchical = median) |>                        # keep M4's middle value for each fund
  mutate(on_its_own   = apply(draws_nopool[, b1_cols], 2, median) * 40,   # middle value when the fund is estimated alone
         type_average = type_slopes$median[match(bond_type, type_slopes$bond_type)],   # the fund's bond-type average from M4
         bond_type    = factor(bond_type, levels = type_names),           # set the order the bond types appear in
         etf          = factor(etf, levels = rev(etf_levels)))            # set the order the funds appear in

# Show the table and save it for the slides
shrink
write_csv(shrink, "output/tables/shrinkage_table.csv")

# Draw an arrow for each fund, from its slope on its own to its slope in M4
fig_shrink <- ggplot(shrink, aes(y = etf)) +
  geom_vline(xintercept = 0, colour = "grey80") +                                       # thin line at no effect
  geom_vline(aes(xintercept = type_average), linetype = "dashed", colour = "grey40") +   # dashed line at the bond type's average
  geom_segment(aes(x = on_its_own, xend = hierarchical, yend = etf),
               arrow = arrow(length = unit(0.15, "cm"))) +                    # arrow showing how far M4 moved the slope
  geom_point(aes(x = on_its_own), shape = 21, fill = "white", size = 3) +     # hollow dot: the fund estimated on its own
  geom_point(aes(x = hierarchical), size = 3) +                               # filled dot: the fund in M4
  facet_grid(bond_type ~ ., scales = "free_y", space = "free_y") +            # one panel per bond type, stacked
  labs(x = "Change in premium (%) for a 40-point rise in the VIX", y = NULL,
       title = "Shrinkage: M4 pulls each fund's slope towards its bond type's average",
       subtitle = "Hollow dot: fund on its own. Filled dot: M4. Dashed line: bond-type average") +   # add labels and a title
  theme_minimal()                                                             # use a clean, plain style

# Show the figure and save it for the slides
fig_shrink
ggsave("output/figures/fig8_shrinkage.png", fig_shrink, width = 8, height = 4.5, dpi = 300)

# Crisis weeks check ----
# Keep only the LQD and HYG weeks where the VIX was at least 20 points above its usual level
stress_weeks <- train |> filter(etf %in% c("LQD", "HYG"), vix_c >= 20)
# For each fund, count those weeks, how many had a discount bigger than 2%, and what share that is
stress_weeks |> group_by(etf) |> summarise(weeks = n(), big_discounts = sum(premium < -2), share = mean(premium < -2))
# List those weeks with their date, centred VIX and premium, oldest first
stress_weeks |> select(etf, date, vix_c, premium) |> arrange(etf, date)

# Save the list of crisis weeks as a table file so we can put it on a slide
write_csv(stress_weeks |> select(etf, date, vix_c, premium) |> arrange(etf, date), "output/tables/crisis_weeks.csv")

# Prior predictive simulation ----
# Fix the random number generator so we get the same numbers every time we run this
set.seed(417)
# Set how many fake datasets to make from the priors
n_sims <- 1000
# Find which fund (1 to 7) each training week belongs to
etf_row <- match(train$etf, etf_levels)
# Find which bond type (1 to 3) each training week belongs to
type_row <- etf_type[etf_row]
# Make an empty grid for the fake premiums: one row per fake dataset, one column per training week
fake <- matrix(NA, nrow = n_sims, ncol = nrow(train))
# Repeat everything inside the curly brackets once for each fake dataset
for (s in 1:n_sims) {
  # Draw a starting level for each of the 7 funds: Normal, mean 0, SD 0.5
  b0 <- rnorm(7, mean = 0, sd = 0.5)
  # Draw an average VIX slope for each of the 3 bond types: Normal, mean 0, SD 0.1
  mu_b <- rnorm(3, mean = 0, sd = 0.1)
  # Draw how far funds of the same type can differ: anywhere between 0 and 0.1
  sigma_b <- runif(1, min = 0, max = 0.1)
  # Draw each fund's VIX slope around its own bond type's average
  b1 <- rnorm(7, mean = mu_b[etf_type], sd = sigma_b)
  # Draw how much of last week's premium carries over, for each bond type: Normal, mean 0.5, SD 0.25
  b2 <- rnorm(3, mean = 0.5, sd = 0.25)
  # Draw the noise precision from Gamma(1, 0.25) and turn it into a typical noise size
  sigma <- 1 / sqrt(rgamma(1, shape = 1, rate = 0.25))
  # Draw how fat the tails are from Gamma(2, 0.1)
  nu <- rgamma(1, shape = 2, rate = 0.1)
  # Work out each week's expected premium from the drawn values, the real VIX and last week's real premium
  mu <- b0[etf_row] + b1[etf_row] * train$vix_c + b2[type_row] * train$lag_premium
  # Add random Student-t noise to get fake premiums, and store them in row s of the grid
  fake[s, ] <- mu + sigma * rt(nrow(train), df = nu)
}

# Show the 1st, 5th, 50th, 95th and 99th percentiles of the real premiums
quantile(train$premium, c(0.01, 0.05, 0.5, 0.95, 0.99))
# Show the same percentiles for all the fake premiums together
quantile(fake, c(0.01, 0.05, 0.5, 0.95, 0.99))
# Show the share of fake premiums bigger than 10% either way
mean(abs(fake) > 10)
# Work out the typical size of the real premiums (the middle value, ignoring plus or minus signs)
real_typical <- median(abs(train$premium))
# Work out the same typical size for each fake dataset
fake_typical <- apply(abs(fake), 1, median)
# Show the real typical size
real_typical
# Show the middle 95% range of the typical size across the 1,000 fake datasets
quantile(fake_typical, c(0.025, 0.5, 0.975))

# Put the 1,000 fake typical sizes into a table so ggplot can use them
prior_check <- tibble(typical = fake_typical)
# Start a chart of the fake typical sizes
fig7 <- ggplot(prior_check, aes(x = typical)) +
  # Draw a histogram of the fake typical sizes using 40 bars
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  # Draw a red vertical line at the real data's typical size
  geom_vline(xintercept = real_typical, colour = "red", linewidth = 1) +
  # Use a log scale on the x axis so small and large sizes both show clearly
  scale_x_log10() +
  # Add the title and axis labels
  labs(title = "Prior predictive check (M4)",
       subtitle = "Typical premium size in 1,000 datasets simulated from the priors; red line = real data",
       x = "Typical premium size, % (median of |premium|, log scale)",
       y = "Number of fake datasets") +
  # Use a clean white background
  theme_minimal()
# Show the chart
fig7
# Save the chart as a picture file for the slides
ggsave("output/figures/fig7_prior_predictive.png", fig7, width = 8, height = 4.5, dpi = 300)

# Put the real and fake percentiles side by side in one table
prior_pred_table <- tibble(percentile = c("1%", "5%", "50%", "95%", "99%"),
                           real = quantile(train$premium, c(0.01, 0.05, 0.5, 0.95, 0.99)),
                           prior_predictive = quantile(fake, c(0.01, 0.05, 0.5, 0.95, 0.99)))
# Save that table as a file for the slides
write_csv(prior_pred_table, "output/tables/prior_predictive_table.csv")

# Prior vs posterior plot ----
# Load M4's saved draws as one big table: one row per draw, one column per parameter
d_m4 <- as.matrix(readRDS("output/m4_samples.rds"))
# List the six quantities to plot: panel name, M4 column, what to multiply by, the prior's mean and SD on that scale, and the plotting range
fig9_info <- tibble(
  label      = c("VIX effect: Treasury", "VIX effect: Inv. grade", "VIX effect: High yield",
                 "Persistence: Treasury", "Persistence: Inv. grade", "Persistence: High yield"),
  column     = c("mu_b[1]", "mu_b[2]", "mu_b[3]", "b2[1]", "b2[2]", "b2[3]"),
  multiply   = c(40, 40, 40, 1, 1, 1),
  prior_mean = c(0, 0, 0, 0.5, 0.5, 0.5),
  prior_sd   = c(0.1 * 40, 0.1 * 40, 0.1 * 40, 0.25, 0.25, 0.25),
  x_from     = c(-1, -1, -1, -0.4, -0.4, -0.4),
  x_to       = c(0.6, 0.6, 0.6, 1.4, 1.4, 1.4))
# Make two empty lists to collect the prior curves and the posterior draws
prior_list <- list()
post_list <- list()
# Repeat the steps inside the curly brackets for each of the six quantities
for (r in 1:nrow(fig9_info)) {
  # Make 400 evenly spaced points across this quantity's plotting range
  x <- seq(fig9_info$x_from[r], fig9_info$x_to[r], length.out = 400)
  # Work out the height of the prior curve at each of those points
  prior_list[[r]] <- tibble(label = fig9_info$label[r], x = x,
                            density = dnorm(x, fig9_info$prior_mean[r], fig9_info$prior_sd[r]))
  # Take this quantity's posterior draws from M4, multiplied onto the right scale
  # Take this quantity's posterior draws from M4, multiplied onto the right scale
  values <- d_m4[, fig9_info$column[r]] * fig9_info$multiply[r]
  # Show what share of draws falls outside the plotting range
  print(mean(values < fig9_info$x_from[r] | values > fig9_info$x_to[r]))
  # Keep only the draws inside the plotting range, so a few extreme draws don't stretch the x axis
  post_list[[r]] <- tibble(label = fig9_info$label[r],
                           value = values[values >= fig9_info$x_from[r] & values <= fig9_info$x_to[r]])
}
# Stack the six prior curves into one table
prior_curves <- bind_rows(prior_list)
# Stack the six sets of posterior draws into one table
post_draws <- bind_rows(post_list)
# Fix the panel order so it matches the list above
prior_curves$label <- factor(prior_curves$label, levels = fig9_info$label)
post_draws$label <- factor(post_draws$label, levels = fig9_info$label)

# Start the chart
fig9 <- ggplot() +
  # Draw each posterior as a shaded curve
  geom_density(data = post_draws, aes(x = value), fill = "steelblue", alpha = 0.5, colour = "steelblue4") +
  # Draw each prior as a dashed line
  geom_line(data = prior_curves, aes(x = x, y = density), linetype = "dashed", colour = "grey30") +
  # Draw a thin grey line at zero
  geom_vline(xintercept = 0, colour = "grey80") +
  # Give each quantity its own panel with its own scales, in 2 rows
  facet_wrap(~ label, nrow = 2, scales = "free") +
  # Add the title and axis labels
  labs(title = "Prior vs posterior (M4)",
       subtitle = "Dashed = prior, shaded = posterior. VIX effects are in % per 40 VIX points",
       x = "Value", y = "Density") +
  # Use a clean white background
  theme_minimal()
# Show the chart
fig9
# Save the chart as a picture file for the slides
ggsave("output/figures/fig9_prior_vs_posterior.png", fig9, width = 9, height = 5, dpi = 300)

# For each quantity, show the posterior's spread, the prior's spread, and how many times narrower the posterior is
post_draws |> group_by(label) |> summarise(posterior_sd = sd(value)) |>
  mutate(prior_sd = fig9_info$prior_sd, times_narrower = prior_sd / posterior_sd)
