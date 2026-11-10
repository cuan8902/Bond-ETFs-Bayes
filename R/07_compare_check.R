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
