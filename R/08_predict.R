# 08_predict.R: posterior predictive checks, the 2022 out-of-sample test and the crisis scenario for M4
# Only LOADS saved files (m4_samples.rds and the weekly data); it refits nothing.

# Turn on the packages for tables, plots and MCMC output
library(tidyverse)
library(coda)

set.seed(417)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Setup ----

# Load the weekly table and split it into the training years and 2022
weekly <- readRDS("data/processed/etf_weekly.rds")
train <- weekly |> filter(set == "train")
test  <- weekly |> filter(set == "test")

# Number the bond types and funds the same way as in 06c_fit_m4.R
type_levels <- c("Treasury", "Investment grade", "High yield")
etf_levels  <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY")
etf_type    <- c(1, 1, 1, 2, 2, 3, 3)

# A discount beyond 2% means a premium below -2 (premium is in %)
disc_cut <- -2
summary(train$premium)   # check: values should look like percentages, e.g. -5 to 2

# Load M4's draws, stack the 3 chains, and keep a few thousand evenly spaced draws
m4 <- as.matrix(readRDS("output/m4_samples.rds"))
S <- 4000
m4 <- m4[round(seq(1, nrow(m4), length.out = S)), ]

# Pull out each parameter as a matrix with one row per draw
b0 <- m4[, paste0("b0[", 1:7, "]")]   # each fund's intercept
b1 <- m4[, paste0("b1[", 1:7, "]")]   # each fund's VIX slope
b2 <- m4[, paste0("b2[", 1:3, "]")]   # each bond type's persistence
sigma <- if ("sigma" %in% colnames(m4)) m4[, "sigma"] else 1 / sqrt(m4[, "tau"])   # scale of the Student-t errors
nu <- m4[, "nu"]                      # degrees of freedom

# Simulate premiums from M4 for a set of weeks: one row per draw, one column per week
simulate_m4 <- function(etf, type, vix_c, lag_premium) {
  n <- length(etf)
  mu <- b0[, etf] + b1[, etf] * matrix(vix_c, S, n, byrow = TRUE) +
        b2[, type] * matrix(lag_premium, S, n, byrow = TRUE)          # predicted premium for every draw and week
  mu + sigma * matrix(rt(S * n, df = rep(nu, n)), S, n)               # add Student-t noise using each draw's nu
}

# ---- 1. Posterior predictive checks ----

train_etf  <- as.integer(factor(train$etf, levels = etf_levels))
train_type <- as.integer(factor(train$type, levels = type_levels))

# Simulate a full fake training dataset for every draw, using the real VIX and real lagged premium
fake <- simulate_m4(train_etf, train_type, train$vix_c, train$lag_premium)

# Work out the two statistics for each fake dataset and for the real data
sim_share <- rowMeans(fake < disc_cut)            # share of weeks with a discount beyond 2%
sim_worst <- apply(fake, 1, min)                  # the worst discount
real_share <- mean(train$premium < disc_cut)
real_worst <- min(train$premium)

# Medians and percentiles (not means, because nu is so low), plus where the real value sits
ppc_row <- function(stat, real, sim) {
  tibble(statistic = stat,
         real = round(real, 3),
         sim_5 = round(quantile(sim, 0.05), 3),
         sim_median = round(median(sim), 3),
         sim_95 = round(quantile(sim, 0.95), 3),
         real_percentile = round(100 * mean(sim <= real), 1))   # % of fake datasets at or below the real value
}
ppc_table <- bind_rows(
  ppc_row("Share of weeks with discount beyond 2%", real_share, sim_share),
  ppc_row("Worst discount (%)", real_worst, sim_worst))
print(ppc_table)
write_csv(ppc_table, "output/tables/ppc_table.csv")

# Plot the simulated values with the real value as a red line
ppc_plot <- bind_rows(
  tibble(statistic = "Share of weeks with discount beyond 2%", value = sim_share, real = real_share),
  tibble(statistic = "Worst discount (%)", value = sim_worst, real = real_worst)) |>
  group_by(statistic) |>
  filter(value >= quantile(value, 0.01), value <= quantile(value, 0.99)) |>   # trim the most extreme 1% each side so the plot is readable
  ungroup()
p10 <- ggplot(ppc_plot, aes(value)) +
  geom_histogram(bins = 50, fill = "grey70", colour = "white") +
  geom_vline(aes(xintercept = real), colour = "red", linewidth = 1) +
  facet_wrap(~ statistic, scales = "free") +
  labs(title = "Posterior predictive checks for M4",
       subtitle = "Grey: 4,000 fake datasets from M4. Red: the real data.",
       x = NULL, y = "Number of fake datasets") +
  theme_minimal()
ggsave("output/figures/fig10_ppc.png", p10, width = 9, height = 4, dpi = 300, bg = "white")

# ---- 2. 2022 out-of-sample test ----

# Drop any 2022 weeks without a previous week's premium
test <- test |> filter(!is.na(lag_premium), !is.na(premium))
test_etf  <- as.integer(factor(test$etf, levels = etf_levels))
test_type <- as.integer(factor(test$type, levels = type_levels))

# Predict each week from that week's VIX (already centred, so not re-centred) and last week's ACTUAL premium
pred <- simulate_m4(test_etf, test_type, test$vix_c, test$lag_premium)
test <- test |>
  mutate(pred_median = apply(pred, 2, median),
         lower90 = apply(pred, 2, quantile, 0.05),
         upper90 = apply(pred, 2, quantile, 0.95),
         inside90 = premium >= lower90 & premium <= upper90,
         abs_err_m4 = abs(premium - pred_median),
         abs_err_naive = abs(premium - lag_premium))      # naive benchmark: same as last week

# Coverage and average absolute error, overall and for each bond type
test_summary <- function(d) {
  summarise(d,
            weeks = n(),
            coverage_90 = round(100 * mean(inside90), 1),
            mae_m4 = round(mean(abs_err_m4), 3),
            mae_naive = round(mean(abs_err_naive), 3))
}
test_2022 <- bind_rows(
  test |> test_summary() |> mutate(type = "All", .before = 1),
  test |> group_by(type) |> test_summary() |> ungroup())
print(test_2022)
write_csv(test_2022, "output/tables/test_2022.csv")

# Plot the actual 2022 premiums against M4's 90% intervals, one panel per fund
test <- test |>
  group_by(etf) |> mutate(week = row_number()) |> ungroup() |>
  mutate(etf = factor(etf, levels = etf_levels))
x_var <- if ("date" %in% names(test)) "date" else "week"
p11 <- ggplot(test, aes(.data[[x_var]])) +
  geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = pred_median), colour = "steelblue") +
  geom_point(aes(y = premium, colour = inside90), size = 0.8) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "red"), name = "Inside 90%") +
  facet_wrap(~ etf, scales = "free_y") +
  labs(title = "2022 out-of-sample test",
       subtitle = "Blue band: M4's 90% predictive interval. Points: actual premium.",
       x = NULL, y = "Premium (%)") +
  theme_minimal()
ggsave("output/figures/fig11_test_2022.png", p11, width = 10, height = 6, dpi = 300, bg = "white")

# ---- 3. Crisis scenario ----

# VIX held at 60 for 4 weeks, centred on the training-period mean VIX
vix_mean <- train$vix[1] - train$vix_c[1]
vix_c_crisis <- 60 - vix_mean

# Start from a zero premium and feed each week's simulated premium into the next
crisis_path <- function(etf) {
  type <- etf_type[etf]
  p <- matrix(NA, S, 4)
  prev <- rep(0, S)
  for (w in 1:4) {
    p[, w] <- b0[, etf] + b1[, etf] * vix_c_crisis + b2[, type] * prev + sigma * rt(S, df = nu)
    prev <- p[, w]
  }
  p
}
crisis <- list(LQD = crisis_path(4), HYG = crisis_path(6))

# P(discount beyond 2%) after 4 weeks, plus the 5%, 50% and 95% percentiles
crisis_scenario <- imap_dfr(crisis, function(p, fund) {
  wk4 <- p[, 4]
  tibble(fund = fund,
         type = ifelse(fund == "LQD", "Investment grade", "High yield"),
         prob_discount_week4 = round(mean(wk4 < disc_cut), 3),               # chance of a discount beyond 2% in week 4
         prob_discount_any = round(mean(rowSums(p < disc_cut) > 0), 3),      # chance of one at least once in the 4 weeks
         p5 = round(quantile(wk4, 0.05), 3),
         p50 = round(median(wk4), 3),
         p95 = round(quantile(wk4, 0.95), 3))
})
print(crisis_scenario)
write_csv(crisis_scenario, "output/tables/crisis_scenario.csv")

# Plot the median path and 5-95% band over the 4 weeks for each fund
crisis_plot <- imap_dfr(crisis, function(p, fund) {
  tibble(fund = fund, week = 0:4,
         p5 = c(0, apply(p, 2, quantile, 0.05)),
         p50 = c(0, apply(p, 2, median)),
         p95 = c(0, apply(p, 2, quantile, 0.95)))
})
p12 <- ggplot(crisis_plot, aes(week)) +
  geom_ribbon(aes(ymin = p5, ymax = p95), fill = "firebrick", alpha = 0.25) +
  geom_line(aes(y = p50), colour = "firebrick", linewidth = 1) +
  geom_hline(yintercept = disc_cut, linetype = "dashed") +
  facet_wrap(~ fund) +
  labs(title = "Crisis scenario: VIX held at 60 for 4 weeks",
       subtitle = "Line: median premium. Band: 5% to 95%. Dashed: 2% discount.",
       x = "Week", y = "Premium (%)") +
  theme_minimal()
ggsave("output/figures/fig12_crisis.png", p12, width = 8, height = 4, dpi = 300, bg = "white")

