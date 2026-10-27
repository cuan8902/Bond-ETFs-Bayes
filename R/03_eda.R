# 03_eda.R: exploratory figures

# Turn on the tidyverse tools for tables and plots
library(tidyverse)

# Load the finished weekly table
model_data <- readRDS("data/processed/etf_weekly.rds") |>
  mutate(type = factor(type, levels = c("Treasury", "Investment grade", "High yield")))   # set the order the bond types appear in

# Plot each fund's weekly premium over time, with one panel per bond type
fig_premium_time <- ggplot(model_data, aes(x = date, y = premium, colour = etf)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2020-03-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.15) +                       # shade March 2020 in grey
  geom_hline(yintercept = 0, linetype = "dashed") +                       # draw a dashed line at zero premium
  geom_line() +                                                           # draw each fund's premium as a line
  facet_wrap(~ type, ncol = 1) +                                          # stack one panel per bond type
  labs(x = NULL, y = "Premium to NAV (%)", colour = "ETF",
       title = "Weekly bond ETF premiums, 2019–2022") +                   # add the axis label, legend title and title
  theme_minimal()                                                         # use a clean, plain style

# Show the figure in the Plots pane
fig_premium_time

# Save the figure as an image for the slides
ggsave("output/figures/fig1_premium_time.png", fig_premium_time,
       width = 9, height = 7, dpi = 300)

# Show GOVT's three highest weekly premiums and their dates
model_data |> filter(etf == "GOVT") |> slice_max(premium, n = 3) |> select(etf, date, close, nav, premium)

# Plot each week's premium against that week's VIX, with one panel per bond type
fig_premium_vix <- ggplot(model_data, aes(x = vix, y = premium)) +
  geom_hline(yintercept = 0, linetype = "dashed") +            # draw a dashed line at zero premium
  geom_point(alpha = 0.4) +                                    # draw one partly see-through dot per fund-week
  geom_smooth(method = "lm", se = FALSE) +                     # add a straight best-fit line through the dots
  facet_wrap(~ type, nrow = 1) +                               # put one panel per bond type side by side
  labs(x = "VIX (index points)", y = "Premium to NAV (%)",
       title = "Weekly premium against market stress, by bond type") +   # add axis labels and a title
  theme_minimal()                                              # use a clean, plain style

# Show the figure in the Plots pane
fig_premium_vix

# Save the figure as an image for the slides
ggsave("output/figures/fig2_premium_vix.png", fig_premium_vix,
       width = 10, height = 4, dpi = 300)

# Label each week as before the Fed's announcement, just after it, or a normal week
model_data <- model_data |>
  mutate(period = case_when(
    date >= as.Date("2020-02-24") & date <  as.Date("2020-03-23") ~ "Crisis, before Fed",
    date >= as.Date("2020-03-23") & date <= as.Date("2020-04-30") ~ "Crisis, after Fed",
    TRUE ~ "Other weeks"))

# Plot each week's premium against that week's VIX, with one panel per bond type
fig_premium_vix <- ggplot(model_data, aes(x = vix, y = premium)) +
  geom_hline(yintercept = 0, linetype = "dashed") +                  # draw a dashed line at zero premium
  geom_point(aes(colour = period), alpha = 0.6) +                    # draw one dot per fund-week, coloured by period
  geom_smooth(method = "lm", se = FALSE, colour = "black") +         # add a straight best-fit line through all the dots
  scale_colour_manual(values = c("Crisis, before Fed" = "red",
                                 "Crisis, after Fed"  = "blue",
                                 "Other weeks"        = "grey60")) + # choose the colour for each period
  facet_wrap(~ type, nrow = 1) +                                     # put one panel per bond type side by side
  labs(x = "VIX (index points)", y = "Premium to NAV (%)", colour = NULL,
       title = "Weekly premium against market stress, by bond type") +   # add axis labels and a title
  theme_minimal() +                                                  # use a clean, plain style
  scale_x_continuous(breaks = c(20, 40, 60)) +                       # show only 20, 40 and 60 on the VIX axis
  theme(legend.position = "bottom")                                  # move the colour key below the chart                                                   # use a clean, plain style

# Show the figure in the Plots pane
fig_premium_vix

# Save the figure as an image for the slides
ggsave("output/figures/fig2_premium_vix.png", fig_premium_vix,
       width = 10, height = 4, dpi = 300)


# Count how many fund-weeks fall at each premium size, with one panel per bond type
fig_hist <- ggplot(model_data, aes(x = premium)) +
  geom_histogram(bins = 60) +                                   # sort the premiums into 60 equal-width bars and count each
  facet_wrap(~ type, nrow = 1, scales = "free_y") +             # put one panel per bond type side by side
  labs(x = "Premium to NAV (%)", y = "Number of fund-weeks",
       title = "How weekly premiums are spread out") +          # add axis labels and a title
  theme_minimal()                                               # use a clean, plain style

# Compare the premiums with what a Normal distribution would give
fig_qq <- ggplot(model_data, aes(sample = premium)) +
  stat_qq(alpha = 0.5) +                                        # plot the sorted premiums against matching Normal values
  stat_qq_line() +                                              # draw the line the dots would follow if premiums were Normal
  facet_wrap(~ type, nrow = 1) +                                # put one panel per bond type side by side
  labs(x = "Normal values", y = "Premium to NAV (%)",
       title = "Normal QQ plot of weekly premiums") +           # add axis labels and a title
  theme_minimal()                                               # use a clean, plain style

# Show both figures in the Plots pane
fig_hist
fig_qq

# Save both figures as images for the slides
ggsave("output/figures/fig3a_premium_histogram.png", fig_hist, width = 10, height = 4, dpi = 300)
ggsave("output/figures/fig3b_premium_qq.png", fig_qq, width = 10, height = 4, dpi = 300)

# Work out how closely each fund's premium matches its own premium 1 to 10 weeks earlier
acf_data <- model_data |>
  arrange(etf, date) |>                                           # sort by fund, then by date
  group_by(etf, type) |>                                          # work on each fund separately
  reframe(lag = 1:10,                                             # look 1 to 10 weeks back
          acf = acf(premium, lag.max = 10, plot = FALSE)$acf[2:11])   # the correlation for each number of weeks back

# Plot the correlations for each fund, with one panel per bond type
fig_acf <- ggplot(acf_data, aes(x = lag, y = acf, colour = etf)) +
  geom_hline(yintercept = 0) +                                    # draw a line at zero correlation
  geom_line() +                                                   # join each fund's correlations with a line
  geom_point() +                                                  # mark each correlation with a dot
  facet_wrap(~ type, nrow = 1) +                                  # put one panel per bond type side by side
  scale_x_continuous(breaks = 1:10) +                             # label every week from 1 to 10
  labs(x = "Weeks back", y = "Correlation with this week's premium", colour = "ETF",
       title = "How long premiums last") +                        # add axis labels and a title
  theme_minimal() +                                               # use a clean, plain style
  theme(legend.position = "bottom")                               # move the colour key below the chart

# Show the figure in the Plots pane
fig_acf

# Save the figure as an image for the slides
ggsave("output/figures/fig4_premium_acf.png", fig_acf, width = 10, height = 4, dpi = 300)

# Show the correlation between each fund's premium and its premium the week before, with and without the crisis weeks
model_data |>
  group_by(etf) |>                                                            # work on each fund separately
  summarise(all_weeks = cor(premium, lag_premium),                            # correlation using every week
            no_crisis = cor(premium[period == "Other weeks"],
                            lag_premium[period == "Other weeks"]))           # correlation leaving out February to April 2020
