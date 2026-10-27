# 02_clean.R: build the weekly dataset

# Turn on the tidyverse tools for working with tables
library(tidyverse)

# Load the saved tables of daily closing prices, NAVs and VIX values
prices <- readRDS("data/raw/prices_yahoo.rds")
navs   <- readRDS("data/raw/navs_ishares.rds")
vix    <- readRDS("data/raw/vix_yahoo.rds") |>
  filter(!is.na(vix))                          # keep only the days that have a VIX value

# Put each day's price and NAV side by side and work out the premium
daily <- prices |>
  inner_join(navs, by = c("etf", "date")) |>   # match rows with the same fund and the same date
  mutate(premium = 100 * (close - nav) / nav)  # add a column: how far the price is above or below the NAV, in %

# Checks
summary(daily$premium)                                        # show the smallest, largest and typical premium
daily |> filter(etf == "LQD", date == as.Date("2019-01-02"))  # show LQD's row for 2 January 2019
daily |> group_by(etf) |> slice_min(premium, n = 1)           # show each fund's lowest premium and the day it happened

# Look more closely at the unusual results
daily |> slice_max(premium, n = 5)                            # show the five highest premiums across all funds
daily |> filter(etf == "GOVT",
                date >= as.Date("2021-08-20"),
                date <= as.Date("2021-09-02"))                # show GOVT's rows from 20 August to 2 September 2021
daily |> filter(etf == "HYG", month(date) == 3, year(date) == 2020) |>
  print(n = 25)                                               # show all of HYG's rows for March 2020

# Keep only the last trading day of each week for each fund
weekly <- daily |>
  mutate(week = floor_date(date, "week", week_start = 1)) |>   # label each day with the Monday its week starts on
  group_by(etf, week) |>                                         # treat each fund's week as its own group
  slice_max(date, n = 1) |>                                      # keep the latest day in each group
  ungroup()                                                      # turn the groups back into one ordinary table

# Count how many weeks each fund has
count(weekly, etf)

# Show GOVT's weekly rows for August 2021
weekly |> filter(etf == "GOVT", month(date) == 8, year(date) == 2021)

# Show HYG's weekly rows for early April 2020
weekly |> filter(etf == "HYG", date >= as.Date("2020-04-01"), date <= as.Date("2020-04-20"))

# Add the VIX, last week's premium and bond type, then split into training and test data
model_data <- weekly |>
  left_join(vix, by = "date") |>                    # add that day's VIX value to each row
  arrange(etf, date) |>                             # sort by fund, then by date from oldest to newest
  group_by(etf) |>                                  # work on each fund separately
  mutate(lag_premium = lag(premium)) |>             # add a column holding the fund's premium from the week before
  ungroup() |>                                      # turn the groups back into one ordinary table
  filter(year(date) >= 2019) |>                     # remove the December 2018 weeks
  mutate(type = case_when(                          # add a column naming each fund's bond type
    etf %in% c("GOVT", "IEF", "TLT") ~ "Treasury",
    etf %in% c("LQD", "IGIB")       ~ "Investment grade",
    etf %in% c("HYG", "USHY")       ~ "High yield"),
    set = if_else(year(date) <= 2021, "train", "test"))   # label 2019–2021 as training data and 2022 as test data

# Subtract the average training-period VIX from every VIX value
model_data <- model_data |>
  mutate(vix_c = vix - mean(vix[set == "train"], na.rm = TRUE))

# Count each fund's weeks in the training and test data
count(model_data, etf, set)

# Count the missing values in each column
colSums(is.na(model_data))

# Save the finished weekly table for the modelling scripts
saveRDS(model_data, "data/processed/etf_weekly.rds")

# Show the average of the centred VIX over the training years
model_data |> filter(set == "train") |> summarise(mean(vix_c))

# Add up each fund's dividends for each week
week_divs <- daily |>
  mutate(week = floor_date(date, "week", week_start = 1)) |>   # label each day with the Monday its week starts on
  group_by(etf, week) |>                                         # treat each fund's week as its own group
  summarise(div_week = sum(dividend), .groups = "drop")          # add up the dividends paid that week

# Add next week's NAV change, with any dividend paid counted back in
model_data <- model_data |>
  left_join(week_divs, by = c("etf", "week")) |>                 # add each week's total dividend to its row
  arrange(etf, date) |>                                          # sort by fund, then by date
  group_by(etf) |>                                               # work on each fund separately
  mutate(nav_change_next = 100 * (lead(nav) + lead(div_week) - nav) / nav) |>   # % change from this week's NAV to next week's, plus next week's dividend
  ungroup()                                                      # turn the groups back into one ordinary table

# Show the smallest, largest and typical weekly NAV changes
summary(model_data$nav_change_next)

# Show the three biggest weekly NAV falls
model_data |> slice_min(nav_change_next, n = 3) |> select(etf, date, nav_change_next)

# Show the three biggest weekly NAV rises
model_data |> slice_max(nav_change_next, n = 3) |> select(etf, date, nav_change_next)

# Save the finished weekly table for the modelling scripts
saveRDS(model_data, "data/processed/etf_weekly.rds")

