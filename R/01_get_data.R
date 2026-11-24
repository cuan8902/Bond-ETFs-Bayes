# 01_get_data.R: download the VIX and ETF prices
library(quantmod)
library(tidyverse)
#using data from yahoo 
getSymbols("^VIX", src = "yahoo", from = "2019-01-01", to = "2023-01-01") #sends request to yahoo for ^VIX(ticker) history
getSymbols("LQD", src = "yahoo", from = "2019-01-01", to = "2023-01-01") #sends request to yahoo for LQD(ticker) history

plot(Cl(VIX))   # "Cl" is quantmode function taking the closing price
head(LQD) # returns the price of lqd for first 6 trading days


read_nav <- function(ticker) {
  read_csv(paste0("data/raw/nav_", ticker, ".csv"), #reads file data/raw/nav_LQD.csv"
           col_types = cols(.default = col_character())) |> # tells to read as a text (|>) is the pipe, takes reulst on left and feeds into next line
    set_names(c("date", "nav", "dividend", "shares")) |> # renaming the column in excel sheets by position
    transmute(etf      = ticker, #creates coloumns listed below  eg etf = ticker adds a columns , LQD on every row to tell funds apart once 7 added 
              date     = mdy(date), # turns text like "jan 02 2019 into real date R can sort flter and match with yahoo
              nav      = parse_number(nav), #turns string to integer
              dividend = replace_na(parse_number(dividend, na = "--"),0)) # turns string to integer or "--" into 0
          
}

lqd_nav <- read_nav("LQD")
lqd_nav |> filter(date == as.Date("2019-01-02")) #|> pull(nav) is optional after to see the real nav valueat that time #checking the code read the file corrrectly 

# Read all seven iShares funds into one table
tickers <- c("GOVT", "IEF", "TLT", "LQD", "IGIB", "HYG", "USHY") #giving each ticker the variable "ticker" so it can be ran through read_nav fucntion
navs <- map(tickers, read_nav) |> list_rbind() #runs read_nav for each ticker and then stacks the seven tables on top of eavh other

count(navs, etf)   # rows per fund

# Daily closing prices for all seven funds
get_close <- function(ticker) {
  x <- getSymbols(ticker, src = "yahoo", from = "2018-12-01", to = "2023-01-01",
                  auto.assign = FALSE)
  tibble(etf = ticker, date = index(x), close = as.numeric(Cl(x)))
}

prices <- map(tickers, get_close) |> list_rbind()
count(prices, etf)   # about 1,030 rows each

saveRDS(prices, "data/raw/prices_yahoo.rds")   # keep a copy in case Yahoo changes

# VIX daily closes, from December 2018 to match the prices
vix <- getSymbols("^VIX", src = "yahoo", from = "2018-12-01", to = "2023-01-01",
                  auto.assign = FALSE)
vix <- tibble(date = index(vix), vix = as.numeric(Cl(vix)))
saveRDS(vix, "data/raw/vix_yahoo.rds")

saveRDS(navs, "data/raw/navs_ishares.rds")
