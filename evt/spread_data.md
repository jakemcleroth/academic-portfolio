Data Analysis
================
Jake Mc Leroth
2025-11-30

``` r
# EVT QQ Plots
library(ReIns)

# Read in data
library(readxl)

# Modifying data
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
library(ISOweek)

# Plotting
library(ggplot2)

# Stats libraries
library(tseries)
```

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

``` r
library(urca)
library(moments)

# EVT Libraries
library(ismev)
```

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.

``` r
library(extRemes)
```

    ## Loading required package: Lmoments

    ## Loading required package: distillery

    ## 
    ## Attaching package: 'extRemes'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     qqnorm, qqplot

``` r
library(fExtremes)
```

    ## 
    ## Attaching package: 'fExtremes'

    ## The following objects are masked from 'package:ReIns':
    ## 
    ##     dgpd, pgpd, qgpd, rgpd, VaR

``` r
# Read in data
data_raw <- read_xlsx("/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/EVT/project/data/raw_data/firstrand_data.xlsx") %>%
  rename(
    date       = `Date`,
    last_price = `Last Price`,
    ask_price  = `Ask Price`,
    bid_price  = `Bid Price`
  )

# Make is chronological
data_raw <- data_raw %>% arrange(date)
```

``` r
# First Issue: Missing values
colSums(is.na(data_raw)) # Output number of missing values per column
```

    ##       date last_price  ask_price  bid_price 
    ##          0          6          5          6

``` r
data_raw[rowSums(is.na(data_raw)) > 0, ] # Output missing values
```

    ## # A tibble: 11 × 4
    ##    date                last_price ask_price bid_price
    ##    <dttm>                   <dbl>     <dbl>     <dbl>
    ##  1 1999-02-22 00:00:00       571.       NA        NA 
    ##  2 1999-06-25 00:00:00       535.       NA        NA 
    ##  3 1999-12-27 00:00:00        NA       698.       NA 
    ##  4 2000-03-21 00:00:00        NA       621.      617.
    ##  5 2000-04-21 00:00:00        NA       612.      609.
    ##  6 2000-04-24 00:00:00        NA       612.      609.
    ##  7 2000-04-27 00:00:00        NA       646.      642.
    ##  8 2000-05-01 00:00:00        NA       646.      638.
    ##  9 2008-10-23 00:00:00      1012.       NA        NA 
    ## 10 2011-06-17 00:00:00      1938        NA        NA 
    ## 11 2016-01-21 00:00:00      4007        NA        NA

``` r
data <- data_raw %>%
  drop_na(bid_price, ask_price) # Drop missing values
```

``` r
# Check Monotinicity (Ask must always be larger than bid)
sum(data$ask_price < data$bid_price, na.rm = TRUE) # Number of Errors
```

    ## [1] 1

``` r
data %>% filter(ask_price < bid_price) # Output the error
```

    ## # A tibble: 1 × 4
    ##   date                last_price ask_price bid_price
    ##   <dttm>                   <dbl>     <dbl>     <dbl>
    ## 1 1999-05-27 00:00:00       479.      479.      482.

``` r
data <- data %>% filter(ask_price >= bid_price) # Remove error
```

``` r
# Check that last price is always between bid and ask
sum(data$last_price < data$bid_price | data$last_price > data$ask_price) # output number of errors
```

    ## [1] NA

``` r
data %>% filter(last_price < bid_price | last_price > ask_price) # Output errors
```

    ## # A tibble: 3 × 4
    ##   date                last_price ask_price bid_price
    ##   <dttm>                   <dbl>     <dbl>     <dbl>
    ## 1 1999-03-19 00:00:00       522.      30.7      28.2
    ## 2 2001-09-17 00:00:00       601.     592.      589. 
    ## 3 2015-03-30 00:00:00      5592     5600      5600

``` r
data <- data %>% 
  filter(!(last_price < bid_price | last_price > ask_price)) # Remove errors
```

``` r
# Duplicates (should be 0)
sum(duplicated(data))
```

    ## [1] 0

``` r
# Checks date ordering (should be FALSE) 
any(diff(data$date) < 0)
```

    ## [1] FALSE

``` r
# 1. Cleaned quotes (what you already have, up to `data`)
quotes <- data %>%
  mutate(
    spread = ask_price - bid_price,
    spread_rel = (ask_price - bid_price)/(ask_price + bid_price) * 100 * 2
  )

# 2. Consecutive-log-change data frame
log_spread_consec <- quotes %>%
  mutate(
    log_spread_ret = log(spread_rel/lag(spread_rel))*100
  ) %>%
  filter(is.finite(log_spread_ret)) %>%    # kills first row + zeros/NaNs
  select(date, log_spread_ret)

log_spread_daily <- quotes %>%
  mutate(
    day_gap          = as.numeric(date - lag(date)),
    log_spread_daily = if_else(
      day_gap == 1 &
      spread_rel > 0 &
      lag(spread_rel) > 0,
      log(spread_rel/lag(spread_rel))*100,
      NA_real_
    )
  ) %>%
  filter(is.finite(log_spread_daily)) %>%
  select(date, log_spread_daily)
```

``` r
data_spreads <- data %>%
  mutate(
    spread = ask_price - bid_price,
    spread_quoted    = ((ask_price - bid_price) / (ask_price + bid_price)) * 2 * 100,
    log_spread_ret = log(spread_quoted/lag(spread_quoted)) * 100
  )

# Drop first row
data_spreads <- data_spreads %>%
  slice(-(which(is.na(log_spread_ret))[1]))

sum(!is.finite(data_spreads$log_spread_ret)) # number of infinite returns
```

    ## [1] 6

``` r
# drop undefined / infinite log-returns
data_spreads <- data_spreads %>%
  filter(is.finite(log_spread_ret)) # Check this
```
