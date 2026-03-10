Bid-Ask Spreads
================
Jake Mc Leroth
2025-09-16

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

# Data Preprocessing

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
data <- data_raw
# data <- data_raw %>% drop_na() # Drop missing values
```

``` r
# Check Monotinicity (Ask must always be larger than bid)
sum(data$ask_price < data$bid_price, na.rm = TRUE) # Number of Errors
```

    ## [1] 1

``` r
data %>% filter(ask_price <= bid_price) # Output the error
```

    ## # A tibble: 5 × 4
    ##   date                last_price ask_price bid_price
    ##   <dttm>                   <dbl>     <dbl>     <dbl>
    ## 1 1998-08-20 00:00:00       638.      638.      638.
    ## 2 1999-02-10 00:00:00       579.      579.      579.
    ## 3 1999-05-27 00:00:00       479.      479.      482.
    ## 4 2000-08-14 00:00:00       587.      587.      587.
    ## 5 2015-03-30 00:00:00      5592      5600      5600

``` r
# Instead of removing rows, set bad quotes to NA
bad_idx <- data$ask_price <= data$bid_price
data$ask_price[bad_idx] <- NA
data$bid_price[bad_idx] <- NA
```

``` r
# Check that last price is always between bid and ask
sum(data$last_price < data$bid_price | data$last_price > data$ask_price, na.rm = TRUE) # output number of errors
```

    ## [1] 2

``` r
# Output the problematic rows
data %>% filter(last_price < bid_price | last_price > ask_price)
```

    ## # A tibble: 2 × 4
    ##   date                last_price ask_price bid_price
    ##   <dttm>                   <dbl>     <dbl>     <dbl>
    ## 1 1999-03-19 00:00:00       522.      30.7      28.2
    ## 2 2001-09-17 00:00:00       601.     592.      589.

``` r
# Set only the first problematic last_price to NA
bad_index <- which(data$last_price < data$bid_price | data$last_price > data$ask_price)[1]
data$last_price[bad_index] <- NA
data$ask_price[bad_index] <- NA
data$bid_price[bad_index] <- NA
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
data_spreads <- data %>%
  mutate(
    spread        = ask_price - bid_price,
    spread_rel = ((ask_price - bid_price) / (ask_price + bid_price)) * 2 * 100,
    spread_dif = (spread_rel - lag(spread_rel)),
    spread_log = log(spread_rel/lag(spread_rel)) * 100
  )

# Drop all NA / ±Inf log-returns (including the first lagged one)
sum(!is.finite(data_spreads$spread_log)) # number of non-finite returns
```

    ## [1] 25

``` r
data_spreads <- data_spreads %>%
  filter(is.finite(spread_log))
```

``` r
ggplot(data_spreads[6083:nrow(data_spreads),],aes(x=date)) +
  geom_line(aes(y=last_price, colour="Last"), linewidth=0.4) +
  geom_line(aes(y=bid_price,  colour="Bid"),  linewidth=0.4) +
  geom_line(aes(y=ask_price,  colour="Ask"),  linewidth=0.4) +
  labs(
    title  = "FirstRand LTD Price with Bid–Ask Quotes",
    x      = "Date",
    y      = "Price",
    colour = "Series"
  ) +
  theme_minimal()
```

![](evt_project_v1-copy_files/figure-gfm/prices%20plots-1.png)<!-- -->

# Exploratory Data Analysis

``` r
print("Spreads: ")
```

    ## [1] "Spreads: "

``` r
summary(data_spreads$spread)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   0.800   2.000   4.000   6.705   8.070 109.000

``` r
print("Quoted Spreads (%): ")
```

    ## [1] "Quoted Spreads (%): "

``` r
summary(data_spreads$spread_rel)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.01153 0.08304 0.20942 0.32696 0.43557 9.23130

``` r
print("Change in spreads: ")
```

    ## [1] "Change in spreads: "

``` r
summary(data_spreads$spread_dif)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -7.704493 -0.149533 -0.000055  0.000280  0.151226  8.653712

``` r
print("Log Return (%): ")
```

    ## [1] "Log Return (%): "

``` r
summary(data_spreads$spread_log)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -428.7857  -77.5297   -0.1608    0.0155   79.4589  468.8319

## Spreads

``` r
# # Plot 1: Absolute spread
ggplot(data_spreads, aes(x = date, y = spread)) +
  geom_line(color = "steelblue") +
  labs(title = "Bid-Ask Spread Over Time", x = "Time", y = "Spread") +
  theme_minimal()
```

![](evt_project_v1-copy_files/figure-gfm/absolute%20spread%20time%20series%20plots-1.png)<!-- -->

``` r
acf(data_spreads$spread, main = "ACF of Bid-Ask Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/absolute%20spread%20time%20series%20plots-2.png)<!-- -->

``` r
acf(data_spreads$spread, type = "partial", main = "PACF of Bid-Ask Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/absolute%20spread%20time%20series%20plots-3.png)<!-- -->

``` r
adf.test(data_spreads$spread)
```

    ## Warning in adf.test(data_spreads$spread): p-value smaller than printed p-value

    ## 
    ##  Augmented Dickey-Fuller Test
    ## 
    ## data:  data_spreads$spread
    ## Dickey-Fuller = -12.709, Lag order = 18, p-value = 0.01
    ## alternative hypothesis: stationary

``` r
kpss.test(data_spreads$spread, null = "Level")
```

    ## Warning in kpss.test(data_spreads$spread, null = "Level"): p-value smaller than
    ## printed p-value

    ## 
    ##  KPSS Test for Level Stationarity
    ## 
    ## data:  data_spreads$spread
    ## KPSS Level = 11.753, Truncation lag parameter = 11, p-value = 0.01

## Relative Spreads

``` r
ggplot(data_spreads, aes(x = date, y = spread_rel)) +
  geom_line(color = "darkgreen") +
  labs(title = "Relative Spread Over Time", x = "Time", y = "Spread") +
  theme_minimal()
```

![](evt_project_v1-copy_files/figure-gfm/Quoted%20spreads%20time%20series%20plots-1.png)<!-- -->

``` r
acf(data_spreads$spread_rel, main = "ACF of Relative Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/Quoted%20spreads%20time%20series%20plots-2.png)<!-- -->

``` r
acf(data_spreads$spread_rel, type = "partial", main = "PACF of Relative Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/Quoted%20spreads%20time%20series%20plots-3.png)<!-- -->

``` r
adf.test(data_spreads$spread_rel)
```

    ## Warning in adf.test(data_spreads$spread_rel): p-value smaller than printed
    ## p-value

    ## 
    ##  Augmented Dickey-Fuller Test
    ## 
    ## data:  data_spreads$spread_rel
    ## Dickey-Fuller = -13.505, Lag order = 18, p-value = 0.01
    ## alternative hypothesis: stationary

``` r
kpss.test(data_spreads$spread_rel, null = "Level")
```

    ## Warning in kpss.test(data_spreads$spread_rel, null = "Level"): p-value smaller
    ## than printed p-value

    ## 
    ##  KPSS Test for Level Stationarity
    ## 
    ## data:  data_spreads$spread_rel
    ## KPSS Level = 33.236, Truncation lag parameter = 11, p-value = 0.01

## Log change of spreads

``` r
# Time series line plot
ggplot(data_spreads, aes(x = date, y = spread_log)) +
  geom_line(color = "firebrick") +
  labs(
    title = "Change in Log of Relative Spreads Over Time",
    x = "Time",
    y = "Spread (%)"
  ) +
  theme_minimal()
```

![](evt_project_v1-copy_files/figure-gfm/log%20return%20spread%20time%20series%20plots-1.png)<!-- -->

``` r
# ACF and PACF
acf(data_spreads$spread_log, main = "ACF of Change in Log of Relative Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/log%20return%20spread%20time%20series%20plots-2.png)<!-- -->

``` r
acf(data_spreads$spread_log, type = "partial", main = "PACF of Change in Log of Relative Spreads")
```

![](evt_project_v1-copy_files/figure-gfm/log%20return%20spread%20time%20series%20plots-3.png)<!-- -->

``` r
adf.test(data_spreads$spread_log)
```

    ## Warning in adf.test(data_spreads$spread_log): p-value smaller than printed
    ## p-value

    ## 
    ##  Augmented Dickey-Fuller Test
    ## 
    ## data:  data_spreads$spread_log
    ## Dickey-Fuller = -31.287, Lag order = 18, p-value = 0.01
    ## alternative hypothesis: stationary

``` r
kpss.test(data_spreads$spread_log, null = "Level")
```

    ## Warning in kpss.test(data_spreads$spread_log, null = "Level"): p-value greater
    ## than printed p-value

    ## 
    ##  KPSS Test for Level Stationarity
    ## 
    ## data:  data_spreads$spread_log
    ## KPSS Level = 0.0037026, Truncation lag parameter = 11, p-value = 0.1

# EVT Modelling

``` r
block_size <- 22
n_blocks <- floor(length(data_spreads$spread_log) / block_size)

# Trim to an exact multiple of the block size
spread_log_trimmed <- tail(data_spreads$spread_log, n_blocks * block_size)


cat("MLE: \n")
```

    ## MLE:

``` r
gev <- gevFit(x = spread_log_trimmed, block = block_size, type = "mle", title = "GEV fit daily returns", description = "Maximum liklihood method is used")
mle_estimates <- gev@fit$par.ests
names(mle_estimates) <- c("gamma", "mu", "sigma")
mle_vcv <- gev@fit$varcov
mle_se <- gev@fit$par.ses
names(mle_se) <- c("gamma", "mu", "sigma")
colnames(mle_vcv) <- c("gamma", "mu", "sigma")
rownames(mle_vcv) <- c("gamma", "mu", "sigma")
print(mle_estimates)
```

    ##       gamma          mu       sigma 
    ##  -0.1633239 210.5369087  61.3828081

``` r
mle_estimates_lower <- mle_estimates - qnorm(0.975)*mle_se
mle_estimates_upper <- mle_estimates + qnorm(0.975)*mle_se
gev_ci <- cbind(mle_estimates_lower,mle_estimates_upper)
colnames(gev_ci) <- c("Lower", "Upper")
rownames(gev_ci) <- c("gamma", "mu", "sigma")
gev_ci
```

    ##             Lower      Upper
    ## gamma  -0.2365169  -0.090131
    ## mu    202.8966432 218.177174
    ## sigma  55.9412881  66.824328

``` r
cat("\n")
```

``` r
cat("SE: \n")
```

    ## SE:

``` r
mle_se
```

    ##      gamma         mu      sigma 
    ## 0.03734402 3.89816626 2.77633673

``` r
# Reshape into a matrix: each column is a block
x_mat <- matrix(spread_log_trimmed, nrow = block_size, ncol = n_blocks)
# Block maxima
block_maxima <- apply(x_mat, 2, max)

start_vals <- list(
  location = mle_estimates["mu"],
  scale    = mle_estimates["sigma"],
  shape    = mle_estimates["gamma"]
)
# extRemes
gev_ext <- fevd(
  x      = block_maxima,
  type   = "GEV",
  method = "MLE",
  initial = start_vals,
  units = "Log Return (%)", 
  period.basis = "month"
)

plot.fevd(x = gev_ext, type = "qq", main = "GEV QQ Plot")
```

![](evt_project_v1-copy_files/figure-gfm/QQ-1.png)<!-- -->

``` r
plot.fevd(x = gev_ext, type ="rl", main = "Return Level Plot")
```

![](evt_project_v1-copy_files/figure-gfm/QQ-2.png)<!-- -->

``` r
ci(gev_ext, type = "parameter", which.par =  3, xrange = c(-0.235, -0.07),method = "proflik", verbose = TRUE, nint=1000,
           main = "Profile Likelihood of Gamma")
```

    ## 
    ##  Preparing to calculate  95 % CI for  shape parameter 
    ## 
    ##  Model is   fixed 
    ## 
    ##  Using Profile Likelihood Method.
    ## 
    ##  Calculating profile likelihood.  This may take a few moments.

![](evt_project_v1-copy_files/figure-gfm/Profile%20Likelihood-1.png)<!-- -->

    ## 
    ##  Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value =  -1716.409

    ## fevd(x = block_maxima, type = "GEV", method = "MLE", initial = start_vals, 
    ##     units = "Log Return (%)", period.basis = "month")
    ## 
    ## [1] "Profile Likelihood"
    ## 
    ## [1] "shape: -0.163"
    ## 
    ## [1] "95% Confidence Interval: (-0.2259, -0.0816)"

``` r
# ismev
cat("ismev: \n")
```

    ## ismev:

``` r
gev_ismev <- gev.fit(block_maxima)
```

    ## $conv
    ## [1] 0
    ## 
    ## $nllh
    ## [1] 1714.488
    ## 
    ## $mle
    ## [1] 210.5271980  61.3870908  -0.1633306
    ## 
    ## $se
    ## [1] 3.898278 2.777075 0.037348

``` r
gev.diag(gev_ismev)
```

![](evt_project_v1-copy_files/figure-gfm/sanity%20check%20for%20parameter%20estimates-1.png)<!-- -->

``` r
gev_q <- function(p, mu, sigma, gamma){
  if (gamma == 0){
    mu - sigma * log(-log(p))
  } else {
    mu + (sigma/gamma) * ((-log(p))^(-gamma) - 1)
  }
}

gev_kappa <- function(p, gamma, sigma) {
  t  <- -log(1 - p)              # t = -log(1-p)
  a  <- t^(-gamma)               # (-log(1-p))^(-gamma)

  dq_dsigma <- (1 / gamma) * (a - 1)

  dq_dgamma <- - (sigma / gamma^2) * (a - 1) -
               (sigma / gamma)    *  a * log(t)

  dq_dmu <- 1

  # Return in order matching varcov: (gamma, mu, sigma)
  c(gamma = dq_dgamma,
    mu    = dq_dmu,
    sigma = dq_dsigma)
}

gev_CI <- function(p_tail, est, vcv){
  p_q <- 1 - p_tail  # CDF prob
  
  z_975 <- qnorm(0.975)
  
  qhat <- gev_q(
    p     = p_q,                    # quantile prob
    mu    = est[["mu"]],
    sigma = est[["sigma"]],
    gamma = est[["gamma"]]
  )

  kappa <- gev_kappa(
    p     = p_tail,                 # *** tail prob here! ***
    gamma = est[["gamma"]],
    sigma = est[["sigma"]]
  )

  var_q <- as.numeric(t(kappa) %*% vcv %*% kappa)
  se_q  <- sqrt(var_q)

  ci <- c(
    lower = qhat - z_975 * se_q,
    upper = qhat + z_975 * se_q
  )

  list(qhat = qhat, se = se_q, CI = ci)
}

# Compute both:
res_mle_001  <- gev_CI(0.01,  mle_estimates, mle_vcv)
res_mle_0001 <- gev_CI(0.001, mle_estimates, mle_vcv)
res_mle_00001 <- gev_CI(0.0001, mle_estimates, mle_vcv)

cat("p=0.01: \n \n")
```

    ## p=0.01: 
    ## 

``` r
cat("Estimate: ")
```

    ## Estimate:

``` r
res_mle_001$qhat
```

    ## [1] 409.0731

``` r
cat("CI: \n")
```

    ## CI:

``` r
res_mle_001$CI
```

    ##    lower    upper 
    ## 384.3318 433.8144

``` r
cat("\n")
```

``` r
cat("p=0.001: \n \n")
```

    ## p=0.001: 
    ## 

``` r
cat("Estimate: ")
```

    ## Estimate:

``` r
res_mle_0001$qhat
```

    ## [1] 464.7361

``` r
cat("CI: \n")
```

    ## CI:

``` r
res_mle_0001$CI
```

    ##    lower    upper 
    ## 420.6852 508.7869

``` r
cat("p=0.0001: \n \n")
```

    ## p=0.0001: 
    ## 

``` r
cat("Estimate: ")
```

    ## Estimate:

``` r
res_mle_00001$qhat
```

    ## [1] 502.8681

``` r
cat("CI: \n")
```

    ## CI:

``` r
res_mle_00001$CI
```

    ##    lower    upper 
    ## 438.8997 566.8365

``` r
ci(gev_ext, type="Return.level", method = "proflik",
           xrange = c(384,460), verbose = TRUE, nint=1000,
           main = "Profile Likelihood", return.period = 100, xlab = "100-month return level (p = 0.01)")
```

    ## 
    ##  Preparing to calculate  95 % CI for  100-month return level 
    ## 
    ##  Model is   fixed 
    ## 
    ##  Using Profile Likelihood Method.
    ## 
    ##  Calculating profile likelihood.  This may take a few moments.

![](evt_project_v1-copy_files/figure-gfm/Quantile%20Profile%20Likelihood-1.png)<!-- -->

    ## 
    ##  Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value =  -1716.409

    ## fevd(x = block_maxima, type = "GEV", method = "MLE", initial = start_vals, 
    ##     units = "Log Return (%)", period.basis = "month")
    ## 
    ## [1] "Profile Likelihood"
    ## 
    ## [1] "100-month return level: 409.086"
    ## 
    ## [1] "95% Confidence Interval: (390.0215, 441.8595)"

``` r
ci(gev_ext, type="Return.level", method = "proflik",
           xrange = c(420,550), verbose = TRUE, nint=1000,
           main = "Profile Likelihood", return.period = 1000, xlab = "1000-month return level (p = 0.001)")
```

    ## 
    ##  Preparing to calculate  95 % CI for  1000-month return level 
    ## 
    ##  Model is   fixed 
    ## 
    ##  Using Profile Likelihood Method.
    ## 
    ##  Calculating profile likelihood.  This may take a few moments.

![](evt_project_v1-copy_files/figure-gfm/Quantile%20Profile%20Likelihood-2.png)<!-- -->

    ## 
    ##  Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value =  -1716.409

    ## fevd(x = block_maxima, type = "GEV", method = "MLE", initial = start_vals, 
    ##     units = "Log Return (%)", period.basis = "month")
    ## 
    ## [1] "Profile Likelihood"
    ## 
    ## [1] "1000-month return level: 464.76"
    ## 
    ## [1] "95% Confidence Interval: (433.7068, 526.7414)"

``` r
ci(gev_ext, type="Return.level", method = "proflik",
           xrange = c(450,625), verbose = TRUE, nint=1000,
           main = "Profile Likelihood", return.period = 10000, xlab = "10000-month return level (p = 0.0001)")
```

    ## 
    ##  Preparing to calculate  95 % CI for  10000-month return level 
    ## 
    ##  Model is   fixed 
    ## 
    ##  Using Profile Likelihood Method.
    ## 
    ##  Calculating profile likelihood.  This may take a few moments.

![](evt_project_v1-copy_files/figure-gfm/Quantile%20Profile%20Likelihood-3.png)<!-- -->

    ## 
    ##  Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value =  -1716.409

    ## fevd(x = block_maxima, type = "GEV", method = "MLE", initial = start_vals, 
    ##     units = "Log Return (%)", period.basis = "month")
    ## 
    ## [1] "Profile Likelihood"
    ## 
    ## [1] "10000-month return level: 502.902"
    ## 
    ## [1] "95% Confidence Interval: (459.7253, 598.1753)"

``` r
cat("GEV Upper Endpoint: \n")
```

    ## GEV Upper Endpoint:

``` r
as.numeric(mle_estimates["mu"] - mle_estimates["sigma"]/mle_estimates["gamma"])
```

    ## [1] 586.3716

``` r
theta_blocks <- function(x, u, r){
  x <- as.numeric(x)
  n <- length(x)
  k <- floor(n / r)
  
  if(k <= 0) stop("Block length r is too large for the sample size n.")
  
  # numerator: proportion of blocks whose max exceeds u
  exceed_blocks <- logical(k)
  for(j in seq_len(k)){
    idx_start <- (j - 1) * r + 1
    idx_end   <- j * r
    block_max <- max(x[idx_start:idx_end], na.rm = TRUE)
    exceed_blocks[j] <- block_max > u
  }
  num <- mean(exceed_blocks)  # k^{-1} sum 1{M > u}
  
  # denominator: r * (n^{-1} sum 1{X_i > u})
  tail_prob_hat <- mean(x > u, na.rm = TRUE)
  denom <- r * tail_prob_hat
  
  theta_hat <- num / denom
  theta_hat
}

theta_blocks_sliding <- function(x, u, r){
  x <- as.numeric(x)
  n <- length(x)
  
  if(r >= n) stop("For sliding blocks, r must be < n.")
  
  n_blocks <- n - r + 1
  
  exceed_blocks <- logical(n_blocks)
  for(i in seq_len(n_blocks)){
    block_max <- max(x[i:(i + r - 1)], na.rm = TRUE)
    exceed_blocks[i] <- block_max > u
  }
  num <- mean(exceed_blocks)  # (n-r+1)^{-1} sum 1{M > u}
  
  tail_prob_hat <- mean(x > u, na.rm = TRUE)
  denom <- r * tail_prob_hat
  
  theta_hat <- num / denom
  theta_hat
}


# Threshold probabilities
p_vals <- c(0.98, 0.99, 0.995, 0.999)

# Block size
r <- block_size

# Storage
results <- data.frame(
  p = p_vals,
  u = NA,
  theta_blocks = NA,
  theta_sliding = NA
)

for(i in seq_along(p_vals)){
  
  p <- p_vals[i]
  u <- quantile(spread_log_trimmed, probs = p)
  
  theta_B <- theta_blocks(spread_log_trimmed, u, r)
  theta_B_sl <- theta_blocks_sliding(spread_log_trimmed, u, r)
  
  results$u[i] <- u
  results$theta_blocks[i] <- theta_B
  results$theta_sliding[i] <- theta_B_sl
}

print(results)
```

    ##       p        u theta_blocks theta_sliding
    ## 1 0.980 270.1981    0.7481481     0.7319128
    ## 2 0.990 298.6159    0.8676471     0.8334823
    ## 3 0.995 324.2127    0.9411765     0.8918194
    ## 4 0.999 378.7496    1.0000000     0.8858803

``` r
theta <- 0.88

gev_q_stationary <- function(p_tail, n, theta, mu, sigma, gamma){
  # effective argument inside the log:
  t <- -n * theta * log(1 - p_tail)   # t = -nθ log(1-p)
  
  if (gamma == 0){
    mu - sigma * log(t)
  } else {
    mu + (sigma / gamma) * (t^(-gamma) - 1)
  }
}

theta <- 0.88
n_block <- 22   # block size, same as in your GEV fit

gev_kappa_stationary <- function(p_tail, n, theta, gamma, sigma) {
  # t = - n * theta * log(1 - p_tail)
  t <- -n * theta * log(1 - p_tail)
  a <- t^(-gamma)
  
  dq_dsigma <- (1 / gamma) * (a - 1)
  
  dq_dgamma <- - (sigma / gamma^2) * (a - 1) -
               (sigma / gamma)    *  a * log(t)
  
  dq_dmu <- 1
  
  # Return in order matching varcov: (gamma, mu, sigma)
  c(gamma = dq_dgamma,
    mu    = dq_dmu,
    sigma = dq_dsigma)
}

gev_CI_stationary <- function(p_tail, n, theta, est, vcv){
  z_975 <- qnorm(0.975)
  
  # point estimate
  qhat <- gev_q_stationary(
    p_tail = p_tail,
    n      = n,
    theta  = theta,
    mu     = est[["mu"]],
    sigma  = est[["sigma"]],
    gamma  = est[["gamma"]]
  )
  
  # gradient
  kappa <- gev_kappa_stationary(
    p_tail = p_tail,
    n      = n,
    theta  = theta,
    gamma  = est[["gamma"]],
    sigma  = est[["sigma"]]
  )
  
  # delta-method variance
  var_q <- as.numeric(t(kappa) %*% vcv %*% kappa)
  se_q  <- sqrt(var_q)
  
  ci <- c(
    lower = qhat - z_975 * se_q,
    upper = qhat + z_975 * se_q
  )
  
  list(qhat = qhat, se = se_q, CI = ci)
}
res_stat_001   <- gev_CI_stationary(0.01,   n_block, theta, mle_estimates, mle_vcv)
res_stat_0001  <- gev_CI_stationary(0.001,  n_block, theta, mle_estimates, mle_vcv)
res_stat_00001 <- gev_CI_stationary(0.0001, n_block, theta, mle_estimates, mle_vcv)

cat("p_tail = 0.01\n")
```

    ## p_tail = 0.01

``` r
res_stat_001$qhat
```

    ## [1] 298.7062

``` r
res_stat_001$CI
```

    ##    lower    upper 
    ## 288.6671 308.7453

``` r
cat("\np_tail = 0.001\n")
```

    ## 
    ## p_tail = 0.001

``` r
res_stat_0001$qhat
```

    ## [1] 389.019

``` r
res_stat_0001$CI
```

    ##    lower    upper 
    ## 368.9764 409.0615

``` r
cat("\np_tail = 0.0001\n")
```

    ## 
    ## p_tail = 0.0001

``` r
res_stat_00001$qhat
```

    ## [1] 450.8879

``` r
res_stat_00001$CI
```

    ##    lower    upper 
    ## 412.5958 489.1800

## Peaks Over Threshold

``` r
MeanExcess(data_spreads$spread_log[data_spreads$spread_log>0])
```

![](evt_project_v1-copy_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
gpd_reins <- GPDmle(data = data_spreads$spread_log[data_spreads$spread_log>0], plot = T)
```

![](evt_project_v1-copy_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
u = 230

log_ret_spreads_declustered <- decluster(x = data_spreads$spread_log, threshold = u, method = "runs", r = 5)
log_ret_spreads_declustered
```

    ## 
    ##  data_spreads$spread_log  declustered via runs  declustering.
    ## 
    ##  Estimated extremal index (intervals estimate) =  0.7894754 
    ## 
    ##  Number of clusters =  209 
    ## 
    ##  Run length =  5

``` r
gpd <- gpdFit(log_ret_spreads_declustered, u = u, type = "mle")
cat("MLE: \n")
```

    ## MLE:

``` r
gpd@fit$par.ests
```

    ##         xi       beta 
    ## -0.1958732 60.6623559

``` r
# Extract MLEs
par_hat <- gpd@fit$par.ests      # xi and beta (scale)
names(par_hat) <- c("gamma", "sigma")

# Extract standard errors
se_hat <- gpd@fit$par.ses
names(se_hat) <- c("gamma", "sigma")

# 95% Wald CI
z <- qnorm(0.975)

ci_lower <- par_hat - z * se_hat
ci_upper <- par_hat + z * se_hat

gpd_ci <- cbind(
  Estimate = par_hat,
  SE = se_hat,
  Lower = ci_lower,
  Upper = ci_upper
)

gpd_ci
```

    ##         Estimate         SE      Lower      Upper
    ## gamma -0.1958732 0.04712672 -0.2882399 -0.1035065
    ## sigma 60.6623559 5.01145032 50.8400938 70.4846181

``` r
# extRemes
start_vals <- list(
  scale    = gpd@fit$par.ests["beta"],
  shape    = gpd@fit$par.ests["xi"]
)

gpd_ext <- fevd(x = log_ret_spreads_declustered, threshold = u, type = "GP", method = "MLE", initial = start_vals, units = "Log Return (%)", period.basis = "month" )
```

``` r
plot.fevd(x = gpd_ext, type = "qq", main = "GEV QQ Plot")
```

![](evt_project_v1-copy_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot.fevd(x = gpd_ext, type ="rl", main = "Return Level Plot")
```

![](evt_project_v1-copy_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
gpd_q <- function(p_tail, sigma, gamma){
  if (gamma == 0){
    -sigma * log(p_tail)
  } else {
    (sigma / gamma) * (p_tail^(-gamma) - 1)
  }
}

gamma_hat <- as.numeric(gpd@fit$par.ests[1])
sigma_hat <- as.numeric(gpd@fit$par.ests[2])


cat("Quantile estimate for p=0.01: \n")
```

    ## Quantile estimate for p=0.01:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.01 <- gpd_q(0.01,sigma_hat, gamma_hat)
q_0.01 + u
```

    ## [1] 414.0419

``` r
cat("Quantile estimate for p=0.001: \n")
```

    ## Quantile estimate for p=0.001:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.001 <- gpd_q(0.001,sigma_hat, gamma_hat)
q_0.001 + u
```

    ## [1] 459.6589

``` r
cat("Quantile estimate for p=0.0001: \n")
```

    ## Quantile estimate for p=0.0001:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.0001 <- gpd_q(0.0001,sigma_hat, gamma_hat)
q_0.0001 + u
```

    ## [1] 488.7161

``` r
gpd_p_underlying <- function(p, u, sigma, xi, zeta_u, theta){
  # m = 1/p return period in time units
  m <- 1 / p
  
  arg <- m * zeta_u * theta   # = θ ζu / p
  if(arg <= 0) stop("θ * ζu / p must be > 0")
  
  if(xi == 0){
    # Gumbel case
    u + sigma * log(arg)
  } else {
    # General GPD case
    u + (sigma / xi) * (arg^xi - 1)
  }
}

zeta_u <- mean(data_spreads$spread_log > u, na.rm = TRUE)
theta <- 0.7894754
cat("Quantile estimate for p=0.01: \n")
```

    ## Quantile estimate for p=0.01:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.01_f <- gpd_p_underlying(0.01, u = u, sigma = par_hat["sigma"], xi = par_hat["gamma"], zeta_u = zeta_u, theta = theta)
q_0.01_f 
```

    ##    sigma 
    ## 292.3936

``` r
cat("Quantile estimate for p=0.001: \n")
```

    ## Quantile estimate for p=0.001:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.001_f <- gpd_p_underlying(0.001, u = u, sigma = par_hat["sigma"], xi = par_hat["gamma"], zeta_u = zeta_u, theta = theta)
q_0.001_f  
```

    ##    sigma 
    ## 382.1712

``` r
cat("Quantile estimate for p=0.0001: \n")
```

    ## Quantile estimate for p=0.0001:

``` r
cat("MLE:")
```

    ## MLE:

``` r
q_0.0001_f <- gpd_p_underlying(0.0001, u = u, sigma = par_hat["sigma"], xi = par_hat["gamma"], zeta_u = zeta_u, theta = theta)
q_0.0001_f 
```

    ##    sigma 
    ## 439.3579

``` r
cat("End point: \n")
```

    ## End point:

``` r
as.numeric(u - par_hat["sigma"]/par_hat["gamma"])
```

    ## [1] 539.7022
