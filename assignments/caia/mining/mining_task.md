Impala Analysis
================
Jake Mc Leroth
2025-03-05

``` r
#install.packages("tseries")
#install.packages("urca")
#install.packages("readxl")
library(tseries)
```

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

``` r
library(urca)
library(readxl)
library(tsfeatures)
library(fpp3)
```

    ## Registered S3 method overwritten by 'tsibble':
    ##   method               from 
    ##   as_tibble.grouped_df dplyr

    ## ── Attaching packages ──────────────────────────────────────────── fpp3 1.0.1 ──

    ## ✔ tibble      3.2.1     ✔ tsibble     1.1.6
    ## ✔ dplyr       1.1.4     ✔ tsibbledata 0.4.1
    ## ✔ tidyr       1.3.1     ✔ feasts      0.4.1
    ## ✔ lubridate   1.9.3     ✔ fable       0.4.1
    ## ✔ ggplot2     3.5.1

    ## ── Conflicts ───────────────────────────────────────────────── fpp3_conflicts ──
    ## ✖ lubridate::date()       masks base::date()
    ## ✖ dplyr::filter()         masks stats::filter()
    ## ✖ tsibble::intersect()    masks base::intersect()
    ## ✖ tsibble::interval()     masks lubridate::interval()
    ## ✖ dplyr::lag()            masks stats::lag()
    ## ✖ tsibble::setdiff()      masks base::setdiff()
    ## ✖ tsibble::union()        masks base::union()
    ## ✖ feasts::unitroot_kpss() masks tsfeatures::unitroot_kpss()
    ## ✖ feasts::unitroot_pp()   masks tsfeatures::unitroot_pp()

We investigate the relationship between the Impala Platinum share price
and the average price of Platinum, Palladium, and Rhodium. Looking at
the excel spreadsheet, the correlation between the monthly returns is
0.37, which is relatively high. Looking further, the year on year
returns have a correlation of 90%. The strong relationship between the
two is evidence supporting the idea that the value of the commodities
that a firm sells drives the firms performance. Evidence against this
idea is that over the period of 2014 to 2024, Impala had a return of
-24%, whereas the average between Platinum, Palladium, and Rhodium had a
return of 105%. We therefore take a more robust approach to determening
if there is an underlying relationship. This is done by testing for
cointegration, which do by using the Engle-Granger and the KPSS test.
Should the two processes be cointegrated, that means there is a long
term relationship between them, where the spread of the prices is a mean
reverting/stationary process. Thus, divergences of the series should be
short term and will be followed by convergences later.

``` r
data <- read_excel("CAIA.xlsx", sheet = "IMPALA PLATINUM")
impala <- data$last_price
pgm_ave <- data$pgm

impala_ts <- ts(impala)
pgm_ave_ts <- ts(pgm_ave)
ln_pgm_ave_tes <- ts(log(pgm_ave))
```

``` r
# Step 1: Run OLS regression
model <- lm(impala_ts ~ pgm_ave_ts)

# Step 2: Extract residuals
residuals <- residuals(model)

# Step 3: ADF test on residuals
summary(ur.df(residuals, type = "none", lags = 24, selectlags = "AIC"))
```

    ## 
    ## ############################################### 
    ## # Augmented Dickey-Fuller Test Unit Root Test # 
    ## ############################################### 
    ## 
    ## Test regression none 
    ## 
    ## 
    ## Call:
    ## lm(formula = z.diff ~ z.lag.1 - 1 + z.diff.lag)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -52.239  -9.706  -2.234   9.736  48.465 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## z.lag.1     -0.18218    0.06773  -2.690  0.00832 **
    ## z.diff.lag1 -0.09345    0.10033  -0.931  0.35377   
    ## z.diff.lag2 -0.14318    0.09733  -1.471  0.14427   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.13 on 104 degrees of freedom
    ## Multiple R-squared:  0.1336, Adjusted R-squared:  0.1086 
    ## F-statistic: 5.345 on 3 and 104 DF,  p-value: 0.001831
    ## 
    ## 
    ## Value of test-statistic is: -2.6901 
    ## 
    ## Critical values for test statistics: 
    ##       1pct  5pct 10pct
    ## tau1 -2.58 -1.95 -1.62

``` r
summary(ur.df(residuals, type = "none", lags = 24, selectlags = "BIC"))
```

    ## 
    ## ############################################### 
    ## # Augmented Dickey-Fuller Test Unit Root Test # 
    ## ############################################### 
    ## 
    ## Test regression none 
    ## 
    ## 
    ## Call:
    ## lm(formula = z.diff ~ z.lag.1 - 1 + z.diff.lag)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -48.997 -10.522  -2.905   8.923  48.940 
    ## 
    ## Coefficients:
    ##            Estimate Std. Error t value Pr(>|t|)   
    ## z.lag.1    -0.21226    0.06492  -3.269  0.00146 **
    ## z.diff.lag -0.05483    0.09737  -0.563  0.57455   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.22 on 105 degrees of freedom
    ## Multiple R-squared:  0.1156, Adjusted R-squared:  0.0987 
    ## F-statistic: 6.859 on 2 and 105 DF,  p-value: 0.001586
    ## 
    ## 
    ## Value of test-statistic is: -3.2694 
    ## 
    ## Critical values for test statistics: 
    ##       1pct  5pct 10pct
    ## tau1 -2.58 -1.95 -1.62

``` r
# number of difference to make residuals stationary
unitroot_ndiffs(residuals)
```

    ## ndiffs 
    ##      0

According to the KPSS and Engle-Granger tests, the residuals are
stationary, meaning a linear combination of the share price and the
commodities price is stationary i.e. the prices are cointegrated. This
means that there is a long term relationship between the two time
series, where there is short term deviations around a mean and these
spreads will revert back to the mean over time. This supports the notion
that the price of the underlying commodities is what, in the long run,
drives the performance of the commodity producing firm.
