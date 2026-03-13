###### Libraries ######
library(readxl)
library(dplyr)

###### Functions ######
ewma_beta <- function(index_returns, stock_returns, lambda) {
  # This function assumes the first value of returns is NA
  
  # Check returns are the same length 
  if (length(index_returns) != length(stock_returns)) {
    stop("Index returns and stock returns have different lengths")
  }
  
  # Initialize variables for EWMA calc
  len <- length(index_returns)
  ewma_var_index <- numeric(len)
  ewma_covar <- numeric(len)
  
  # Use overall variance as starting values
  ewma_var_index[1] <- var(index_returns[-1])
  ewma_covar[1] <- cov(stock_returns[-1], index_returns[-1])
  
  # Calculate EWMA variances & covariance
  for (i in 2:len) {
    ewma_var_index[i] <- lambda * ewma_var_index[i - 1] + (1 - lambda) * index_returns[i]^2
    ewma_covar[i] <- lambda * ewma_covar[i-1] + (1 - lambda) * stock_returns[i] * index_returns[i]
  }
  
  # Calculate beta
  ewma_beta <- ewma_covar[len]/ewma_var_index[len]
  
  return(c(ewma_var_index[len], ewma_beta))
}

ols_beta <- function(index_returns, stock_returns) {
  return(cov(index_returns, stock_returns, use = "pairwise.complete.obs")/var(index_returns, na.rm = T))
}

###### Data ######
data_daily <- read_excel("data/q1.xlsx", sheet = "Daily")
data_weekly <- read_excel("data/q1.xlsx", sheet = "Weekly")

data_daily <- data_daily %>% mutate(
  index_returns = log(Index) - lag(log(Index), 1),
  stock_returns = log(Stock) - lag(log(Stock), 1)
)

data_weekly <- data_weekly %>% mutate(
  index_returns = log(Index) - lag(log(Index), 1),
  stock_returns = log(Stock) - lag(log(Stock), 1)
)


###### Betas ######
# a) OLS estimation using weekly data since 31 December 2001 (all data)
beta_a <- ols_beta(index_returns = data_weekly$index_returns, stock_returns = data_weekly$stock_returns)
# b) OLS estimation using weekly data since 28 December 2006
beta_b <- ols_beta(index_returns = data_weekly[data_weekly$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$index_returns, stock_returns = data_weekly[data_weekly$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$stock_returns)
# OLS estimation using daily data since 31 December 2001 (all data)
beta_c <- ols_beta(index_returns = data_daily$index_returns, stock_returns = data_daily$stock_returns)
# OLS estimation using daily data since 28 December 2006
beta_d <- ols_beta(index_returns = data_daily[data_daily$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$index_returns, stock_returns = data_daily[data_daily$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$stock_returns)

# EWMA estimation using weekly data with a smoothing constant of 0.95 (all data)
e <- ewma_beta(index_returns = data_weekly$index_returns, stock_returns = data_weekly$stock_returns, lambda = 0.95)
beta_e <- e[2]
# EWMA estimation using weekly data with a smoothing constant of 0.9 (all data)
f <- ewma_beta(index_returns = data_weekly$index_returns, stock_returns = data_weekly$stock_returns, lambda = 0.9)
beta_f <- f[2]
# EWMA estimation using daily data with a smoothing constant of 0.95 (all data)
g <- ewma_beta(index_returns = data_daily$index_returns, stock_returns = data_daily$stock_returns, lambda = 0.95)
beta_g <- g[2]
# EWMA estimation using daily data with a smoothing constant of 0.9 (all data)
h <- ewma_beta(index_returns = data_daily$index_returns, stock_returns = data_daily$stock_returns, lambda = 0.9)
beta_h <- h[2]

###### Yearly Vol ######
vol_a <- sd(data_weekly$index_returns,na.rm = T)*sqrt(52)
vol_b <- sd(data_weekly[data_weekly$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$index_returns, na.rm = T)*sqrt(52)
vol_c <- sd(data_daily$index_returns,na.rm = T)*sqrt(250)
vol_d <- sd(data_daily[data_daily$Date >= as.POSIXct("08/01/2016", format = "%d/%m/%Y", tz = "UTC"), ]$index_returns, na.rm = T)*sqrt(250)

vol_e <- sqrt(e[1]*52)
vol_f <- sqrt(f[1]*52)
vol_g <- sqrt(g[1]*250)
vol_h <- sqrt(h[1]*250)

###### VaR ######
var_a <- -qnorm(0.01)*beta_a*vol_a*sqrt(10/250)
var_b <- -qnorm(0.01)*beta_b*vol_b*sqrt(10/250)
var_c <- -qnorm(0.01)*beta_c*vol_c*sqrt(10/250)
var_d <- -qnorm(0.01)*beta_d*vol_d*sqrt(10/250)
var_e <- -qnorm(0.01)*beta_e*vol_e*sqrt(10/250)
var_f <- -qnorm(0.01)*beta_f*vol_f*sqrt(10/250)
var_g <- -qnorm(0.01)*beta_g*vol_g*sqrt(10/250)
var_h <- -qnorm(0.01)*beta_h*vol_h*sqrt(10/250)

###### Results Output ######
### Beta estimates ###
beta_a
beta_b
beta_c
beta_d
beta_e
beta_f
beta_g
beta_h

### Index Vol Estimate ###
vol_a
vol_b
vol_c
vol_d
vol_e
vol_f
vol_g
vol_h

### VaR % ###
var_a
var_b
var_c
var_d
var_e
var_f
var_g
var_h

### Rand VaR ###
var_a*1000000
var_b*1000000
var_c*1000000
var_d*1000000
var_e*1000000
var_f*1000000
var_g*1000000
var_h*1000000


