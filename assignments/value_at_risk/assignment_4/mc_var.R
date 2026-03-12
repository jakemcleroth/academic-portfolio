################## Libraries ##################
library(readxl)
library(dplyr)



################## Functions ##################

bsm_futures_price <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("F, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  d_2 <- d_1 - sigma * sqrt(tau)
  price <- omega * exp(-r * tau) * (f * pnorm(omega * d_1) - k * pnorm(omega * d_2))
  
  return(price)
}

bsm_spot_price <- function(s, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 (call) or -1 (put)")
  }
  
  if (!all(length(s) == length(r), length(r) == length(sigma))) {
    stop("s, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(s / k) + (r + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
  d_2 <- d_1 - sigma * sqrt(tau)
  
  price <- omega * (s * pnorm(omega * d_1) - k * exp(-r * tau) * pnorm(omega * d_2))
  
  return(price)
}

calculate_implied_vol <- function(f, k, r, t_days, market_price, omega) {
  T_ <- t_days / 365
  
  implied_vol_function <- function(sigma) {
    bsm_futures_price(f, k, r, sigma, t = 0, T_, omega) - market_price
  }
  
  result <- uniroot(implied_vol_function, c(0.0001, 5), tol = 1e-8)$root
  return(result)
}

discounted_pnl <- function(h, r, bsm_h, bsm_t, sign) {
  if (sign != 1 && sign != -1) {
    stop("sign must be 1 or -1")
  }
  return(sign*(exp(-r*(h/365))*bsm_h - bsm_t))
}

value_at_risk <- function(pnl, pv, alpha, h) {
  return(-pv*quantile(pnl, probs = alpha, na.rm = T)*sqrt(h))
}

etl <- function(pnl, pv, alpha, h) {
  threshold <- quantile(pnl, probs = alpha, na.rm = T)
  excedences <- pnl[pnl<threshold]
  return(-pv*mean(excedences, na.rm=T)*sqrt(h))
}

################## Data ##################
data <- read_excel("sa_data.xlsx")

data <- data %>%
  mutate(
    # 1 day returns
    log_Discount_1d = log(Discount) - lag(log(Discount), 1),
    log_Index_1d    = log(Index)    - lag(log(Index), 1),
    log_Vol_1d      = log(Vol)      - lag(log(Vol), 1)
  )



################## IV.5.11 ##################
### Parameters ###
index_price <- as.numeric(data[nrow(data),2]) - 10
point_value <- 250
risk_free_rate <- as.numeric(data[nrow(data),4])
maturity_days <- 30
strike <- as.numeric(data[nrow(data),2])
call_price <- 1700
put_price <- 1700
imp_vol_put <- calculate_implied_vol(index_price, strike, risk_free_rate, maturity_days, put_price, -1) 
imp_vol_call <- calculate_implied_vol(index_price, strike, risk_free_rate, maturity_days, call_price, 1) 
mean_ret <- mean(data$log_Index_1d, na.rm = T)
stdev_ret <- sd(data$log_Index_1d, na.rm = T)
n_sims <- 1000000

### Sim ###
sim <- data.frame(index_ret_t = numeric(n_sims), index_ret_norm = numeric(n_sims))
set.seed(123)
sim$index_ret_t <- rt(n = n_sims, df = 6)*stdev_ret + mean_ret
sim$index_ret_norm <-  rnorm(n = n_sims, mean = mean_ret, sd = stdev_ret)



### Getting vol equation ###
reg <-  lm(log_Vol_1d ~ 0 + log_Index_1d + I(log_Index_1d^2), data = data)
beta_1 <- as.numeric(reg$coefficients[1])
beta_2 <-  as.numeric(reg$coefficients[2])
se <- sd(reg$residuals)*rnorm(n_sims)


### Calculating Vol ###
sim$vol_ret_t_quad <- beta_1*sim$index_ret_t + beta_2*(sim$index_ret_t)^2 + se
sim$vol_ret_t_lin <- beta_1*sim$index_ret_t + se
sim$vol_ret_norm_quad <- beta_1*sim$index_ret_norm + beta_2*(sim$index_ret_norm)^2 + se
sim$vol_ret_norm_lin <- beta_1*sim$index_ret_norm + se


### Simulated Prices and Vol ###
sim$index_t <- index_price*exp(sim$index_ret_t)
sim$index_norm <- index_price*exp(sim$index_ret_norm)

sim$call_vol_t_lin <- imp_vol_call*exp(sim$vol_ret_t_lin)
sim$call_vol_t_quad <- imp_vol_call*exp(sim$vol_ret_t_quad)
sim$call_vol_norm_lin <- imp_vol_call*exp(sim$vol_ret_norm_lin)
sim$call_vol_norm_quad <- imp_vol_call*exp(sim$vol_ret_norm_quad)

sim$put_vol_t_lin <- imp_vol_put*exp(sim$vol_ret_t_lin)
sim$put_vol_t_quad <- imp_vol_put*exp(sim$vol_ret_t_quad)
sim$put_vol_norm_lin <- imp_vol_put*exp(sim$vol_ret_norm_lin)
sim$put_vol_norm_quad <- imp_vol_put*exp(sim$vol_ret_norm_quad)

### Simulated Options Prices ###
sim$call_price_t_lin <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$call_vol_t_lin, t = 1/365, T_ = maturity_days/365, omega = 1)
sim$call_price_t_quad <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$call_vol_t_quad, t = 1/365, T_ = maturity_days/365, omega = 1)
sim$call_price_norm_lin  <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$call_vol_norm_lin, t = 1/365, T_ = maturity_days/365, omega = 1)
sim$call_price_norm_quad  <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$call_vol_norm_quad, t = 1/365, T_ = maturity_days/365, omega = 1)

sim$put_price_t_lin <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$put_vol_t_lin, t = 1/365, T_ = maturity_days/365, omega = -1)
sim$put_price_t_quad <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$put_vol_t_quad, t = 1/365, T_ = maturity_days/365, omega = -1)
sim$put_price_norm_lin  <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$put_vol_norm_lin, t = 1/365, T_ = maturity_days/365, omega = -1)
sim$put_price_norm_quad  <- bsm_futures_price(sim$index_t, strike, r = rep(risk_free_rate, n_sims), sigma = sim$put_vol_norm_quad, t = 1/365, T_ = maturity_days/365, omega = -1)

### PnL Distributions ###
sim$long_call_pnl_t_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_t_lin, bsm_t = call_price, sign = 1)
sim$short_call_pnl_t_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_t_lin, bsm_t = call_price, sign = -1)
sim$long_call_pnl_t_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_t_quad, bsm_t = call_price, sign = 1)
sim$short_call_pnl_t_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_t_quad, bsm_t = call_price, sign = -1)
sim$long_call_pnl_norm_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_norm_lin, bsm_t = call_price, sign = 1)
sim$short_call_pnl_norm_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_norm_lin, bsm_t = call_price, sign = -1)
sim$long_call_pnl_norm_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_norm_quad, bsm_t = call_price, sign = 1)
sim$short_call_pnl_norm_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$call_price_norm_quad, bsm_t = call_price, sign = -1)

sim$long_put_pnl_t_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_t_lin, bsm_t = put_price, sign = 1)
sim$short_put_pnl_t_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_t_lin, bsm_t = put_price, sign = -1)
sim$long_put_pnl_t_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_t_quad, bsm_t = put_price, sign = 1)
sim$short_put_pnl_t_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_t_quad, bsm_t = put_price, sign = -1)
sim$long_put_pnl_norm_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_norm_lin, bsm_t = put_price, sign = 1)
sim$short_put_pnl_norm_lin <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_norm_lin, bsm_t = put_price, sign = -1)
sim$long_put_pnl_norm_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_norm_quad, bsm_t = put_price, sign = 1)
sim$short_put_pnl_norm_quad <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = sim$put_price_norm_quad, bsm_t = put_price, sign = -1)


### VaR ###
# student t quadratic vol
value_at_risk(sim$long_call_pnl_t_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_call_pnl_t_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$long_put_pnl_t_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_put_pnl_t_quad, pv = point_value, alpha = 0.01, h = 1)

# student t lin vol
value_at_risk(sim$long_call_pnl_t_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_call_pnl_t_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$long_put_pnl_t_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_put_pnl_t_lin, pv = point_value, alpha = 0.01, h = 1)

# normal quadratic vol
value_at_risk(sim$long_call_pnl_norm_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_call_pnl_norm_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$long_put_pnl_norm_quad, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_put_pnl_norm_quad, pv = point_value, alpha = 0.01, h = 1)

# normal lin vol
value_at_risk(sim$long_call_pnl_norm_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_call_pnl_norm_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$long_put_pnl_norm_lin, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(sim$short_put_pnl_norm_lin, pv = point_value, alpha = 0.01, h = 1)

### Remove variables for next question ###
rm(index_price, point_value, risk_free_rate, maturity_days, strike, call_price, put_price,
   imp_vol_put, imp_vol_call, mean_ret, stdev_ret, n_sims, sim, beta_1, beta_2, se, reg)






################## IV.5.14 ##################
data <- read_excel("sa_data.xlsx")

data <- data %>%
  mutate(
    # 1 day returns
    log_Discount_1d = log(Discount) - lag(log(Discount), 1),
    log_Index_1d    = log(Index)    - lag(log(Index), 1),
    log_Vol_1d      = log(Vol)      - lag(log(Vol), 1)
  )

data <- data %>% filter_all(all_vars(is.finite(.)))
### Parameters ###
index_price <- as.numeric(data[nrow(data),2])
point_value <- 250
risk_free_rate <- as.numeric(data[nrow(data),4])
maturity_days <- 90
strike <- as.numeric(data[nrow(data),2]) - 15
imp_vol <- 0.2
call_price <- bsm_spot_price(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, 1)
n_sims <- 1000000
h <- 10


### mean, variance/covariance & Cholesky matrix ###
mu_1d <- as.vector(colMeans(data[2:nrow(data),c(5,6,7)]))
mu_10d <- mu_1d*h
names(mu_1d) <- c("log_Discount_1d", "log_Index_1d", "log_Vol_1d")
names(mu_10d) <- c("log_Discount_10d", "log_Index_10d", "log_Vol_10d")

vc_1d <- cov(data[2:nrow(data),c(5,6,7)])
vc_10d <- vc_1d*h
colnames(vc_10d) <- c("log_Discount_10d", "log_Index_10d", "log_Vol_10d")
rownames(vc_10d) <- c("log_Discount_10d", "log_Index_10d", "log_Vol_10d")

chol_1d <- chol(vc_1d)
chol_10d <- chol(vc_10d)


### 10-day Sim ###
set.seed(123)
ret_sims_10d <- matrix(rnorm(3*n_sims), ncol = 3)%*%chol_10d + matrix(rep(as.vector(mu_10d),n_sims), nrow = n_sims, ncol = 3, byrow = TRUE)
value_sims_10d <- exp(ret_sims_10d)%*%diag(c(risk_free_rate, index_price, imp_vol))
price_sims_10d <- bsm_spot_price(value_sims_10d[,2], strike, value_sims_10d[,1], value_sims_10d[,3], 10/365, maturity_days/365, 1)
pnl_sims_10d <- discounted_pnl(h, risk_free_rate, price_sims_10d, call_price, 1)
value_at_risk(pnl_sims_10d, 250, 0.01, 1) # set h = 1 since we are already using 10 day returns


### 1-day Sim ###
ret_sims_1d_1 <- matrix(rnorm(3*n_sims), ncol = 3)%*%chol_1d + matrix(rep(as.vector(mu_1d),n_sims), nrow = n_sims, ncol = 3, byrow = TRUE)
ret_sims_1d_cumulative <- ret_sims_1d_1

for (i in 2:10) {
  ret_sims_1d_cumulative <- ret_sims_1d_cumulative + matrix(rnorm(3*n_sims), ncol = 3)%*%chol_1d + matrix(rep(as.vector(mu_1d),n_sims), nrow = n_sims, ncol = 3, byrow = TRUE)
}
value_sims_10d_multi <- exp(ret_sims_1d_cumulative)%*%diag(c(risk_free_rate, index_price, imp_vol))
price_sims_10d_multi <- bsm_spot_price(value_sims_10d_multi[,2], strike, value_sims_10d_multi[,1], value_sims_10d_multi[,3], 10/365, maturity_days/365, 1)
pnl_sims_10d_multi <- discounted_pnl(h, risk_free_rate, price_sims_10d_multi, call_price, 1)
value_at_risk(pnl_sims_10d_multi, 250, 0.01, 1) # set h = 1 since we are already using 10 day returns



### Remove variables for next question ###
rm(index_price, point_value, risk_free_rate, maturity_days, strike, imp_vol, call_price, n_sims, h,
   mu_1d, mu_10d, vc_1d, vc_10d, chol_1d, chol_10d,
   ret_sims_10d, value_sims_10d, price_sims_10d, pnl_sims_10d,
   ret_sims_1d_1, ret_sims_1d_cumulative, value_sims_10d_multi, price_sims_10d_multi, pnl_sims_10d_multi, i, data)



################## IV.5.20 ##################

### Parameters ###
h <- 10
risk_free_rate <- 0.05
### Index values ###
value_ftse100 <- 6220.80
value_sp500 <- 1418.30
value_dax <- 6596.92
value_vftse <- 0.1280
value_vix <- 0.1156
value_vdax <- 0.1385

### Point Values ###
point_value_ftse100 <- 10     
point_value_sp500 <- 125.00    
point_value_dax <- 3.75      

### Exchange Rates ###
exch_us  = 0.5
exch_eur = 0.75

### VC Matrix ###
cv_mat_1d <- matrix(c(
  0.000134745, 6.34148e-05, 0.000146388, -0.000490705, -0.00024177, -0.000331826,
  6.34148e-05, 0.000133116, 0.000118051, -0.000259339, -0.000456837, -0.000243348,
  0.000146388, 0.000118051, 0.000279472, -0.000591391, -0.000414085, -0.00055151,
  -0.000490705, -0.000259339, -0.000591391, 0.003444646, 0.001468086, 0.002039951,
  -0.00024177, -0.000456837, -0.000414085, 0.001468086, 0.00288307, 0.001327802,
  -0.000331826, -0.000243348, -0.00055151, 0.002039951, 0.001327802, 0.002058965
), nrow = 6, byrow = TRUE)
rownames(cv_mat_1d) <- c("FTSE 100", "S&P 500", "DAX 30", "Vftse", "Vix", "Vdax")
colnames(cv_mat_1d) <- c("FTSE 100", "S&P 500", "DAX 30", "Vftse", "Vix", "Vdax")
cv_mat_10d <- cv_mat_1d*10
chol_10d <- chol(cv_mat_10d)

### Greeks ###
# Delta
delta_ftse <- -0.5
delta_sp <- -0.2
delta_dax <- 0.7

# Gamma
gamma_ftse <- -0.005
gamma_sp <- -0.001
gamma_dax <- 0.004

# Vega
vega_ftse <- -150
vega_sp <- -100
vega_dax <- 200

### Simulate values ###
n_sims <- 1000000
ret_sims_10d <- matrix(rnorm(6*n_sims), ncol = 6)%*%chol_10d

delta_pnl <- exp(-risk_free_rate*h/365)*(delta_ftse*ret_sims_10d[,1]*point_value_ftse100*value_ftse100 + delta_sp*ret_sims_10d[,2]*point_value_sp500*value_sp500 + delta_dax*ret_sims_10d[,3]*point_value_dax*value_dax)
value_at_risk(delta_pnl, 1, 0.01, 1)
delta_gamma_pnl <- delta_pnl + 0.5*exp(-risk_free_rate*h/365)*(gamma_ftse*(ret_sims_10d[,1]^2)*point_value_ftse100*(value_ftse100)^2 + gamma_sp*(ret_sims_10d[,2])^2*point_value_sp500*(value_sp500)^2 + gamma_dax*(ret_sims_10d[,3])^2*point_value_dax*(value_dax)^2)
value_at_risk(delta_gamma_pnl, 1, 0.01, 1)
delta_gamma_vega_pnl <- delta_gamma_pnl + exp(-risk_free_rate*h/365)*(vega_ftse*(ret_sims_10d[,4])*point_value_ftse100 + vega_sp*(ret_sims_10d[,5])*point_value_sp500 + vega_dax*(ret_sims_10d[,6])^2*point_value_dax)
value_at_risk(delta_gamma_vega_pnl, 1, 0.01, 1)
