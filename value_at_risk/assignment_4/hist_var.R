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

bsm_futures_delta <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("F, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  
  delta <- omega * exp(-r * tau) * pnorm(omega * d_1)
  
  return(delta)
}

bsm_futures_gamma <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  
  gamma <- exp(-r * tau) * dnorm(d_1) / (f * sigma * sqrt(tau))
  
  return(gamma)
}

bsm_futures_vega <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  
  vega <- exp(-r * tau) * f * dnorm(d_1) * sqrt(tau)
  
  return(vega)
}

etl <- function(pnl, pv, alpha, h) {
  threshold <- quantile(pnl, probs = alpha, na.rm = T)
  excedences <- pnl[pnl<threshold]
  return(-pv*mean(excedences, na.rm=T)*sqrt(h))
}

bsm_futures_theta <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  d_2 <- d_1 - sigma * sqrt(tau)
  
  theta <- exp(-r * tau) * (
    - (f * dnorm(d_1) * sigma) / (2 * sqrt(tau)) + 
      omega * r * f * pnorm(omega * d_1) -
      omega * r * k * pnorm(omega * d_2)
  )
  
  return(theta)
}

bsm_futures_rho <- function(f, k, r, sigma, t, T_, omega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  d_2 <- d_1 - sigma * sqrt(tau)
  
  rho <- -tau * exp(-r * tau) * (omega * f * pnorm(omega * d_1) - omega * k * pnorm(omega * d_2))
  
  return(rho)
}


bsm_futures_volga <- function(f, k, r, sigma, t, T_, omega, vega) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  d_2 <- d_1 - sigma * sqrt(tau)
  
  volga <- vega * d_1 * d_2 / sigma
  
  return(volga)
}

bsm_futures_vanna <- function(f, k, r, sigma, t, T_, omega, volga) {
  if (omega != 1 && omega != -1) {
    stop("omega must be 1 or -1")
  }
  
  if (!all(length(f) == length(r), length(r) == length(sigma))) {
    stop("f, r, and sigma must be vectors of the same length")
  }
  
  tau <- T_ - t
  if (tau <= 0) stop("T_ must be greater than t")
  
  d_1 <- (log(f / k) + 0.5 * sigma^2 * tau) / (sigma * sqrt(tau))
  
  vanna <- - volga / (f * sqrt(tau)*d_1 )
  
  return(vanna)
}


################## Data ##################
data <- read_excel("/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_4/data/sa_data.xlsx")

data <- data %>%
  mutate(
    # 1 day returns
    log_Discount_1d = log(Discount) - lag(log(Discount), 1),
    log_Index_1d    = log(Index)    - lag(log(Index), 1),
    log_Vol_1d      = log(Vol)      - lag(log(Vol), 1),
    
    # 10 day overlapping returns
    log_Discount_10d = log(Discount) - lag(log(Discount), 10),
    log_Index_10d    = log(Index)    - lag(log(Index), 10),
    log_Vol_10d      = log(Vol)      - lag(log(Vol), 10)
  )


################## IV.5.3 ##################
### Parameters ###
index_price <- as.numeric(data[nrow(data),2])
point_value <- 250
volatility <- 0.20          
risk_free_rate <- 0.00
maturity_days <- 30
strike <- as.numeric(data[nrow(data),2])

### Option Prices ###
call_price <- bsm_futures_price(f = index_price,k = strike,r = risk_free_rate,sigma = volatility, t = 0, T_ = maturity_days/365, omega = 1)
put_price <- bsm_futures_price(f = index_price,k = strike,r = risk_free_rate,sigma = volatility, t = 0, T_ = maturity_days/365, omega = -1)

### Historical Simulation ###
data$discount_hist_sim_1d <- risk_free_rate*exp(data$log_Discount_1d)
data$discount_hist_sim_10d <- risk_free_rate*exp(data$log_Discount_10d)
data$index_hist_sim_1d <- index_price*exp(data$log_Index_1d)
data$index_hist_sim_10d <- index_price*exp(data$log_Index_10d)
data$vol_hist_sim_1d <- volatility*exp(data$log_Vol_1d)
data$vol_hist_sim_10d <- volatility*exp(data$log_Vol_10d)


### Pricing Using Simulated Values ###
data$call_prices_1d <- bsm_futures_price(data$index_hist_sim_1d, k = strike, r = data$discount_hist_sim_1d, sigma = data$vol_hist_sim_1d, t = 1/365, T_ = 30/365, omega = 1)
data$put_prices_1d <- bsm_futures_price(data$index_hist_sim_1d, k = strike, r = data$discount_hist_sim_1d, sigma = data$vol_hist_sim_1d, t = 1/365, T_ = 30/365, omega = -1)

data$call_prices_10d <- bsm_futures_price(data$index_hist_sim_10d, k = strike, r = data$discount_hist_sim_10d, sigma = data$vol_hist_sim_10d, t = 10/365, T_ = 30/365, omega = 1)
data$put_prices_10d <- bsm_futures_price(data$index_hist_sim_10d, k = strike, r = data$discount_hist_sim_10d, sigma = data$vol_hist_sim_10d, t = 10/365, T_ = 30/365, omega = -1)

### PnL Distributions ###
data$call_pnl_1d <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = data$call_prices_1d, bsm_t = call_price, sign = 1)
data$call_pnl_10d <- discounted_pnl(h = 10, r = risk_free_rate, bsm_h = data$call_prices_10d, bsm_t = call_price, sign = 1)
data$put_pnl_1d <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = data$put_prices_1d, bsm_t = put_price, sign = 1)
data$put_pnl_10d <- discounted_pnl(h = 10, r = risk_free_rate, bsm_h = data$put_prices_10d, bsm_t = put_price, sign = 1)

### 10-day dynamic var ###
value_at_risk(pnl = data$call_pnl_1d, pv = point_value, alpha = 0.01, h = 10)
value_at_risk(pnl = -data$call_pnl_1d, pv = point_value, alpha = 0.01, h = 10) # Short

value_at_risk(pnl = data$put_pnl_1d, pv = point_value, alpha = 0.01, h = 10)
value_at_risk(pnl = -data$put_pnl_1d, pv = point_value, alpha = 0.01, h = 10) # Short

### 10-day static var ###
value_at_risk(pnl = data$call_pnl_10d, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(pnl = -data$call_pnl_10d, pv = point_value, alpha = 0.01, h = 1) # Short

value_at_risk(pnl = data$put_pnl_10d, pv = point_value, alpha = 0.01, h = 1)
value_at_risk(pnl = -data$put_pnl_10d, pv = point_value, alpha = 0.01, h = 1) # Short

### Removing columns for next question ###
rm(index_price,point_value,volatility,risk_free_rate,maturity_days,strike,call_price,put_price)








################## IV.5.4 ##################
data <- data[,c(1:7)]
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

### Greeks ###
delta_put <- bsm_futures_delta(index_price, strike, risk_free_rate, imp_vol_put, 0, maturity_days/365, -1)
delta_call <- bsm_futures_delta(index_price, strike, risk_free_rate, imp_vol_call, 0, maturity_days/365, 1)

### Historical Simulation ###
data$discount_hist_sim_1d <- risk_free_rate*exp(data$log_Discount_1d)
data$index_hist_sim_1d <- index_price*exp(data$log_Index_1d)
data$put_vol_hist_sim_1d <- imp_vol_put*exp(data$log_Vol_1d)
data$call_vol_hist_sim_1d <- imp_vol_call*exp(data$log_Vol_1d)

### Pricing ###
data$put_prices_1d <- bsm_futures_price(data$index_hist_sim_1d, k = strike, r = data$discount_hist_sim_1d, sigma = data$put_vol_hist_sim_1d, t = 1/365, T_ = 30/365, omega = -1)
data$call_prices_1d <- bsm_futures_price(data$index_hist_sim_1d, k = strike, r = data$discount_hist_sim_1d, sigma = data$call_vol_hist_sim_1d, t = 1/365, T_ = 30/365, omega = 1)

## PnL Distributions ###
data$long_put_pnl_1d <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = data$put_prices_1d, bsm_t = put_price, sign = 1)
data$short_put_pnl_1d <- -data$long_put_pnl_1d
data$long_call_pnl_1d <- discounted_pnl(h = 1, r = risk_free_rate, bsm_h = data$call_prices_1d, bsm_t = call_price, sign = 1)
data$short_call_pnl_1d <- -data$long_call_pnl_1d 

data$delta_hedged_short_put_pnl_1d <- data$short_put_pnl_1d + delta_put*(exp(-risk_free_rate/365)*data$index_hist_sim_1d - index_price)
data$delta_hedged_long_put_pnl_1d <- data$long_put_pnl_1d - delta_put*(exp(-risk_free_rate/365)*data$index_hist_sim_1d - index_price)
data$delta_hedged_short_call_pnl_1d <- data$short_call_pnl_1d + delta_call*(exp(-risk_free_rate/365)*data$index_hist_sim_1d - index_price)
data$delta_hedged_long_call_pnl_1d <- data$long_call_pnl_1d - delta_call*(exp(-risk_free_rate/365)*data$index_hist_sim_1d - index_price)

### 1-day VaR short put###
value_at_risk(pnl = data$short_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
value_at_risk(pnl = data$delta_hedged_short_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day ETL short put ###
etl(pnl = data$short_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
etl(pnl = data$delta_hedged_short_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day VaR long put ###
value_at_risk(pnl = data$long_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
value_at_risk(pnl = data$delta_hedged_long_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day ETL long put ###
etl(pnl = data$long_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
etl(pnl = data$delta_hedged_long_put_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day VaR short call ###
value_at_risk(pnl = data$short_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
value_at_risk(pnl = data$delta_hedged_short_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day ETL short call ###
etl(pnl = data$short_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
etl(pnl = data$delta_hedged_short_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day VaR long call ###
value_at_risk(pnl = data$long_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
value_at_risk(pnl = data$delta_hedged_long_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### 1-day ETL long call ###
etl(pnl = data$long_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Unhedged
etl(pnl = data$delta_hedged_long_call_pnl_1d, pv = 250, alpha = 0.01, h = 1) # Hedged

### Removing columns for next question ###
rm(index_price,point_value,risk_free_rate,imp_vol_put,imp_vol_call,maturity_days,strike,call_price,put_price, delta_call, delta_put)






################## IV.5.6 ##################
data <- data[,c(1:7)]
### Parameters ###
option_1 <- 1800
strike_1 <- as.numeric(data[nrow(data),2]) 
omega_1 <- 1
position_1 <- -1
position_1_weight <- -1

option_2 <- 1700
strike_2 <- as.numeric(data[nrow(data),2]) - 15
omega_2 <- -1
position_2 <- 1

option_3 <- 1850
strike_3 <- as.numeric(data[nrow(data),2]) + 15
omega_3 <- 1
position_3 <- 1

index_price <- as.numeric(data[nrow(data),2]) - 10
point_value <- 250
risk_free_rate <- as.numeric(data[nrow(data),4])
maturity_1 <- 60
maturity_2 <- 30
maturity_3 <- 90

### Implied Vol ###
imp_vol_1 <-  calculate_implied_vol(index_price, strike_1, risk_free_rate, maturity_1, option_1, omega_1)
imp_vol_2 <-  calculate_implied_vol(index_price, strike_2, risk_free_rate, maturity_2, option_2, omega_2)
imp_vol_3 <-  calculate_implied_vol(index_price, strike_3, risk_free_rate, maturity_3, option_3, omega_3)

### Greeks ###
delta_1 <- position_1*bsm_futures_delta(index_price, strike_1, risk_free_rate, imp_vol_1, 0, maturity_1/365, omega_1)
delta_2 <- position_2*bsm_futures_delta(index_price, strike_2, risk_free_rate, imp_vol_2, 0, maturity_2/365, omega_2)
delta_3 <- position_3*bsm_futures_delta(index_price, strike_3, risk_free_rate, imp_vol_3, 0, maturity_3/365, omega_3)

gamma_1 <- position_1*bsm_futures_gamma(index_price, strike_1, risk_free_rate, imp_vol_1, 0, maturity_1/365, omega_1)
gamma_2 <- position_2*bsm_futures_gamma(index_price, strike_2, risk_free_rate, imp_vol_2, 0, maturity_2/365, omega_2)
gamma_3 <- position_3*bsm_futures_gamma(index_price, strike_3, risk_free_rate, imp_vol_3, 0, maturity_3/365, omega_3)

vega_1 <- position_1*bsm_futures_vega(index_price, strike_1, risk_free_rate, imp_vol_1, 0, maturity_1/365, omega_1)
vega_2 <- position_2*bsm_futures_vega(index_price, strike_2, risk_free_rate, imp_vol_2, 0, maturity_2/365, omega_2)
vega_3 <- position_3*bsm_futures_vega(index_price, strike_3, risk_free_rate, imp_vol_3, 0, maturity_3/365, omega_3)

### Hedging ###
position_2_weight <- -(gamma_1*vega_3 - gamma_3*vega_1)/(gamma_2*vega_3-gamma_3*vega_2) 
position_3_weight <- -(-gamma_1*vega_2+gamma_2*vega_1)/(gamma_2*vega_3-gamma_3*vega_2) 
portfolio_delta <- delta_1+position_2_weight*delta_2+position_3_weight*delta_3


### Historical Simulation ###
data$discount_hist_sim_1d <- risk_free_rate*exp(data$log_Discount_1d)
data$index_hist_sim_1d <- index_price*exp(data$log_Index_1d)
data$vol_hist_sim_1d_1 <- imp_vol_1*exp(data$log_Vol_1d)
data$vol_hist_sim_1d_2 <- imp_vol_2*exp(data$log_Vol_1d)
data$vol_hist_sim_1d_3 <- imp_vol_3*exp(data$log_Vol_1d)

### Pricing ###
data$prices_1d_1 <- bsm_futures_price(data$index_hist_sim_1d, k = strike_1, r = data$discount_hist_sim_1d, sigma = data$vol_hist_sim_1d_1, t = 1/365, T_ = maturity_1/365, omega = omega_1)
data$prices_1d_2 <- bsm_futures_price(data$index_hist_sim_1d, k = strike_2, r = data$discount_hist_sim_1d, sigma = data$vol_hist_sim_1d_2, t = 1/365, T_ = maturity_2/365, omega = omega_2)
data$prices_1d_3 <- bsm_futures_price(data$index_hist_sim_1d, k = strike_3, r = data$discount_hist_sim_1d, sigma = data$vol_hist_sim_1d_3, t = 1/365, T_ = maturity_3/365, omega = omega_3)

### PnL ###
port_pv <- position_1_weight*option_1 + position_2_weight*option_2 + position_3_weight*option_3
data$gamma_vega_hedged_port_pv <- position_1_weight*data$prices_1d_1 + position_2_weight*data$prices_1d_2  + position_3_weight*data$prices_1d_3
data$gamma_vega_hedged_pnl <- discounted_pnl(1, risk_free_rate, data$gamma_vega_hedged_port_pv, port_pv, 1)
data$delta_gamma_vega_hedged_pnl <- data$gamma_vega_hedged_pnl - portfolio_delta*discounted_pnl(1, risk_free_rate, data$index_hist_sim_1d, index_price, 1)


### 10-day VaR ###
value_at_risk(data$delta_gamma_vega_hedged_pnl, point_value, 0.01, 10)

### 10-day ETL ### 
etl(data$delta_gamma_vega_hedged_pnl, point_value, 0.01, 10)

### Removing columns for next question ###
rm(
  option_1, strike_1, omega_1, position_1, position_1_weight,
  option_2, strike_2, omega_2, position_2,
  option_3, strike_3, omega_3, position_3,
  index_price, point_value, risk_free_rate,
  maturity_1, maturity_2, maturity_3,
  imp_vol_1, imp_vol_2, imp_vol_3,
  delta_1, delta_2, delta_3,
  gamma_1, gamma_2, gamma_3,
  vega_1, vega_2, vega_3,
  position_2_weight, position_3_weight,
  portfolio_delta, port_pv
)



################## IV.5.7 ##################
data <- data[,c(1:4)]
data <- data %>%
  mutate(
    # 1 day returns
    log_Discount_1d = log(Discount) - lag(log(Discount), 1),
    log_Index_1d    = log(Index)    - lag(log(Index), 1),
    log_Vol_1d      = log(Vol)      - lag(log(Vol), 1),
    
    # 10 day overlapping returns
    log_Discount_10d = log(Discount) - lag(log(Discount), 10),
    log_Index_10d    = log(Index)    - lag(log(Index), 10),
    log_Vol_10d      = log(Vol)      - lag(log(Vol), 10)
  )

### Parameters ###
index_price <-  as.numeric(data[nrow(data),2]) - 10
point_value <- 250
risk_free_rate <- as.numeric(data[nrow(data),4])
maturity_days <- 30
strike <-  as.numeric(data[nrow(data),2])
put_price <- 1700
imp_vol <- calculate_implied_vol(index_price, strike, risk_free_rate, maturity_days, put_price, -1) 
omega <- -1 # change to 1 for long

### Greeks (Long Position) ###
delta <- bsm_futures_delta(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega)
gamma <- bsm_futures_gamma(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega)
vega <- bsm_futures_vega(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega)/100 # Excel divides by 100
theta <- bsm_futures_theta(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega)/365 # Excel divides by 365
volga <- bsm_futures_volga(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega, vega)
vanna <- bsm_futures_vanna(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega, volga)
rho <- bsm_futures_rho(index_price, strike, risk_free_rate, imp_vol, 0, maturity_days/365, omega)/100 # Excel divides by 100


### Historical Simulation ###
data$discount_hist_sim_1d <- risk_free_rate*exp(data$log_Discount_1d)
data$discount_hist_sim_10d <- risk_free_rate*exp(data$log_Discount_10d)
data$index_hist_sim_1d <- index_price*exp(data$log_Index_1d)
data$index_hist_sim_10d <- index_price*exp(data$log_Index_10d)
data$vol_hist_sim_1d <- imp_vol*exp(data$log_Vol_1d)
data$vol_hist_sim_10d <- imp_vol*exp(data$log_Vol_10d)

### changes ###
data <- data %>%
  mutate(
    # 1-day changes
    dr_1d  = discount_hist_sim_1d - risk_free_rate,
    ds_1d  = index_hist_sim_1d - index_price,
    dv_1d = vol_hist_sim_1d - imp_vol,
    
    # 10-day changes
    dr_10d  = discount_hist_sim_10d - risk_free_rate,
    ds_10d  = index_hist_sim_10d - index_price,
    dv_10d = vol_hist_sim_10d - imp_vol
  )

### PnL ###
data <- data %>%
  mutate(
    # 1-day
    pnl_d_1d     = exp(-risk_free_rate / 365) * delta * omega * ds_1d,
    pnl_dg_1d    = pnl_d_1d + 0.5 * exp(-risk_free_rate / 365) * gamma * ds_1d^2 * omega,
    pnl_dgv_1d   = pnl_dg_1d + exp(-risk_free_rate / 365) * vega * 100 * dv_1d * omega,
    pnl_dgvt_1d  = pnl_dgv_1d + exp(-risk_free_rate / 365) * theta * 365 * (1/365) * omega,
    pnl_all_1d   = pnl_dgvt_1d + exp(-risk_free_rate / 365) * omega *(
      rho * dr_1d * 100 +
        0.5 * volga * (dv_1d )^2  +
        vanna * ds_1d * dv_1d  
    ),
    
    # 10-day
    pnl_d_10d    = exp(-risk_free_rate * 10 / 365) * delta * omega * ds_10d,
    pnl_dg_10d   = pnl_d_10d + 0.5 * exp(-risk_free_rate * 10 / 365) * gamma * ds_10d^2 * omega,
    pnl_dgv_10d  = pnl_dg_10d + exp(-risk_free_rate * 10 / 365) * vega * 100 * dv_10d * omega,
    pnl_dgvt_10d = pnl_dgv_10d +  exp(-risk_free_rate * 10 / 365) * theta * 365 * (10/365) * omega,
    pnl_all_10d  = pnl_dgvt_10d + exp(-risk_free_rate * 10 / 365) * omega*(
      rho * dr_10d * 100 +
        0.5 * volga * (dv_10d)^2   +
        vanna * ds_10d * dv_10d 
    )
  )



### VaR ###
# Delta
value_at_risk(data$pnl_d_1d, point_value, 0.01, 10)
value_at_risk(data$pnl_d_10d, point_value, 0.01, 1)
# Delta Gamma
value_at_risk(data$pnl_dg_1d, point_value, 0.01, 10)
value_at_risk(data$pnl_dg_10d, point_value, 0.01, 1)
# Delta Gamma Vega
value_at_risk(data$pnl_dgv_1d, point_value, 0.01, 10)
value_at_risk(data$pnl_dgv_10d, point_value, 0.01, 1)
# Delta Gamma Vega Theta
value_at_risk(data$pnl_dgvt_1d, point_value, 0.01, 10)
value_at_risk(data$pnl_dgvt_10d, point_value, 0.01, 1)
# All
value_at_risk(data$pnl_all_1d, point_value, 0.01, 10)
value_at_risk(data$pnl_all_10d, point_value, 0.01, 1)

### For 1 day calc only ###
# Delta
value_at_risk(data$pnl_d_1d, point_value, 0.01, 1)
# Delta Gamma
value_at_risk(data$pnl_dg_1d, point_value, 0.01, 1)
# Delta Gamma Vega
value_at_risk(data$pnl_dgv_1d, point_value, 0.01, 1)
# Delta Gamma Vega Theta
value_at_risk(data$pnl_dgvt_1d, point_value, 0.01, 1)
# All
value_at_risk(data$pnl_all_1d, point_value, 0.01, 1)















