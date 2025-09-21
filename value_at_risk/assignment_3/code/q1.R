# Question 1: 
# 1. Repeat example IV.4.7 for a 5-day VaR with:
#   (a) A constant volatility of 20%
#   (b) The GARCH parameters as in the textbook.

###### GLOBAL VARIABLES ######
n <- 1000000
h <- 5
levels <- c(0.001, 0.01, 0.05, 0.1)
norm_sims <- matrix(rnorm(n*h),nrow = n, ncol = h)


###### CONSTANT VOLATILITY ######
### PARAMETERS ###
vol <- 0.25
daily_vol <- vol/sqrt(250)

# daily log-returns
returns_cv <- norm_sims*daily_vol

# h-day returns
returns_cv_hday <- rowSums(returns_cv) # Distribution of returns

# Value at Risk
VaR_cv <- -quantile(returns_cv_hday, probs = levels)




###### A-GARCH  #######
### PARAMETERS ###
omega <- 4*10^-6
alpha <- 0.06
lambda <- 0.01
beta <- 0.9
current_ret <- 0.1 # Change to -0.1 and run again from here to get initial negative shock
vol_unc_daily <- sqrt((omega + (lambda^2)*alpha)/(1-(alpha+beta)))



### VARIABLES ###
returns_agarch <- matrix(nrow = n, ncol = h)
vol_agarch <- matrix(nrow = n, ncol = h)


### SIMULATION ###
vol_agarch[,1] <- sqrt(omega + alpha*(current_ret - lambda)^2 + beta*(vol_unc_daily^2))
returns_agarch[,1] <- vol_agarch[,1]*norm_sims[,1]

for (i in 2:h) {
  vol_agarch[,i] <-sqrt(omega + alpha*(returns_agarch[,i-1] - lambda)^2 + beta*(vol_agarch[,i-1]^2))
  returns_agarch[,i] <- vol_agarch[,i]*norm_sims[,i]
}

returns_agarch_hday <- rowSums(returns_agarch)
VaR_agarch <- -quantile(returns_agarch_hday, probs = levels)


###### OUTPUT VALUES ######
VaR_cv
VaR_agarch














