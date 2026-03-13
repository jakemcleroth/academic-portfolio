# 2. Repeat example IV.4.8 with your own values for the parameters and the correlation matrix.
###### VARIABLES ######
h <- 10 # No. of days
n <- 1000000 # No. of sims

# Excess expected returns
mu <- c(0.06, 0.04, 0.03, 0.02, 0.01)
h_day_mu <- mu * h / 250

# Risk factor vol (Annual)
sigma <- c(0.15, 0.20, 0.12, 0.10, 0.18)

# Risk factor sensitivities
theta <- c(0.8, 0.4, 0.3, 0.1, -0.1)

# Correlations
corr_matrix <- matrix(c(
  1.00,  0.30,  0.20,  0.10,  0.15,
  0.30,  1.00,  0.25,  0.20,  0.10,
  0.20,  0.25,  1.00,  0.30,  0.25,
  0.10,  0.20,  0.30,  1.00,  0.35,
  0.15,  0.10,  0.25,  0.35,  1.00
), nrow = 5, byrow = TRUE)

# Covariance 
cov_matrix <- diag(sigma) %*% corr_matrix %*% diag(sigma) # Annual
h_day_cov_matrix <- cov_matrix*10/250 # h-day

# h-day Cholesky matrix
h_day_chol_matrix <- chol(h_day_cov_matrix) # upper triangular



###### SIMULATION ######
norm_sims <- matrix(rnorm(n*5), nrow = n, ncol = 5)
h_day_factor_returns <- norm_sims%*%h_day_chol_matrix + matrix(rep(h_day_mu,n), nrow = n, ncol = 5, byrow = TRUE)

# Check covariance of sim
cov_matrix
cov(h_day_factor_returns)*25

# VaR
h_day_portfolio_returns <- h_day_factor_returns%*%theta
h_day_VaR <- -quantile(h_day_portfolio_returns, probs = 0.01)
h_day_VaR












