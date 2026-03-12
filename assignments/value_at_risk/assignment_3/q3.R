# 3. Repeat example IV.4.13 with your own values for the parameters and discuss your results.

###### VARIABLES ######
h <- 10 # No. of days
n <- 1000000 # No. of sims
crash_prob <- 0.05
weights <- c(1/3, 1/3, 1/3)

# Excess expected returns

mu_crash <- c(-0.4, -0.55, -0.65)
mu_ord <- c(0.07, 0.12, 0.1)
h_day_mu_crash <- mu_crash * h / 250
h_day_mu_ord <- mu_ord * h / 250

# Risk factor vol (Annual)
sigma_crash <- c(0.5, 0.6, 0.55)
sigma_ord <- c(0.2, 0.25, 0.22)

# Correlations
corr_matrix_crash <- matrix(c(
  1.00, 0.85, 0.80,
  0.85, 1.00, 0.75,
  0.80, 0.75, 1.00
), nrow = 3, byrow = TRUE)

corr_matrix_ord <- matrix(c(
  1.00, 0.40, 0.60,
  0.40, 1.00, 0.30,
  0.60, 0.30, 1.00
), nrow = 3, byrow = TRUE)

# Covariance 

cov_matrix_crash <- diag(sigma_crash) %*% corr_matrix_crash %*% diag(sigma_crash) # Annual
h_day_cov_matrix_crash <- cov_matrix_crash*10/250 # h-day
cov_matrix_ord <- diag(sigma_ord) %*% corr_matrix_ord %*% diag(sigma_ord) # Annual
h_day_cov_matrix_ord <- cov_matrix_ord*10/250 # h-day



# normal values to compare to
mu <- crash_prob * mu_crash + (1 - crash_prob) * mu_ord
h_day_mu <- mu * h / 250
cov_matrix <- crash_prob * cov_matrix_crash + (1 - crash_prob) * cov_matrix_ord
h_day_cov_matrix <- cov_matrix*10/250 # h-day


# h-day Cholesky matrix
h_day_chol_matrix <- chol(h_day_cov_matrix)
h_day_chol_matrix_crash <- chol(h_day_cov_matrix_crash) # upper triangular
h_day_chol_matrix_ord <- chol(h_day_cov_matrix_ord) # upper triangular


###### SIMULATION ######
norm_sims <- matrix(rnorm(n*3), nrow = n, ncol = 3)
crash <- rbinom(n = n, size = 1, prob = crash_prob)
h_day_stock_returns <-  norm_sims%*%h_day_chol_matrix + matrix(rep(h_day_mu,n), nrow = n, ncol = 3, byrow = TRUE)
h_day_stock_returns_crash <- norm_sims%*%h_day_chol_matrix_crash + matrix(rep(h_day_mu_crash,n), nrow = n, ncol = 3, byrow = TRUE)
h_day_stock_returns_ord <- norm_sims%*%h_day_chol_matrix_ord + matrix(rep(h_day_mu_ord,n), nrow = n, ncol = 3, byrow = TRUE)

# Check covariance of sim
cov_matrix
cov(h_day_stock_returns)*25

cov_matrix_crash
cov(h_day_stock_returns_crash)*25

cov_matrix_ord
cov(h_day_stock_returns_ord)*25

# VaR
h_day_portfolio_returns <- h_day_stock_returns%*%weights
h_day_portfolio_returns_crash <- h_day_stock_returns_crash%*%weights
h_day_portfolio_returns_ord <- h_day_stock_returns_ord%*%weights
h_day_portfolio_returns_mixture <- h_day_portfolio_returns_crash*crash + h_day_portfolio_returns_ord*(1-crash)
h_day_VaR <- -quantile(h_day_portfolio_returns, probs = c(0.001, 0.01, 0.05))
h_day_VaR_mixture <- -quantile(h_day_portfolio_returns_mixture, probs = c(0.001, 0.01, 0.05))
h_day_VaR
h_day_VaR_mixture



