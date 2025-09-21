#install.packages("rugarch")
# install.packages("readxl")
# install.packages("xts")
library(rugarch)
library(readxl)
library(xts)

########## Question 3 ##########
### Get Return Data ###
alsi <- as.data.frame(read_excel("/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_2/data/data_assignment_2.xlsx", sheet = "ALLSHARE"))
alsi$returns <- c(NA, diff(log(alsi$price)))
returns <- alsi$returns[-1]
dates <- alsi$date[-1]


### Fit GARCH(1,1) ###
init_var <- var(returns) # Initial variance (sigma_0) is square of mean return
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),  # GARCH(1,1)
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE)    # zero mean
)
fit = ugarchfit(spec = spec, data = returns, out.sample = 0, fit.control = list(rec.init = init_var)) #rec.init sets initial variance (sigma_0)
print(fit)
params_garch <- coef(fit)
mu_garch <- params_garch[[1]]
omega_garch <- params_garch[[2]]
alpha_garch <- params_garch[[3]]
beta_garch <- params_garch[[4]]
unconditional_vol_garch <- sqrt(uncvariance(fit)*250)
alsi$vol_garch <- c(NA,as.numeric(sigma(fit)*sqrt(250)))

### GARCH Adjusted Returns ###
alsi$ret_garch <- alsi$returns*alsi$vol_garch[nrow(alsi)]/alsi$vol_garch


### Fit A-GARCH(1,1) ###
agarch_loglik <- function(params, returns) {
  # Parameters
  mu <- params[1]  # Mean of returns
  omega <- params[2]  # Constant term in variance equation
  alpha <- params[3]  # Coefficient for shifted squared residual term
  lambda <- params[4]  # Shift parameter
  beta <- params[5]  # GARCH term
  
  n <- length(returns)
  sigma2 <- numeric(n)  # Conditional variance
  epsilon <- numeric(n)  # Residuals
  
  # Initialize conditional variance
  sigma2_0 <- var(returns) 
  r_0 <- mean(returns)
  epsilon_0 <- r_0 - mu
  sigma2[1] <- omega + alpha * (epsilon_0 - lambda)^2 + beta * sigma2_0
  epsilon[1] <- returns[1] - mu
  # Compute residuals and conditional variance
  for (t in 2:n) {
    epsilon[t] <- returns[t] - mu
    sigma2[t] <- omega + alpha * (epsilon[t-1] - lambda)^2 + beta * sigma2[t-1]
  }
  
  # Log-likelihood (assuming normal distribution)
  loglik <- -0.5 * sum(log(sigma2) + (epsilon^2) / sigma2)
  
  return(-loglik)  # Return negative log-likelihood for minimization
}

# Initial parameter guesses
params_init <- c(
  mu = mean(returns),  # Mean of returns
  omega = var(returns) * 0.1,  # Small fraction of unconditional variance
  alpha = 0.1,  # ARCH term
  lambda = 0.02,  # Shift parameter
  beta = 0.8  # GARCH term
)

# Minimize negative log-likelihood
fit <- optim(
  par = params_init,
  fn = agarch_loglik,
  returns = returns,
  method = "L-BFGS-B",
  lower = c(-Inf, 1e-6, 1e-6, -Inf, 1e-6),  # Ensure omega, alpha, beta > 0
  upper = c(Inf, Inf, Inf, Inf, 1 - 1e-6)  # Ensure alpha + beta < 1
)


print(fit$par)

# Parameters
mu_agarch <- fit$par[1]  
omega_agarch <- fit$par[2]  
alpha_agarch <- fit$par[3]  
lambda <- fit$par[4]  
beta_agarch <- fit$par[5]  

# Compute unconditional volatility
unconditional_vol_agarch <- as.numeric(sqrt(250*((omega_agarch + lambda^2*alpha_agarch)/ (1 - alpha_agarch - beta_agarch))))

# Initialize conditional volatility
n <- length(returns)
sigma2 <- numeric(n)  # Conditional vol
epsilon <- numeric(n)  # Residuals

# Initialize conditional variance
sigma2_0 <- var(returns)  # Use unconditional variance as initial value
r_0 <- mean(returns)
epsilon_0 <- r_0 - mu_agarch
sigma2[1] <- omega_agarch + alpha_agarch * (epsilon_0 - lambda)^2 + beta_agarch * sigma2_0
epsilon[1] <- returns[1] - mu_agarch

# Compute conditional volatility through time
for (t in 2:n) {
  epsilon[t] <- returns[t] - mu_agarch
  sigma2[t] <- omega_agarch + alpha_agarch * (epsilon[t-1] - lambda)^2 + beta_agarch * sigma2[t-1]
}

# Convert conditional variance to conditional volatility
conditional_vol <- sqrt(sigma2*250)
alsi$vol_agarch <- c(NA, conditional_vol)  # Add NA for the first row to align with dates





### A-GARCH Adjusted Returns ###
alsi$ret_agarch <- alsi$returns*alsi$vol_agarch[nrow(alsi)]/alsi$vol_agarch

### Plots ###
plot(x = alsi$date[-1],y = alsi$vol_agarch[-1], type = "l", xlab = "Date", ylab = "Volatility",main = "Conditional Volatility", col = "red", lwd = 2)
lines(x = alsi$date[-1], y = alsi$vol_garch[-1], col = "blue", lwd = 2)
legend("topright",  legend = c("AGARCH Conditional Volatility", "GARCH Conditional Volatility"), col = c("red", "blue"), lwd = 2, bty = "n")

# Calculate the range for the y-axis
y_min <- min(c(alsi$returns[-1], alsi$ret_garch[-1]), na.rm = TRUE) 
y_max <- max(c(alsi$returns[-1], alsi$ret_garch[-1]), na.rm = TRUE)  
y_buffer <- 0.2 * (y_max - y_min)  
y_lim <- c(y_min - y_buffer, y_max + y_buffer)

# Plot
plot(x = alsi$date[-1], y = alsi$returns[-1], type = "l", col = "green", lwd = 2, xlab = "Date", ylab = "Returns", main = "Returns vs Adjusted Returns", ylim = y_lim)
lines(x = alsi$date[-1], y = alsi$ret_garch[-1], col = "blue", lwd = 2)
legend("topright",  legend = c("Returns", "GARCH Adjusted Returns"), col = c("green", "blue"), lwd = 2, bty = "n")

# Calculate the range for the y-axis
y_min <- min(c(alsi$returns[-1], alsi$ret_agarch[-1]), na.rm = TRUE) 
y_max <- max(c(alsi$returns[-1], alsi$ret_agarch[-1]), na.rm = TRUE)  
y_buffer <- 0.2 * (y_max - y_min)  
y_lim <- c(y_min - y_buffer, y_max + y_buffer)

# Plot
plot(x = alsi$date[-1], y = alsi$returns[-1], type = "l", col = "green", lwd = 2, xlab = "Date", ylab = "Returns", main = "Returns vs Adjusted Returns", ylim = y_lim)
lines(x = alsi$date[-1], y = alsi$ret_agarch[-1], col = "red", lwd = 2)
legend("topright",  legend = c("Returns", "AGARCH Adjusted Returns"), col = c("green", "red"), lwd = 2, bty = "n")



### Find Historical VaR for adjusted returns ###
var_alphas <- c(0.001, 0.01, 0.05, 0.1)
hist_var <- numeric(length = 4)
garch_var <- numeric(length = 4)
agarch_var <- numeric(length = 4)
for (i in 1:length(var_alphas)) {
  hist_var[i] <- -quantile(alsi$returns[-1], probs = var_alphas[i])
  garch_var[i] <- -quantile(alsi$ret_garch[-1], probs = var_alphas[i])
  agarch_var[i] <- -quantile(alsi$ret_agarch[-1], probs = var_alphas[i])
}
hist_var
garch_var
agarch_var













########## Question 4 ##########
### Initial Values ###
sigma_0 <- 0.1/sqrt(250) # choose from: alsi$vol_agarch[nrow(alsi)]/sqrt(250) or 0.1/sqrt(250)
r_0 <- alsi$returns[nrow(alsi)]
sigma_1 <- sqrt(omega_agarch+alpha_agarch*(r_0-lambda)^2+beta_agarch*sigma_0^2)


eps_values <- sqrt(250)*alsi$returns[-1]/alsi$vol_agarch[-1]
set.seed(123)
bootstrap <- replicate(10, sample(eps_values, size = 1000, replace = TRUE))
boot_values <- as.data.frame(bootstrap)
colnames(boot_values) <- paste0("Bootstrap_", 1:10)


ret_1 <- boot_values[,1]*sigma_1
fhs_values <- as.data.frame(rep(sigma_1,1000))
colnames(fhs_values) <- c("sigma_1")
fhs_values$ret_1 <- fhs_values[,1]*boot_values[,1]

fhs_values$sigma_2 <- NA
fhs_values$ret_2 <- NA
fhs_values$sigma_3 <- NA
fhs_values$ret_3 <- NA
fhs_values$sigma_4 <- NA
fhs_values$ret_4 <- NA
fhs_values$sigma_5 <- NA
fhs_values$ret_5 <- NA
fhs_values$sigma_6 <- NA
fhs_values$ret_6 <- NA
fhs_values$sigma_7 <- NA
fhs_values$ret_7 <- NA
fhs_values$sigma_8 <- NA
fhs_values$ret_8 <- NA
fhs_values$sigma_9 <- NA
fhs_values$ret_9 <- NA
fhs_values$sigma_10 <- NA
fhs_values$ret_10 <- NA


for (i in 2:10) {
  fhs_values[[paste0("sigma_",i)]] <- sqrt(omega_agarch+alpha_agarch*(fhs_values[[paste0("ret_",i-1)]]-lambda)^2+beta_agarch*fhs_values[[paste0("sigma_",i-1)]]^2)
  fhs_values[[paste0("ret_",i)]] <- fhs_values[[paste0("sigma_",i)]]*boot_values[,i]
}

ret_columns <- grep("^ret_", colnames(fhs_values), value = TRUE)

fhs_values$ret_sum <- rowSums(fhs_values[, ret_columns])


agarch_10_day_VaR <-  sqrt(10)*(0.1/sqrt(250))*agarch_var/(alsi$vol_agarch[nrow(alsi)]/sqrt(250))   # sqrt(10)*agarch_var or sqrt(10)*(0.1/sqrt(250))*agarch_var/(alsi$vol_agarch[nrow(alsi)]/sqrt(250)) 
fhs_10_day_VaR <- numeric(length = 4)
for (i in 1:4) {
  fhs_10_day_VaR[i] <- -quantile(fhs_values$ret_sum, probs = var_alphas[i])
}

agarch_10_day_VaR
fhs_10_day_VaR



