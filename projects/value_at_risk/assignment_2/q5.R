# 3.6
# 1) unadjusted VaR (10 day)
# 2) using adjusted x tilda (3.15)
# 3) adjust portfolio returns
# 4) adjust stocks ignoring cor
library(readxl)
# Data
data <- as.data.frame(read_excel("data_assignment_2.xlsx", sheet = "ALL"))
data <- data[data$Date<= as.POSIXct("01/04/2020", format = "%d/%m/%Y", tz = "UTC"), ] # Uncomment for Pre-covid
data$ret_firstrand <- c(NA, diff(log(data$firstrand)))
data$ret_sappi <- c(NA, diff(log(data$sappi)))

# Plot Citigroup prices (left y-axis)
plot(data$Date, data$firstrand, type = "l", col = "blue", lwd = 2, xlab = "Date", ylab = "Firstrand Price", main = "Firstrand and Sappi Price Series")
par(new = TRUE)  # Allow overlay on the same plot
plot(data$Date, data$sappi, type = "l", col = "red", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(side = 4)  # Add right y-axis
mtext("Sappi Price", side = 4, line = 3)  # Label for right y-axis
legend("topright", legend = c("Firstrand", "Sappi"), col = c("blue", "red"), lwd = 2, bty = "n")

portfolio_weights <- c(0.3, 0.7)
lambda <- 0.94
alpha_levels <- c(0.001, 0.01, 0.05, 0.1)

# Get EWMA
data$ret_port <- data$ret_firstrand*portfolio_weights[1] + data$ret_sappi*portfolio_weights[2]
data$ewma_var_firstrand <- NA
data$ewma_vol_firstrand <- NA
data$ewma_var_firstrand[1] <- var(data$ret_firstrand[-1])
data$ewma_var_sappi <- NA
data$ewma_vol_sappi <- NA
data$ewma_var_sappi[1] <- var(data$ret_sappi[-1])
data$ewma_var_port <- NA
data$ewma_vol_port <- NA
data$ewma_var_port[1] <- var(data$ret_port[-1])
data$ewma_covar <- NA
data$ewma_covar[1] <- cov(data$ret_firstrand[-1], data$ret_sappi[-1])

for (i in 2:nrow(data)) {
  data$ewma_var_firstrand[i] <- lambda * data$ewma_var_firstrand[i - 1] + (1 - lambda) * data$ret_firstrand[i]^2
  data$ewma_var_sappi[i] <- lambda * data$ewma_var_sappi[i - 1] + (1 - lambda) * data$ret_sappi[i]^2
  data$ewma_var_port[i] <- lambda * data$ewma_var_port[i - 1] + (1 - lambda) * data$ret_port[i]^2
  data$ewma_covar[i] <- lambda * data$ewma_covar[i-1] + (1 - lambda) * data$ret_firstrand[i] * data$ret_sappi[i]
}

data$ewma_vol_firstrand[-1] <- sqrt(250*data$ewma_var_firstrand[-1])
data$ewma_vol_sappi[-1] <- sqrt(250*data$ewma_var_sappi[-1])
data$ewma_vol_port[-1] <- sqrt(250*data$ewma_var_port[-1])

# plot EWMA
plot(x = data$Date, y = data$ewma_vol_sappi, type = "l", col = "red", lwd = 2, xlab = "Date", ylab = "Volatility", main = "EWMA Volatility")
lines(x = data$Date, y = data$ewma_vol_firstrand,type = "l", col = "blue", lwd = 2)
legend("topright",  legend = c("Sappi", "FirsRand"), col = c("red", "blue"), lwd = 2, bty = "n")

# Unadjausted VaR
unadjusted_var <- numeric(length = 4)
for (i in 1:4){
  unadjusted_var[i] <- -quantile(data$ret_port[-1], probs = alpha_levels[i])*sqrt(10)
}
unadjusted_var

# Adjusting returns
data$vol_adj_ret_firstrand <- data$ret_firstrand*data$ewma_vol_firstrand[nrow(data)]/data$ewma_vol_firstrand
data$vol_adj_ret_sappi <- data$ret_sappi*data$ewma_vol_sappi[nrow(data)]/data$ewma_vol_sappi
data$vol_adj_ret_port <- data$ret_port*data$ewma_vol_port[nrow(data)]/data$ewma_vol_port

# VaR from adjusting portfolio returns (method b)
adjusted_port_return_var <- numeric(length = 4)
for (i in 1:4){
  adjusted_port_return_var[i] <- -quantile(data$vol_adj_ret_port[-1], probs = alpha_levels[i])*sqrt(10)
}
adjusted_port_return_var

# VaR from adjusting stock returns without accounting for correlation (method c)
data$stock_vol_adj_ret_port <- portfolio_weights[1]*data$vol_adj_ret_firstrand + portfolio_weights[2]*data$vol_adj_ret_sappi
adjusted_stock_ret_var <- numeric(length = 4)

for (i in 1:4) {
  adjusted_stock_ret_var[i] <- -quantile(data$stock_vol_adj_ret_port[-1], probs = alpha_levels[i])*sqrt(10)
}
adjusted_stock_ret_var

# VaR from adjusting stock returns by accounting for correlation (method a)
# We can use firstrand adjusted returns and then adjust sappi specifically, then follow standard procedure
data$correl_adj_ret_sappi <- data$ret_firstrand*(data$ewma_covar[nrow(data)]*sqrt(data$ewma_var_sappi-data$ewma_covar^2/data$ewma_var_firstrand)/sqrt(data$ewma_var_firstrand[nrow(data)])-data$ewma_covar*sqrt(data$ewma_var_sappi[nrow(data)]-data$ewma_covar[nrow(data)]^2/data$ewma_var_firstrand[nrow(data)])/sqrt(data$ewma_var_firstrand))/sqrt(data$ewma_var_firstrand*data$ewma_var_sappi-data$ewma_covar^2) + data$ret_sappi*sqrt(data$ewma_var_sappi[nrow(data)]-data$ewma_covar[nrow(data)]^2/data$ewma_var_firstrand[nrow(data)])/sqrt(data$ewma_var_sappi-data$ewma_covar^2/data$ewma_var_firstrand)                                             
data$stock_correl_adj_ret_port <- portfolio_weights[1]*data$vol_adj_ret_firstrand + portfolio_weights[2]*data$correl_adj_ret_sappi
adjusted_correl_stock_ret_var <- numeric(length = 4)
for (i in 1:4) {
  adjusted_correl_stock_ret_var[i] <- -quantile(data$stock_correl_adj_ret_port[-1], probs = alpha_levels[i])*sqrt(10)
}
adjusted_correl_stock_ret_var




# final VaR table
var_table <- data.frame(
  Alpha = alpha_levels,
  Unadjusted = unadjusted_var,
  `Method a` = adjusted_correl_stock_ret_var,
  `Method b` = adjusted_port_return_var,
  `Method c `= adjusted_stock_ret_var
)

var_table



# To Do: 3.7
# 1) Estimate and plot betas
# 2) Use most recent beta to decompose returns into systematic and specific
# 3) Calculate EWMA vol of systematic and specific returns
# 4) Use EWMA to adjust systematic and specific returns
# 5) Estimate total, systematic and specific VaR using returns and adjusted returns for T = 30 Oct 2006 and T = 21 April 2008
# 6) Estimate again using linear VaR

# Market EWMA
data$ret_alsi <- c(NA, diff(log(data$alsi)))
data$ewma_var_alsi <- NA
data$ewma_var_alsi[1] <- var(data$ret_alsi[-1])
data$ewma_covar_alsi_firstrand <- NA
data$ewma_covar_alsi_firstrand[1] <- cov(data$ret_alsi[-1], data$ret_firstrand[-1])
data$ewma_covar_alsi_sappi <- NA
data$ewma_covar_alsi_sappi[1] <- cov(data$ret_alsi[-1], data$ret_sappi[-1])

for (i in 2:nrow(data)) {
  data$ewma_var_alsi[i] <- lambda * data$ewma_var_alsi[i - 1] + (1 - lambda) * data$ret_alsi[i]^2
  data$ewma_covar_alsi_firstrand[i] <- lambda * data$ewma_covar_alsi_firstrand[i-1] + (1 - lambda) * data$ret_alsi[i] * data$ret_firstrand[i]
  data$ewma_covar_alsi_sappi[i] <- lambda * data$ewma_covar_alsi_sappi[i-1] + (1 - lambda) * data$ret_alsi[i] * data$ret_sappi[i]
}

# Calculate Beta's
data$beta_firstrand <- c(NA, data$ewma_covar_alsi_firstrand[-1]/data$ewma_var_alsi[-1])
data$beta_sappi <- c(NA, data$ewma_covar_alsi_sappi[-1]/data$ewma_var_alsi[-1])

# Plot Beta's
plot(data$Date, data$beta_sappi, type = "l", col = "red", lwd = 2, xlab = "Date", ylab = "Beta", main = "Firstrand and Sappi Beta Estimates")
lines(data$Date, data$beta_firstrand, type = "l", col = "blue")
legend("topright",  legend = c("Sappi", "Firstrand"), col = c("red", "blue"), lwd = 2, bty = "n")

firstrand_beta <- data$beta_firstrand[nrow(data)] 
sappi_beta <- data$beta_sappi[nrow(data)]

port_beta <- portfolio_weights[1]*firstrand_beta + portfolio_weights[2]*sappi_beta
data$ret_systematic <- port_beta*data$ret_alsi
data$ret_specific <- data$ret_port-data$ret_systematic
data$ewma_var_sys <- NA
data$ewma_var_sys[1] <- var(data$ret_systematic[-1])
data$ewma_var_spec <- NA
data$ewma_var_spec[1] <- var(data$ret_specific[-1])

# EWMA for systematic and specific returns
for (i in 2:nrow(data)) {
  data$ewma_var_sys[i] <- lambda * data$ewma_var_sys[i - 1] + (1 - lambda) * data$ret_systematic[i]^2
  data$ewma_var_spec[i] <- lambda * data$ewma_var_spec[i - 1] + (1 - lambda) * data$ret_specific[i]^2
}

data$ewma_vol_sys <- sqrt(data$ewma_var_sys*250)
data$ewma_vol_spec <- sqrt(data$ewma_var_spec*250)

data$vol_adj_ret_sys <- data$ret_systematic*data$ewma_vol_sys[nrow(data)]/data$ewma_vol_sys
data$vol_adj_ret_spec <- data$ret_specific*data$ewma_vol_spec[nrow(data)]/data$ewma_vol_spec

# Plot returns
plot(data$Date, data$ret_systematic, type = "l", col = "green", lwd = 2, xlab = "Date", ylab = "Return", main = "Adjusted vs Unadjusted Systematic Returns")
lines(data$Date, data$vol_adj_ret_sys, type = "l", col = "orange")
legend("topright",  legend = c("Unadjusted", "Adjusted"), col = c("green", "orange"), lwd = 2, bty = "n")

plot(data$Date, data$ret_specific, type = "l", col = "green", lwd = 2, xlab = "Date", ylab = "Return", main = "Adjusted vs Unadjusted Specific Returns")
lines(data$Date, data$vol_adj_ret_spec, type = "l", col = "orange")
legend("topright",  legend = c("Unadjusted", "Adjusted"), col = c("green", "orange"), lwd = 2, bty = "n")

# Calculate Value at Risk
unadjusted_systematic_var <- numeric(length = 4)
unadjusted_specific_var <- numeric(length = 4)
adjusted_systematic_var <- numeric(length = 4)
adjusted_specific_var <- numeric(length = 4)
for (i in 1:4){
  unadjusted_systematic_var[i] <- -quantile(data$ret_systematic[-1], probs = alpha_levels[i])*sqrt(10)
  unadjusted_specific_var[i] <- -quantile(data$ret_specific[-1], probs = alpha_levels[i])*sqrt(10)
  adjusted_systematic_var[i] <- -quantile(data$vol_adj_ret_sys[-1], probs = alpha_levels[i])*sqrt(10)
  adjusted_specific_var[i] <- -quantile(data$vol_adj_ret_spec[-1], probs = alpha_levels[i])*sqrt(10)
}
unadjusted_systematic_var
unadjusted_specific_var
adjusted_systematic_var
adjusted_specific_var

var_table <- data.frame(
  `Historical VaR` = c("Unadjusted", "adjusted"),
  Total = c(unadjusted_var[2], adjusted_port_return_var[2]),
  Systematic = c(unadjusted_systematic_var[2], adjusted_systematic_var[2]),
  Specific = c(unadjusted_specific_var[2], adjusted_specific_var[2])
)
t(var_table)

total_normal_var <- sqrt(10)*sd(data$ret_port[-1])*qnorm(0.99)
adjusted_total_normal_var <- sqrt(10)*sd(data$vol_adj_ret_port[-1])*qnorm(0.99)
systematic_normal_var <- sqrt(10)*sd(data$ret_systematic[-1])*qnorm(0.99)
adjusted_systematic_normal_var <- sqrt(10)*sd(data$vol_adj_ret_sys[-1])*qnorm(0.99)
specific_normal_var <- sqrt(10)*sd(data$ret_specific[-1])*qnorm(0.99)
adjusted_specific_normal_var <- sqrt(10)*sd(data$vol_adj_ret_spec[-1])*qnorm(0.99)

var_table <- data.frame(
  `Normal VaR` = c("Unadjusted", "adjusted"),
  Total = c(total_normal_var, adjusted_total_normal_var),
  Systematic = c(systematic_normal_var, adjusted_systematic_normal_var),
  Specific = c(specific_normal_var, adjusted_specific_normal_var)
)
t(var_table)



