library(readxl)


### Data ###
alsi <- as.data.frame(read_excel("/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_2/data/data_assignment_2.xlsx", sheet = "ALLSHARE"))
alsi$returns <- c(NA, diff(log(alsi$price)))
plot(x = alsi$date, y = alsi$price, type = "l", main = "JSE All Share Index", ylab = "Prices", xlab = "Date")
plot(x = alsi$date[-1],y = alsi$returns[-1], type = "l", , main = "JSE All Share Index Returns", ylab = "Log Returns", xlab = "Date")


# Calculate VaR
alsi$hist500_VaR <- NA
alsi$norm500_VaR <- NA
alsi$hist2000_VaR <- NA
alsi$norm2000_VaR <- NA
alsi$dif500 <- NA
alsi$dif2000 <- NA

inv_norm <- qnorm(0.99, 0, 1)

for (i in 501:nrow(alsi)){
  temp_returns <- alsi$returns[(i-500+1):i]
  alsi$hist500_VaR[i] <- -quantile(temp_returns, probs = 0.01)
  alsi$norm500_VaR[i] <- sd(temp_returns)*inv_norm - mean(temp_returns)
  alsi$dif500[i] <- alsi$hist500_VaR[i] - alsi$norm500_VaR[i]
}

for (i in 2001:nrow(alsi)){
  temp_returns <- alsi$returns[(i-2000+1):i]
  alsi$hist2000_VaR[i] <- -quantile(temp_returns, probs = 0.01)
  alsi$norm2000_VaR[i] <- sd(temp_returns)*inv_norm - mean(temp_returns)
  alsi$dif2000[i] <- alsi$hist2000_VaR[i] - alsi$norm2000_VaR[i]
}

# Find appropriate y-axis limits for both VaR series
y_limits <- range(c(alsi$hist500_VaR, alsi$hist2000_VaR), na.rm = TRUE)

# Plot Hist VaR
plot(x = alsi$date, y = alsi$hist500_VaR, type = "l", main = "JSE All Share Historical VaR", ylab = "VaR", xlab = "Date", col = "blue", ylim = y_limits)  
lines(x = alsi$date, y = alsi$hist2000_VaR, col = "red")
legend("topright", legend = c("500-day VaR", "2000-day VaR"), col = c("blue", "red"), lty = 1)

# Find appropriate y-axis limits for both VaR series
y_limits <- range(c(alsi$norm500_VaR, alsi$norm2000_VaR), na.rm = TRUE)

# Plot normal VaR
plot(x = alsi$date, y = alsi$norm500_VaR, type = "l", main = "JSE All Share Normal Linear VaR", ylab = "VaR", xlab = "Date", col = "blue", ylim = y_limits)
lines(x = alsi$date, y = alsi$norm2000_VaR, col = "red")
legend("topright", legend = c("500-day VaR", "2000-day VaR"), col = c("blue", "red"), lty = 1)


y_limits <- range(c(alsi$dif500, alsi$dif2000), na.rm = TRUE)

# Plot difference
plot(x = alsi$date, y = alsi$dif500, type = "l", main = "JSE All Share VaR Differences", ylab = "VaR", xlab = "Date", col = "blue", ylim = y_limits)  
lines(x = alsi$date, y = alsi$dif2000, col = "red")
legend("topright", legend = c("500-day Difference", "2000-day Difference"), col = c("blue", "red"), lty = 1)













