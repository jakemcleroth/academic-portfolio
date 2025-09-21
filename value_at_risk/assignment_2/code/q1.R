library(readxl)
library(dplyr)
library(purrr)

# Data
file_path <- "/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_2/data/data_assignment_2.xlsx"
data <- read_xlsx(file_path, sheet = "ALLSHARE")
data$returns <- c(NA, diff(log(data$price)))
plot(data$date, data$price, main = "JSE ALSI", xlab = "Date", ylab = "Price", type = "l", col = "red")
# data <- data[data$date >= as.POSIXct("01/01/2005", format = "%d/%m/%Y", tz = "UTC"), ] # Uncomment for second time period
# data <- data[data$date >= as.POSIXct("01/01/2015", format = "%d/%m/%Y", tz = "UTC"), ] # Uncomment for third time period


log_returns <- na.omit(data$returns)
h_values <- c(1,2,4,6,8,10,12,14,16,18,20)
h_day_returns <- data.frame(h = h_values, h_day_returns = I(vector("list", length(h_values))))  

for (i in 1:length(h_values)){
  h_day_returns$h_day_returns[[i]] <- tibble(index = seq_along(log_returns), log_returns) %>%
    mutate(group = (index - 1) %/% h_values[i]) %>%  
    group_by(group) %>%
    summarise(sum_returns = sum(log_returns), .groups = "drop") %>%
    pull(sum_returns)
}

str(h_day_returns)
alpha_values <- c(0.001, 0.01, 0.05, 0.10)  

# Compute all quantiles and store as separate columns
h_day_returns <- h_day_returns %>%
  mutate(
    quantile_0.1 = sapply(h_day_returns, function(x) quantile(x, probs = 0.001, na.rm = TRUE)),
    quantile_1   = sapply(h_day_returns, function(x) quantile(x, probs = 0.01, na.rm = TRUE)),
    quantile_5   = sapply(h_day_returns, function(x) quantile(x, probs = 0.05, na.rm = TRUE)),
    quantile_10  = sapply(h_day_returns, function(x) quantile(x, probs = 0.10, na.rm = TRUE))
  )
str(h_day_returns)

# Extract 1-day quantiles
quantile_0.1_h1 <- h_day_returns$quantile_0.1[h_day_returns$h == 1]
quantile_1_h1   <- h_day_returns$quantile_1[h_day_returns$h == 1]
quantile_5_h1   <- h_day_returns$quantile_5[h_day_returns$h == 1]
quantile_10_h1  <- h_day_returns$quantile_10[h_day_returns$h == 1]

# Compute log differences
h_day_returns <- h_day_returns %>%
  mutate(
    log_q_0.1 = log(quantile_0.1/quantile_0.1_h1),
    log_q_1   = log(quantile_1/quantile_1_h1),
    log_q_5   = log(quantile_5/quantile_5_h1),
    log_q_10  = log(quantile_10/quantile_10_h1)
  )

print(h_day_returns)

scale_exp_0.1 <- var(h_day_returns$log_q_0.1)/cov(h_day_returns$log_q_0.1,log(h_day_returns$h))
scale_exp_1 <- var(h_day_returns$log_q_1)/cov(h_day_returns$log_q_1,log(h_day_returns$h))
scale_exp_5 <- var(h_day_returns$log_q_5)/cov(h_day_returns$log_q_5,log(h_day_returns$h))
scale_exp_10 <- var(h_day_returns$log_q_10)/cov(h_day_returns$log_q_10,log(h_day_returns$h))

scale_exp_0.1
scale_exp_1
scale_exp_5
scale_exp_10

# Change title for each time period
plot(h_day_returns$log_q_0.1, log(h_day_returns$h), main = "ALSI log-log plot from 1995 onwards vs 0.1% quantile ratio", xlab = "log(Quantile Ratio)", ylab = "log(h)")
plot(h_day_returns$log_q_1, log(h_day_returns$h), main = "ALSI log-log plot from 1995 onwards vs 1% quantile ratio", xlab = "log(Quantile Ratio)", ylab = "log(h)")
plot(h_day_returns$log_q_5, log(h_day_returns$h), main = "ALSI log-log plot from 1995 onwards vs 5% quantile ratio", xlab = "log(Quantile Ratio)", ylab = "log(h)")
plot(h_day_returns$log_q_10, log(h_day_returns$h), main = "ALSI log-log plot from 1995 onwards vs 10% quantile ratio", xlab = "log(Quantile Ratio)", ylab = "log(h)")

