library("readxl")
library("ggplot2")

### GET DATA
file_path <- "/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_1/raw_data/q5_data_formated.xlsx"
all_share <- read_xlsx(file_path, sheet = "ALLSHARE")
top_40 <- read_xlsx(file_path, sheet = "TOP40")
indi_25 <- read_xlsx(file_path, sheet = "INDUSTRIAL25")

# ALL SHARE RETURNS
all_share$returns <- c(NA, diff(log(all_share$`JSE ALL SHARE`)))
all_share_returns <- na.omit(all_share$returns)

# TOP 40 RETURNS
top_40$returns <- c(NA, diff(log(top_40$`JSE TOP 40`)))
top_40_returns <- na.omit(top_40$returns)

# INDUSTRIAL 25 RETURNS
indi_25$returns <- c(NA, diff(log(indi_25$`JSE INDUSTRIAL 25`)))
indi_25_returns <- na.omit(indi_25$returns)

### a)
### PLOT RETURNS
# ALL SHARE
ggplot(all_share, aes(x = Date, y = returns)) + 
  geom_line(color = "#1f77b4", size = 1) +
  labs(title = "JSE ALL SHARE RETURNS", x = "Date", y = "Log Returns") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_line(color = "gray", size = 0.2)) +
  theme(panel.grid.minor = element_line(color = "gray", size = 0.1))

# ALL SHARE HISTOGRAM
ggplot(all_share, aes(x = returns)) +
  geom_histogram(binwidth = 0.01, fill = "#1f77b4", color = "black", alpha = 0.7) +
  labs(title = "DISTRIBUTION OF ALL SHARE RETURNS", x = "Log Returns", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# TOP 40
ggplot(top_40, aes(x = Date, y = returns)) + 
  geom_line(color = "#ff7f0e", size = 1) +
  labs(title = "JSE TOP 40 RETURNS", x = "Date", y = "Log Returns") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_line(color = "gray", size = 0.2)) +
  theme(panel.grid.minor = element_line(color = "gray", size = 0.1))

# TOP 40 HISTOGRAM
ggplot(top_40, aes(x = returns)) +
  geom_histogram(binwidth = 0.01, fill = "#ff7f0e", color = "black", alpha = 0.7) +
  labs(title = "DISTRIBUTION OF TOP 40 RETURNS", x = "Log Returns", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# INDUSTRIAL 25
ggplot(indi_25, aes(x = Date, y = returns)) + 
  geom_line(color = "#2ca02c", size = 1) +
  labs(title = "JSE INDUSTRIAL 25 RETURNS", x = "Date", y = "Log Returns") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_line(color = "gray", size = 0.2)) +
  theme(panel.grid.minor = element_line(color = "gray", size = 0.1))

# INDUSTRIAL 25 HISTOGRAM
ggplot(indi_25, aes(x = returns)) +
  geom_histogram(binwidth = 0.01, fill = "#2ca02c", color = "black", alpha = 0.7) +
  labs(title = "DISTRIBUTION OF INDUSTRIAL 25 RETURNS", x = "Log Returns", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

### QQ PLOTS

# ALL SHARE
ggplot(data.frame(returns = all_share_returns), aes(sample = returns)) + 
  stat_qq() + 
  stat_qq_line() + 
  labs(title = "QQ Plot for JSE ALL SHARE RETURNS") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# TOP 40
ggplot(data.frame(returns = top_40_returns), aes(sample = returns)) + 
  stat_qq() + 
  stat_qq_line() + 
  labs(title = "QQ Plot for JSE TOP 40 RETURNS") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# INDUSTRIAL 25
ggplot(data.frame(returns = indi_25_returns), aes(sample = returns)) + 
  stat_qq() + 
  stat_qq_line() + 
  labs(title = "QQ Plot for JSE INDUSTRIAL 25 RETURNS") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


### TESTING FOR NORMALITY
# Kolmogorov-Smirnov Test
ks.test(all_share_returns, "pnorm", mean = mean(all_share_returns), sd = sd(all_share_returns))
ks.test(top_40_returns, "pnorm", mean = mean(top_40_returns), sd = sd(top_40_returns))
ks.test(indi_25_returns, "pnorm", mean = mean(indi_25_returns), sd = sd(indi_25_returns))

# Shapiro-Wilk Test
shapiro.test(all_share_returns)
shapiro.test(top_40_returns)
shapiro.test(indi_25_returns)

# Anderson-Darling Test
# install.packages("nortest")
library("nortest")
ad.test(all_share_returns)
ad.test(top_40_returns)
ad.test(indi_25_returns)

### b)
### TESTING INDEPENDENCE
# ACF PLOTS
acf(all_share_returns, lag.max =10, main = "Autocorrelation of ALL SHARE Returns")
acf(top_40_returns, lag.max = 10, main = "Autocorrelation of TOP 40 Returns")
acf(indi_25_returns, lag.max = 10, main = "Autocorrelation of INDUSTRIAL 25 Returns")

# Ljung-Box Test
Box.test(all_share_returns, lag = 4, type = "Ljung-Box")
Box.test(top_40_returns, lag = 4, type = "Ljung-Box")
Box.test(indi_25_returns, lag = 4, type = "Ljung-Box")





