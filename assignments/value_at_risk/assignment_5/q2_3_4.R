###### Libraries ######
library(readxl)
library(dplyr)

###### Data ######
data <- read_excel("data/q2_3_4.xlsx")
data$Vol <- data$Vol*100

###### Question 2 ######
### Params ###
n <- 2000
point <- 100
alpha <- 0.01

### Returns ###
data <- data %>% mutate(
  index_returns = log(Index) - lag(log(Index), 1)
)

### Backtest ###
data$VaR <- NA
data$PnL <- NA
data$Indicator_1 <- NA
data$Indicator_2 <- NA
data$Indicator_3 <- NA
for (i in 252:2251) {
  data$VaR[i] <- qnorm(alpha)*sd(data$index_returns[(i-250):(i-1)])*data$Index[i-1]*point
  data$PnL[i] <- (data$Index[i] - data$Index[i-1])*point
  data$Indicator_1[i] <- ifelse(data$VaR[i]>data$PnL[i],1,0)
  
}

for (i in 253:2251) {
  data$Indicator_2[i] <- ifelse(data$Indicator_1[i-1]==data$Indicator_1[i], ifelse(data$Indicator_1[i]==1,1,0),0)
  data$Indicator_3[i] <- ifelse(data$Indicator_1[i-1]==data$Indicator_1[i], ifelse(data$Indicator_1[i]==0,1,0),0)
}

### Independence Test ###
expected_exceptions = alpha*n
n_1 <- sum(data$Indicator_1, na.rm = T)
n_0 <- n-n_1

n_00 <- sum(data$Indicator_3, na.rm = T)
n_11 <- sum(data$Indicator_2, na.rm = T)
n_01 <- n_1-n_11
n_10 <- n_0-n_00
p_exp <- alpha
p_obs <- n_1/n
p_01 <- n_01/(n_00+n_01)
p_11 <- n_11/(n_10+n_11)

numerator <- (p_obs^n_1) * ((1 - p_obs)^n_0)
denominator <- (p_01^n_01) * ((1 - p_01)^n_00) * (p_11^n_11) * ((1 - p_11)^n_10)
ln_lr <- ifelse(n_11>0, log(numerator / denominator), 0)
test_stat <- -2*ln_lr

if (qchisq(0.1,1,lower.tail = F)<test_stat) {
  print("Reject at 0.1")
} else {
  print("Fail to Reject at 0.1")
}

if (qchisq(0.05,1,lower.tail = F)<test_stat) {
  print("Reject at 0.05")
} else {
  print("Fail to Reject at 0.05")
}

if (qchisq(0.01,1,lower.tail = F)<test_stat) {
  print("Reject at 0.01")
} else {
  print("Fail to Reject at 0.01")
}


###### Question 3 ######

### Vol difference ###
data <- data %>% mutate(
  delta_vol = Vol - lag(Vol, 1)
  
)

### Lagged vol difference ###
data <- data %>% mutate(
  lagged_delta_vol = lag(delta_vol,1)
)

### Fit regression ###
model_full <- lm(Indicator_1 ~ lagged_delta_vol, data = data, na.action = na.omit)
summary(model_full)
data$null_model <- data$Indicator_1-alpha

# RSS under H0
RSS_restricted <- sum(data$null_model^2, na.rm = T)

# RSS under full model
RSS_full <- sum(resid(model_full)^2)

# degrees of freedom
df1 <- 2  
df2 <- df.residual(model_full)

# F-statistic
F_stat <- ((RSS_restricted - RSS_full) / df1) / (RSS_full / df2)

# p-value
p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)

F_stat
p_value

if (p_value<0.05) {
  print("Reject")
} else {
  print("Fail to Reject")
}


###### Question 4 ######
rm(list = ls())


data <- read_excel("data/q2_3_4.xlsx")
data$Vol <- data$Vol*100
data <- data %>% mutate(
  index_returns = log(Index) - lag(log(Index), 1)
)

### params ###
lambda <- 0.94
len <- nrow(data)
alpha <- 0.01
point <- 100
n <- 2000
test_prob <- 0.1
chi_sqr_crit <- qchisq(test_prob,1,lower.tail = F)

### EWMA Variance ###
data$ewma_var_index <- NA
data$ewma_var_index[1] <- var(data$index_returns[-1])
data$ewma_var_index[3] <- lambda * data$ewma_var_index[1] + (1 - lambda) * data$index_returns[2]^2

for (i in 4:len) {
  data$ewma_var_index[i] <- lambda * data$ewma_var_index[i - 1] + (1 - lambda) * data$index_returns[i-1]^2
}

### Backtest ###
data$VaR <- NA
data$PnL <- NA
data$Indicator_1 <- NA
data$Indicator_2 <- NA
data$Indicator_3 <- NA
for (i in 252:2251) {
  data$VaR[i] <- qnorm(alpha)*sqrt(data$ewma_var_index[i])*data$Index[i-1]*point
  data$PnL[i] <- (data$Index[i] - data$Index[i-1])*point
  data$Indicator_1[i] <- ifelse(data$VaR[i]>data$PnL[i],1,0)
  
}

for (i in 253:2251) {
  data$Indicator_2[i] <- ifelse(data$Indicator_1[i-1]==data$Indicator_1[i], ifelse(data$Indicator_1[i]==1,1,0),0)
  data$Indicator_3[i] <- ifelse(data$Indicator_1[i-1]==data$Indicator_1[i], ifelse(data$Indicator_1[i]==0,1,0),0)
}

### Tests ###
expected_exceptions = alpha*n
n_1 <- sum(data$Indicator_1, na.rm = T)
n_0 <- n-n_1

n_00 <- sum(data$Indicator_3, na.rm = T)
n_11 <- sum(data$Indicator_2, na.rm = T)
n_01 <- n_1-n_11
n_10 <- n_0-n_00
p_exp <- alpha
p_obs <- n_1/n
p_01 <- n_01/(n_00+n_01)
p_11 <- n_11/(n_10+n_11)

# Unconditional Coverage
numerator_uc <- (p_exp^n_1) * ((1 - p_exp)^n_0)
denominator_uc <- (p_obs^n_1) * ((1 - p_obs)^n_0)
ln_lr_uc <- log(numerator_uc / denominator_uc)
test_stat_unc <- -2 * ln_lr_uc
if (chi_sqr_crit<test_stat_unc) {
  print("Reject")
} else {
  print("Fail to Reject")
}

# Independence
numerator <- (p_obs^n_1) * ((1 - p_obs)^n_0)
denominator <- (p_01^n_01) * ((1 - p_01)^n_00) * (p_11^n_11) * ((1 - p_11)^n_10)
ln_lr <- ifelse(n_11>0, log(numerator / denominator), 0)
test_stat_ind <- -2*ln_lr
if (chi_sqr_crit<test_stat_ind) {
  print("Reject")
} else {
  print("Fail to Reject")
}

# Conditional Coverage
test_stat_con <- -2*(ln_lr+ln_lr_uc)
if (qchisq(test_prob,2,lower.tail = F)<test_stat_con) {
  print("Reject")
} else {
  print("Fail to Reject")
}

