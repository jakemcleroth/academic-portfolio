Variance Reduction Techniques for Pricing Asian Options
================
Nicole Lyne & Jake Mc Leroth
3 May 2024

- [1 Introduction](#1-introduction)
- [2 Financial Background](#2-financial-background)
  - [2.1 Overview of Options](#21-overview-of-options)
  - [2.2 Geometric Brownian Motion](#22-geometric-brownian-motion)
  - [2.3 Options Pricing](#23-options-pricing)
  - [2.4 Asian Options](#24-asian-options)
    - [2.4.1 Background](#241-background)
- [3 Monte Carlo Simulation for Pricing Asian
  Options](#3-monte-carlo-simulation-for-pricing-asian-options)
  - [3.1 Monte Carlo Integration](#31-monte-carlo-integration)
  - [3.2 Pricing Asian Options](#32-pricing-asian-options)
    - [3.2.1 Direct Method](#321-direct-method)
    - [3.2.2 Implementation of Method to Price Asian Call
      Option](#322-implementation-of-method-to-price-asian-call-option)
    - [3.2.3 Variance of Asian Call Option
      Price](#323-variance-of-asian-call-option-price)
- [4 Variance Reduction](#4-variance-reduction)
  - [4.1 Antithetic Variables Method](#41-antithetic-variables-method)
    - [4.1.1 Description](#411-description)
    - [4.1.2 Implementation of Method to Price Asian Call
      Option](#412-implementation-of-method-to-price-asian-call-option)
    - [4.1.3 Variance of Asian Call Option
      Price](#413-variance-of-asian-call-option-price)
  - [4.2 n-Simplex Antithetic Variate
    Method](#42-n-simplex-antithetic-variate-method)
    - [4.2.1 The Simplex](#421-the-simplex)
    - [4.2.2 n-Simplex Method](#422-n-simplex-method)
    - [4.2.3 Implementation of Method to Price Asian Call
      Option](#423-implementation-of-method-to-price-asian-call-option)
- [5 Analyis](#5-analyis)
  - [5.1 Results](#51-results)
  - [5.2 Analysis of results](#52-analysis-of-results)
- [6 Conclusion](#6-conclusion)
- [7 References](#7-references)

<style>
  table {
    border-collapse: collapse;
    width: 100%;
    border-top: 1px solid black; /* Added border-top */
  }
  th, td {
    border: 1px solid black;
    padding: 8px;
    text-align: center; 
  }
  th {
    border-top: 1px solid black;
    border-bottom: 1px solid black;
    background-color: lightgray; /* Added background color for headers */
  }
</style>

# 1 Introduction

The project aims to price Asian call options using Monte Carlo
simulation, comparing direct, antithetic variables, and n-simplex
methods. Asian options, offering unique advantages in hedging and low
liquidity scenarios, are path-dependent options whose payoff relies on
the average of the underlying asset. Monte Carlo simulation, a method
for estimating complex quantities, will be used to simulate future stock
prices under geometric Brownian motion. By implementing variance
reduction techniques like antithetic variables and n-simplex methods,
the project aims to improve pricing accuracy and efficiency. Overall, it
explores diverse methods for pricing Asian options which is particularly
of interest in the field of Financial Risk Management.

# 2 Financial Background

## 2.1 Overview of Options

An option is a contract between a buyer - the holder of the contract -
and a seller, which gives the holder the right, but not the obligation,
to buy or sell the underlying asset by a certain date at a certain
price, known as the strike price. A call option is the right to buy the
underlying asset, while a put option is the right to sell the underlying
asset. An in-the-money call option is where the strike price is less
than the current asset price. An out-the-money call option is when the
strike price is more than the current asset price. The expiration date
or maturity date is the date agreed upon in the contract and this is
when the contract ends (Hull 2017). The most basic type of option is a
European call option on stocks, where the holder has the right to buy
the underlying stock only on the expiration date, and not any time
before or after. There are many types of options, each with their own
unique characteristics and payoff structure, catering to different risk
preferences and market conditions. The most common options are European
and American options(Hull 2017). Asian options, another type of option,
are the options priced in this project. More specifically, the prices of
arithmetic average Asian call options are estimated.

The payoff is given by (Zhang 2009):

![\phi(S) = (S(T) - K)^+](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi%28S%29%20%3D%20%28S%28T%29%20-%20K%29%5E%2B "\phi(S) = (S(T) - K)^+")

for a European call option

and by

![\phi(S) = (K- S(T))^+](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi%28S%29%20%3D%20%28K-%20S%28T%29%29%5E%2B "\phi(S) = (K- S(T))^+")

for a European put option, where -
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
is the expiration date, -
![S(T)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%28T%29 "S(T)")
is the spot price of the underlying asset at expiration T, -
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
is the strike price of the option, -
![(x)^+](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x%29%5E%2B "(x)^+")
represents the maximum of
![x \\\text{and} \\0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%5C%20%5Ctext%7Band%7D%20%5C%200 "x \ \text{and} \ 0").

## 2.2 Geometric Brownian Motion

To model the behaviour of stock prices, the most commonly used model is
geometric Brownian motion. It can be formulated as the following
stochastic differential equation (SDE) as seen in Hull (2017):

![dS_t=\mu S_t dt + \sigma S_t d W_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;dS_t%3D%5Cmu%20S_t%20dt%20%2B%20%5Csigma%20S_t%20d%20W_t "dS_t=\mu S_t dt + \sigma S_t d W_t")

where: -
![dS](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;dS "dS")
is the instantaneous change in the stock price, -
![dW](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;dW "dW")
is the Wiener process
(![\triangle W_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctriangle%20W_t "\triangle W_t")
is a random variable from the normal distribution
![\text{N}(0,\triangle t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BN%7D%280%2C%5Ctriangle%20t%29 "\text{N}(0,\triangle t)")), -
\$\$ represents the volatility of the stock price
(![\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma "\sigma")
is assumed to be constant over time), -
![\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
represents the growth rate of the stock price when there is no risk
(also referred to as the ‘drift’), -
![\mu dt](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20dt "\mu dt")
represents the deterministic return within the time interval.

with the discrete time version of the SDE as:

![\triangle S_t = \mu S_t \triangle t + \sigma S_t \triangle W_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctriangle%20S_t%20%3D%20%5Cmu%20S_t%20%5Ctriangle%20t%20%2B%20%5Csigma%20S_t%20%5Ctriangle%20W_t "\triangle S_t = \mu S_t \triangle t + \sigma S_t \triangle W_t")

The solution to the discrete SDE is as follows: where
![Z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z "Z")
is a standard normal random variable.

## 2.3 Options Pricing

The main idea behind options pricing, or derivative pricing in general,
is known as risk neutral valuation. The idea is that the price of an
option should be equal to the present value of the expected payoff.
According to Hull (2017), this can be done in the following steps: 1.
Assume expected rate of return is the risk free rate. 2. Calculate the
expected payoff of the option. For example, a European call option:

![E\[\text{payoff}\] = E\[(S(T) - K)^+\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5B%5Ctext%7Bpayoff%7D%5D%20%3D%20E%5B%28S%28T%29%20-%20K%29%5E%2B%5D "E[\text{payoff}] = E[(S(T) - K)^+]")

3.  Discount the expected payoff at the risk free rate. For the European
    call option:

![C = e^{-rt}E\[\text{payoff}\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%20%3D%20e%5E%7B-rt%7DE%5B%5Ctext%7Bpayoff%7D%5D "C = e^{-rt}E[\text{payoff}]")

This European call option can be priced using the Black-Scholes-Merton
framework with a closed form solution. However, there is no closed form
solution for pricing arithmetic Asian options. Hence we need to use
simulation.

## 2.4 Asian Options

### 2.4.1 Background

An Asian option is similar to a European option since the holder can
only exercise on the date of maturity. However, the payoff depends on
the average of the stock price instead of just the price on the date of
maturity. There are many ways to calculate the average, leading to a
variety of options within the Asian class. In this project, the discrete
arithmetic average is used. Asian options were introduced in the late
1970s in the oil market and now they are traded in the
over-the-counter(OTC) market. They are path-dependent because the payoff
depends on the average of the asset price, hence the full path is
required to evaluate the payoff. Since the sum of lognormal random
variables is not log normal, the average does not have a known
distribution. As a result, there is no closed form pricing formula
(Nielsen 2001).

Discrete arithmetic average Asian options have the following payoff
structure:

![\Phi(S) = \left(\bar{S}  - K \right)^{+} \quad \text{for a call}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%28S%29%20%3D%20%5Cleft%28%5Cbar%7BS%7D%20%20-%20K%20%5Cright%29%5E%7B%2B%7D%20%5Cquad%20%5Ctext%7Bfor%20a%20call%7D "\Phi(S) = \left(\bar{S}  - K \right)^{+} \quad \text{for a call}")

or

![\Phi(S) = \left( K - \bar{S}  \right)^{+} \quad \text{for a put}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%28S%29%20%3D%20%5Cleft%28%20K%20-%20%5Cbar%7BS%7D%20%20%5Cright%29%5E%7B%2B%7D%20%5Cquad%20%5Ctext%7Bfor%20a%20put%7D "\Phi(S) = \left( K - \bar{S}  \right)^{+} \quad \text{for a put}")

where

![\text{where} \\\bar{S} = \frac{1}{m} \sum\_{i=1}^{m} S(t_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7Bwhere%7D%20%5C%20%5Cbar%7BS%7D%20%3D%20%5Cfrac%7B1%7D%7Bm%7D%20%5Csum_%7Bi%3D1%7D%5E%7Bm%7D%20S%28t_i%29 "\text{where} \ \bar{S} = \frac{1}{m} \sum_{i=1}^{m} S(t_i)")

and

For the remainder of this project, the discrete arithmetic average Asian
option will be referred to as an Asian option.

# 3 Monte Carlo Simulation for Pricing Asian Options

Monte Carlo methods are statistical methods that use random sampling to
solve problems or estimate quantities that may be difficult or
impossible to solve analytically. After World War II (1940s), Monte
Carlo methods were developed, although the idea of random sampling came
about as early as 1777 (Rizzo 2019).

## 3.1 Monte Carlo Integration

Consider a function g(x). Suppose that we would like to determine
![\int\_{a}^{b} g(x) \\ dx](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cint_%7Ba%7D%5E%7Bb%7D%20g%28x%29%20%5C%2C%20dx "\int_{a}^{b} g(x) \, dx"),
provided the integral exists. If
![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X")
is a random variable with probability density function
![f(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29 "f(x)"),
then the expected value of the random variable
![Y = g(X)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20g%28X%29 "Y = g(X)")
is given by:

![E\[g(X)\] = \int\_{-\infty}^{\infty} g(x)f(x) \\ dx](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5Bg%28X%29%5D%20%3D%20%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D%20g%28x%29f%28x%29%20%5C%2C%20dx "E[g(X)] = \int_{-\infty}^{\infty} g(x)f(x) \, dx")

When a random sample is available from the distribution of
![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X"),
an unbiased estimator of
![E\[g(X)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5Bg%28X%29%5D "E[g(X)]")
is simply the sample mean (Rizzo 2019).

Suppose
![f(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29 "f(x)")
is a probability density function, such that
![f(x) \geq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%5Cgeq%200 "f(x) \geq 0")
for all
![x \in \mathbb{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%5Cin%20%5Cmathbb%7BR%7D "x \in \mathbb{R}")
and
![\int\_{A} f(x) \\ dx = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cint_%7BA%7D%20f%28x%29%20%5C%2C%20dx%20%3D%201 "\int_{A} f(x) \, dx = 1").
In order to estimate the integral
![\theta = \int\_{A} g(x)f(x) \\ dx](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20%5Cint_%7BA%7D%20g%28x%29f%28x%29%20%5C%2C%20dx "\theta = \int_{A} g(x)f(x) \, dx"),
generate a random sample
![x_1, x_2, ..., x_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%2C%20x_2%2C%20...%2C%20x_n "x_1, x_2, ..., x_n")
from the distribution of
![f(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29 "f(x)")
and calculate the sample mean

![\hat{\theta} = \frac{1}{n} \sum\_{j=1}^{n} g(x_j)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20g%28x_j%29 "\hat{\theta} = \frac{1}{n} \sum_{j=1}^{n} g(x_j)")

The core principle of Monte Carlo integration relies on the strong law
of large numbers:

![\hat{\theta} = \frac{1}{n} \sum\_{j=1}^{n} g(x_j) \rightarrow  E\[g(X)\]\quad a.s.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20g%28x_j%29%20%5Crightarrow%20%20E%5Bg%28X%29%5D%5Cquad%20a.s. "\hat{\theta} = \frac{1}{n} \sum_{j=1}^{n} g(x_j) \rightarrow  E[g(X)]\quad a.s.")

![\hat{\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D "\hat{\theta}")
as defined above is the Monte Carlo estimator of
![E\[g(X)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5Bg%28X%29%5D "E[g(X)]")
(Rizzo 2019).

## 3.2 Pricing Asian Options

### 3.2.1 Direct Method

We estimate
![\theta = E\[g(X)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5Bg%28X%29%5D "\theta = E[g(X)]")
by
![\hat{\theta}\_{DIR} = \frac{1}{n} \sum\_{j=1}^{n} g(X_j)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D_%7BDIR%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20g%28X_j%29 "\hat{\theta}_{DIR} = \frac{1}{n} \sum_{j=1}^{n} g(X_j)")
where
![X_1, X_2, ..., X_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_1%2C%20X_2%2C%20...%2C%20X_n "X_1, X_2, ..., X_n")
are independent and identically distributed. The direct method refers to
the basic Monte Carlo simulation without any variance reduction
techniques. Random samples are generated directly from the distribution
of interest.

In the case of pricing Asian call options, the direct method involves
generating random samples of the underlying asset’s price at various
future time points based on a stochastic model, such as geometric
Brownian motion. These simulated price paths are then used to calculate
the payoff of the option at the expiration date. The option’s price is
estimated by discounting back the average of the payoffs. Recall the
payoff of a Asian call option is (Zhang 2009):

![\Phi(S) = \left( \bar{S} - K \right)^{+}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%28S%29%20%3D%20%5Cleft%28%20%5Cbar%7BS%7D%20-%20K%20%5Cright%29%5E%7B%2B%7D "\Phi(S) = \left( \bar{S} - K \right)^{+}")

By risk-neutral valuation, the price of the Asian call option can be
expressed as:

Monte Carlo simulation can be applied to price Asian options.
Substituting
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
from above with the payoff of the Asian call option, the following is
obtained:

![\frac{1}{n} \sum\_{j=1}^{n} \left( \bar{s} -K\right)^{+}\_j \rightarrow E\left\[\left( \bar{S} -K\right)^{+}\right\] \quad a.s.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j%20%5Crightarrow%20E%5Cleft%5B%5Cleft%28%20%5Cbar%7BS%7D%20-K%5Cright%29%5E%7B%2B%7D%5Cright%5D%20%5Cquad%20a.s. "\frac{1}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j \rightarrow E\left[\left( \bar{S} -K\right)^{+}\right] \quad a.s.")

where: -
![\left( \bar{s} -K\right)^{+}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j "\left( \bar{s} -K\right)^{+}_j")
is the payoff of the Asian call option from a simulated stock path, - n
is the number of Monte Carlo simulations. Hence, the estimated price of
the call option can be expressed as follows:

![\hat{C} = \frac{e^{-rt}}{n} \sum\_{j=1}^{n} \left( \bar{s} -K\right)^{+}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D%20%3D%20%5Cfrac%7Be%5E%7B-rt%7D%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j "\hat{C} = \frac{e^{-rt}}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j")

### 3.2.2 Implementation of Method to Price Asian Call Option

To simulate stock paths, we make use of the geometric Brownian motion
explained earlier. The simulation goes as follows: 1. Generate m random
variable observations,
![z_1, z_2, ...,z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_1%2C%20z_2%2C%20...%2Cz_m "z_1, z_2, ...,z_m"),
from the standard normal distribution. 2. Simulate the path of the stock
at m time points using:

![s_t = s\_{t-1}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z_t \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_t%20%3D%20s_%7Bt-1%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20z_t%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_t = s_{t-1}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z_t \sqrt{\triangle t}\right)}")

3.  Repeat n times to generate n paths.

10 simulated paths are shown below:
![](stoch_sim_project_files/figure-gfm/Asian%20Call%20Option%2010-1.png)<!-- -->

In the simulation above, it can be seen that the stock price path begins
at 100 for each simulation since this is the stock price at time 0. The
number of time steps chosen is 250 since there are roughly 250 trading
days in 1 year. Each “day’s” price (i.e. time step 1, 2, 3 etc.) is
determined by the price the day before and is modelled according to
geometric Brownian motion. This creates one path. This process is
repeated n times to create a wide variety of different outcomes. In this
case, ten stock price paths are simulated, so the simulation has
produced ten potential trajectories of the stock price, based on random
fluctuations.

The Asian call option is then priced with the following steps: 1.
Simulate a stock path with geometric Brownian motion. 2. Sum up the
prices on each of the m days and divide by m to obtain the average price
of the asset over the entire term \$ \_{i=1}^{m} s(t_i) = {s}![.
3. Subtract K from the average stock price over the entire term (](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;.%0A3.%20Subtract%20K%20from%20the%20average%20stock%20price%20over%20the%20entire%20term%20%28 ".
3. Subtract K from the average stock price over the entire term (")-K![). Then the payoff of the Asian call option is max](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%29.%20Then%20the%20payoff%20of%20the%20Asian%20call%20option%20is%20max "). Then the payoff of the Asian call option is max")({s}-K,0)\$.
4. Repeat steps 1-3 n times to obtain n different payoffs. 5. Take the
average of all n payoffs and discount it with
![e^{-rt}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;e%5E%7B-rt%7D "e^{-rt}")
where t is the expiration date.

The more paths that are generated, the more accuracte the price will be.
Shown below is an example of simulating 10000 stock paths.
![](stoch_sim_project_files/figure-gfm/Asian%20Call%20Option%2010000-1.png)<!-- -->

### 3.2.3 Variance of Asian Call Option Price

![\theta = E\[\text{payoff}\] = E\left\[\left( \bar{S} -K\right)^{+}\_j\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5B%5Ctext%7Bpayoff%7D%5D%20%3D%20E%5Cleft%5B%5Cleft%28%20%5Cbar%7BS%7D%20-K%5Cright%29%5E%7B%2B%7D_j%5Cright%5D "\theta = E[\text{payoff}] = E\left[\left( \bar{S} -K\right)^{+}_j\right]")

![\hat{\theta}\_{DIR} = \frac{1}{n} \sum\_{j=1}^{n} \left( \bar{s} -K\right)^{+}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D_%7BDIR%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bn%7D%20%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j "\hat{\theta}_{DIR} = \frac{1}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j")

![\therefore](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctherefore "\therefore")
the variance follows as

\$\$ the variance of the price of the Asian call options follows as:

# 4 Variance Reduction

Functions such as
![E\[g(X)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5Bg%28X%29%5D "E[g(X)]")
can be estimated using Monte Carlo integration. There are different
methods that can be applied to reduce the variance of the sample mean
estimator of
![\theta = E\[g(X)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5Bg%28X%29%5D "\theta = E[g(X)]")
(Rizzo 2019). Some of these methods include antithetic variables,
control variables, importance sampling, Latin hypercube sampling and
Quasi-monte carlo methods.

## 4.1 Antithetic Variables Method

### 4.1.1 Description

Suppose there are two identically distributed random variables
![Y_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_1 "Y_1")
and
![Y_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_2 "Y_2").
If
![Y_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_1 "Y_1")
and
![Y_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_2 "Y_2")
are independent, then:

![Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2))](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Var%28%5Cfrac%7BY_1%20%2B%20Y_2%7D%7B2%7D%29%20%3D%20%5Cfrac%7B1%7D%7B4%7D%28Var%28Y_1%29%20%2B%20Var%28Y_2%29%29 "Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2))")

otherwise, in general,

![Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2) + 2Cov(Y_1, Y_2))](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Var%28%5Cfrac%7BY_1%20%2B%20Y_2%7D%7B2%7D%29%20%3D%20%5Cfrac%7B1%7D%7B4%7D%28Var%28Y_1%29%20%2B%20Var%28Y_2%29%20%2B%202Cov%28Y_1%2C%20Y_2%29%29 "Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2) + 2Cov(Y_1, Y_2))")

From the above, it can be seen that if
![Y_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_1 "Y_1")
and
![Y_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_2 "Y_2")
are negatively correlated, then the variance of
![\left(\frac{Y_1 + Y_2}{2} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7BY_1%20%2B%20Y_2%7D%7B2%7D%20%5Cright%29 "\left(\frac{Y_1 + Y_2}{2} \right)")
is smaller than when
![Y_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_1 "Y_1")
and
![Y_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_2 "Y_2")
are independent (Rizzo 2019).

So, for the normal case, suppose that (Zhang 2009):

![\theta = E\[Y\] = E\[g(X)\] \quad \text{where} \\X \sim N(0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5BY%5D%20%3D%20E%5Bg%28X%29%5D%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5C%20X%20%5Csim%20N%280%2C1%29 "\theta = E[Y] = E[g(X)] \quad \text{where} \ X \sim N(0,1)")

where the direct Monte Carlo estimate is:

![\hat{\theta}\_{DIR} = \frac{1}{n} \sum\_{i=1}^{n} g(X_i) \quad \text{where}  \\X_i \overset{\text{i.i.d.}}{\sim} N(0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D_%7BDIR%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%20g%28X_i%29%20%5Cquad%20%5Ctext%7Bwhere%7D%20%20%5C%20X_i%20%5Coverset%7B%5Ctext%7Bi.i.d.%7D%7D%7B%5Csim%7D%20N%280%2C1%29 "\hat{\theta}_{DIR} = \frac{1}{n} \sum_{i=1}^{n} g(X_i) \quad \text{where}  \ X_i \overset{\text{i.i.d.}}{\sim} N(0,1)")

and the antithetic variable estimate is:

![\hat{\theta}\_{AV} = \frac{1}{\frac{n}{2}} \sum\_{i=1}^{\frac{n}{2}} \frac{g(X_i)+g(-X_i)}{2} \quad \text{where}  \\X_i \overset{\text{i.i.d.}}{\sim} N(0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D_%7BAV%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Csum_%7Bi%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Cfrac%7Bg%28X_i%29%2Bg%28-X_i%29%7D%7B2%7D%20%5Cquad%20%5Ctext%7Bwhere%7D%20%20%5C%20X_i%20%5Coverset%7B%5Ctext%7Bi.i.d.%7D%7D%7B%5Csim%7D%20N%280%2C1%29 "\hat{\theta}_{AV} = \frac{1}{\frac{n}{2}} \sum_{i=1}^{\frac{n}{2}} \frac{g(X_i)+g(-X_i)}{2} \quad \text{where}  \ X_i \overset{\text{i.i.d.}}{\sim} N(0,1)")

where
![X_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_i "X_i")
and
![-X_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-X_i "-X_i")
are the antithetic variates. The two antithetic variates are negatively
correlated. Therefore, if the function
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
is monotone increasing or decreasing, variance reduction can be achieved
by using this method (Zhang 2009).

### 4.1.2 Implementation of Method to Price Asian Call Option

The following expected value needs to be estimated for Asian call
options:

![\theta = E\[\text{payoff}\] = E\left\[\left(\bar{S}  -K\right)^{+}\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5B%5Ctext%7Bpayoff%7D%5D%20%3D%20E%5Cleft%5B%5Cleft%28%5Cbar%7BS%7D%20%20-K%5Cright%29%5E%7B%2B%7D%5Cright%5D "\theta = E[\text{payoff}] = E\left[\left(\bar{S}  -K\right)^{+}\right]")

1.  Generate m random standard normal observations
    (![z_1, z_2, ... , z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_1%2C%20z_2%2C%20...%20%2C%20z_m "z_1, z_2, ... , z_m")).
2.  Create another m standard normal observations that are negatively
    correlated with
    ![z_1, z_2,...,z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_1%2C%20z_2%2C...%2Cz_m "z_1, z_2,...,z_m")
    (![-z_1, -z_2, ... , -z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-z_1%2C%20-z_2%2C%20...%20%2C%20-z_m "-z_1, -z_2, ... , -z_m")).
3.  Generate a stock price path using
    ![z_1, z_2, ... , z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_1%2C%20z_2%2C%20...%20%2C%20z_m "z_1, z_2, ... , z_m")
    and another path using
    ![-z_1, -z_2, ... , -z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-z_1%2C%20-z_2%2C%20...%20%2C%20-z_m "-z_1, -z_2, ... , -z_m").
    Remember the model underlying the stock price behaviour:

![s\_{t_i} = s\_{t\_{i-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z\_{t\_{i}} \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_%7Bt_i%7D%20%3D%20s_%7Bt_%7Bi-1%7D%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20z_%7Bt_%7Bi%7D%7D%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_{t_i} = s_{t_{i-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z_{t_{i}} \sqrt{\triangle t}\right)}")

![s\_{t_i}^{'} = s\_{t\_{i-1}}^{'}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (-z\_{t\_{i}}) \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_%7Bt_i%7D%5E%7B%27%7D%20%3D%20s_%7Bt_%7Bi-1%7D%7D%5E%7B%27%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20%28-z_%7Bt_%7Bi%7D%7D%29%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_{t_i}^{'} = s_{t_{i-1}}^{'}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (-z_{t_{i}}) \sqrt{\triangle t}\right)}")

This will generate stock prices for m time steps. 4. Calculate the
payoff from each path (there will be 2) by summing up the m stock prices
and dividing by m and then subtracting K from this result:

![\text{Payoff} = \left(\bar{S} - K \right)^{+} \quad \text{for one path}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BPayoff%7D%20%3D%20%5Cleft%28%5Cbar%7BS%7D%20-%20K%20%5Cright%29%5E%7B%2B%7D%20%5Cquad%20%5Ctext%7Bfor%20one%20path%7D "\text{Payoff} = \left(\bar{S} - K \right)^{+} \quad \text{for one path}")

and

![\text{Payoff} =\left( \bar{S}' - K\right)^{+} \quad \text{for the other path}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BPayoff%7D%20%3D%5Cleft%28%20%5Cbar%7BS%7D%27%20-%20K%5Cright%29%5E%7B%2B%7D%20%5Cquad%20%5Ctext%7Bfor%20the%20other%20path%7D "\text{Payoff} =\left( \bar{S}' - K\right)^{+} \quad \text{for the other path}")

5.  Now use Monte Carlo simulation to repeat this process
    ![\frac{n}{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Bn%7D%7B2%7D "\frac{n}{2}")
    times to obtain n simulated paths representing n different outcomes
    and calculating n different payoffs.

![\hat{\theta} = \frac{1}{\frac{n}{2}} \sum\_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}\_j+\left( \bar{s}^{'} -K\right)^{+}\_j}{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Cfrac%7B%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j%2B%5Cleft%28%20%5Cbar%7Bs%7D%5E%7B%27%7D%20-K%5Cright%29%5E%7B%2B%7D_j%7D%7B2%7D "\hat{\theta} = \frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}")

6.  The price of a Asian call option is calculated by discounting the
    payoff:

![\hat{C}\_{AV} = e^{-rt}\left\[\frac{1}{\frac{n}{2}} \sum\_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}\_j+\left( \bar{s}^{'} -K\right)^{+}\_j}{2}\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7BAV%7D%20%3D%20e%5E%7B-rt%7D%5Cleft%5B%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Cfrac%7B%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j%2B%5Cleft%28%20%5Cbar%7Bs%7D%5E%7B%27%7D%20-K%5Cright%29%5E%7B%2B%7D_j%7D%7B2%7D%5Cright%5D "\hat{C}_{AV} = e^{-rt}\left[\frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}\right]")

### 4.1.3 Variance of Asian Call Option Price

![\theta = E\[\text{payoff}\]  = E\left\[\left( \bar{S} -K\right)^{+}\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5B%5Ctext%7Bpayoff%7D%5D%20%20%3D%20E%5Cleft%5B%5Cleft%28%20%5Cbar%7BS%7D%20-K%5Cright%29%5E%7B%2B%7D%5Cright%5D "\theta = E[\text{payoff}]  = E\left[\left( \bar{S} -K\right)^{+}\right]")

![\hat{\theta}\_{AV} = \frac{1}{\frac{n}{2}} \sum\_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}\_j+\left( \bar{s}^{'} -K\right)^{+}\_j}{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D_%7BAV%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7B2%7D%7D%20%5Cfrac%7B%5Cleft%28%20%5Cbar%7Bs%7D%20-K%5Cright%29%5E%7B%2B%7D_j%2B%5Cleft%28%20%5Cbar%7Bs%7D%5E%7B%27%7D%20-K%5Cright%29%5E%7B%2B%7D_j%7D%7B2%7D "\hat{\theta}_{AV} = \frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}")

![\therefore](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctherefore "\therefore")
the variance follows as (Rizzo 2019): \`\`{=tex} where Y is the payoff
using standard normal random variables
(![Z_1, Z_2, ... , Z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z_1%2C%20Z_2%2C%20...%20%2C%20Z_m "Z_1, Z_2, ... , Z_m"))
and Y’ is the payoff using standard normal random variables that are
negatively correlated to
(![Z_1, Z_2, ... , Z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z_1%2C%20Z_2%2C%20...%20%2C%20Z_m "Z_1, Z_2, ... , Z_m"))
i.e. (![-Z_1, -Z_2, ... , -Z_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-Z_1%2C%20-Z_2%2C%20...%20%2C%20-Z_m "-Z_1, -Z_2, ... , -Z_m")).
Since the payoff function is monotone increasing, this method will
reduce variance.

![\therefore](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctherefore "\therefore")
the variance of the price of the Asian call option is:

where Y and Y’ are defined as above.

## 4.2 n-Simplex Antithetic Variate Method

### 4.2.1 The Simplex

A simplex is an object in geometry that generalises the idea of a
triangle or tetrahedron to D-dimensions. In D dimensions, a simplex is
determined by D+1 points. These are called the vertices of the simplex.
A simplex in D dimensions is the convex hull of D+1 points
![\mathbf{v}\_1, \mathbf{v}\_2,...,\mathbf{v}\_{D+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1%2C%20%5Cmathbf%7Bv%7D_2%2C...%2C%5Cmathbf%7Bv%7D_%7BD%2B1%7D "\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}").
In 1D, this would be a straight line, 2D a triangle, 3D a tetrahedron
and 4D a Pentatope. A regular D-simplex is a simplex where the euclidean
distance between each of the D+1 points is equal and the angles between
each face is the same. A regular D-simplex centred at the origin of a
unit D-dimensional hypersphere is of interest as it has specific
properties that will be shown later, which allow us to achieve variance
reduction. From here on, this will be referred to as the regular
D-simplex.The idea of the regular simplex is used to extend antithetic
variables to multiple dimensions for random variables with symmetric
probability distributions (Devriendt and Van Mieghem 2019). An example
of a regular 2-simplex is shown below.

![](stoch_sim_project_files/figure-gfm/Unit%20Circle-1.png)<!-- -->

In two dimensions, as above, a simplex is determined by three points.
These are points on a unit circle which correspond to the vertices of an
inscribed equilateral triangle. Any three points can be chosen on the
unit circle as long as these points form an equilateral triangle. One
such example is:
![\left(1,0\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%281%2C0%5Cright%29 "\left(1,0\right)"),
![\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B2%7D%2C%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%5Cright%29 "\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)")
and
![\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B2%7D%2C%20-%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%5Cright%29 "\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)")
(Park and Choe 2016). Confirm these points are correct by calculating
the Euclidean distance between each point and ensuring it is the same
for each pair of points.

### 4.2.2 n-Simplex Method

Let
![\mathbf{v}\_1, \mathbf{v}\_2,...,\mathbf{v}\_{D+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1%2C%20%5Cmathbf%7Bv%7D_2%2C...%2C%5Cmathbf%7Bv%7D_%7BD%2B1%7D "\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}")
be the vertices of a regular D-simplex. The vertices have the following
properties:

![\lVert \mathbf{v}\_i \rVert ^2 = 1 \quad \forall i=1,2, ..., D+1 \quad (\text{by definition})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5ClVert%20%5Cmathbf%7Bv%7D_i%20%5CrVert%20%5E2%20%3D%201%20%5Cquad%20%5Cforall%20i%3D1%2C2%2C%20...%2C%20D%2B1%20%5Cquad%20%28%5Ctext%7Bby%20definition%7D%29 "\lVert \mathbf{v}_i \rVert ^2 = 1 \quad \forall i=1,2, ..., D+1 \quad (\text{by definition})")

and

<center>

![\mathbf{v}\_i ' \mathbf{v}\_j = -\frac{1}{D} \quad \text{for} \\i\neq j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i%20%27%20%5Cmathbf%7Bv%7D_j%20%3D%20-%5Cfrac%7B1%7D%7BD%7D%20%5Cquad%20%5Ctext%7Bfor%7D%20%5C%20i%5Cneq%20j "\mathbf{v}_i ' \mathbf{v}_j = -\frac{1}{D} \quad \text{for} \ i\neq j")
(see Parks and Wills (2002))

</center>

Let \$={Z}\_1, {Z}\_2,…,{Z}\_D \$ be a vector of D standard normal
random variables, then:

![\mathbf{v}\_i'\mathbf{Z} \sim N(0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i%27%5Cmathbf%7BZ%7D%20%5Csim%20N%280%2C1%29 "\mathbf{v}_i'\mathbf{Z} \sim N(0,1)")

and

![\text{Cov}(\mathbf{v}\_i'\mathbf{Z}, \mathbf{v}\_j'\mathbf{Z})= -\frac{1}{D}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BCov%7D%28%5Cmathbf%7Bv%7D_i%27%5Cmathbf%7BZ%7D%2C%20%5Cmathbf%7Bv%7D_j%27%5Cmathbf%7BZ%7D%29%3D%20-%5Cfrac%7B1%7D%7BD%7D "\text{Cov}(\mathbf{v}_i'\mathbf{Z}, \mathbf{v}_j'\mathbf{Z})= -\frac{1}{D}")

The D-simplex average of
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
is: \$\$

S() = \_{i=1}^{D+1}g(\_i’) \$\$

![\therefore](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctherefore "\therefore")
the Monte Carlo estimator of
![E\[S(\mathbf{Z})\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;E%5BS%28%5Cmathbf%7BZ%7D%29%5D "E[S(\mathbf{Z})]")
is:

![\frac{1}{\frac{n}{D+1}}\sum\_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}\_j)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%20S%28%5Cmathbf%7BZ%7D_j%29 "\frac{1}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)")

It is now shown that using this estimator will reduce variance from the
direct method: It is shown in theory that two monotone increasing
functions where the variables in each function are negatively correlated
has a covariance of less than or equal to zero. Hence, Hence, a variance
reduction is achieved with a higher dimesnions resulting in a higher
reduction. It follows that the case for which the dimension is one is
simply the antithetic variables method.

### 4.2.3 Implementation of Method to Price Asian Call Option

Consider the implementation of the n-simplex method to price Asian call
options. Let D be the number of dimensions, therefore D+1 points are
used to determine the simplex: Remember, the following expected value
needs to be estimated for Asian call options:

![\theta = E\[\text{payoff}\] = E\left\[\left(\bar{S}  -K\right)^{+}\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20E%5B%5Ctext%7Bpayoff%7D%5D%20%3D%20E%5Cleft%5B%5Cleft%28%5Cbar%7BS%7D%20%20-K%5Cright%29%5E%7B%2B%7D%5Cright%5D "\theta = E[\text{payoff}] = E\left[\left(\bar{S}  -K\right)^{+}\right]")

1.  Compute D+1 vertices of the regular D-simplex:
    ![\mathbf{v}\_1, \mathbf{v}\_2,...,\mathbf{v}\_{D+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1%2C%20%5Cmathbf%7Bv%7D_2%2C...%2C%5Cmathbf%7Bv%7D_%7BD%2B1%7D "\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}").
2.  Generate D vectors,
    ![\mathbf{z} = \mathbf{z}\_1,\mathbf{z}\_2,...,\mathbf{z}\_D](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%20%3D%20%5Cmathbf%7Bz%7D_1%2C%5Cmathbf%7Bz%7D_2%2C...%2C%5Cmathbf%7Bz%7D_D "\mathbf{z} = \mathbf{z}_1,\mathbf{z}_2,...,\mathbf{z}_D"),
    each containing m observations of standard normal random variables.

![\mathbf{z}^{(1)}=z_1^{(1)}, z_2^{(1)}, ..., z_m^{(1)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%5E%7B%281%29%7D%3Dz_1%5E%7B%281%29%7D%2C%20z_2%5E%7B%281%29%7D%2C%20...%2C%20z_m%5E%7B%281%29%7D "\mathbf{z}^{(1)}=z_1^{(1)}, z_2^{(1)}, ..., z_m^{(1)}")

![\mathbf{z}^{(2)}=z_1^{(2)}, z_2^{(2)}, ..., z_m^{(2)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%5E%7B%282%29%7D%3Dz_1%5E%7B%282%29%7D%2C%20z_2%5E%7B%282%29%7D%2C%20...%2C%20z_m%5E%7B%282%29%7D "\mathbf{z}^{(2)}=z_1^{(2)}, z_2^{(2)}, ..., z_m^{(2)}")

<center>

. . .

</center>

![\mathbf{z}^{(D)}=z_1^{(D)}, z_2^{(D)}, ..., z_m^{(D)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%5E%7B%28D%29%7D%3Dz_1%5E%7B%28D%29%7D%2C%20z_2%5E%7B%28D%29%7D%2C%20...%2C%20z_m%5E%7B%28D%29%7D "\mathbf{z}^{(D)}=z_1^{(D)}, z_2^{(D)}, ..., z_m^{(D)}")

3.  Generate D+1 stock price paths. Remember the model underlying the
    stock price behaviour:

![s\_{t_j}^{(1)} = s\_{t\_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x\_{ij}) \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_%7Bt_j%7D%5E%7B%281%29%7D%20%3D%20s_%7Bt_%7Bj-1%7D%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20%28x_%7Bij%7D%29%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_{t_j}^{(1)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}")

![s\_{t_j}^{(2)} = s\_{t\_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x\_{ij}) \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_%7Bt_j%7D%5E%7B%282%29%7D%20%3D%20s_%7Bt_%7Bj-1%7D%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20%28x_%7Bij%7D%29%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_{t_j}^{(2)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}")

<center>

. . .

</center>

![s\_{t_j}^{(D+1)} = s\_{t\_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x\_{ij}) \sqrt{\triangle t}\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s_%7Bt_j%7D%5E%7B%28D%2B1%29%7D%20%3D%20s_%7Bt_%7Bj-1%7D%7De%5E%7B%5Cleft%28%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%29%5Ctriangle%20t%20%2B%20%5Csigma%20%28x_%7Bij%7D%29%20%5Csqrt%7B%5Ctriangle%20t%7D%5Cright%29%7D "s_{t_j}^{(D+1)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}")

for
![j = 1,...,m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j%20%3D%201%2C...%2Cm "j = 1,...,m")
and
![i = 1,2,...,D+1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%3D%201%2C2%2C...%2CD%2B1 "i = 1,2,...,D+1"),
where
![x\_{ij}=\[ \mathbf{z}\_j^{(1)},\mathbf{z}\_j^{(2)},...,\mathbf{z}\_j^{(D)} \]' \mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bij%7D%3D%5B%20%5Cmathbf%7Bz%7D_j%5E%7B%281%29%7D%2C%5Cmathbf%7Bz%7D_j%5E%7B%282%29%7D%2C...%2C%5Cmathbf%7Bz%7D_j%5E%7B%28D%29%7D%20%5D%27%20%5Cmathbf%7Bv%7D_i "x_{ij}=[ \mathbf{z}_j^{(1)},\mathbf{z}_j^{(2)},...,\mathbf{z}_j^{(D)} ]' \mathbf{v}_i")
This generates the stock price path for m time steps. 4. Calculate the
payoff from each path (there are D+1 paths) by summing up the m stock
prices and dividing by m and then subtracting K from this result:

![\text{Payoff} = (\bar{S} - K)^+ \\\text{for D+1 paths}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BPayoff%7D%20%3D%20%28%5Cbar%7BS%7D%20-%20K%29%5E%2B%20%5C%20%5Ctext%7Bfor%20D%2B1%20paths%7D "\text{Payoff} = (\bar{S} - K)^+ \ \text{for D+1 paths}")

5.  Then calculate the average of the the D+1 payoffs:

![S(\mathbf{Z}) = \frac{1}{D+1} \sum\_{i=1}^{D+1}\text{payoff}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%28%5Cmathbf%7BZ%7D%29%20%3D%20%5Cfrac%7B1%7D%7BD%2B1%7D%20%5Csum_%7Bi%3D1%7D%5E%7BD%2B1%7D%5Ctext%7Bpayoff%7D_i "S(\mathbf{Z}) = \frac{1}{D+1} \sum_{i=1}^{D+1}\text{payoff}_i")

6.  Repeat steps 2 to 5
    ![\frac{n}{D+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Bn%7D%7BD%2B1%7D "\frac{n}{D+1}")
    times.
7.  Calculate the mean of the
    ![\frac{n}{D+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Bn%7D%7BD%2B1%7D "\frac{n}{D+1}")
    average payoffs.

![\frac{1}{\frac{n}{D+1}}\sum\_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}\_j)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%7D%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%20S%28%5Cmathbf%7BZ%7D_j%29 "\frac{1}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)")

8.  Discount this using the risk-free rate to obtain the price of the
    Asian call option.

![\frac{e^{-rt}}{\frac{n}{D+1}}\sum\_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}\_j)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Be%5E%7B-rt%7D%7D%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%5Csum_%7Bj%3D1%7D%5E%7B%5Cfrac%7Bn%7D%7BD%2B1%7D%7D%20S%28%5Cmathbf%7BZ%7D_j%29 "\frac{e^{-rt}}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)")

# 5 Analyis

## 5.1 Results

The methods were compared using the following constants:

|  S  |  r  | ![\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma "\sigma") |  t  |  m  |
|:---:|:---:|:---------------------------------------------------------------------------------------------------------:|:---:|:---:|
| 100 | 5%  |                                                    20%                                                    |  1  | 250 |

where: - S is the initial stock price. - r is the risk-free rate of
interest. -
![\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma "\sigma")
is the volatility of the underlying asset. - t is the expiration date. -
m is the number of of time steps.

The price and standard error for each method was obtained and compared.
All methods were implemented using the algorithms that have been
described. The n-simplex method was done in 2,3 and 4 dimensions. For
the 2-Simplex, the points
![\left(1,0\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%281%2C0%5Cright%29 "\left(1,0\right)"),
![\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B2%7D%2C%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%5Cright%29 "\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)")
and
![\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B2%7D%2C%20-%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%5Cright%29 "\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)")
were used. For the 3-Simplex, the points
![\left(1,0,0\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%281%2C0%2C0%5Cright%29 "\left(1,0,0\right)"),
![\left(-\frac{1}{3}, \frac{\sqrt{5}}{3}, \frac{\sqrt{3}}{3}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B3%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B3%7D%2C%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B3%7D%5Cright%29 "\left(-\frac{1}{3}, \frac{\sqrt{5}}{3}, \frac{\sqrt{3}}{3}\right)"),
![\left(-\frac{1}{3}, \frac{-3-\sqrt{5}}{6},\frac{\sqrt{15}- \sqrt{3}}{6}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B3%7D%2C%20%5Cfrac%7B-3-%5Csqrt%7B5%7D%7D%7B6%7D%2C%5Cfrac%7B%5Csqrt%7B15%7D-%20%5Csqrt%7B3%7D%7D%7B6%7D%5Cright%29 "\left(-\frac{1}{3}, \frac{-3-\sqrt{5}}{6},\frac{\sqrt{15}- \sqrt{3}}{6}\right)")
and
![\left(-\frac{1}{3}, \frac{3-\sqrt{5}}{6},\frac{-\sqrt{15}- \sqrt{3}}{6}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B3%7D%2C%20%5Cfrac%7B3-%5Csqrt%7B5%7D%7D%7B6%7D%2C%5Cfrac%7B-%5Csqrt%7B15%7D-%20%5Csqrt%7B3%7D%7D%7B6%7D%5Cright%29 "\left(-\frac{1}{3}, \frac{3-\sqrt{5}}{6},\frac{-\sqrt{15}- \sqrt{3}}{6}\right)").
For the 4-Simplex, the points
![\left(1,0,0,0\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%281%2C0%2C0%2C0%5Cright%29 "\left(1,0,0,0\right)"),
![\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%5Cright%29 "\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)"),
![\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20-%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%5Cright%29 "\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)"),
![\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B4%7D%2C%20-%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20-%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%5Cright%29 "\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)"),
![\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28-%5Cfrac%7B1%7D%7B4%7D%2C%20-%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20-%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%2C%20%5Cfrac%7B%5Csqrt%7B5%7D%7D%7B4%7D%5Cright%29 "\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)")
were used.

The results are shown in the following table:

**Table 1: Price and Standard Error from the Respective Methods**

| K       | n      | Direct Price | Direct Standard Error | Antithetic Price | Antithetic Standard Error | 2-Simplex Price | 2-Simplex Standard Error | 3-Simplex Price | 3-Simplex Standard Error | 4-Simplex Price | 4-Simplex Standard Error |
|:--------|:-------|--------------|:----------------------|:-----------------|:--------------------------|:----------------|:-------------------------|:----------------|:-------------------------|:----------------|:-------------------------|
| **90**  | 1000   | 12.5533239   | 0.3268601             | 12.7055017       | 0.0877486                 | 12.6131679      | 0.0633758                | 12.7477091      | 0.0505096                | 12.5424661      | 0.0409757                |
|         | 10000  | 12.7679562   | 0.1043739             | 12.5559164       | 0.024916                  | 12.6059742      | 0.0184413                | 12.5912949      | 0.0154244                | 12.5827063      | 0.0134648                |
|         | 100000 | 12.614105    | 0.0331502             | 12.6098847       | 0.0081651                 | 12.5981434      | 0.0059089                | 12.6006781      | 0.0049104                | 12.6151473      | 0.0042368                |
| **100** | 1000   | 5.6471417    | 0.2537597             | 5.7741284        | 0.1208521                 | 5.5907299       | 0.085626                 | 5.8907499       | 0.0814799                | 5.6898131       | 0.0547328                |
|         | 10000  | 5.7702871    | 0.0810186             | 5.7137139        | 0.0386967                 | 5.7523675       | 0.0281161                | 5.8060386       | 0.0230676                | 5.7521011       | 0.0194926                |
|         | 100000 | 5.7962818    | 0.0252594             | 5.7766894        | 0.0123769                 | 5.7753943       | 0.0087051                | 5.7850407       | 0.0071326                | 5.7781959       | 0.0061595                |
| **110** | 1000   | 2.3031334    | 0.1688262             | 1.9106531        | 0.0970899                 | 1.9785303       | 0.0745425                | 1.8750571       | 0.0630819                | 1.8457502       | 0.0494986                |
|         | 10000  | 1.9469827    | 0.0490162             | 2.0264369        | 0.0316543                 | 1.9929129       | 0.0236604                | 2.0263966       | 0.0196898                | 2.0566903       | 0.0169199                |
|         | 100000 | 2.0031709    | 0.0154814             | 2.0072226        | 0.0100539                 | 1.9993602       | 0.007527                 | 1.9995149       | 0.00622                  | 2.0103173       | 0.0054524                |

The variance reduction achieved from each method is summarised as
follows:

**Table 2: Variance Reduction Compared to the Direct Method**

|         |        | **Antithetic**         | **2-Simplex**          | **3-Simplex**          | **4-Simplex**          |
|---------|--------|:-----------------------|:-----------------------|:-----------------------|:-----------------------|
| **K**   | **n**  | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** |
| **90**  | 1000   | 0.9279296              | 0.9624056              | 0.9761206              | 0.9842845              |
|         | 10000  | 0.9430132              | 0.9687823              | 0.978161               | 0.9833575              |
|         | 100000 | 0.9393329              | 0.9682283              | 0.9780586              | 0.9836652              |
| **100** | 1000   | 0.7731895              | 0.8861414              | 0.8969007              | 0.9534789              |
|         | 10000  | 0.7718719              | 0.8795687              | 0.9189343              | 0.9421144              |
|         | 100000 | 0.7599087              | 0.8812328              | 0.9202642              | 0.940537               |
| **110** | 1000   | 0.6692737              | 0.8050476              | 0.8603856              | 0.9140381              |
|         | 10000  | 0.582953               | 0.7669957              | 0.838637               | 0.8808436              |
|         | 100000 | 0.5782544              | 0.7636135              | 0.8385778              | 0.875963               |

## 5.2 Analysis of results

The output from the tables above show that significant variance
reduction was achieved for all the variance reduction methods. The
standard error decreases as the number of repetitions increase, as
expected.The n-simplex method also outperforms the antithetic variable
method. It is seen that as the number of dimensions for the n-simplex
increases, the variance reduction increases. This agrees with the theory
of the variance reduction shown in 4.2.2 . It is interesting to note
that the variance reduction for all methods decreases as the strike
price increases. This could be due to the fact that more of the paths
with out-the-money call options will provide a payoff of zero compared
to in-the-money call options. The covariance between any random variable
and zero is always zero, hence a zero payoff does not provide negative
covariance. This means more payoffs of zero results in lower variance
reduction.

# 6 Conclusion

This project illustrates the usefulness of using Monte Carlo simulation
to price Asian options. Furthermore, we see how variance reduction
methods can improve the accuracy of the prices. The n-simplex is shown
to be a useful variance reduction technique that outperforms antithetic
variables in terms of variance reduction. However, it is slightly more
complicated and harder to implement.

# 7 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-SimplexDevriendt" class="csl-entry">

Devriendt, Karel, and Piet Van Mieghem. 2019. “The Simplex Geometry of
Graphs.” *Journal of Complex Networks* 7 (4): 469–90.

</div>

<div id="ref-HullTextbook" class="csl-entry">

Hull, John C. 2017. *Options, Futures, and Other Derivatives*. Global
edition. Pearson Education.

</div>

<div id="ref-AboutAsianOptions" class="csl-entry">

Nielsen, Lars B. 2001. “Pricing Asian Options.” *Aarhus Universitet,
Institut for Matematiske Fag, Afdeling for Nationaløkonomi*.

</div>

<div id="ref-SimplexMethod" class="csl-entry">

Park, Jong Jun, and Geon Ho Choe. 2016. “A New Variance Reduction Method
for Option Pricing Based on Sampling the Vertices of a Simplex.”
*Quantitative Finance* 16 (8): 1165–73.

</div>

<div id="ref-SimplexProperties" class="csl-entry">

Parks, Harold R, and Dean C Wills. 2002. “An Elementary Calculation of
the Dihedral Angle of the Regular n-Simplex.” *The American Mathematical
Monthly* 109 (8): 756–58.

</div>

<div id="ref-StochSimTextbook" class="csl-entry">

Rizzo, Maria L. 2019. *Statistical Computing with r*. Chapman; Hall/CRC.

</div>

<div id="ref-PricingAsianOptions" class="csl-entry">

Zhang, Hongbin. 2009. *Pricing Asian Options Using Monte Carlo Methods*.

</div>

</div>
