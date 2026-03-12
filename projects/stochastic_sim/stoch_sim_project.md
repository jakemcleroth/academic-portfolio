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

The payoff is given by (Zhang 2009): $$ \phi(S) = (S(T) - K)^+ $$ for a
European call option and by $$ \phi(S) = (K- S(T))^+ $$ for a European
put option, where <br> - \$ T \$ is the expiration date,<br> - \$ S(T)
\$ is the spot price of the underlying asset at expiration T, <br> - \$
K \$ is the strike price of the option, <br> - \$ (x)^+ \$ represents
the maximum of \$ x    0\$. <br>

## 2.2 Geometric Brownian Motion

To model the behaviour of stock prices, the most commonly used model is
geometric Brownian motion. It can be formulated as the following
stochastic differential equation (SDE) as seen in Hull (2017): $$
dS_t=\mu S_t dt + \sigma S_t d W_t
$$ where:<br> - $dS$ is the instantaneous change in the stock
price,<br> - $dW$ is the Wiener process ($\triangle W_t$ is a random
variable from the normal distribution $\text{N}(0,\triangle t)$), <br> -
\$\$ represents the volatility of the stock price ($\sigma$ is assumed
to be constant over time), <br> - $\mu$ represents the growth rate of
the stock price when there is no risk (also referred to as the ‘drift’),
<br> - $\mu dt$ represents the deterministic return within the time
interval.<br>

with the discrete time version of the SDE as: $$
\triangle S_t = \mu S_t \triangle t + \sigma S_t \triangle W_t
$$

The solution to the discrete SDE is as follows: where $Z$ is a standard
normal random variable.

## 2.3 Options Pricing

The main idea behind options pricing, or derivative pricing in general,
is known as risk neutral valuation. The idea is that the price of an
option should be equal to the present value of the expected payoff.
According to Hull (2017), this can be done in the following
steps:<br> 1. Assume expected rate of return is the risk free rate. <br>
2. Calculate the expected payoff of the option. For example, a European
call option: $$E[\text{payoff}] = E[(S(T) - K)^+]$$ 3. Discount the
expected payoff at the risk free rate. For the European call option:
$$ C = e^{-rt}E[\text{payoff}]$$ This European call option can be priced
using the Black-Scholes-Merton framework with a closed form solution.
However, there is no closed form solution for pricing arithmetic Asian
options. Hence we need to use simulation.

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
structure: $$
\Phi(S) = \left(\bar{S}  - K \right)^{+} \quad \text{for a call}$$ or
$$\Phi(S) = \left( K - \bar{S}  \right)^{+} \quad \text{for a put}$$
where <br>
$$\text{where} \ \bar{S} = \frac{1}{m} \sum_{i=1}^{m} S(t_i)$$

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
$\int_{a}^{b} g(x) \, dx$, provided the integral exists. If \$ X \$ is a
random variable with probability density function \$ f(x) \$, then the
expected value of the random variable \$ Y = g(X) \$ is given by:
$$E[g(X)] = \int_{-\infty}^{\infty} g(x)f(x) \, dx$$ When a random
sample is available from the distribution of $X$, an unbiased estimator
of \$ E$$g(X)$$ \$ is simply the sample mean (Rizzo 2019). <br> Suppose
$f(x)$ is a probability density function, such that $f(x) \geq 0$ for
all $x \in \mathbb{R}$ and $\int_{A} f(x) \, dx = 1$. In order to
estimate the integral $\theta = \int_{A} g(x)f(x) \, dx$, generate a
random sample $x_1, x_2, ..., x_n$ from the distribution of $f(x)$ and
calculate the sample mean
$$\hat{\theta} = \frac{1}{n} \sum_{j=1}^{n} g(x_j)$$

The core principle of Monte Carlo integration relies on the strong law
of large numbers: $$
\hat{\theta} = \frac{1}{n} \sum_{j=1}^{n} g(x_j) \rightarrow  E[g(X)]\quad a.s.
$$ $\hat{\theta}$ as defined above is the Monte Carlo estimator of
$E[g(X)]$ (Rizzo 2019).

## 3.2 Pricing Asian Options

### 3.2.1 Direct Method

We estimate $\theta = E[g(X)]$ by
$\hat{\theta}_{DIR} = \frac{1}{n} \sum_{j=1}^{n} g(X_j)$ where
$X_1, X_2, ..., X_n$ are independent and identically distributed. The
direct method refers to the basic Monte Carlo simulation without any
variance reduction techniques. Random samples are generated directly
from the distribution of interest.

In the case of pricing Asian call options, the direct method involves
generating random samples of the underlying asset’s price at various
future time points based on a stochastic model, such as geometric
Brownian motion. These simulated price paths are then used to calculate
the payoff of the option at the expiration date. The option’s price is
estimated by discounting back the average of the payoffs. <br> Recall
the payoff of a Asian call option is (Zhang 2009):
$$\Phi(S) = \left( \bar{S} - K \right)^{+}$$

By risk-neutral valuation, the price of the Asian call option can be
expressed as:

Monte Carlo simulation can be applied to price Asian options.
Substituting $g$ from above with the payoff of the Asian call option,
the following is obtained: $$
\frac{1}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j \rightarrow E\left[\left( \bar{S} -K\right)^{+}\right] \quad a.s.
$$ where:<br> - $\left( \bar{s} -K\right)^{+}_j$ is the payoff of the
Asian call option from a simulated stock path,<br> - n is the number of
Monte Carlo simulations. <br> Hence, the estimated price of the call
option can be expressed as follows:
$$\hat{C} = \frac{e^{-rt}}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j$$

### 3.2.2 Implementation of Method to Price Asian Call Option

To simulate stock paths, we make use of the geometric Brownian motion
explained earlier. The simulation goes as follows:<br> 1. Generate m
random variable observations, $z_1, z_2, ...,z_m$, from the standard
normal distribution. <br> 2. Simulate the path of the stock at m time
points using: $$
s_t = s_{t-1}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z_t \sqrt{\triangle t}\right)}
$$ 3. Repeat n times to generate n paths.

10 simulated paths are shown below:<br>
![](stoch_sim_project_files/figure-gfm/Asian%20Call%20Option%2010-1.png)<!-- -->
<br> In the simulation above, it can be seen that the stock price path
begins at 100 for each simulation since this is the stock price at time
0. The number of time steps chosen is 250 since there are roughly 250
trading days in 1 year. Each “day’s” price (i.e. time step 1, 2, 3 etc.)
is determined by the price the day before and is modelled according to
geometric Brownian motion. This creates one path. This process is
repeated n times to create a wide variety of different outcomes. In this
case, ten stock price paths are simulated, so the simulation has
produced ten potential trajectories of the stock price, based on random
fluctuations. <br>

The Asian call option is then priced with the following steps: <br> 1.
Simulate a stock path with geometric Brownian motion. 2. Sum up the
prices on each of the m days and divide by m to obtain the average price
of the asset over the entire term \$ \_{i=1}^{m} s(t_i) = {s}$. <br>
3. Subtract K from the average stock price over the entire term ($-K$). Then the payoff of the Asian call option is max$({s}-K,0)\$.
<br> 4. Repeat steps 1-3 n times to obtain n different payoffs.<br> 5.
Take the average of all n payoffs and discount it with $e^{-rt}$ where t
is the expiration date. <br> <br> The more paths that are generated, the
more accuracte the price will be. Shown below is an example of
simulating 10000 stock paths.
![](stoch_sim_project_files/figure-gfm/Asian%20Call%20Option%2010000-1.png)<!-- -->

### 3.2.3 Variance of Asian Call Option Price

$$
\theta = E[\text{payoff}] = E\left[\left( \bar{S} -K\right)^{+}_j\right]
$$ $$
\hat{\theta}_{DIR} = \frac{1}{n} \sum_{j=1}^{n} \left( \bar{s} -K\right)^{+}_j
$$ $\therefore$ the variance follows as

\$\$ the variance of the price of the Asian call options follows as:

# 4 Variance Reduction

Functions such as $E[g(X)]$ can be estimated using Monte Carlo
integration. There are different methods that can be applied to reduce
the variance of the sample mean estimator of $\theta = E[g(X)]$ (Rizzo
2019). Some of these methods include antithetic variables, control
variables, importance sampling, Latin hypercube sampling and Quasi-monte
carlo methods.

## 4.1 Antithetic Variables Method

### 4.1.1 Description

Suppose there are two identically distributed random variables $Y_1$ and
$Y_2$. If $Y_1$ and $Y_2$ are independent, then:
$$Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2))$$
otherwise, in general,
$$Var(\frac{Y_1 + Y_2}{2}) = \frac{1}{4}(Var(Y_1) + Var(Y_2) + 2Cov(Y_1, Y_2))$$
From the above, it can be seen that if $Y_1$ and $Y_2$ are negatively
correlated, then the variance of $\left(\frac{Y_1 + Y_2}{2} \right)$ is
smaller than when $Y_1$ and $Y_2$ are independent (Rizzo 2019).

So, for the normal case, suppose that (Zhang 2009): $$
\theta = E[Y] = E[g(X)] \quad \text{where} \ X \sim N(0,1)
$$ where the direct Monte Carlo estimate is: $$
\hat{\theta}_{DIR} = \frac{1}{n} \sum_{i=1}^{n} g(X_i) \quad \text{where}  \ X_i \overset{\text{i.i.d.}}{\sim} N(0,1)
$$ and the antithetic variable estimate is: $$
\hat{\theta}_{AV} = \frac{1}{\frac{n}{2}} \sum_{i=1}^{\frac{n}{2}} \frac{g(X_i)+g(-X_i)}{2} \quad \text{where}  \ X_i \overset{\text{i.i.d.}}{\sim} N(0,1)
$$ where $X_i$ and $-X_i$ are the antithetic variates. The two
antithetic variates are negatively correlated. Therefore, if the
function $g$ is monotone increasing or decreasing, variance reduction
can be achieved by using this method (Zhang 2009).

### 4.1.2 Implementation of Method to Price Asian Call Option

The following expected value needs to be estimated for Asian call
options: $$
\theta = E[\text{payoff}] = E\left[\left(\bar{S}  -K\right)^{+}\right]
$$ 1. Generate m random standard normal observations
($z_1, z_2, ... , z_m$). <br> 2. Create another m standard normal
observations that are negatively correlated with $z_1, z_2,...,z_m$
($-z_1, -z_2, ... , -z_m$). <br> 3. Generate a stock price path using
$z_1, z_2, ... , z_m$ and another path using $-z_1, -z_2, ... , -z_m$.
Remember the model underlying the stock price behaviour: $$
s_{t_i} = s_{t_{i-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma z_{t_{i}} \sqrt{\triangle t}\right)}
$$ $$
s_{t_i}^{'} = s_{t_{i-1}}^{'}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (-z_{t_{i}}) \sqrt{\triangle t}\right)}
$$ This will generate stock prices for m time steps. <br> 4. Calculate
the payoff from each path (there will be 2) by summing up the m stock
prices and dividing by m and then subtracting K from this result: $$ 
\text{Payoff} = \left(\bar{S} - K \right)^{+} \quad \text{for one path}
$$ and $$
\text{Payoff} =\left( \bar{S}' - K\right)^{+} \quad \text{for the other path} 
$$ 5. Now use Monte Carlo simulation to repeat this process
$\frac{n}{2}$ times to obtain n simulated paths representing n different
outcomes and calculating n different payoffs. $$
\hat{\theta} = \frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}
$$ 6. The price of a Asian call option is calculated by discounting the
payoff: $$
\hat{C}_{AV} = e^{-rt}\left[\frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}\right]
$$

### 4.1.3 Variance of Asian Call Option Price

$$
\theta = E[\text{payoff}]  = E\left[\left( \bar{S} -K\right)^{+}\right]
$$ $$
\hat{\theta}_{AV} = \frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}
$$ $\therefore$ the variance follows as (Rizzo 2019): \`\`{=tex} where Y
is the payoff using standard normal random variables
($Z_1, Z_2, ... , Z_m$) and Y’ is the payoff using standard normal
random variables that are negatively correlated to
($Z_1, Z_2, ... , Z_m$) i.e. ($-Z_1, -Z_2, ... , -Z_m$). Since the
payoff function is monotone increasing, this method will reduce
variance.

$\therefore$ the variance of the price of the Asian call option is:

where Y and Y’ are defined as above.

## 4.2 n-Simplex Antithetic Variate Method

### 4.2.1 The Simplex

A simplex is an object in geometry that generalises the idea of a
triangle or tetrahedron to D-dimensions. In D dimensions, a simplex is
determined by D+1 points. These are called the vertices of the simplex.
A simplex in D dimensions is the convex hull of D+1 points
$\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}$. In 1D, this would be
a straight line, 2D a triangle, 3D a tetrahedron and 4D a Pentatope. A
regular D-simplex is a simplex where the euclidean distance between each
of the D+1 points is equal and the angles between each face is the same.
A regular D-simplex centred at the origin of a unit D-dimensional
hypersphere is of interest as it has specific properties that will be
shown later, which allow us to achieve variance reduction. From here on,
this will be referred to as the regular D-simplex.The idea of the
regular simplex is used to extend antithetic variables to multiple
dimensions for random variables with symmetric probability distributions
(Devriendt and Van Mieghem 2019). An example of a regular 2-simplex is
shown below.

![](stoch_sim_project_files/figure-gfm/Unit%20Circle-1.png)<!-- --> <br>
In two dimensions, as above, a simplex is determined by three points.
These are points on a unit circle which correspond to the vertices of an
inscribed equilateral triangle. Any three points can be chosen on the
unit circle as long as these points form an equilateral triangle. One
such example is: $\left(1,0\right)$,
$\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)$ and
$\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)$ (Park and Choe 2016).
Confirm these points are correct by calculating the Euclidean distance
between each point and ensuring it is the same for each pair of points.

### 4.2.2 n-Simplex Method

Let $\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}$ be the vertices of
a regular D-simplex. The vertices have the following properties: $$
\lVert \mathbf{v}_i \rVert ^2 = 1 \quad \forall i=1,2, ..., D+1 \quad (\text{by definition})
$$ and

<center>

$\mathbf{v}_i ' \mathbf{v}_j = -\frac{1}{D} \quad \text{for} \ i\neq j$
(see Parks and Wills (2002))

</center>

<br> Let \$={Z}\_1, {Z}\_2,…,{Z}\_D \$ be a vector of D standard normal
random variables, then: $$\mathbf{v}_i'\mathbf{Z} \sim N(0,1)$$ and $$
\text{Cov}(\mathbf{v}_i'\mathbf{Z}, \mathbf{v}_j'\mathbf{Z})= -\frac{1}{D}
$$ The D-simplex average of $g$ is: $$
S(\mathbf{Z}) = \frac{1}{D+1} \sum_{i=1}^{D+1}g(\mathbf{v}_i'\mathbf{Z})
$$ $\therefore$ the Monte Carlo estimator of $E[S(\mathbf{Z})]$ is: $$
\frac{1}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)
$$ It is now shown that using this estimator will reduce variance from
the direct method: It is shown in theory that two monotone increasing
functions where the variables in each function are negatively correlated
has a covariance of less than or equal to zero. Hence, Hence, a variance
reduction is achieved with a higher dimesnions resulting in a higher
reduction. It follows that the case for which the dimension is one is
simply the antithetic variables method.

### 4.2.3 Implementation of Method to Price Asian Call Option

Consider the implementation of the n-simplex method to price Asian call
options. Let D be the number of dimensions, therefore D+1 points are
used to determine the simplex: <br> Remember, the following expected
value needs to be estimated for Asian call options: $$
\theta = E[\text{payoff}] = E\left[\left(\bar{S}  -K\right)^{+}\right]
$$ 1. Compute D+1 vertices of the regular D-simplex:
$\mathbf{v}_1, \mathbf{v}_2,...,\mathbf{v}_{D+1}$. 2. Generate D
vectors, $\mathbf{z} = \mathbf{z}_1,\mathbf{z}_2,...,\mathbf{z}_D$, each
containing m observations of standard normal random variables.
$$\mathbf{z}^{(1)}=z_1^{(1)}, z_2^{(1)}, ..., z_m^{(1)}
$$ $$\mathbf{z}^{(2)}=z_1^{(2)}, z_2^{(2)}, ..., z_m^{(2)}
$$

<center>

.<br> .<br> . <br>

</center>

$$
\mathbf{z}^{(D)}=z_1^{(D)}, z_2^{(D)}, ..., z_m^{(D)}
$$ 3. Generate D+1 stock price paths. Remember the model underlying the
stock price behaviour: $$
s_{t_j}^{(1)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}
$$ $$
s_{t_j}^{(2)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}
$$

<center>

.<br> .<br> . <br>

</center>

$$
s_{t_j}^{(D+1)} = s_{t_{j-1}}e^{\left((r-\frac{\sigma^2}{2})\triangle t + \sigma (x_{ij}) \sqrt{\triangle t}\right)}
$$ for $j = 1,...,m$ and $i = 1,2,...,D+1$, <br> where
$x_{ij}=[ \mathbf{z}_j^{(1)},\mathbf{z}_j^{(2)},...,\mathbf{z}_j^{(D)} ]' \mathbf{v}_i$<br>
This generates the stock price path for m time steps.<br> 4. Calculate
the payoff from each path (there are D+1 paths) by summing up the m
stock prices and dividing by m and then subtracting K from this result:
$$
\text{Payoff} = (\bar{S} - K)^+ \ \text{for D+1 paths}
$$ 5. Then calculate the average of the the D+1 payoffs: $$
S(\mathbf{Z}) = \frac{1}{D+1} \sum_{i=1}^{D+1}\text{payoff}_i
$$ 6. Repeat steps 2 to 5 $\frac{n}{D+1}$ times. <br> 7. Calculate the
mean of the $\frac{n}{D+1}$ average payoffs. <br> $$
\frac{1}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)
$$ 8. Discount this using the risk-free rate to obtain the price of the
Asian call option. <br> $$
\frac{e^{-rt}}{\frac{n}{D+1}}\sum_{j=1}^{\frac{n}{D+1}} S(\mathbf{Z}_j)
$$

# 5 Analyis

## 5.1 Results

The methods were compared using the following constants:

|  S  |  r  | $\sigma$ |  t  |  m  |
|:---:|:---:|:--------:|:---:|:---:|
| 100 | 5%  |   20%    |  1  | 250 |

where: <br> - S is the initial stock price.<br> - r is the risk-free
rate of interest. <br> - $\sigma$ is the volatility of the underlying
asset.<br> - t is the expiration date. <br> - m is the number of of time
steps.<br>

The price and standard error for each method was obtained and compared.
All methods were implemented using the algorithms that have been
described. The n-simplex method was done in 2,3 and 4 dimensions. For
the 2-Simplex, the points $\left(1,0\right)$,
$\left(-\frac{1}{2},\frac{\sqrt{3}}{2}\right)$ and
$\left(-\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)$ were used. For the
3-Simplex, the points $\left(1,0,0\right)$,
$\left(-\frac{1}{3}, \frac{\sqrt{5}}{3}, \frac{\sqrt{3}}{3}\right)$,
$\left(-\frac{1}{3}, \frac{-3-\sqrt{5}}{6},\frac{\sqrt{15}- \sqrt{3}}{6}\right)$
and
$\left(-\frac{1}{3}, \frac{3-\sqrt{5}}{6},\frac{-\sqrt{15}- \sqrt{3}}{6}\right)$.
For the 4-Simplex, the points $\left(1,0,0,0\right)$,
$\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)$,
$\left(-\frac{1}{4}, \frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)$,
$\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}\right)$,
$\left(-\frac{1}{4}, -\frac{\sqrt{5}}{4}, -\frac{\sqrt{5}}{4}, \frac{\sqrt{5}}{4}\right)$
were used.

The results are shown in the following table:

**Table 1: Price and Standard Error from the Respective Methods**

| K       | n      | Direct Price | Direct Standard Error | Antithetic Price | Antithetic Standard Error | 2-Simplex Price | 2-Simplex Standard Error | 3-Simplex Price | 3-Simplex Standard Error | 4-Simplex Price | 4-Simplex Standard Error |
|:--------|:-------|--------------|:----------------------|:-----------------|:--------------------------|:----------------|:-------------------------|:----------------|:-------------------------|:----------------|:-------------------------|
| **90**  | 1000   | 12.307079    | 0.3353537             | 12.6754041       | 0.0815546                 | 12.6555906      | 0.0608964                | 12.8263892      | 0.0534304                | 12.7834754      | 0.0430906                |
|         | 10000  | 12.8170867   | 0.1044421             | 12.6129053       | 0.0256425                 | 12.6268057      | 0.0187483                | 12.602466       | 0.0153161                | 12.6158688      | 0.0135036                |
|         | 100000 | 12.6541917   | 0.0331243             | 12.6225864       | 0.0081949                 | 12.6186997      | 0.0059988                | 12.619477       | 0.0049033                | 12.6149372      | 0.0042143                |
| **100** | 1000   | 5.5641811    | 0.2497639             | 5.8213534        | 0.1307974                 | 5.635155        | 0.0888653                | 5.9343749       | 0.0726026                | 5.860884        | 0.0636978                |
|         | 10000  | 5.8269498    | 0.0803125             | 5.8187024        | 0.0394167                 | 5.7521974       | 0.0274039                | 5.804665        | 0.0225513                | 5.7507879       | 0.0191825                |
|         | 100000 | 5.7862373    | 0.0253422             | 5.7916642        | 0.0124984                 | 5.7563086       | 0.0086521                | 5.8076507       | 0.0071744                | 5.7638815       | 0.0061304                |
| **110** | 1000   | 2.1593544    | 0.1729704             | 2.1153115        | 0.1037591                 | 2.3037734       | 0.0809331                | 1.9295741       | 0.0667774                | 2.1483104       | 0.0637792                |
|         | 10000  | 1.9945776    | 0.0484776             | 2.0688757        | 0.0316793                 | 2.0054124       | 0.0237079                | 1.96631         | 0.0195329                | 1.9926313       | 0.0171746                |
|         | 100000 | 1.98476      | 0.0154589             | 2.0084893        | 0.009999                  | 2.0146016       | 0.0075024                | 1.9771459       | 0.006164                 | 2.0128765       | 0.0054824                |

The variance reduction achieved from each method is summarised as
follows:

**Table 2: Variance Reduction Compared to the Direct Method**

|         |        | **Antithetic**         | **2-Simplex**          | **3-Simplex**          | **4-Simplex**          |
|---------|--------|:-----------------------|:-----------------------|:-----------------------|:-----------------------|
| **K**   | **n**  | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** |
| **90**  | 1000   | 0.9408588              | 0.9670256              | 0.9746154              | 0.9834896              |
|         | 10000  | 0.9397205              | 0.9677766              | 0.9784946              | 0.9832833              |
|         | 100000 | 0.9387946              | 0.9672034              | 0.9780878              | 0.9838132              |
| **100** | 1000   | 0.7257551              | 0.8734084              | 0.9155023              | 0.9349587              |
|         | 10000  | 0.7591237              | 0.8835717              | 0.9211545              | 0.9429515              |
|         | 100000 | 0.7567671              | 0.8834395              | 0.919853               | 0.9414827              |
| **110** | 1000   | 0.6401604              | 0.7810682              | 0.8509554              | 0.864039               |
|         | 10000  | 0.57296                | 0.7608309              | 0.8376499              | 0.8744868              |
|         | 100000 | 0.5816382              | 0.7644745              | 0.8410135              | 0.8742274              |

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
