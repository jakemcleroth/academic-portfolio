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

The payoff is given by (Zhang 2009): $$ `\phi`{=tex}(S) =
(S(T) - K)\^+ $$ for a European call option and by $$ `\phi`{=tex}(S) =
(K- S(T))\^+ $$ for a European put option, where <br> - ( T ) is the
expiration date,<br> - ( S(T) ) is the spot price of the underlying
asset at expiration T, <br> - ( K ) is the strike price of the option,
<br> - ( (x)^+ ) represents the maximum of ( x    0). <br>

## 2.2 Geometric Brownian Motion

To model the behaviour of stock prices, the most commonly used model is
geometric Brownian motion. It can be formulated as the following
stochastic differential equation (SDE) as seen in Hull (2017): $$
dS_t=`\mu `{=tex}S_t dt + `\sigma `{=tex}S_t d W_t $$ where:<br> - (dS)
is the instantaneous change in the stock price,<br> - (dW) is the Wiener
process ((W_t) is a random variable from the normal distribution
((0,t))), <br> - () represents the volatility of the stock price (() is
assumed to be constant over time), <br> - () represents the growth rate
of the stock price when there is no risk (also referred to as the
‘drift’), <br> - (dt) represents the deterministic return within the
time interval.<br>

with the discrete time version of the SDE as: $$ `\triangle `{=tex}S_t =
`\mu `{=tex}S_t `\triangle `{=tex}t + `\sigma `{=tex}S_t
`\triangle `{=tex}W_t $$

The solution to the discrete SDE is as follows: where (Z) is a standard
normal random variable.

## 2.3 Options Pricing

The main idea behind options pricing, or derivative pricing in general,
is known as risk neutral valuation. The idea is that the price of an
option should be equal to the present value of the expected payoff.
According to Hull (2017), this can be done in the following
steps:<br> 1. Assume expected rate of return is the risk free rate. <br>
2. Calculate the expected payoff of the option. For example, a European
call option: $$E\[`\text{payoff}`{=tex}$$ = E$$(S(T) - K)\^+$$\] 3.
Discount the expected payoff at the risk free rate. For the European
call option: $$ C =
e\^{-rt}E\[`\text{payoff}`{=tex}$$\] This European call option can be
priced using the Black-Scholes-Merton framework with a closed form
solution. However, there is no closed form solution for pricing
arithmetic Asian options. Hence we need to use simulation.

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
structure: $$ `\Phi`{=tex}(S) = `\left`{=tex}(`\bar`{=tex}{S} - K
`\right`{=tex})\^{+} `\quad `{=tex}`\text{for a call}`{=tex}$$ or
$$`\Phi`{=tex}(S) = `\left`{=tex}( K - `\bar`{=tex}{S}
`\right`{=tex})\^{+} `\quad `{=tex}`\text{for a put}`{=tex}$$ where <br>
$$`\text{where}`{=tex}  `\bar`{=tex}{S} =
`\frac{1}{m}`{=tex} `\sum`{=tex}\_{i=1}\^{m} S(t_i)$$

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
(*{a}^{b} g(x) , dx), provided the integral exists. If ( X ) is a random
variable with probability density function ( f(x) ), then the expected
value of the random variable ( Y = g(X) ) is given by: $$E\[g(X)$$ =
*{-}^{} g(x)f(x) , dx\] When a random sample is available from the
distribution of (X), an unbiased estimator of ( E$$g(X)$$ ) is simply
the sample mean (Rizzo 2019). <br> Suppose (f(x)) is a probability
density function, such that (f(x) ) for all (x ) and (*{A} f(x) , dx =
1). In order to estimate the integral (= *{A} g(x)f(x) , dx), generate a
random sample (x_1, x_2, …, x_n) from the distribution of (f(x)) and
calculate the sample mean
$$`\hat{\theta}`{=tex} = `\frac{1}{n}`{=tex} `\sum`{=tex}\_{j=1}\^{n}
g(x_j)$$

The core principle of Monte Carlo integration relies on the strong law
of large numbers: $$ `\hat{\theta}`{=tex} = `\frac{1}{n}`{=tex}
`\sum`{=tex}\_{j=1}\^{n} g(x_j)
`\rightarrow  `{=tex}E\[g(X)$$a.s. \] () as defined above is the Monte
Carlo estimator of (E$$g(X)$$) (Rizzo 2019).

## 3.2 Pricing Asian Options

### 3.2.1 Direct Method

We estimate (= E$$g(X)$$) by (*{DIR} = *{j=1}^{n} g(X_j)) where (X_1,
X_2, …, X_n) are independent and identically distributed. The direct
method refers to the basic Monte Carlo simulation without any variance
reduction techniques. Random samples are generated directly from the
distribution of interest.

In the case of pricing Asian call options, the direct method involves
generating random samples of the underlying asset’s price at various
future time points based on a stochastic model, such as geometric
Brownian motion. These simulated price paths are then used to calculate
the payoff of the option at the expiration date. The option’s price is
estimated by discounting back the average of the payoffs. <br> Recall
the payoff of a Asian call option is (Zhang 2009):
$$`\Phi`{=tex}(S) = `\left`{=tex}( `\bar`{=tex}{S} - K
`\right`{=tex})\^{+}$$

By risk-neutral valuation, the price of the Asian call option can be
expressed as:

Monte Carlo simulation can be applied to price Asian options.
Substituting (g) from above with the payoff of the Asian call option,
the following is obtained: $$ `\frac{1}{n}`{=tex}
`\sum`{=tex}\_{j=1}\^{n} `\left`{=tex}( `\bar`{=tex}{s}
-K`\right`{=tex})\^{+}\_j
`\rightarrow `{=tex}E`\left[\left( \bar{S} -K\right)^{+}\right] `{=tex}`\quad `{=tex}a.s.
$$ where:<br> - (( {s}
-K)^{+}\*j) is the payoff of the Asian call option from a simulated stock path,<br> - n is the number of Monte Carlo simulations. <br> Hence, the estimated price of the call option can be expressed as follows: $$`\hat{C}`{=tex}\ =\ `\frac{e^{-rt}}{n}`{=tex}\ `\sum`{=tex}*{j=1}^{n}
`\left`{=tex}( `\bar`{=tex}{s} -K`\right`{=tex})\^{+}\_j$$

### 3.2.2 Implementation of Method to Price Asian Call Option

To simulate stock paths, we make use of the geometric Brownian motion
explained earlier. The simulation goes as follows:<br> 1. Generate m
random variable observations, (z_1, z_2, …,z_m), from the standard
normal distribution. <br> 2. Simulate the path of the stock at m time
points using: $$ s_t =
s\_{t-1}e\^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}z_t `\sqrt{\triangle t}`{=tex}`\right`{=tex})} $$ 3.
Repeat n times to generate n paths.

10 simulated paths are shown below:<br>
![](Stoch-Sim-718-Project_files/figure-gfm/Asian%20Call%20Option%2010-1.png)<!-- -->
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
of the asset over the entire term ( \_{i=1}^{m} s(t_i) = {s}). <br> 3.
Subtract K from the average stock price over the entire term ((-K)).
Then the payoff of the Asian call option is max(({s}-K,0)). <br> 4.
Repeat steps 1-3 n times to obtain n different payoffs.<br> 5. Take the
average of all n payoffs and discount it with (e^{-rt}) where t is the
expiration date. <br> <br> The more paths that are generated, the more
accuracte the price will be. Shown below is an example of simulating
10000 stock paths.
![](Stoch-Sim-718-Project_files/figure-gfm/Asian%20Call%20Option%2010000-1.png)<!-- -->

### 3.2.3 Variance of Asian Call Option Price

$$ `\theta `{=tex}= E\[`\text{payoff}`{=tex}$$ = E\] $$
`\hat{\theta}`{=tex}*{DIR} = `\frac{1}{n}`{=tex} `\sum`{=tex}*{j=1}\^{n}
`\left`{=tex}( `\bar`{=tex}{s} -K`\right`{=tex})\^{+}\_j $$ () the
variance follows as

() the variance of the price of the Asian call options follows as:

# 4 Variance Reduction

Functions such as (E$$g(X)$$) can be estimated using Monte Carlo
integration. There are different methods that can be applied to reduce
the variance of the sample mean estimator of (= E$$g(X)$$) (Rizzo 2019).
Some of these methods include antithetic variables, control variables,
importance sampling, Latin hypercube sampling and Quasi-monte carlo
methods.

## 4.1 Antithetic Variables Method

### 4.1.1 Description

Suppose there are two identically distributed random variables (Y_1) and
(Y_2). If (Y_1) and (Y_2) are independent, then:
$$Var(`\frac{Y_1 + Y_2}{2}`{=tex}) = `\frac{1}{4}`{=tex}(Var(Y_1) +
Var(Y_2))$$ otherwise, in general, $$Var(`\frac{Y_1 + Y_2}{2}`{=tex}) =
`\frac{1}{4}`{=tex}(Var(Y_1) + Var(Y_2) + 2Cov(Y_1, Y_2))$$ From the
above, it can be seen that if (Y_1) and (Y_2) are negatively correlated,
then the variance of (( )) is smaller than when (Y_1) and (Y_2) are
independent (Rizzo 2019).

So, for the normal case, suppose that (Zhang 2009): $$
`\theta `{=tex}= E\[Y$$ = E$$g(X)$$  X N(0,1) \] where the direct Monte
Carlo estimate is: $$
`\hat{\theta}`{=tex}*{DIR} = `\frac{1}{n}`{=tex} `\sum`{=tex}*{i=1}\^{n}
g(X_i) `\quad `{=tex}`\text{where}`{=tex}  X_i
`\overset{\text{i.i.d.}}{\sim}`{=tex} N(0,1) $$ and the antithetic
variable estimate is: $$ `\hat{\theta}`{=tex}*{AV} =
`\frac{1}{\frac{n}{2}}`{=tex} `\sum`{=tex}*{i=1}\^{`\frac{n}{2}`{=tex}}
`\frac{g(X_i)+g(-X_i)}{2}`{=tex} `\quad `{=tex}`\text{where}`{=tex}  X_i
`\overset{\text{i.i.d.}}{\sim}`{=tex} N(0,1) $$ where (X_i) and (-X_i)
are the antithetic variates. The two antithetic variates are negatively
correlated. Therefore, if the function (g) is monotone increasing or
decreasing, variance reduction can be achieved by using this method
(Zhang 2009).

### 4.1.2 Implementation of Method to Price Asian Call Option

The following expected value needs to be estimated for Asian call
options: $$ `\theta `{=tex}= E\[`\text{payoff}`{=tex}$$ = E\] 1.
Generate m random standard normal observations ((z_1, z_2, … , z_m)).
<br> 2. Create another m standard normal observations that are
negatively correlated with (z_1, z_2,…,z_m) ((-z_1, -z_2, … , -z_m)).
<br> 3. Generate a stock price path using (z_1, z_2, … , z_m) and
another path using (-z_1, -z_2, … , -z_m). Remember the model underlying
the stock price behaviour: $$ s\_{t_i} =
s\_{t\_{i-1}}e\^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}z\_{t\_{i}} `\sqrt{\triangle t}`{=tex}`\right`{=tex})} $$
$$ s\_{t_i}\^{'} =
s\_{t\_{i-1}}^{'}e^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}(-z\_{t\_{i}}) `\sqrt{\triangle t}`{=tex}`\right`{=tex})}
$$ This will generate stock prices for m time steps. <br> 4. Calculate
the payoff from each path (there will be 2) by summing up the m stock
prices and dividing by m and then subtracting K from this result:
$$ `\text{Payoff}`{=tex} = `\left`{=tex}(`\bar`{=tex}{S} - K
`\right`{=tex})\^{+} `\quad `{=tex}`\text{for one path}`{=tex} $$ and $$
`\text{Payoff}`{=tex} =`\left`{=tex}( `\bar`{=tex}{S}' -
K`\right`{=tex})\^{+} `\quad `{=tex}`\text{for the other path}`{=tex} $$
5. Now use Monte Carlo simulation to repeat this process () times to
obtain n simulated paths representing n different outcomes and
calculating n different payoffs. $$
`\hat{\theta}`{=tex} = `\frac{1}{\frac{n}{2}}`{=tex}
`\sum`{=tex}*{j=1}\^{`\frac{n}{2}`{=tex}}
`\frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}`{=tex}
$$ 6. The price of a Asian call option is calculated by discounting the
payoff: $$ `\hat{C}`{=tex}*{AV} =
e\^{-rt}`\left[\frac{1}{\frac{n}{2}} \sum_{j=1}^{\frac{n}{2}} \frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}\right]`{=tex}$$

### 4.1.3 Variance of Asian Call Option Price

$$ `\theta `{=tex}= E\[`\text{payoff}`{=tex}$$ = E\] $$
`\hat{\theta}`{=tex}*{AV} = `\frac{1}{\frac{n}{2}}`{=tex}
`\sum`{=tex}*{j=1}\^{`\frac{n}{2}`{=tex}}
`\frac{\left( \bar{s} -K\right)^{+}_j+\left( \bar{s}^{'} -K\right)^{+}_j}{2}`{=tex}
$$ () the variance follows as (Rizzo 2019): \`\`{=tex} where Y is the
payoff using standard normal random variables ((Z_1, Z_2, … , Z_m)) and
Y’ is the payoff using standard normal random variables that are
negatively correlated to ((Z_1, Z_2, … , Z_m)) i.e. ((-Z_1, -Z_2, … ,
-Z_m)). Since the payoff function is monotone increasing, this method
will reduce variance.

() the variance of the price of the Asian call option is:

where Y and Y’ are defined as above.

## 4.2 n-Simplex Antithetic Variate Method

### 4.2.1 The Simplex

A simplex is an object in geometry that generalises the idea of a
triangle or tetrahedron to D-dimensions. In D dimensions, a simplex is
determined by D+1 points. These are called the vertices of the simplex.
A simplex in D dimensions is the convex hull of D+1 points (\_1,
*2,…,*{D+1}). In 1D, this would be a straight line, 2D a triangle, 3D a
tetrahedron and 4D a Pentatope. A regular D-simplex is a simplex where
the euclidean distance between each of the D+1 points is equal and the
angles between each face is the same. A regular D-simplex centred at the
origin of a unit D-dimensional hypersphere is of interest as it has
specific properties that will be shown later, which allow us to achieve
variance reduction. From here on, this will be referred to as the
regular D-simplex.The idea of the regular simplex is used to extend
antithetic variables to multiple dimensions for random variables with
symmetric probability distributions (Devriendt and Van Mieghem 2019). An
example of a regular 2-simplex is shown below.

![](Stoch-Sim-718-Project_files/figure-gfm/Unit%20Circle-1.png)<!-- -->
<br> In two dimensions, as above, a simplex is determined by three
points. These are points on a unit circle which correspond to the
vertices of an inscribed equilateral triangle. Any three points can be
chosen on the unit circle as long as these points form an equilateral
triangle. One such example is: ((1,0)), ((-,)) and ((-, -)) (Park and
Choe 2016). Confirm these points are correct by calculating the
Euclidean distance between each point and ensuring it is the same for
each pair of points.

### 4.2.2 n-Simplex Method

Let (\_1, *2,…,*{D+1}) be the vertices of a regular D-simplex. The
vertices have the following properties: $$
`\lVert `{=tex}`\mathbf{v}`{=tex}\_i `\rVert `{=tex}\^2 = 1
`\quad `{=tex}`\forall `{=tex}i=1,2, ..., D+1
`\quad `{=tex}(`\text{by definition}`{=tex}) $$ and

<center>

(\_i ’ \_j = -  ij) (see Parks and Wills (2002))

</center>

<br> Let (={Z}\_1, {Z}\_2,…,{Z}\_D ) be a vector of D standard normal
random variables, then:
$$`\mathbf{v}`{=tex}\_i'`\mathbf{Z}`{=tex} `\sim `{=tex}N(0,1)$$ and $$
`\text{Cov}`{=tex}(`\mathbf{v}`{=tex}\_i'`\mathbf{Z}`{=tex},
`\mathbf{v}`{=tex}*j'`\mathbf{Z}`{=tex})= -`\frac{1}{D}`{=tex} $$ The
D-simplex average of (g) is: $$ S(`\mathbf{Z}`{=tex}) =
`\frac{1}{D+1}`{=tex}
`\sum`{=tex}*{i=1}^{D+1}g(`\mathbf{v}`{=tex}*i'`\mathbf{Z}`{=tex})\ $$ () the Monte Carlo estimator of (E$$S(`\mathbf{Z}`{=tex})$$) is: $$\ `\frac{1}{\frac{n}{D+1}}`{=tex}`\sum`{=tex}*{j=1}^{`\frac{n}{D+1}`{=tex}}
S(`\mathbf{Z}`{=tex}\_j) $$ It is now shown that using this estimator
will reduce variance from the direct method: It is shown in theory that
two monotone increasing functions where the variables in each function
are negatively correlated has a covariance of less than or equal to
zero. Hence, Hence, a variance reduction is achieved with a higher
dimesnions resulting in a higher reduction. It follows that the case for
which the dimension is one is simply the antithetic variables method.

### 4.2.3 Implementation of Method to Price Asian Call Option

Consider the implementation of the n-simplex method to price Asian call
options. Let D be the number of dimensions, therefore D+1 points are
used to determine the simplex: <br> Remember, the following expected
value needs to be estimated for Asian call options: $$
`\theta `{=tex}= E\[`\text{payoff}`{=tex}$$ = E\] 1. Compute D+1
vertices of the regular D-simplex: (\_1, *2,…,*{D+1}). 2. Generate D
vectors, ( = \_1,\_2,…,\_D), each containing m observations of standard
normal random variables.
$$`\mathbf{z}`{=tex}^{(1)}=z_1^{(1)}, z_2\^{(1)}, ..., z_m\^{(1)} $$
$$`\mathbf{z}`{=tex}^{(2)}=z_1^{(2)}, z_2\^{(2)}, ..., z_m\^{(2)} $$

<center>

.<br> .<br> . <br>

</center>

$$ `\mathbf{z}`{=tex}^{(D)}=z_1^{(D)}, z_2\^{(D)}, ..., z_m\^{(D)} $$ 3.
Generate D+1 stock price paths. Remember the model underlying the stock
price behaviour: $$ s\_{t_j}\^{(1)} =
s\_{t\_{j-1}}e\^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}(x\_{ij}) `\sqrt{\triangle t}`{=tex}`\right`{=tex})} $$
$$ s\_{t_j}\^{(2)} =
s\_{t\_{j-1}}e\^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}(x\_{ij}) `\sqrt{\triangle t}`{=tex}`\right`{=tex})} $$

<center>

.<br> .<br> . <br>

</center>

$$ s\_{t_j}\^{(D+1)} =
s\_{t\_{j-1}}e\^{`\left`{=tex}((r-`\frac{\sigma^2}{2}`{=tex})`\triangle `{=tex}t +
`\sigma `{=tex}(x\_{ij}) `\sqrt{\triangle t}`{=tex}`\right`{=tex})} $$
for (j = 1,…,m) and (i = 1,2,…,D+1), <br> where (x\_{ij}=$$
`\mathbf{z}`{=tex}\_j\^{(1)},`\mathbf{z}`{=tex}\_j\^{(2)},...,`\mathbf{z}`{=tex}\_j\^{(D)}
$$’ \*i)<br> This generates the stock price path for m time steps.<br>
4. Calculate the payoff from each path (there are D+1 paths) by summing
up the m stock prices and dividing by m and then subtracting K from this
result: $$ `\text{Payoff}`{=tex} =
(`\bar`{=tex}{S} - K)\^+  `\text{for D+1 paths}`{=tex} $$ 5. Then
calculate the average of the the D+1 payoffs: $$ S(`\mathbf{Z}`{=tex}) =
`\frac{1}{D+1}`{=tex}
`\sum`{=tex}*{i=1}^{D+1}`\text{payoff}`{=tex}*i\ $$ 6. Repeat steps 2 to 5 () times. <br> 7. Calculate the mean of the () average payoffs. <br> $$\ `\frac{1}{\frac{n}{D+1}}`{=tex}`\sum`{=tex}*{j=1}^{`\frac{n}{D+1}`{=tex}}
S(`\mathbf{Z}`{=tex}*j) $$ 8. Discount this using the risk-free rate to
obtain the price of the Asian call option. <br> $$
`\frac{e^{-rt}}{\frac{n}{D+1}}`{=tex}`\sum`{=tex}*{j=1}\^{`\frac{n}{D+1}`{=tex}}
S(`\mathbf{Z}`{=tex}\_j) $$

# 5 Analyis

## 5.1 Results

The methods were compared using the following constants:

|  S  |  r  | ()  |  t  |  m  |
|:---:|:---:|:---:|:---:|:---:|
| 100 | 5%  | 20% |  1  | 250 |

where: <br> - S is the initial stock price.<br> - r is the risk-free
rate of interest. <br> - () is the volatility of the underlying
asset.<br> - t is the expiration date. <br> - m is the number of of time
steps.<br>

The price and standard error for each method was obtained and compared.
All methods were implemented using the algorithms that have been
described. The n-simplex method was done in 2,3 and 4 dimensions. For
the 2-Simplex, the points ((1,0)), ((-,)) and ((-, -)) were used. For
the 3-Simplex, the points ((1,0,0)), ((-, , )), ((-, ,)) and ((-, ,)).
For the 4-Simplex, the points ((1,0,0,0)), ((-, , , )), ((-, , , -)),
((-, -, , -)), ((-, -, -, )) were used.

The results are shown in the following table:

**Table 1: Price and Standard Error from the Respective Methods**

| K       | n      | Direct Price | Direct Standard Error | Antithetic Price | Antithetic Standard Error | 2-Simplex Price | 2-Simplex Standard Error | 3-Simplex Price | 3-Simplex Standard Error | 4-Simplex Price | 4-Simplex Standard Error |
|:--------|:-------|--------------|:----------------------|:-----------------|:--------------------------|:----------------|:-------------------------|:----------------|:-------------------------|:----------------|:-------------------------|
| **90**  | 1000   | 12.269133    | 0.3091909             | 12.6079316       | 0.0850343                 | 12.6408824      | 0.0622441                | 12.5803482      | 0.0523175                | 12.6601075      | 0.0371124                |
|         | 10000  | 12.5021309   | 0.1038632             | 12.6900288       | 0.0278656                 | 12.5949879      | 0.0184647                | 12.6427349      | 0.0158526                | 12.6517939      | 0.0138116                |
|         | 100000 | 12.5584617   | 0.0329272             | 12.6094516       | 0.008136                  | 12.6141421      | 0.0059168                | 12.6023288      | 0.0048241                | 12.6153872      | 0.0042582                |
| **100** | 1000   | 5.8634422    | 0.251366              | 5.5131091        | 0.1196008                 | 5.7965957       | 0.0935846                | 5.7417604       | 0.0689959                | 5.8748634       | 0.0621291                |
|         | 10000  | 5.8225393    | 0.0801917             | 5.6431412        | 0.038279                  | 5.7966381       | 0.0272829                | 5.7667744       | 0.0222191                | 5.7800728       | 0.0189877                |
|         | 100000 | 5.7155697    | 0.025128              | 5.7959844        | 0.0124327                 | 5.7998836       | 0.0087769                | 5.7855059       | 0.0071562                | 5.7634189       | 0.0061283                |
| **110** | 1000   | 1.8188844    | 0.1454878             | 2.0139693        | 0.1017538                 | 1.9145937       | 0.0772157                | 2.0004032       | 0.0582981                | 1.8452674       | 0.0499755                |
|         | 10000  | 2.0543881    | 0.0503011             | 2.0426364        | 0.0323151                 | 1.9884411       | 0.0234972                | 1.9871272       | 0.0196743                | 1.9676632       | 0.0167926                |
|         | 100000 | 2.0049696    | 0.0155143             | 2.0126702        | 0.01004                   | 2.010395        | 0.0074904                | 1.9938864       | 0.0062662                | 2.0146504       | 0.0054838                |

The variance reduction achieved from each method is summarised as
follows:

**Table 2: Variance Reduction Compared to the Direct Method**

|         |        | **Antithetic**         | **2-Simplex**          | **3-Simplex**          | **4-Simplex**          |
|---------|--------|:-----------------------|:-----------------------|:-----------------------|:-----------------------|
| **K**   | **n**  | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** | **Variance Reduction** |
| **90**  | 1000   | 0.924363               | 0.9594731              | 0.9713688              | 0.9855926              |
|         | 10000  | 0.9280197              | 0.9683946              | 0.9767042              | 0.9823166              |
|         | 100000 | 0.9389466              | 0.9677098              | 0.9785356              | 0.983276               |
| **100** | 1000   | 0.7736113              | 0.8613896              | 0.9246587              | 0.938909               |
|         | 10000  | 0.7721433              | 0.8842499              | 0.9232298              | 0.9439356              |
|         | 100000 | 0.7552003              | 0.8779984              | 0.9188945              | 0.9405219              |
| **110** | 1000   | 0.5108437              | 0.7183185              | 0.8394331              | 0.8820059              |
|         | 10000  | 0.5872786              | 0.7817886              | 0.8470168              | 0.8885501              |
|         | 100000 | 0.5812037              | 0.7668986              | 0.8368676              | 0.875059               |

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
