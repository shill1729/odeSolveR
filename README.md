
# odeSolveR

<!-- badges: start -->
<!-- badges: end -->

A basic ODE solver using **[Euler scheme](https://en.wikipedia.org/wiki/Euler_method)** and fourth-order **[Runge-Kutta scheme](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)** implementations to numerically solve arbitrary first order one dimensional ODEs and systems of first-order ODEs. The default scheme used is the latter RK4 scheme. Also included is the basic SIR model and an extended SIR model called the SIDARTHE model from mathematical epidemiology. Please see the original  **[paper](https://arxiv.org/abs/2003.09861)** by the COVID19 IRCCS San Matteo Pavia Task Force et al. for mathematical details and discussion of the SIDARTHE model that is numerically solved via Euler schemes in this simple ODE-solver.

The documentation should be sufficient for setting up custom problems, however, we also provide below a series of examples from the 1970s to the 2000s demonstrating how to use the generic solver. The examples are pulled from S. Strogatz's excellent *Nonlinear Dynamics and Chaos* (2015) with the original authors listed.

# Table of contents
1. [Installation](#installation)

2. [One dimensional ODEs](#one-dimensional-odes)

3. [The SIR model](#the-sir-model)

4. [Arbitrary first-order ODE systems](#arbitrary-first-order-ode-systems)

    a. [Replicating the hard-coded SIR example](#replicating-the-hard-coded-sir-example)
    
    b. [The Lorenz system](#the-lorenz-system)
    
    c. [Vasquez and Redner political opinion dynamics](#vasquez-and-redner-political-opinion-dynamics)
    
    d. [Jordan and Smith Keynesian Cross model](#jordan-and-smith-keynesian-cross-model)
    
    e. [Eigen and Schuster Hypercycle equation](#eigen-and-schuster-hypercycle-equation)
    
    f. [Haken two-mode laser model](#haken-two-mode-laser-model)
    
    g. [Reversible System example](#reversible-system-example)

5. [SIDARTHE examples](#sidarthe-examples)

## Installation

You can install the latest version of odeSolveR on GitHub with:

``` r
devtools::install_github("shill1729/odeSolveR")
```

## One dimensional ODEs
```r
library(odeSolveR)
# The ODE is given in the form y'(t)=f(t, y) with initial y(t_0)=y_0.
# Here we solve y'=cos(t), y(0)=0. The solution is y=sin(t), of course.
f <- function(t, y) cos(t)
y <- ode1(f = f, y0 = 0, tn = 10, n = 5000)
plot(y)
```

## The SIR model

Running an Euler scheme on the SIR model from mathematical epidemiology.

``` r
library(odeSolveR)
# Solution to SIR model with 2 days between contacts, 3 days between recovery
# 90 days out, 
sir.dat <- sir(beta = 1/2, gamma = 1/3, s0 = 50, i0 = 1, tn = 90, n = 1000)
```

## Arbitrary first-order ODE systems

We can solve arbitrary systems using special syntax. As of now, no error handling is implemented so systems with irregularities, singularities, or blow ups will not be handled with any grace...

### Replicating the hard-coded SIR example

This will replicate the above SIR example (without all the plots...) that is implemented as a hard-coded ODE.

```r
library(odeSolveR)
# Numerical parameters
t0 <- 0
tn <- 90
n <- 1000

# Model parameters
parameters <- list(beta = 1/2, gamma = 1/3)

# TODO: need a less hacky fix
list2env(parameters, envir = .GlobalEnv)

# The functions describing rate of change of each state variable in the ODE system
f <- list(function(t, s, i, r) -beta*s*i, function(t, s, i, r) beta*i*s-gamma*i, function(t, s, i, r) gamma*i)

# Initial conditions
IC <- list(s = 1-0.001, i = 0.001, r = 0)
# Solve the system via RK4 iterations
sol <- ode(f, IC, parameters, tn = 90)

# Plotting trajectories
par(mfrow = c(1, 1))
plot(sol$time, sol$s, type = "l", main = "Fraction of populations", ylim = c(0, 1))
lines(sol$time, sol$i, col = "red")
lines(sol$time, sol$r, col = "blue")
legend(x = "topright", legend = c("susceptible", "infected", "recovered"), col = c("black", "red", "blue"), lty = 1)
```
![SIR](examplePlots/sir.jpeg)

### The Lorenz system

As another example of entering in our own ODEs, we can solve the Lorenz system.

```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 40
n <- 10000

# Model parameters
parameters <- list(sigma = 10, rho = 28, beta = 8/3)

# TODO: need a less hacky fix to load in parameters so they do not be referenced as parameters$parm in function bodies
list2env(parameters, envir = .GlobalEnv)

# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, x, y, z) sigma*(y-x),
          function(t, x, y, z) x*(rho-z)-y,
          function(t, x, y, z) x*y-beta*z
)

# Initial conditions: names and order must match variables appearing after "t" in function definitions above!
IC <- list(x = 1, y = 1, z = 1)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

# Plot 2 dimension phase space of x-z variables
plot(sol$x, sol$z, type = "l", main = "x-z Phase plane")

# Plot 3 dimensional phase portrait
plot3D::lines3D(x = sol$x, y = sol$y, z = sol$z, main = "3 dimensional phase portrait")
```
![Lorenz](examplePlots/lorenz.jpeg)


### Vasquez and Redner political opinion dynamics

A highly simplified model of political opinion dynamics due to Vasquez and Redner (2004 p. 8489) consisting of a population of rightists (x), leftists (y), and centrists (z). The dynamics folows the assumptions that leftists and rightists never talk to each other but when either talks to a centrist, one of them convinces the other to change his or her political opinion with the winner depending on the sign of a parameter (r). If this parameter is positive, extremists always win and persuade the centrist to move to that end of the spectrum. Otherwise, the centrist pulls the others towards the center.
```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 90
n <- 10000

# Model parameters
parameters <- list(r = -0.1)

# TODO: need a less hacky fix to load in parameters so they do not be referenced as parameters$parm in function bodies
list2env(parameters, envir = .GlobalEnv)

# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, x, y, z) r*x*z,
          function(t, x, y, z) r*y*z,
          function(t, x, y, z) -r*x*z-r*y*z
)

# Initial conditions
IC <- list(x = 0.5, y = 0.4, z = 0.10)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

par(mfrow = c(2, 1))
# State trajectories over time
plot(sol$time, sol$x, type = "l", col = "red", ylim = c(0, 1), main = "% of parties over time")
lines(sol$time, sol$y, col = "blue")
lines(sol$time, sol$z, col = "black")

# 3 Dimensional phase portrait
plot3D::lines3D(x = sol$x, y = sol$y, z = sol$z, theta = 120, main = "Phase portrait")
```
![Politics](examplePlots/politics.jpeg)

### Jordan and Smith Keynesian Cross model

Adapted from Exercise 2.24 in Jordan and Smith (1987). A simple model for a national economy is given by three variables: income (I), rate of consumer spending (C) and rate of government spending (G), with two parameters alpha and beta both greater than unity. The parameters here correspond to a semi-oscillating trajectory that takes a long time to reach equillibrium.

```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 300
n <- 5000

# Model parameters
parameters <- list(alpha = 1.01, beta = 1.01, k = 0.2)

# TODO: need a less hacky fix to load in parameters so they do not be referenced as parameters$parm in function bodies
list2env(parameters, envir = .GlobalEnv)

# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, I, C, G) I-alpha*C,
          function(t, I, C, G) beta*(I-C-G),
          function(t, I, C, G) 0 # Or use k*(I-alpha*C) for G=G0+k*I linear growth of expenditures in terms of income
          # or use 2*k*I*(I-alpha*C) for quadratic growth of expenditures in terms of income
)

# Initial conditions
IC <- list(I = 300, C = 100, G = 300*k)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

par(mfrow = c(2, 1))
# State trajectories over time
plot(sol$time, sol$I, type = "l", col = "blue")
lines(sol$time, sol$C, col = "red")
lines(sol$time, sol$G, col = "green")

# Plot 2 dimension phase space of x-z variables
plot(sol$I, sol$C, type = "l")
```
![Economy](examplePlots/economy.jpeg)

### Eigen and Schuster Hypercycle equation

In Eigen and Schuster's (1978) model of pre-biotic evolution, a group of RNA molecules or other self-replicating chemical units are imagined to catalyze each other's reproduction in a closed feedback loop, with one molecule serving as the catalyst for the next. The simplest scheme considered by Eigen and Schuster's is implemented here, in dimensionless form as the hypercycle equations, for three molecules. For the dynamics of the hypercycle equation under larger values of the number of molecules we refer interested readers to Hofbauer and Sigmund (1998).

```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 45
n <- 5000

parameters <- list()

# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, x1, x2, x3) x1*(x3-x1*x3-x2*x1-x3*x2),
          function(t, x1, x2, x3) x2*(x1-x1*x3-x2*x1-x3*x2),
          function(t, x1, x2, x3) x3*(x2-x1*x3-x2*x1-x3*x2)
          
)

# Initial conditions
IC <- list(x1 = 0.3, x2 = 0.6, x3 = 0.1)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

par(mfrow = c(2, 1))
# State trajectories over time
plot(sol$time, sol$x1, type = "l", col = "blue", ylim = c(0, 1), main = "Relative frequency of molecules")
lines(sol$time, sol$x2, col = "red")
lines(sol$time, sol$x3, col = "black")


# 3 Dimensional phase portrait
plot3D::lines3D(x = sol$x1, y = sol$x2, z = sol$x3, main = "Phase portrait")
```
![Hypercycle](examplePlots/hypercycle.jpeg)

### Haken two-mode laser model

According to Haken (1983, p. 129), a two-mode laser produces two different kinds of photons. The number of each kind of photon is modeled via a first order ODE system, implemented below. The rate of change of the number of photons takes a basic "gain-loss" form dependent on the number of excited photons.

```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 45
n <- 5000

# Model parameters
parameters <- list(G1 = 0.5, G2 = 0.4, k1 = 0.3, k2 = 0.2, alpha1 = 0.9, alpha2 = 0.02, N0 = 100)


# TODO: need a less hacky fix to load in parameters so they do not be referenced as parameters$parm in function bodies
list2env(parameters, envir = .GlobalEnv)


# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, n1, n2) G1*(N0-alpha1*n1-alpha2*n2)-k1*n1,
          function(t, n1, n2) G2*(N0-alpha1*n1-alpha2*n2)-k2*n2
)

# Initial conditions
IC <- list(n1 = 0, n2 = 0)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

par(mfrow = c(2, 1))
# State trajectories over time
plot(sol$time, sol$n1, type = "l", col = "blue", ylim = c(min(sol$n1, sol$n2), max(sol$n1, sol$n2)), main = "Number of photons")
lines(sol$time, sol$n2, col = "red")

# 2 Dimensional phase portrait
plot(sol$n1, sol$n2, main = "Phase portrait", type = "l")
```
![Laser](examplePlots/laser.jpeg)

### Reversible System example

Here is an example of a reversible system sometimes referred to as "wallpaper" (Strogatz 2015).

```r
library(odeSolveR)

# Numerical parameters
t0 <- 0
tn <- 100
n <- 5000

parameters <- list()

# Create the functions that describe the rates of change for each state variable in the ODE system
f <- list(function(t, x, y) sin(y),
          function(t, x, y) sin(x)
)

# Initial conditions
IC <- list(x = 1, y = 0.01)
# Solve the system using RK4
sol <- ode(f, IC, parameters, tn = tn, n = n)

par(mfrow = c(2, 1))
# State trajectories over time
plot(sol$time, sol$x, type = "l", col = "blue", ylim = c(min(sol$x, sol$y), max(sol$x, sol$y)), main = "Trajectory over time")
lines(sol$time, sol$y, col = "red")

# 2 Dimensional phase portrait
plot(sol$x, sol$y, main = "Phase portrait", type = "l")
```
![Wallpaper](examplePlots/wallpaper.jpeg)

## SIDARTHE examples

Running an Euler scheme to solve the SIDARTHE model developed by the COVID19 IRCCS San Matteo Pavia TaskForce et al.

### Parameters from original paper, solved using an Euler scheme

```r
library(odeSolveR)
parameters <- list(alpha = 0.57,
                beta = 0.011,
                delta = 0.011,
                gamma = 0.456,
                epsilon = 0.171,
                theta = 0.371,
                zeta = 0.125,
                eta = 0.125,
                mu = 0.012,
                nu = 0.027,
                tau = 0.003,
                lambda = 0.034,
                rho = 0.034, 
                kappa = 0.017,
                sigma = 0.017,
                xi = 0.017)
initial_conditions <- list(i0 = 200/60000000, 
                           d0 = 20/60000000, 
                           a0 = 1/60000000, 
                           r0 = 2/60000000, 
                           thr0 = 0,
                           h0 = 0,
                           e0 = 0)
initial_conditions$s0 <- 1-sum(unlist(initial_conditions))
sir.dat <- sidarthe(parameters, initial_conditions, tn = 300)
```
### SIDARTHE Example 2: manually fitted model parameters to US data, solved using RK4

Below contains an ad-hoc fit to US state-wide data available from the **[NY Times covid-19 dataset](https://github.com/nytimes/covid-19-data)** on GitHub. A model is produced with all 8 of the state variables graphed plus three additional graphs fitting two of the model curves to data and visualizing the death-rate.

```r
library(odeSolveR)
# Choose a terminal date to forecast until
terminalDate <- "2020-08-30"
pop <- 329.45*10^6 # US 2019 population
# Reading NY Times github csv of cumulative cases: replace 'filepath' with the path to the downloaded csv
# Something like '.../Downloads/covid-19-data-master/covid-19-data-master/us-states.csv'
dat <- read.csv(file = filepath)

# Total US cases and deaths
x <- aggregate(x = dat$cases, by = list(dat$date), FUN = sum)$x 
y <- aggregate(x = dat$deaths, by = list(dat$date), FUN = sum)$x
us_dat <- data.frame(date = as.Date(unique(dat$date)), cases = x, deaths = y) 

# Length of data-set in days
N <- length(us_dat$date)

# Forecast dates
extended_dates <- seq(as.Date(us_dat$date[N]), as.Date(terminalDate), by = 1)
extended_dates <- extended_dates[-1] # Remove last day with data
# Number of days to forecast out
Nf <- length(extended_dates)
# Total length of time of data-set + forecast period
tn <- as.numeric(extended_dates[Nf]-us_dat$date[1])
# Number of nodes: ensures 0.1 step
n <- tn*10

# SIDARTHE model parameters
parameters <- list(alpha = 0.248,
                   beta = 0.05,
                   delta = 0.05,
                   gamma = 0.2,
                   epsilon = 0.0019,
                   theta = 0.411,
                   zeta = 0.125,
                   eta = 0.125,
                   mu = 0.1,
                   nu = 0.1,
                   tau = 0.00015,
                   lambda = 0.01,
                   rho = 0.01, 
                   kappa = 0.02,
                   sigma = 0.06,
                   xi = 0.05)
list2env(x = parameters, envir = .GlobalEnv)
# Initial conditions: S is set to 1-sum(...) afterwards, do not let the 0/pop throw you off...
initial_conditions <- list(S = 0/pop,
                           I = 450/pop, 
                           D = 1/pop, 
                           A = 1/pop, 
                           R = 1/pop, 
                           Th = 0,
                           H = 0,
                           E = 0)
initial_conditions$S <- 1-sum(unlist(initial_conditions))

# Computing the basic reproduction ratio
r1 <- epsilon+zeta+lambda
r2 <- eta + rho
r3 <- theta + mu + kappa
r4 <- nu + xi
repo <- alpha/r1 + beta*epsilon/(r1*r2) + gamma*zeta/(r1*r3) + delta*eta*epsilon/(r1*r2*r4) + delta*zeta*theta/(r1*r3*r4)

# Setting the ODE functions
f <- list(function(t, S, I, D, A, R, Th, H, E) -S*(alpha*I+beta*D+gamma*A+delta*R),
          function(t, S, I, D, A, R, Th, H, E) S*(alpha*I+beta*D +gamma*A +delta*R)-(epsilon+zeta+lambda)*I,
          function(t, S, I, D, A, R, Th, H, E) epsilon*I-(eta+rho)*D,
          function(t, S, I, D, A, R, Th, H, E) zeta*I-(theta+mu+kappa)*A,
          function(t, S, I, D, A, R, Th, H, E) eta*D+theta*A-(nu+xi)*R,
          function(t, S, I, D, A, R, Th, H, E) mu*A+nu*R-(sigma+tau)*Th,
          function(t, S, I, D, A, R, Th, H, E) lambda*I+rho*D+kappa*A+xi*R+sigma*Th,
          function(t, S, I, D, A, R, Th, H, E) tau*Th
)
# Solving the ODE system with RK4
sol <- ode(f = f, IC = initial_conditions, parameters = parameters, t0 = 0, tn = tn, n = n)


# Plotting model output and fits to data
par(mfrow = c(2, 2))
plot(sol$time, sol$S, type = "l", main = "Fractions of population", ylim = c(0, 1), xlab = "days", ylab = "% of population")
for(i in 3:9)
{
  lines(sol[, 1], sol[, i], col = i-1)
}
legend(x = "topright", 
       legend = c("Susceptible", "Infected", "Diagnosed", "Ailing", "Recognized", "Threatened", "Healed", "Extinct"), 
       lty = 1, col = c(1:10), cex = 0.5)

# Daily cases
plot(c(us_dat$date, extended_dates), c(us_dat$cases, rep(NA, Nf)), ylim = c(0, max(us_dat$cases, sol$D*pop)), xlab = "Date", ylab = "# of cases", main = "Daily cumulative US cases")
lines(c(us_dat$date, extended_dates), sol$D[sol$time%%1==0]*pop, col = "green")

# Daily deaths
plot(c(us_dat$date, extended_dates), c(us_dat$deaths, rep(NA, Nf)), ylim = c(0, max(us_dat$deaths, sol$E*pop)), xlab = "Date", ylab = "# of deaths", main = "Daily cumulative US deaths")
lines(c(us_dat$date, extended_dates), sol$E[sol$time%%1==0]*pop, col = "grey")

# Death-rate
plot(us_dat$date, us_dat$deaths/us_dat$cases, main = "Death rate")

# Printing the reproductive ratio
print(paste("R0 = ", repo))
```
### SIDARTHE Example 3: manually fitted model parameters to New York state data, solved using RK4
...and another example for just a single state, e.g., New York.

```r
library(odeSolveR)
# Choose a terminal date to forecast until
terminalDate <- "2020-08-30"
pop <- 19.5*10^6 # NY 2019 population
# Reading NY Times github csv of cumulative cases: replace 'filepath' with the path to the downloaded csv
# Something like '.../Downloads/covid-19-data-master/covid-19-data-master/us-states.csv'
dat <- read.csv(file = filepath)

# Get out just New York cases
ny_dat <- dat[which(dat$state == "New York"),]
ny_dat$date <- as.Date(ny_dat$date)

# Length of data-set in days
N <- length(ny_dat$date)

# Forecast dates
extended_dates <- seq(as.Date(ny_dat$date[N]), as.Date(terminalDate), by = 1)
extended_dates <- extended_dates[-1] # Remove last day with data
# Number of days to forecast out
Nf <- length(extended_dates)
# Total length of time of data-set + forecast period
tn <- as.numeric(extended_dates[Nf]-ny_dat$date[1])
# Number of nodes: ensures 0.1 step
n <- tn*10

# SIDARTHE model parameters
parameters <- list(alpha = 0.248,
                   beta = 0.05,
                   delta = 0.05,
                   gamma = 0.2,
                   epsilon = 0.0115,
                   theta = 0.5,
                   zeta = 0.125,
                   eta = 0.125,
                   mu = 0.1,
                   nu = 0.1,
                   tau = 0.00063,
                   lambda = 0.009,
                   rho = 0.005, 
                   kappa = 0.02,
                   sigma = 0.06,
                   xi = 0.05)
# Load in parameters. A bit hacky...
list2env(x = parameters, envir = .GlobalEnv)

# Initial conditions: S is set to 1-sum(...) afterwards, do not let the 0/pop throw you off...
initial_conditions <- list(S = 0/pop,
                           I = 22000/pop, 
                           D = 1/pop, 
                           A = 1/pop, 
                           R = 1/pop, 
                           Th = 0,
                           H = 0,
                           E = 0)
initial_conditions$S <- 1-sum(unlist(initial_conditions))

# Computing the basic reproduction ratio
r1 <- epsilon+zeta+lambda
r2 <- eta + rho
r3 <- theta + mu + kappa
r4 <- nu + xi
repo <- alpha/r1 + beta*epsilon/(r1*r2) + gamma*zeta/(r1*r3) + delta*eta*epsilon/(r1*r2*r4) + delta*zeta*theta/(r1*r3*r4)

# Setting the ODE functions
f <- list(function(t, S, I, D, A, R, Th, H, E) -S*(alpha*I+beta*D+gamma*A+delta*R),
          function(t, S, I, D, A, R, Th, H, E) S*(alpha*I+beta*D +gamma*A +delta*R)-(epsilon+zeta+lambda)*I,
          function(t, S, I, D, A, R, Th, H, E) epsilon*I-(eta+rho)*D,
          function(t, S, I, D, A, R, Th, H, E) zeta*I-(theta+mu+kappa)*A,
          function(t, S, I, D, A, R, Th, H, E) eta*D+theta*A-(nu+xi)*R,
          function(t, S, I, D, A, R, Th, H, E) mu*A+nu*R-(sigma+tau)*Th,
          function(t, S, I, D, A, R, Th, H, E) lambda*I+rho*D+kappa*A+xi*R+sigma*Th,
          function(t, S, I, D, A, R, Th, H, E) tau*Th
)
# Solving the ODE system with RK4
sol <- ode(f = f, IC = initial_conditions, parameters = parameters, t0 = 0, tn = tn, n = n)


# Plotting model output and fits to data
par(mfrow = c(2, 2))
plot(sol$time, sol$S, type = "l", main = "Fractions of population", ylim = c(0, 1), xlab = "days", ylab = "% of population")
for(i in 3:9)
{
  lines(sol[, 1], sol[, i], col = i-1)
}
legend(x = "topright", 
       legend = c("Susceptible", "Infected", "Diagnosed", "Ailing", "Recognized", "Threatened", "Healed", "Extinct"), 
       lty = 1, col = c(1:10), cex = 0.5)

# Daily cases
plot(c(ny_dat$date, extended_dates), c(ny_dat$cases, rep(NA, Nf)), ylim = c(0, max(ny_dat$cases, sol$D*pop)), xlab = "Date", ylab = "# of cases", main = "Daily cumulative NY cases")
lines(c(ny_dat$date, extended_dates), sol$D[sol$time%%1==0]*pop, col = "green")

# Daily deaths
plot(c(ny_dat$date, extended_dates), c(ny_dat$deaths, rep(NA, Nf)), ylim = c(0, max(ny_dat$deaths, sol$E*pop)), xlab = "Date", ylab = "# of deaths", main = "Daily cumulative NY deaths")
lines(c(ny_dat$date, extended_dates), sol$E[sol$time%%1==0]*pop, col = "grey")

# Death-rate
plot(ny_dat$date, ny_dat$deaths/ny_dat$cases, main = "Death rate")

# Printing the reproductive ratio
print(paste("R0 = ", repo))
```

