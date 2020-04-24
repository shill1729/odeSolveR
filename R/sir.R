#' Euler scheme for the SIR model
#'
#' @param beta the mean rate of contact, 1/beta the mean time between contact for a given individual (homogeneous population and mixing)
#' @param gamma the mean rate of recovery, 1/gamma the mean time until recovery.
#' @param s0 the initial susceptible population
#' @param i0 the initial infected population
#' @param r0 the initial recovering population, defaults to zer0
#' @param t0 initial time point
#' @param tn the terminal time point
#' @param n number of sub-intervals in time-grid
#' @param verbose whether to plot the populations over time and print parameter+reproductive ratio
#'
#' @description {A basic Euler scheme is used to numerically solve the ODE system of the standard SIR model.}
#' @return list
#' \itemize{
#' \item parameters, a data.frame containing the contact rate, recovery rate, and r0 the basic reproduction ratio.
#' \item sir, a matrix of the three populations in columns over time, susceptible, infected, recovered}
#' @export sir
#' @examples {sir.dat <- sir(beta = 1/2, gamma = 1/3, s0 = 50, i0 = 1, tn = 90, n = 1000)}
sir <- function(beta, gamma, s0 = 1, i0 = 0.0001, r0 = 0, t0 = 0, tn = 1, n = 1000, verbose = TRUE)
{
  h <- (tn-t0)/n
  s <- matrix(0, n+1)
  I <- matrix(0, n+1)
  r <- matrix(0, n+1)
  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  s[1] <- s0
  I[1] <- i0
  r[1] <- r0
  N <- s[1] + I[1] +r[1]
  for(i in 2:(n+1))
  {
    s[i] <- s[i-1] - beta*s[i-1]*I[i-1]*h/N
    I[i] <- I[i-1] +(beta*I[i-1]*s[i-1]/N -gamma*I[i-1])*h
    r[i] <- r[i-1] + h*gamma*I[i-1]
    N <- s[i] + I[i] +r[i]
  }

  sir.dat <- list(parameters = data.frame(beta = beta, gamma = gamma, r0 = beta/gamma),
                  sir = cbind(time = tt, susceptible = s, infected = I, recovered = r))
  sir.dat$sir <- as.data.frame(sir.dat$sir)
  if(verbose)
  {
    sir <- sir.dat$sir
    graphics::par(mfrow = c(2, 2))
    graphics::plot(sir[, 1], sir[, 2], type = "l", main = "# of Susceptible people")
    graphics::plot(sir[, 1], sir[, 3], type = "l", main = "# of Infected people")
    graphics::plot(sir[, 1], sir[, 4], type = "l", main = "# of Recovered people")
    graphics::plot(sir[, 1], sir[, 3]/(apply(sir[,-1], 1, sum)), type = "l", main = "Infection rate")
    print(sir.dat$parameters)
  }
  return(sir.dat)
}
