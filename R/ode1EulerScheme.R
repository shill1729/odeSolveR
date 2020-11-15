#' Euler scheme for first order ODEs
#'
#' @param dynamics the function of t and y in the explicit ODE \eqn{y'=f(t,y)}
#' @param y0 the initial value
#' @param t0 the initial value input
#' @param tn the terminal value input
#' @param n the number of sub-intervals to use, defaults to 1000
#'
#' @description {One dimensional first order ODE solver.}
#' @details {\code{f} must be a function of time and the state variable \code{y}}
#' @return data.frame of the time grid and solution grid
ode1.EulerScheme <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
{
  h <- (tn-t0)/n
  y <- matrix(0, n+1)
  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  y[1] <- y0
  for(i in 2:(n+1))
  {
    y[i] <- y[i-1] + h*dynamics(tt[i-1], y[i-1])
  }
  return(data.frame(time = tt, state = y))
}
