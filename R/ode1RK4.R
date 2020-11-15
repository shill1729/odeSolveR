#' Runge-Kutta for first order ODE
#'
#' @param dynamics the function of t and y in the explicit ODE \eqn{y'=f(t,y)}
#' @param y0 the initial value
#' @param t0 the initial value input
#' @param tn the terminal value input
#' @param n the number of sub-intervals to use, defaults to 1000
#'
#' @description {One dimensional first order ODE solver via Runge-Kutta fourth-order method.}
#' @details {\code{f} must be a function of time and the state variable \code{y}}
#' @return data.frame of the time grid and solution grid
ode1.RK4 <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
{
  h <- (tn-t0)/n
  y <- matrix(0, n+1)
  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  y[1] <- y0
  for(i in 2:(n+1))
  {
    k1 <- dynamics(tt[i-1], y[i-1])
    k2 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k1)
    k3 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k2)
    k4 <- dynamics(tt[i-1]+h, y[i-1]+h*k3)

    y[i] <- y[i-1] + (h/6)*(k1+2*k2+2*k3+k4)
  }
  return(data.frame(time = tt, state = y))
}

