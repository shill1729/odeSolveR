#' Solve arbitrary first order ODE systems
#'
#' @param dynamics a list of functions with as many variables as IC, and in the same order with the names
#' @param IC a list of variables with the same names as the non-time input as f, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param n number of sub-intervals in time-grid
#' @param engine which scheme to use to solve the ODE; 'RK4' or 'Euler', see details.
#'
#' @description {A solver for arbitrary first-order ODE systems of finite dimension, using either RK4 (default) or an Euler scheme.}
#' @details {The list of functions must have syntax e.g. \code{f(t,x,y,z)} for each element matching the named elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.
#' Further, the engine must be either "RK4" for a fourth-order Runge-Kutta method, or "Euler" for an Euler scheme.}
#' @return data.frame
#' @export ode
ode <- function(dynamics, IC, t0 = 0, tn = 1, n = 1000, engine = "RK4")
{
  if(engine == "RK4")
  {
    # print("Solving ODE using fourth-order Runge-Kutta scheme")
    return(ode.RK4(dynamics, IC, t0, tn, n))
  } else if(engine == "Euler")
  {
    # print("Solving ODE using simple Euler scheme")
    return(ode.EulerScheme(dynamics, IC, t0, tn, n))
  } else
  {
    stop("engine must be either: 'RK4' or 'Euler'")
  }
}
