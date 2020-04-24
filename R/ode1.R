#' Solve for first order ODEs
#'
#' @param f the function of t and y in the explicit ODE \eqn{y'=f(t,y)}
#' @param y0 the initial value
#' @param t0 the initial value input
#' @param tn the terminal value input
#' @param n the number of sub-intervals to use, defaults to 1000
#' @param engine the type of scheme to use, defaults to RK4.
#'
#' @description {One dimensional first order ODE solver using either RK4 or an Euler scheme..}
#' @details {\code{f} must be a function of time and the state variable \code{y}. The engine must either be "RK4" or "Euler"}
#' @return data.frame of the time grid and solution grid
#' @export ode1
ode1 <- function(f, y0, t0 = 0, tn = 1, n = 1000, engine = "RK4")
{
  if(engine == "RK4")
  {
    print("Solving ODE using fourth-order Runge-Kutta scheme")
    return(ode1.RK4(f, y0, t0, tn, n))
  } else if(engine == "Euler")
  {
    print("Solving ODE using simple Euler scheme")
    return(ode1.EulerScheme(f, y0, t0, tn, n))
  } else
  {
    stop("engine must be either: 'RK4' or 'Euler'")
  }
}
