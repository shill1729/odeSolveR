#' Solve arbitrary first order ODE systems
#'
#' @param dynamics a list of functions with as many variables as IC, and in the same order with the names defining the dynamics of the system
#' @param IC a list of variables with the same names as the non-time input as dynamics, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param n number of sub-intervals in time-grid
#'
#' @description {A general Euler scheme implementation for arbitrary first-order ODE systems of finite dimension.}
#' @details {The list of functions each element with arguments \code{(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
#' @return data.frame
ode.EulerScheme <- function(dynamics, IC, t0 = 0, tn = 1, n = 1000)
{
  # Time step size
  h <- (tn-t0)/n
  # Time grid
  tt <- seq(t0, tn, length.out = n+1)
  # Number of state variables
  m <- length(IC)
  # state variables
  states <- list()
  for(i in 1:m)
  {
    states[[i]] <- matrix(0, n+1)
    # Initialize ICs
    states[[i]][1] <- IC[[i]]
  }
  # Over every time-node
  for(j in 2:(n+1))
  {
    # For each state variable
    for(i in 1:m)
    {
      # Gather input
      input <- list()
      input[[1]] <- tt[j-1]
      for(l in 1:m)
      {
        input[[l+1]] <- states[[l]][j-1]
      }
      names(input) <- c("t", names(IC))
      states[[i]][j] <- states[[i]][j-1]+h*do.call(what = dynamics[[i]], args = input)
    }
  }
  states <- data.frame(time = tt, do.call(cbind, states))
  names(states) <- c("time", names(IC))
  return(states)
}
