#' Solve arbitrary first order ODE systems via fourth-order Runge-Kutta scheme
#'
#' @param f a list of functions with as many variables as IC, and in the same order with the names
#' @param IC a list of variables with the same names as the non-time input as f, in the same order with initial values
#' @param parameters a list of auxillary parameters of the model defined globally (it is not actually used in the function call...)
#' @param t0 initial time
#' @param tn ending time
#' @param n number of sub-intervals in time-grid
#'
#' @description {An implementation of the fourth-order Runge-Kutta scheme for arbitrary first-order ODE systems of finite dimension.}
#' @details {The list of functions must have syntax \code{f(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
#' @return data.frame
ode.RK4 <- function(f, IC, parameters, t0 = 0, tn = 1, n = 1000)
{
  # Time step size
  h <- (tn-t0)/n
  # Time grid
  tt <- seq(t0, tn, length.out = n+1)
  # Number of state variables
  m <- length(IC)
  # state variables: list of grids of each state-variable. This will be the returned solution
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
    k1 <- list()
    k2 <- list()
    k3 <- list()
    k4 <- list()
    # For each state variable
    for(i in 1:m)
    {
      # Gather input: we need four evaluations of f(t, ....)
      input1 <- list()
      input2 <- list()
      input3 <- list()
      input4 <- list()

      # Time argument is first, always
      input1[[1]] <- tt[j-1]
      input2[[1]] <- tt[j-1]+0.5*h
      input3[[1]] <- tt[j-1]+0.5*h
      input4[[1]] <- tt[j-1]+h
      for(l in 1:m)
      {
        # The first input just takes the states
        input1[[l+1]] <- states[[l]][j-1]

      }
      names(input1) <- c("t", names(IC))
      # Evaluate k1=f(t, x, y, z,...) for all states
      for(l in 1:m)
      {
        k1[[l]] <- do.call(what = f[[l]], args = input1)
      }
      # Get input ready for k2
      for(l in 1:m)
      {
        # Second input is state[j-1]+0.5*h*k1[l]
        input2[[l+1]] <- states[[l]][j-1]+0.5*h*k1[[l]]

      }
      names(input2) <- c("t", names(IC))
      # Evaluate k2=f(t, x, y, z,...) for all states
      for(l in 1:m)
      {
        k2[[l]] <- do.call(what = f[[l]], args = input2)
      }
      # Get input ready for k3
      for(l in 1:m)
      {
        input3[[l+1]] <- states[[l]][j-1]+0.5*h*k2[[l]]

      }
      names(input3) <- c("t", names(IC))
      # Evaluate k3=f(t, x, y, z,...) for all states
      for(l in 1:m)
      {
        k3[[l]] <- do.call(what = f[[l]], args = input3)
      }
      # Get input ready for k4 the final term
      for(l in 1:m)
      {
        input4[[l+1]] <- states[[l]][j-1]+k3[[l]]
      }
      names(input4) <- c("t", names(IC))
      for(l in 1:m)
      {
        k4[[l]] <- do.call(what = f[[l]], args = input4)
      }

      # The Runge-Kutta method: y[j]=y[j-1]+(h/6)(k1+2k2+2k3+k4), over every state variable and from their respective ODE function
      states[[i]][j] <- states[[i]][j-1]+(h/6)*(k1[[i]]+2*k2[[i]]+2*k3[[i]]+k4[[i]])
    }
  }
  states <- data.frame(time = tt, do.call(cbind, states))
  names(states) <- c("time", names(IC))
  return(states)
}
