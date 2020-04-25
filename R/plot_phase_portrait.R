#' Title
#'
#' @param sol solution object returned from \code{ode}, \code{ode1}, etc.
#' @param states vector of column indexes: must be of length 2 or 3 and contain numbers greater than 2
#'
#' @description {Plots the phase plane or portrait of a 2-dimensional or 3-dimensional solution to a ODE system.
#' Otherwise the user can specify which variables of a higher-dimensional system to plot in a phase-portrait.}
#' @return null
#' @export plot_phase_portrait
plot_phase_portrait <- function(sol, states = NULL)
{

  # Get dimension of state-variables:
  d <- dim(sol)[2]-1
  if(is.null(states))
  {
    # If states are not passed, assume its a 2 or 3-dimensional system and just plot it.
    if(d == 2)
    {
      states <- c(2, 3)
      graphics::plot(sol[, 2], sol[, 3], type = "l", main = "Phase plane", xlab = names(sol)[2], ylab = names(sol)[3])
    } else if(d == 3)
    {
      states <- c(2, 3)
      plot3D::lines3D(sol[, 2], sol[, 3], sol[, 4], type = "l", main = "Phase portrait",)
    }
  } else
  {
    if(1 %in% states)
    {
      stop("Phase diagram does not plot time variable: use columns indexes greater than 1")
    }
    if(length(states) == 2)
    {
      graphics::plot(sol[, states[1]], sol[, states[2]], type = "l", main = "Phase plane", xlab = names(sol)[states[1]], ylab = names(sol)[states[2]])
    } else if(length(states) == 3)
    {
      plot3D::lines3D(sol[, states[1]], sol[, states[2]], sol[, states[3]], type = "l", main = "Phase portrait",)
    } else if(length(states) >3)
    {
      stop("Cannot plot in more than 3-coordinates")
    }
  }

}
