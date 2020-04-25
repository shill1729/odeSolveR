#' Plot state-variables trajectory over time given a numerical solution
#'
#' @param sol solution to plot, object returned from \code{ode}, \code{ode1} etc
#' @param legend_names optional, will otherwise use the names given to \code{sol}
#' @param legend_loc where to place the legend: "topleft", "topright", etc.
#' @param legend_size number between 0 and 1 for scaling the legend box
#'
#' @description {Function for plotting trajectories of the solution to an ODE system.}
#' @return null
#' @export plot_trajectories
plot_trajectories <- function(sol, legend_names = NULL, legend_loc = "topright", legend_size = 0.9)
{
  # Check to see if solution is finite before graphing
  if(!all(apply(sol, 2, is.finite)))
  {
    stop("Solution blows up or is NaN at some point")
  }
  # Get legend names if not passed
  if(is.null(legend_names))
  {
    legend_names <-  c(names(sol[,-1]))
  }
  m <- dim(sol)[2]
  graphics::plot(sol$time, sol[, 2], type = "l", main = "Solution trajectory over time", xlab = "time", ylab = "state", ylim = c(min(sol[,-1]), max(sol[, -1])))
  for(i in 3:m)
  {
    graphics::lines(sol$time, sol[, i], col = i-1)
  }
  graphics::legend(x = legend_loc, legend = legend_names, col = c(1:(m-1)), lty = 1, cex = legend_size)
}
