% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_phase_portrait.R
\name{plot_phase_portrait}
\alias{plot_phase_portrait}
\title{Title}
\usage{
plot_phase_portrait(sol, states = NULL)
}
\arguments{
\item{sol}{solution object returned from \code{ode}, \code{ode1}, etc.}

\item{states}{vector of column indexes: must be of length 2 or 3 and contain numbers greater than 2}
}
\value{
null
}
\description{
{Plots the phase plane or portrait of a 2-dimensional or 3-dimensional solution to a ODE system.
Otherwise the user can specify which variables of a higher-dimensional system to plot in a phase-portrait.}
}
