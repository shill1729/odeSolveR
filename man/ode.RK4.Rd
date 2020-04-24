% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odeRK4.R
\name{ode.RK4}
\alias{ode.RK4}
\title{Solve arbitrary first order ODE systems via fourth-order Runge-Kutta scheme}
\usage{
ode.RK4(f, IC, parameters, t0 = 0, tn = 1, n = 1000)
}
\arguments{
\item{f}{a list of functions with as many variables as IC, and in the same order with the names}

\item{IC}{a list of variables with the same names as the non-time input as f, in the same order with initial values}

\item{parameters}{a list of auxillary parameters of the model defined globally (it is not actually used in the function call...)}

\item{t0}{initial time}

\item{tn}{ending time}

\item{n}{number of sub-intervals in time-grid}
}
\value{
data.frame
}
\description{
{An implementation of the fourth-order Runge-Kutta scheme for arbitrary first-order ODE systems of finite dimension.}
}
\details{
{The list of functions must have syntax \code{f(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
}