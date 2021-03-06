% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ode1EulerScheme.R
\name{ode1.EulerScheme}
\alias{ode1.EulerScheme}
\title{Euler scheme for first order ODEs}
\usage{
ode1.EulerScheme(dynamics, y0, t0 = 0, tn = 1, n = 1000)
}
\arguments{
\item{dynamics}{the function of t and y in the explicit ODE \eqn{y'=f(t,y)}}

\item{y0}{the initial value}

\item{t0}{the initial value input}

\item{tn}{the terminal value input}

\item{n}{the number of sub-intervals to use, defaults to 1000}
}
\value{
data.frame of the time grid and solution grid
}
\description{
{One dimensional first order ODE solver.}
}
\details{
{\code{f} must be a function of time and the state variable \code{y}}
}
