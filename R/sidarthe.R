#' Euler scheme for the SIDARTHE model
#'
#' @param parameters named list of parameters see details
#' @param initial_conditions named list of initial conditions, see details
#' @param t0 initial time point
#' @param tn the terminal time point
#' @param n number of sub-intervals in time-grid
#' @param verbose whether to plot the populations over time and print parameter+reproductive ratio
#'
#' @description {A basic Euler scheme is used to numerically solve the ODE system of the novel SIDARTHE model.}
#' @details {The parameters must be a named list containing values
#' \itemize{
#' \item \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta} "represent the transimission rate
#' (i.e.  the probability of disease transmission in a singlecontact times the average number of contacts
#' per person) due to contacts between a Susceptible subject and an Infected, a Diagnosed, an Ailing,
#' a Recognised subject. Typically, \eqn{\alpha} is larger than \eqn{\gamma} (assuming that people tend to
#' avoid contacts with subjects showing symptoms, even though diagnosis has not been made yet), which in turn
#' is probably larger than \eqn{\beta} and \eqn{\delta} (assuming that subjects who have been diagnosed are
#' properly isolated). These parameters can be modified by social distancing policies (e.g., closing schools,
#' remote working, etc.). The risk of contagion due to Threatened subjects, treated in proper ICUs, is assumed negligible."
#'
#' \item \eqn{\epsilon} and \eqn{\theta} "capture capture  the  probability  rate  of  detection,  relative
#' to  asymptomatic  and  mildly  symptomatic  casesrespectively.  These parameters, also modifiable,
#' reflect the level of attention on the disease and the numberof tests performed over the population. Note that
#' \eqn{\theta} is typically larger than \eqn{\epsilon}, since symptomatic individuals are more likely to get tested."
#'
#' \item \eqn{\zeta} and \eqn{\eta} "denote the probability rate at which an infected subject,
#' respectively not aware and aware of being infected, develops clinically relevant symptoms, and
#' are probably comparable.  These parameters are disease-dependent and hardly modifiable."
#'
#' \item \eqn{\mu} and \eqn{\nu} "respectively denote the rate at which undetected and detected infected
#' subjects develop life-threatening symptoms, and are likely to be comparable if there is no known specific
#' treatment that is effective against the disease, otherwise \eqn{\mu} is likely to be larger.
#' These parameters can be reduced by means of improved therapies and acquisition of immunity against
#' the virus."
#'
#' \item \eqn{\tau} "denotes the mortality rate (for infected subjects with life-threatening symptoms) and
#' can be reduced bymeans of improved therapies and treatments."
#'
#' \item \eqn{\lambda, \kappa, \xi, \rho} and \eqn{\sigma} "denote the rate of recovery for the five
#' classes of infected subjects, and may differ significantly if an appropriate treatment for the
#' disease is known and adopted to diagnosed patients, while are probably comparable otherwise.
#' These parameters can be increased thanks to improved treatments and acquisition of immunity against the
#' virus."
#' }
#' These descriptions of the model parameters are directly taken from the 2020-03-21 paper "A SIDARTHE model of COVID-19 Epidemic in Italy"
#' by the COVID19 IRCCS San Matteo Pavia Task Force et al. to mitigate any misinterpretation of the input into the model.
#'
#' (Note the named list must actually have the greek letter names written out, i.e. \code{alpha}, \code{beta},...)
#'
#' Finally the list \code{initial_conditions} must contain the following, all as fractions of the total population
#' \itemize{
#' \item \code{s0} the initial level of susceptible individuals
#' \item \code{i0} the initial level of infected individuals
#' \item \code{d0} the initial level of diagnosed individuals
#' \item \code{a0} the initial level of ailing individuals
#' \item \code{r0} the initial level of recognized individuals
#' \item \code{thr0} the initial level of threatened individuals
#' \item \code{h0} the initial level of heal individuals
#' \item \code{e0} the initial level of extinct individuals}
#'
#' It is best to specify all but \code{s0} and then set the initial susceptible fraction to 1 less the sum of the rest.
#' }
#' @return list consisting of (1) parameters, a data.frame containing the contact rates, etc and (2) sidarthe, a matrix whose columns consist of the
#' eight (fractions of) populations over time:
#' \itemize{
#' \item suspected
#' \item infected
#' \item diagnosed
#' \item ailing
#' \item recognized
#' \item threatened
#' \item healed
#' \item extinct}
#' and finally (3) the basic reproductive ratio as defined for this extended SIR model as a consequence of Proposition 2 in the original paper.
#' @export sidarthe
#' @references {The Euler scheme implemented by this function is for the model developed in the following paper. It is 8 coupled non-linear ODEs describing dynamics of sub-populations under a pandemic.
#' The descriptions within the novel paper were used to document the initial conditions and parameters of this model in this package.
#' \itemize{\item "A SIDARTHE Model of COVID-19 Epidemic in Italy" by the COVID19 IRCCS San Matteo Pavia Task Force et al. \href{https://arxiv.org/abs/2003.09861}{https://arxiv.org/abs/2003.09861}}}
sidarthe <- function(parameters, initial_conditions, t0 = 0, tn = 1, n = 1000, verbose = TRUE)
{
  # The massive parameter list detailing mean interaction/recovery rates
  alpha <- parameters$alpha
  beta <- parameters$beta
  gamma <- parameters$gamma
  delta <- parameters$delta
  epsilon <- parameters$epsilon
  theta <- parameters$theta
  zeta <- parameters$zeta
  eta <- parameters$eta
  mu <- parameters$mu
  nu <- parameters$nu
  tau <- parameters$tau
  lambda <- parameters$lambda
  kappa <- parameters$kappa
  xi  <- parameters$xi
  rho <- parameters$rho
  sigma <- parameters$sigma

  # The initial condition parameter list
  s0 <- initial_conditions$s0
  i0 <- initial_conditions$i0
  d0 <- initial_conditions$d0
  a0 <- initial_conditions$a0
  r0 <- initial_conditions$r0
  thr0 <- initial_conditions$thr0
  h0 <- initial_conditions$h0
  e0 <- initial_conditions$e0

  # Time step size
  h <- (tn-t0)/n
  # the population sub-samples
  s <- matrix(0, n+1)
  I <- matrix(0, n+1)
  D <- matrix(0, n+1)
  A <- matrix(0, n+1)
  R <- matrix(0, n+1)
  Thr <- matrix(0, n+1)
  H <- matrix(0, n+1)
  E <- matrix(0, n+1)

  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  s[1] <- s0
  I[1] <- i0
  D[1] <- d0
  A[1] <- a0
  R[1] <- r0
  Thr[1] <- thr0
  H[1] <- h0
  E[1] <- e0

  # N <- s[1] + I[1] +D[1] + A[1] + R[1] + Thr[1] + H[1] + E[1]

  # The time-stepping Euler scheme
  for(i in 2:(n+1))
  {
    s[i] <- s[i-1] - h*s[i-1]*(alpha*I[i-1]+beta*D[i-1]+gamma*A[i-1]+delta*R[i-1])
    I[i] <- I[i-1] + h*(s[i-1]*(alpha*I[i-1]+beta*D[i-1]+gamma*A[i-1]+delta*R[i-1])-(epsilon+zeta+lambda)*I[i-1])
    D[i] <- D[i-1] + h*(epsilon*I[i-1]-(eta+rho)*D[i-1])
    A[i] <- A[i-1] + h*(zeta*I[i-1]-(theta+mu+kappa)*A[i-1])
    R[i] <- R[i-1] + h*(eta*D[i-1]+theta*A[i-1]-(nu+xi)*R[i-1])
    Thr[i] <- Thr[i-1] + h*(mu*A[i-1] + nu*R[i-1] - (sigma+tau)*Thr[i-1])
    H[i] <- H[i-1] + h*(lambda*I[i-1] + rho*D[i-1] + kappa*A[i-1] +xi*R[i-1] +sigma*Thr[i-1])
    E[i] <- E[i-1] + h*tau*Thr[i-1]

    # We omit population absolute levels and use fractional levels.
    # N <- s[i] + I[i] +r[i]
  }

  # Computing the basic reproduction ratio
  r1 <- epsilon+zeta+lambda
  r2 <- eta + rho
  r3 <- theta + mu + kappa
  r4 <- nu + xi
  repo <- alpha/r1 + beta*epsilon/(r1*r2) + gamma*zeta/(r1*r3) + delta*eta*epsilon/(r1*r2*r4) + delta*zeta*theta/(r1*r3*r4)

  # The object to return: a list containing parameters, the model output, and the reproduction number
  sir.dat <- list(parameters = as.data.frame(parameters),
                  sidarthe = cbind(time = tt,
                              susceptible = s,
                              infected = I,
                              diagnosed = D,
                              ailing = A,
                              recognized = R,
                              threatened = Thr,
                              healed = H,
                              extinct = E),
                  reproductive_ratio = repo)
  colnames(sir.dat$sidarthe) <- c("time" , "susceptible", "infected", "diagnosed", "ailing", "recognized", "threatened", "healed", "extinct")
  sir.dat$sidarthe <- as.data.frame(sir.dat$sidarthe)
  if(verbose)
  {
    sir <- sir.dat$sidarthe
    # graphics::par(mfrow = c(1, 1))
    graphics::plot(sir[, 1], sir[, 2], type = "l", main = "Fractions of population", ylim = c(0, 1), xlab = "days", ylab = "% of population")
    for(i in 3:9)
    {
      graphics::lines(sir[, 1], sir[, i], col = i-1)
    }
    graphics::legend(x = "topright", legend = c("Susceptible", "Infected", "Diagnosed", "Ailing", "Recognized", "Threatened", "Healed", "Extinct"), lty = 1, col = c(1:10))
    print(sir.dat$parameters)
    print(paste("Reproductive Ratio =", sir.dat$reproductive_ratio))
  }
  return(sir.dat)
}
