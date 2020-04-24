#' Classify equilibrium points of a dynamical system
#'
#' @param A the jacobian matrix evaluated at a given equilibrium point (model dependent)
#' @description {Given the Jacobian matrix evaluated at a fixed point, this function classifies the stability and whether it is a sink or source.}
#' @details {The classification scheme is based off the Poincare diagram in the \eqn{Tr(A), det(A)} plane. The user must know the form of the Jacobian matrix
#' as well as the fixed points and be able to compute the Jacobian at each of the fixed points. These are heavily model dependent.}
#' @return list
#' @export classify_equilibrium
classify_equilibrium <- function(A)
{
  stab <- ""
  type <- ""

  # Trace and determinant and discriminant of characteristic polynomial associated to eigenvalues
  tau <- sum(diag(A))
  Delta <- det(A)
  discriminant <- tau^2-4*Delta
  ev <- eigen(x = A)$values

  # Real distinct eigenvalues
  if(discriminant > 0)
  {
    # Both eigenvalues negative
    if(tau < 0 && Delta > 0)
    {
      stab <- "stable"
      type <- "node"
    } else if(tau > 0 && Delta > 0) # Both eigenvalues positive
    {
      stab <- "unstable"
      type <- "node"
    } else if(Delta < 0) # Eigenvalues of opposite sign
    {
      stab <- "unstable"
      type <- "saddle point"
    } else if(tau > 0 && Delta == 0)
    {
      stab <- "unstable"
      type <- "degenerate line"
    } else if(tau < 0 && Delta == 0)
    {
      stab <- "stable"
      type <= "degenerate line"
    }
  } else if(discriminant == 0)
  {
    if(ev[1] < 0)
    {
      stab <- "stable"
      type <- "degenerate node"
    } else if(ev[1] > 0)
    {
      stab <- "unstable"
      type <- "degenerate node"
    }
  } else if(discriminant < 0)
  {
    if(tau == 0)
    {
      type <- "center"
    } else if(tau != 0)
    {
      if(tau < 0)
      {
        stab <- "stable"
      } else if(tau > 0)
      {
        stab <- "unstable"
      }
      type <- "spiral point"
    }
  }
  classification <- paste(stab, type)
  output <- list(jacobian = A, eigenvalues = ev, trace = tau, det = Delta, discr = discriminant, classification = classification)
  return(output)
}

