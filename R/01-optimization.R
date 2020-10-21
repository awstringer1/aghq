### Optimization ###
# This script contains a function for optimizing
# the un-normalized posterior. The mode and curvature
# are required for AGHQ.

#' Obtain function information necessary for performing quadrature
#'
#' This function computes the two pieces of information needed about
#' the log posterior to do adaptive quadrature: the mode, and the hessian at the mode.
#' It is designed for use within \code{aghq::aghq}, but is exported in case users
#' need to debug the optimization process and documented in case users want to
#' write their own optimizations.
#'
#' @param ff A list with three elements:
#' \itemize{
#' \item{\code{fn}}{: function taking argument \code{theta} and returning a numeric
#' value representing the log-posterior at \code{theta}}
#' \item{\code{gr}}{: function taking argument \code{theta} and returning a numeric
#' vector representing the gradient of the log-posterior at \code{theta}}
#' \item{\code{he}}{: function taking argument \code{theta} and returning a numeric
#' matrix representing the hessian of the log-posterior at \code{theta}}
#' }
#' The user may wish to use \code{numDeriv::grad} or \code{numDeriv::hessian} to
#' obtain these.
#' @param startingvalue Value to start the optimization. \code{ff$fn(startingvalue)},
#' \code{ff$gr(startingvalue)}, and \code{ff$he(startingvalue)} must all return
#' appropriate values without error.
#' @param control A list with elements
#' \itemize{
#' \item{\code{method}: }{optimization method to use:
#' \itemize{
#' \item{'sparse_trust' (default): }{\code{trustOptim::trust.optim}}
#' \item{'sparse': }{\code{trust::trust}}
#' \item{'BFGS': }{\code{optim(...,method = "BFGS")}}
#' }
#' }
#' }
#' @return A list with elements
#' \itemize{
#' \item{\code{solution}: }{the mode of the log posterior}
#' \item{\code{hessian}: }{the hessian of the log posterior at the mode}
#' \item{\code{convergence}: }{specific to the optimizer used, a message indicating whether it converged}
#' }
#'
#' @examples
#' # Poisson/Exponential example
#' logfteta <- function(eta,y) {
#'   sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
#' }
#'
#' y <- rpois(10,5) # Mode should be sum(y) / (10 + 1)
#'
#' objfunc <- function(x) logfteta(x,y)
#' funlist <- list(
#'   fn = objfunc,
#'   gr = function(x) numDeriv::grad(objfunc,x),
#'   he = function(x) numDeriv::hessian(objfunc,x)
#' )
#'
#' optimize_theta(funlist,1.5)
#' optimize_theta(funlist,1.5,control = list(method = "trust"))
#' optimize_theta(funlist,1.5,control = list(method = "BFGS"))
#'
#' @importFrom stats optim
#' @importFrom utils installed.packages
#' @importFrom methods as
#' @export
optimize_theta <- function(ff,startingvalue,control = default_control()) {
  optfunc <- function(x) -1 * ff$fn(x)
  optgrad <- function(x) as.numeric(-1 * ff$gr(x))
  opthess <- function(x) as.matrix(-1 * ff$he(x))

  method <- control$method[1]

  if (method == "sparse_trust") {
    if (!("trustOptim" %in% rownames(installed.packages()))) stop("Method = sparse_trust requires the trustOptim package, but you do not have this package installed.")
    opt <- trustOptim::trust.optim(
      x = startingvalue,
      fn = optfunc,
      gr = optgrad,
      hs = function(x) as(opthess(x),"dgCMatrix"),
      method = "Sparse"
    )
    out <- list(
      solution = opt$solution,
      hessian = opt$hessian,
      convergence = opt$status
    )
  }
  else if (method == "trust") {
    if (!("trustOptim" %in% rownames(installed.packages()))) stop("Method = trust requires the trust package, but you do not have this package installed.")
    funlist <- function(x) {
      list(
        value = optfunc(x),
        gradient = optgrad(x),
        hessian = opthess(x)
      )
    }
    opt <- trust::trust(
      objfun = funlist,
      parinit = startingvalue,
      rinit = 1
    )
    out <- list(
      solution = opt$argument,
      hessian = opt$hessian,
      convergence = opt$converged
    )
  }
  else if (method == "BFGS") {
    opt <- optim(startingvalue,optfunc,optgrad,method = "BFGS")
    out <- list(
      solution = opt$par,
      hessian = opthess(opt$par),
      convergence = opt$convergence
    )
  }
  else {
    stop(paste0("Unknown optimization method: ",method))
  }

  out
}