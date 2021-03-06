### Misc functions ###
# This script contains miscillaneous functions for supporting the functions
# in the AGHQ package. Some, but nort all, are exported.

#' Default control arguments for \code{aghq::aghq()}.
#'
#' Run \code{default_control()} to print the list of valid control parameters
#' and their defaults, and run with named arguments to change the defaults.
#'
#' @param ... You can provide a named value for any control parameter and its
#' value will be set accordingly. See \code{?aghq} and examples here.
#'
#' @return A list of argument values.
#'
#' @details Valid options are:
#' \itemize{
#' \item{\code{method}: }{optimization method to use:
#' \itemize{
#' \item{'sparse_trust' (default): }{\code{trustOptim::trust.optim} with \code{method = 'sparse'}}
#' \item{'SR1' (default): }{\code{trustOptim::trust.optim} with \code{method = 'SR1'}}
#' \item{'trust': }{\code{trust::trust}}
#' \item{'BFGS': }{\code{optim(...,method = "BFGS")}}
#' }
#' Default is 'sparse_trust'.
#' }
#' \item \code{negate}: default \code{FALSE}. Multiply the functions in \code{ff} by \code{-1}?
#' The reason for having this option is for full compatibility with \code{TMB}:
#' while of course \code{TMB} allows you to code up your log-posterior any way you like,
#' all of its excellent features including its automatic Laplace approximation and MCMC
#' sampling with \code{tmbstan} assume you have coded your template to return the
#' **negated** log-posterior. However, by default, \code{aghq} assumes you have
#' provided the log-posterior **without negation**. Set \code{negate = TRUE} if you
#' have provided a template which computes the **negated** log-posterior and its
#' derivatives.
#' \item \code{ndConstruction}: construct a multivariate quadrature rule using a \code{"product"}
#' rule or a \code{"sparse"} grid? Default \code{"product"}. See \code{?mvQuad::createNIGrid()}.
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'polynomial'}
#' option (default) uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' }
#'
#' @examples
#'
#' default_control()
#' default_control(method = "trust")
#' default_control(negate = TRUE)
#'
#' @export
#'
default_control <- function(...) {
  out <- list(
    method = c("sparse_trust","trust","BFGS"),
    negate = FALSE,
    ndConstruction = "product",
    interpolation = 'polynomial'
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}

#' Default control arguments for \code{aghq::marginal_laplace()}.
#'
#' Run \code{default_control_marglaplace()} to print the list of valid control parameters
#' and their defaults, and run with named arguments to change the defaults.
#'
#' @param ... You can provide a named value for any control parameter and its
#' value will be set accordingly. See \code{?marginal_laplace} and examples here.
#'
#' @return A list of argument values.
#'
#' @details Valid options are:
#' \itemize{
#' \item{\code{method}: }{optimization method to use for the \code{theta} optimization:
#' \itemize{
#' \item{'sparse_trust' (default): }{\code{trustOptim::trust.optim}}
#' \item{'sparse': }{\code{trust::trust}}
#' \item{'BFGS': }{\code{optim(...,method = "BFGS")}}
#' }
#' }
#' \item{\code{inner_method}: }{optimization method to use for the \code{W} optimization; same
#' options as for \code{method}. Default \code{inner_method} is 'sparse_trust' and default \code{method} is 'BFGS'.
#' }
#' \item \code{negate}: default \code{FALSE}. Multiply the functions in \code{ff} by \code{-1}?
#' The reason for having this option is for full compatibility with \code{TMB}:
#' while of course \code{TMB} allows you to code up your log-posterior any way you like,
#' all of its excellent features including its automatic Laplace approximation and MCMC
#' sampling with \code{tmbstan} assume you have coded your template to return the
#' **negated** log-posterior. However, by default, \code{aghq} assumes you have
#' provided the log-posterior **without negation**. Set \code{negate = TRUE} if you
#' have provided a template which computes the **negated** log-posterior and its
#' derivatives. **Note** that I don't expect there to be any reason to need this
#' argument for \code{marginal_laplace}; if you are doing a marginal Laplace approximation
#' using the automatic Laplace approximation provided by \code{TMB}, you should
#' check out \code{aghq::marginal_laplace_tmb()}.
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'polynomial'}
#' option (default) uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' }
#'
#' @examples
#'
#' default_control_marglaplace()
#' default_control_marglaplace(method = "trust")
#' default_control_marglaplace(method = "trust",inner_method = "trust")
#' default_control_marglaplace(negate = TRUE)
#'
#' @export
#'
default_control_marglaplace <- function(...) {
  out <- list(
    method = c("BFGS","sparse_trust","trust"),
    inner_method = c("sparse_trust","trust","BFGS"),
    negate = FALSE,
    ndConstruction = "product",
    interpolation = 'polynomial'
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}

#' Default control arguments for \code{aghq::marginal_laplace()}.
#'
#' Run \code{default_control_marglaplace()} to print the list of valid control parameters
#' and their defaults, and run with named arguments to change the defaults.
#'
#' @param ... You can provide a named value for any control parameter and its
#' value will be set accordingly. See \code{?marginal_laplace} and examples here.
#'
#' @return A list of argument values.
#'
#' @details Valid options are:
#' \itemize{
#' \item{\code{method}: }{optimization method to use for the \code{theta} optimization:
#' \itemize{
#' \item{'sparse_trust' (default): }{\code{trustOptim::trust.optim}}
#' \item{'sparse': }{\code{trust::trust}}
#' \item{'BFGS': }{\code{optim(...,method = "BFGS")}}
#' }
#' }
#' \item \code{negate}: {default \code{TRUE}. Assumes that your \code{TMB} function
#' template computes the **negated** log-posterior, which it must if you're using \code{TMB}'s automatic
#' Laplace approximation, which you must be if you're using this function!}.
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'polynomial'}
#' option (default) uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' }
#'
#' @examples
#'
#' default_control_marglaplace()
#' default_control_marglaplace(method = "trust")
#' default_control_marglaplace(method = "trust",inner_method = "trust")
#' default_control_marglaplace(negate = TRUE)
#'
#' @export
#'
default_control_tmb <- function(...) {
  out <- list(
    method = c("BFGS","sparse_trust","trust"),
    negate = TRUE,
    numhessian = TRUE,
    ndConstruction = 'product',
    interpolation = 'polynomial'
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}

# Logsumexp
# From: https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp
# Accessed on: 2021/03/27 10:49AM
logsumexp <- function(l) {
  if (length(l) == 0) return(c())
  if (length(l) == 1) return(l)
  n <- length(l)
  L <- sort(l, decreasing = TRUE)
  S <- rep(L[1], n)
  for (k in 1:(n-1)) {
    S[k+1] <- max(L[k+1], S[k]) + log1p(exp(-abs(L[k+1] - S[k])))
  }
  S[n]
}
logdiffexp <- function(a,b) {
  # computes log(exp(a) - exp(b)) where a > b
  if (length(b) == 0) return(a)
  if (b > a) stop("Error: negative weights outweigh the positive weights for this quadrature rule. This would result in negative density estimates and is obviously incorrect.\n")
  b + log(exp(a-b) - 1)
}
