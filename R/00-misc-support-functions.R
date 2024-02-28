


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
#' \item{'BFGS' (default): }{\code{optim(...,method = "BFGS")}}
#' \item{'sparse_trust': }{\code{trustOptim::trust.optim}}
#' \item{'SR1': }{\code{trustOptim::trust.optim} with \code{method = 'SR1'}}
#' \item{'sparse': }{\code{trust::trust}}
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
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'auto'} option
#' (default) chooses for you and should always work well. The \code{'polynomial'}
#' option uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' \item{numhessian}: logical, default \code{FALSE}. Replace the \code{ff$he} with a numerically-differentiated
#' version, by calling \code{numDeriv::jacobian} on \code{ff$gr}. Used mainly for \code{TMB} with the automatic
#' Laplace approximation, which does not have an automatic Hessian.
#' \item{onlynormconst}: logical, default \code{FALSE}. Skip everything after the calculation of the log integral,
#' and just return the numeric value of the log integral. Saves computation time, and most useful in cases
#' where \code{aghq} is being used as a step in a more complicated procedure.
#' \item{method_summaries}: default \code{'reuse'}, method to use to compute moments and marginals. Choosing
#' \code{'correct'} corresponds to the approximations suggested in the *Stochastic Convergence...* paper,
#' which attain the same rate of convergence as the approximation to the marginal likelihood. See \code{?compute_moment}.
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
    method = c("BFGS","sparse_trust","trust"),
    negate = FALSE,
    ndConstruction = "product",
    interpolation = 'auto',
    numhessian = FALSE,
    onlynormconst = FALSE,
    method_summaries = c('reuse','correct'),
    verbose=FALSE
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
#' \item{'BFGS' (default): }{\code{optim(...,method = "BFGS")}}
#' \item{'sparse_trust': }{\code{trustOptim::trust.optim}}
#' \item{'SR1': }{\code{trustOptim::trust.optim} with \code{method = 'SR1'}}
#' \item{'sparse': }{\code{trust::trust}}
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
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'auto'} option
#' (default) chooses for you and should always work well. The \code{'polynomial'}
#' option uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' \item{numhessian}: logical, default \code{FALSE}. Replace the \code{ff$he} with a numerically-differentiated
#' version, by calling \code{numDeriv::jacobian} on \code{ff$gr}. Used mainly for \code{TMB} with the automatic
#' Laplace approximation, which does not have an automatic Hessian.
#' \item{onlynormconst}: logical, default \code{FALSE}. Skip everything after the calculation of the log integral,
#' and just return the numeric value of the log integral. Saves computation time, and most useful in cases
#' where \code{aghq} is being used as a step in a more complicated procedure.
#' \item{method_summaries}: default \code{'reuse'}, method to use to compute moments and marginals. Choosing
#' \code{'correct'} corresponds to the approximations suggested in the *Stochastic Convergence...* paper,
#' which attain the same rate of convergence as the approximation to the marginal likelihood. See \code{?compute_moment}.
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
    inner_method = c("BFGS","sparse_trust","trust"),
    negate = FALSE,
    ndConstruction = "product",
    interpolation = 'auto',
    numhessian = FALSE,
    onlynormconst = FALSE,
    method_summaries = c('reuse','correct'),
    verbose = FALSE
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}

#' Default control arguments for \code{aghq::marginal_laplace_tmb()}.
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
#' \item{'BFGS' (default): }{\code{optim(...,method = "BFGS")}}
#' \item{'sparse_trust': }{\code{trustOptim::trust.optim}}
#' \item{'SR1': }{\code{trustOptim::trust.optim} with \code{method = 'SR1'}}
#' \item{'sparse': }{\code{trust::trust}}
#' }
#' }
#' \item \code{negate}: {default \code{TRUE}. Assumes that your \code{TMB} function
#' template computes the **negated** log-posterior, which it must if you're using \code{TMB}'s automatic
#' Laplace approximation, which you must be if you're using this function!}.
#' \item \code{interpolation}: how to interpolate the marginal posteriors. The \code{'auto'} option
#' (default) chooses for you and should always work well. The \code{'polynomial'}
#' option uses \code{polynom::poly.calc()} to construct a global polynomial interpolant
#' and has been observed to be unstable as the number of quadrature points gets larger, which
#' is obviously a bad thing. Try \code{'spline'} instead, which uses a cubic B-Spline
#' interpolant from \code{splines::interpSpline()}.
#' \item{numhessian}: logical, default \code{TRUE}. Replace the \code{ff$he} with a numerically-differentiated
#' version, by calling \code{numDeriv::jacobian} on \code{ff$gr}. Used mainly for \code{TMB} with the automatic
#' Laplace approximation, which does not have an automatic Hessian.
#' \item{onlynormconst}: logical, default \code{FALSE}. Skip everything after the calculation of the log integral,
#' and just return the numeric value of the log integral. Saves computation time, and most useful in cases
#' where \code{aghq} is being used as a step in a more complicated procedure.
#' \item{method_summaries}: default \code{'reuse'}, method to use to compute moments and marginals. Choosing
#' \code{'correct'} corresponds to the approximations suggested in the *Stochastic Convergence...* paper,
#' which attain the same rate of convergence as the approximation to the marginal likelihood. See \code{?compute_moment}.
#' }
#'
#' @examples
#'
#' default_control_tmb()
#' default_control_tmb(method = "trust")
#'
#' @export
#'
default_control_tmb <- function(...) {
  out <- list(
    method = c("BFGS","sparse_trust","trust"),
    negate = TRUE,
    numhessian = TRUE,
    ndConstruction = 'product',
    interpolation = 'auto',
    onlynormconst = FALSE,
    method_summaries = c('reuse','correct'),
    verbose = FALSE
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}

#' Validate a control list
#'
#' This function checks that the names and value types for a supplied `control` list
#' are valid and are unlikely to cause further errors within `aghq` and related functions.
#' Users should not have to worry about this and should just use `default_control()` and related
#' constructors.
#'
#' @param control A list.
#' @param type One of `c('aghq','marglapace','tmb')`. The type of control object to validate. Will
#' basically validate against the arguments required by `aghq`, `marginal_laplace`, and `marginal_laplace_tmb`,
#' respectively.
#' @param ... Not used.
#'
#' @return Logical, `TRUE` if the list of control arguments is valid, else `FALSE`.
#'
#' @details To users reading this: just use `default_control()`, `default_control_marglaplace()`, or `default_control_tmb()`
#' as appropriate, to ensure that your control arguments are correct. This function just exists to provide more
#' descriptive error messages in the event that an incompatible list is provided.
#'
#' @examples
#' validate_control(default_control())
#' validate_control(default_control_marglaplace(),type = "marglaplace")
#' validate_control(default_control_tmb(),type = "tmb")
#'
#' @export
#'
validate_control <- function(control,type = c('aghq','marglaplace','tmb'),...) {
  type <- type[1]
  if (!(type%in%c('aghq','marglaplace','tmb'))) stop(paste0("Unrecognized type argument. Should be one of 'aghq','marglaplace','tmb'. You provided: ",type))

  # Choose the default args to check against
  args_to_check <- switch(type,aghq = default_control(),marglaplace = default_control_marglaplace(),tmb = default_control_tmb())

  names_to_check <- names(args_to_check)
  names_provided <- names(control)
  if (!setequal(names_to_check,names_provided)) {
    missing_names <- setdiff(names_to_check,names_provided)
    extra_names <- setdiff(names_provided,names_to_check)
    msg <- paste0("Problem with the provided control list. It should have names: ",paste0(names_to_check,collapse = ", "),". Yours has names: ",paste0(names_provided,collapse = ", "),".\n")
    if (length(missing_names) > 0) msg <- paste0(msg,"You are therefore missing: ",paste0(missing_names,collapse = ", "),".\n")
    if (length(extra_names) > 0) msg <- paste0(msg,"You have provided extra names: ",paste0(extra_names,collapse=", "),".\n")
    msg <- paste0(msg,"Try using default_control(), default_control_marglaplace(), or default_control_tmb() to construct your control list.\n")
    stop(msg)
  }

  types_to_check <- Map(class,args_to_check)
  types_provided <- Map(class,control)
  msg <- paste0("Problem with provided control list. Arguments ",names(types_to_check)," should have types ",types_to_check,".")
  bad <- FALSE

  for (arg in names(types_to_check)) {
    if (types_to_check[[arg]] != types_provided[[arg]]) {
      bad <- TRUE
      msg <- paste0(msg," Argument ",arg,", has type ",types_provided[[arg]],", but should have type ",types_to_check[[arg]])
    }
  }
  if (bad) stop(msg)

  TRUE
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
  # b + log(exp(a-b) - 1)
  b + log(expm1(a-b))
}

logsumexpweights <- function(pp,ww) {
  logdiffexp(
    logsumexp(log(ww[ww>0]) + pp[ww>0]),
    logsumexp(log(-ww[ww<0]) + pp[ww<0])
  )
}

# Splice a vector
splice <- function(v,t,j) {
  # Insert t into v such that if vnew = splice(v,t,j) then vnew[j] == t
  if (j == 1) return(c(t,v))
  n <- length(v)
  if (j == (n+1)) return(c(v,t))
  if (j<=0) stop("j must be >0")
  if (j>(n+1)) stop("j must be <= n+1")
  c(v[1:(j-1)],t,v[j:n])
}

# invert positive def matrix using eigenvalues
# approximate by nearest pos def matrix if necessary
# resulting matrix wont produce errors with mvQuad::rescale
safeInverse = function(H, control=list(), ...) {
  if(max(abs(as.matrix(H) - t(as.matrix(H)))) > sqrt(.Machine$double.eps)) { 
    warning("matrix not symmetric")
    H2 = as.matrix(H)
    H = (H2 + t(H2))/2
  }

  Heigen = eigen(H, symmetric=TRUE)
  if(identical(control$verbose, TRUE)) {
    cat("hessian eigenvalues", paste(Heigen$vaules, collapse=', '), '\n')
  }

  if(!all(Heigen$values>0) ) {
    warning("negative eigenvalues in H, approxmiating with pracma::nearest_spd")
    if(requireNamespace("pracma")) {
      Hfix = pracma::nearest_spd(H)
      Heigen = eigen(Hfix, symmetric=TRUE)
    } else {
      warning("pracma package not available, taking absolute values of eigenvalues")
      Heigen$values = abs(Heigen$values)
    }

  }

  result = Matrix::forceSymmetric(
    Heigen$vectors %*%
      Matrix::Diagonal(x=1/Heigen$values) %*%
      t(Heigen$vectors)
  )

  attributes(result)$logDet = -sum(log(Heigen$values))
  result
}