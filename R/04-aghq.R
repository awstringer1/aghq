## Adaptive Gauss-Hermite Quadrature ##
# This file contains the main package function, aghq,
# and related functions.


#' Adaptive Gauss-Hermite Quadrature
#'
#' Normalize the log-posterior distribution using Adaptive Gauss-Hermite Quadrature.
#' This function takes in a function and its gradient and Hessian, and returns
#' a list of information about the normalized posterior, with methods for summarizing
#' and plotting.
#'
#' When \code{k = 1} the AGHQ method is a Laplace approximation, and you should use
#' the \code{aghq::laplace_approximation} function, since some of the methods for
#' \code{aghq} objects won't work with only one quadrature point. Objects of
#' class \code{laplace} have different methods suited to this case. See \code{?aghq::laplace_approximation}.
#'
#' @inheritParams optimize_theta
#' @inheritParams normalize_logpost
#' @param optresults Optional. A list of the results of the optimization of the log
#' posterior, formatted according to the output of \code{aghq::optimize_theta}. The
#' \code{aghq::aghq} function handles the optimization for you; passing this list
#' overrides this, and is useful for when you know your optimization is too difficult to be
#' handled by general-purpose software. See the software paper for several examples of this.
#' If you're unsure whether this option is needed for your problem then it probably is not.
#' @param basegrid Optional. Provide an object of class \code{NIGrid} from the \code{mvQuad}
#' package, representing the base quadrature rule that will be adapted. This is only
#' for users who want more complete control over the quadrature, and is not necessary
#' if you are fine with the default option which basically corresponds to
#' \code{mvQuad::createNIGrid(length(theta),'GHe',k,'product')}. **Note**: the \code{mvQuad}
#' functions used within \code{aghq} operate on grids in memory, so your \code{basegrid}
#' object will be changed after you run \code{aghq}.
#' @param transformation Optional. Do the quadrature for parameter \code{theta}, but
#' return summaries and plots for parameter \code{g(theta)}.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#'
#' @return An object of class \code{aghq} which is a list containing elements:
#' \itemize{
#' \item{normalized_posterior: }{The output of the \code{normalize_logpost} function, which
#' itself is a list with elements:
#' \itemize{
#' \item{\code{nodesandweights}: }{a dataframe containing the nodes and weights for the adaptive quadrature rule, with the un-normalized and normalized log posterior evaluated at the nodes.}
#' \item{\code{thegrid}: }{a \code{NIGrid} object from the \code{mvQuad} package, see \code{?mvQuad::createNIGrid}.}
#' \item{\code{lognormconst}: }{the actual result of the quadrature: the log of the normalizing constant of the posterior.}
#' }}
#' \item{marginals: }{a list of the same length as \code{startingvalue} of which element \code{j}
#' is the result of calling \code{aghq::marginal_posterior} with that \code{j}. This is
#' a tbl_df/tbl/data.frame containing the normalized log marginal posterior
#' for theta_j evaluated at the original quadrature points. Has columns
#' \code{"thetaj","logpost_normalized","weights"}, where \code{j} is the \code{j} you specified.
#' }
#' \item{optresults: }{information and results from the optimization of the log posterior, the
#' result of calling \code{aghq::optimize_theta}. This a list with elements:
#' \itemize{
#' \item{\code{ff}: }{the function list that was provided}
#' \item{\code{mode}: }{the mode of the log posterior}
#' \item{\code{hessian}: }{the hessian of the log posterior at the mode}
#' \item{\code{convergence}: }{specific to the optimizer used, a message indicating whether it converged}
#' }
#' }
#' \item{control: }{the control parameters passed}
#' }
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#'
#' @family quadrature
#'
#' @export
#'
aghq <- function(ff,k,startingvalue,transformation = default_transformation(),optresults = NULL,basegrid = NULL,control = default_control(),...) {

  validate_control(control)
  validate_transformation(transformation)
  transformation <- make_transformation(transformation)
  # ff <- make_moment_function(ff)
  # validate_moment(ff)

  # If they provided a basegrid, get the k from that. If they also provided a k, compare them and issue a warning
  if (!is.null(basegrid)) {
    if (missing(k)) {
      k <- max(as.numeric(basegrid$level))
    } else {
      k2 <- max(as.numeric(basegrid$level))
      if (k != k2) {
        warning(paste0("You provided a basegrid and a specified number of quadrature points k. You do not need to specify k if you supply a basegrid. Further, they don't match: your grid has k = ",k2,", but you specified k = ",k,". Proceeding with k = ",k2,", from the supplied grid.\n"))
        k <- k2
      }
    }
  }

  # Optimization
  if (is.null(optresults)) utils::capture.output(optresults <- optimize_theta(ff,startingvalue,control,...))

  # Normalization
  normalized_posterior <- normalize_logpost(optresults,k,basegrid = basegrid,ndConstruction = control$ndConstruction,...)

  if (control$onlynormconst) return(normalized_posterior$lognormconst)

  out <- list(
    normalized_posterior = normalized_posterior,
    # marginals = marginals,
    optresults = optresults,
    control = control,
    transformation = transformation
  )
  class(out) <- "aghq"
  # Marginals
  d <- length(startingvalue)
  marginals <- vector(mode = "list",length = d)
  if (control$method_summaries[1] == 'correct') {
    for (j in 1:d) marginals[[j]] <- marginal_posterior.aghq(out,j,method = 'correct')
  } else {
    for (j in 1:d) marginals[[j]] <- marginal_posterior.aghq(out,j,method = 'reuse')
  }
  out$marginals <- marginals
  out
}

#' Summary statistics computed using AGHQ
#'
#' The \code{summary.aghq} method computes means, standard deviations, and
#' quantiles of the transformed parameter.
#' The associated print method
#' prints these along with diagnostic and other information about
#' the quadrature.
#'
#' @param object The return value from \code{aghq::aghq}. Summaries are computed for
#' \code{object$transformation$fromtheta(theta)}.
#' @param ... not used.
#'
#' @return A list of class \code{aghqsummary}, which has a print method. Elements:
#' \itemize{
#' \item{mode: }{the mode of the log posterior}
#' \item{hessian: }{the hessian of the log posterior at the mode}
#' \item{covariance: }{the inverse of the hessian of the log posterior at the mode}
#' \item{cholesky: }{the upper Cholesky triangle of the hessian of the log posterior at the mode}
#' \item{quadpoints: }{the number of quadrature points used in each dimension}
#' \item{dim: }{the dimension of the parameter space}
#' \item{summarytable: }{a table containing the mean, median, mode, standard deviation
#' and quantiles of each transformed parameter, computed according to the posterior normalized
#' using AGHQ}
#' }
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#' # Summarize and automatically call its print() method when called interactively:
#' summary(thequadrature)
#' # or, compute the summary and save for further processing:
#' ss <- summary(thequadrature)
#' str(ss)
#'
#' @export
#'
summary.aghq <- function(object,...) {
  d <- get_param_dim(object)

  thetanames <- colnames(object$normalized_posterior$nodesandweights)[1:d]

  do_correct <- object$control$method_summaries[1] == 'correct'

  # Moments
  if (do_correct) {
    themeans <- compute_moment(object,nn = 1,type = 'raw',method = 'correct')
    thesds <- sqrt(compute_moment(object,nn = 2,type = 'central',method = 'correct'))
  } else {
    themeans <- compute_moment(object$normalized_posterior,object$transformation$fromtheta)
    thesds <- numeric(d)
    for (j in 1:d) thesds[j] <- sqrt(compute_moment(object$normalized_posterior,function(x) (object$transformation$fromtheta(x) - themeans[j])^2)[j])
    names(thesds) <- names(themeans) <- thetanames
  }


  themoments <- cbind(themeans,thesds)
  colnames(themoments) <- c('mean','sd')
  themoments <- as.data.frame(themoments)

  # Quantiles
  thequants <- vector(mode = 'list',length = d)
  for (j in 1:d) thequants[[j]] <- compute_quantiles(object$marginals[[j]],q = c(.025,.5,.975),transformation = object$transformation,interpolation = object$control$interpolation)

  names(thequants) <- thetanames
  thequants <- t(as.data.frame(thequants))
  colnames(thequants)[2] <- 'median'

  # thesummary <- cbind(themoments,thequants,data.frame(mode = object$optresults$mode))
  thesummary <- cbind(themoments,thequants)

  # thesummary <- thesummary[ ,c('mean','median','mode','sd','2.5%','97.5%')]
  thesummary <- thesummary[ ,c('mean','sd','2.5%','median','97.5%')]


  out <- list()
  class(out) <- "aghqsummary"
  out$mode <- object$optresults$mode
  out$hessian <- object$optresults$hessian
  out$lognormconst <- object$normalized_posterior$lognormconst
  out$covariance <- solve(out$hessian)
  # out$cholesky <- chol(out$covariance)
  out$quadpoints <- as.numeric(object$normalized_posterior$grid$level)
  out$dim <- length(out$quadpoints)
  out$summarytable <- thesummary

  out
}

#' Print method for AGHQ objects
#'
#' Pretty print the object-- just gives some basic information and then suggests
#' the user call \code{summary(...)}.
#'
#' @param x An object of class \code{aghq}.
#' @param ... not used.
#'
#' @return Silently prints summary information.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#' thequadrature
#'
#' @export
#'
print.aghq <- function(x,...) {
cat("aghq object. \n
Use summary(...) to see more detailed information,\n
plot(...) to see plots of marginal distributions, and \n
?compute_moment to see a list of useful summary methods.\n ")
}

#' Print method for AGHQ summary objects
#'
#' Print the summary of an \code{aghq} object. Almost always called by invoking
#' \code{summary(...)} interactively in the console.
#'
#' @param x The result of calling \code{summary(...)} on an object of class \code{aghq}.
#' @param ... not used.
#'
#' @return Silently prints summary information.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#' # Summarize and automatically call its print() method when called interactively:
#' summary(thequadrature)
#'
#' @family quadrature
#'
#' @export
#'
print.aghqsummary <- function(x,...) {
  cat("AGHQ on a",x$dim,"dimensional posterior with ",x$quadpoints,"quadrature points\n\n")
  cat("The posterior mode is:",x$mode,"\n\n")
  cat("The log of the normalizing constant/marginal likelihood is:",x$lognormconst,"\n\n")
  # cat("The posterior Hessian at the mode is:\n")
  # print(as.matrix(x$hessian))
  # cat("\n")
  cat("The covariance matrix used for the quadrature is...\n")
  print(as.matrix(x$covariance))
  # cat("\n")
  # cat("...and its Cholesky is:\n")
  # print(as.matrix(x$cholesky))
  cat("\n")
  cat("Here are some moments and quantiles for the transformed parameter:\n\n")
  print(x$summarytable)
  cat("\n")
}

#' Plot method for AGHQ objects
#'
#' Plot the marginal pdf and cdf of the transformed parameter from an \code{aghq} object.
#'
#' @param x The return value of \code{aghq::aghq}. Plots are created for the marginal
#' pdf and cdf of \code{x$transformation$fromtheta(theta)}.
#' @param ... not used.
#'
#' @return Silently plots.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#' plot(thequadrature)
#'
#' @family quadrature
#'
#' @export
#'
plot.aghq <- function(x,...) {
  d <- length(x$marginals)
  # Compute pdf and cdf
  pdfandcdf <- vector(mode = 'list',length = d)
  for (j in 1:d) pdfandcdf[[j]] <- compute_pdf_and_cdf(x$marginals[[j]],interpolation = x$control$interpolation,transformation = x$transformation)

  for (j in 1:d) {
    graphics::plot(pdfandcdf[[j]]$pdf ~ pdfandcdf[[j]]$transparam,
         type = 'l',
         main = paste0("Marg. Post., theta",j),
         xlab = paste0("theta",j),
         ylab = "Density")
    graphics::plot(pdfandcdf[[j]]$cdf ~ pdfandcdf[[j]]$transparam,
         type = 'l',
         main = paste0("Marg. CDF, theta",j),
         xlab = paste0("theta",j),
         ylab = "CDF")
  }
}

#' Laplace Approximation
#'
#' Wrapper function to implement a Laplace approximation to the posterior. A
#' Laplace approximation is AGHQ with \code{k = 1} quadrature points.
#' However, the returned
#' object is of a different class \code{laplace}, and a different summary
#' method is given for it. It is especially useful for high-dimensional problems where
#' the curse of dimensionality renders the use of \code{k > 1} quadrature points
#' infeasible. The summary method reflects the fact that the user may
#' be using this for a high-dimensional problem, and no plot method is given,
#' because there isn't anything
#' interesting to plot.
#'
#' @inheritParams aghq
#'
#' @return An object of class \code{laplace} with summary and plot methods. This
#' is simply a list with elements \code{lognormconst} containing the log of the
#' approximate normalizing constant, and \code{optresults} containing the optimization
#' results formatted the same way as \code{optimize_theta} and \code{aghq}.
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#'
#' @family quadrature
#'
#' @export
#'
laplace_approximation <- function(ff,startingvalue,optresults = NULL,control = default_control(),...) {

  validate_control(control)

  # Negate it if asked
  if (control$negate) {
    ffa <- list(
      fn = function(theta) -1 * ff$fn(theta),
      gr = function(theta) -1 * ff$gr(theta),
      he = function(theta) -1 * ff$he(theta)
    )
  } else {
    ffa <- ff
  }

  if(is.null(optresults)) utils::capture.output(optresults <- optimize_theta(ffa,startingvalue,control,...))
  lognorm <- normalize_logpost(optresults,1,...)
  out <- list(lognormconst = lognorm,optresults = optresults)
  class(out) <- "laplace"
  out
}

#' Print method for AGHQ objects
#'
#' Pretty print the object-- just gives some basic information and then suggests
#' the user call \code{summary(...)}.
#'
#' @param x An object of class \code{aghq}.
#' @param ... not used.
#'
#' @return Silently prints summary information.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thequadrature <- aghq(funlist2d,3,c(0,0))
#' thequadrature
#'
#' @export
#'
print.laplace <- function(x,...) {
cat("laplace object. \n
Use summary(...) to see more detailed information.\n
Further summaries are not available for this class; if you
actually wanted aghq with one quadrature point, use
aghq(...,k = 1). Note that some summaries may not
behave as expected in this case; you are recommended
to use k >= 3 if you want moments, quantiles, etc.\n")
}

#' Summary method for Laplace Approximation objects
#'
#' Summary method for objects of class \code{laplace}. Similar
#' to the method for objects of class \code{aghq}, but assumes the
#' problem is high-dimensional and does not compute or
#' print any large objects or summaries. See \code{summary.aghq} for
#' further information.
#'
#' @param object An object of class \code{laplace}.
#' @param ... not used.
#'
#' @return Silently prints summary information.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thelaplace <- laplace_approximation(funlist2d,c(0,0))
#' # Summarize and automatically call its print() method when called interactively:
#' summary(thelaplace)
#'
#' @family quadrature
#'
#' @export
summary.laplace <- function(object,...) {
  d <- length(object$optresults$mode)

  out <- list()
  class(out) <- "laplacesummary"
  out$mode <- object$optresults$mode
  out$lognormconst <- object$lognormconst
  out$dim <- d

  out
}

#' Print method for laplacesummary objects
#'
#' Print the summary of an \code{laplace} object. Almost always called by invoking
#' \code{summary(...)} interactively in the console.
#'
#' @param x The result of calling \code{summary(...)} on an object of class \code{laplace}.
#' @param ... not used.
#'
#' @return Silently prints summary information.
#'
#' @family quadrature
#'
#' @examples
#'
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#'
#' thelaplace <- laplace_approximation(funlist2d,c(0,0))
#' # Summarize and automatically call its print() method when called interactively:
#' summary(thelaplace)
#'
#' @family quadrature
#'
#' @export
#'
print.laplacesummary <- function(x,...) {
  cat("Laplace approximation for a",x$dim,"dimensional posterior\n\n")
  if (x$dim > 5) {
    cat("The posterior mode is:",c(round(x$mode[1:5],3),"..."),"\n\n")
  } else {
    cat("The posterior mode is:",round(x$mode,3),"\n\n")
  }
  cat("The log of the normalizing constant/marginal likelihood is:",x$lognormconst,"\n\n")
  cat("\n")
}

#' Marginal Laplace approximation
#'
#' Implement the marginal Laplace approximation of Tierney and Kadane (1986) for
#' finding the marginal posterior \code{(theta | Y)} from an unnormalized joint posterior
#' \code{(W,theta,Y)} where \code{W} is high dimensional and \code{theta} is low dimensional.
#' See the \code{AGHQ} software paper for a detailed example, or Stringer et. al. (2020).
#'
#' @param ff A function list similar to that required by \code{aghq}. However, each
#' function now takes arguments \code{W} and \code{theta}. Explicitly, this is a
#' list containing elements:
#' \itemize{
#' \item{\code{fn}}{: function taking arguments \code{W} and \code{theta} and returning a numeric
#' value representing the log-joint posterior at \code{W,theta}}
#' \item{\code{gr}}{: function taking arguments \code{W} and \code{theta} and returning a numeric
#' vector representing the gradient with respect to \code{W} of the log-joint posterior at \code{W,theta}}
#' \item{\code{he}}{: function taking arguments \code{W} and \code{theta} and returning a numeric
#' matrix representing the hessian with respect to \code{W} of the log-joint posterior at \code{W,theta}}
#' }
#' @param startingvalue A list with elements \code{W} and \code{theta}, which are numeric
#' vectors to start the optimizations for each variable. If you're using this method
#' then the log-joint posterior should be concave and these optimizations should not be
#' sensitive to starting values.
#' @param transformation Optional. Do the quadrature for parameter \code{theta}, but
#' return summaries and plots for parameter \code{g(theta)}. This applies to the \code{theta}
#' parameters only, not the \code{W} parameters.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#' @param control A list with elements
#' \itemize{
#' \item{\code{method}: }{optimization method to use for the \code{theta} optimization:
#' \itemize{
#' \item{'sparse_trust' (default): }{\code{trustOptim::trust.optim}}
#' \item{'sparse': }{\code{trust::trust}}
#' \item{'BFGS': }{\code{optim(...,method = "BFGS")}}
#' }
#' }
#' \item{\code{inner_method}: }{optimization method to use for the \code{W} optimization; same
#' options as for \code{method}}
#' }
#' Default \code{inner_method} is 'sparse_trust' and default \code{method} is 'BFGS'.
#'
#' @inheritParams aghq
#'
#' @return If \code{k > 1}, an object of class \code{marginallaplace}, which includes
#' the result of calling \code{aghq::aghq} on
#' the Laplace approximation of \code{(theta|Y)}, .... See software paper for full details.
#' If \code{k = 1}, an object of class \code{laplace} which is the result of calling
#' \code{aghq::laplace_approximation} on
#' the Laplace approximation of \code{(theta|Y)}.
#'
#' @family quadrature
#'
#' @examples
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' objfunc2dmarg <- function(W,theta) objfunc2d(c(W,theta))
#' objfunc2dmarggr <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::grad(fn,W)
#' }
#' objfunc2dmarghe <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::hessian(fn,W)
#' }
#'
#' funlist2dmarg <- list(
#'   fn = objfunc2dmarg,
#'   gr = objfunc2dmarggr,
#'   he = objfunc2dmarghe
#' )
#
# themarginallaplace <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0))
#' @export
#'
marginal_laplace <- function(ff,k,startingvalue,transformation = default_transformation(),optresults = NULL,control = default_control_marglaplace(),...) {

  validate_control(control,type = 'marglaplace')
  validate_transformation(transformation)
  transformation <- make_transformation(transformation)

  # Negate it if asked
  if (control$negate) {
    ffa <- list(
      fn = function(theta) -1 * ff$fn(theta),
      gr = function(theta) -1 * ff$gr(theta),
      he = function(theta) -1 * ff$he(theta)
    )
  } else {
    ffa <- ff
  }
  # Dimension of W space
  Wd <- length(startingvalue$W)

  # Mechanism for saving and reusing starting values
  # and hessians/modes/chols
  thetatable <- t(as.data.frame(startingvalue$theta))
  colnames(thetatable) <- paste0('theta',1:length(startingvalue$theta))
  rownames(thetatable) <- NULL

  Wlist <- list()
  Wlist[[1]] <- startingvalue$W

  Hlist <- list()
  Hlist[[1]] <- -1 * ffa$he(W = startingvalue$W,theta = startingvalue$theta)

  envtouse <- environment()

  get_thetaidx <- function(theta,whichenv = parent.frame()) {
    # Return the closest theta in the thetatable
    thetatable <- whichenv$thetatable
    if (is.null(thetatable)) return(0)
    if (nrow(thetatable) == 0) return(0)

    which.min(apply(thetatable,1,function(x) sum(abs(x - theta))))
  }

  theta_in_table <- function(theta,whichenv = parent.frame()) {
    # logical: is theta in the thetatable?
    thetatable <- whichenv$thetatable
    any(apply(thetatable,1,function(x) all(x == theta)))
  }

  get_Wstart <- function(theta,whichenv = parent.frame()) {
    whichenv$Wlist[[get_thetaidx(theta,whichenv)]]
  }

  get_H <- function(theta,whichenv = parent.frame()) {
    whichenv$Hlist[[get_thetaidx(theta,whichenv)]]
  }

  add_elements <- function(theta,W,H,whichenv = parent.frame()) {
    # Check if that theta already exists, if so, overwrite
    if (theta_in_table(theta,whichenv)) {
      Wlist <- whichenv$Wlist
      Wlist[[get_thetaidx(theta,whichenv)]] <- W
      assign("Wlist",Wlist,envir = whichenv)

      Hlist <- whichenv$Hlist
      Hlist[[get_thetaidx(theta,whichenv)]] <- H
      assign("Hlist",Hlist,envir = whichenv)
    } else {
      newrow <- t(as.data.frame(theta))
      colnames(newrow) <- paste0('theta',1:length(theta))
      rownames(newrow) <- NULL

      thetatable <- rbind(
        whichenv$thetatable,
        newrow
      )
      assign("thetatable",thetatable,envir = whichenv)

      Wlist <- whichenv$Wlist
      ll <- length(Wlist)
      Wlist[[ll+1]] <- W
      assign("Wlist",Wlist,envir = whichenv)

      Hlist <- whichenv$Hlist
      ll <- length(Hlist)
      Hlist[[ll+1]] <- H
      assign("Hlist",Hlist,envir = whichenv)
    }
  }

  # Write a function for fixed theta that computes the laplace approximation
  log_posterior_theta <- function(theta,whichenv = envtouse) {
    # cat('theta = ',theta,'\n')
    ffinner <- list(
      fn = function(W) ffa$fn(W,theta),
      gr = function(W) ffa$gr(W,theta),
      he = function(W) ffa$he(W,theta)
    )

    # If theta is in the table already then create the optresults
    Wstart <- get_Wstart(theta,whichenv = whichenv)
    optresults <- NULL
    # if (theta_in_table(theta,whichenv)) {
    #   H <- get_H(theta,whichenv)
    #   optresults <- list(
    #     ff = ffinner,
    #     mode = Wstart,
    #     hessian = H,
    #     convergence = 0
    #   )
    # }

    utils::capture.output(lap <- laplace_approximation(ffinner,Wstart,optresults = optresults,control = default_control(method = control$inner_method,negate = FALSE)))

    add_elements(theta,lap$optresults$mode,lap$optresults$hessian,whichenv = whichenv)

    as.numeric(lap$lognormconst)
  }

  ## Do the quadrature ##
  # Create the function list
  # TODO: better optimization like they do for GAMs
  ffouter <- list(
    fn = log_posterior_theta,
    gr = function(theta) numDeriv::grad(log_posterior_theta,theta,method = 'simple'),
    he = function(theta) numDeriv::hessian(log_posterior_theta,theta)
  )
  # If requesting an "outer" Laplace approximation, return it
  if (k == 1) return(aghq::laplace_approximation(ffouter,startingvalue$theta,control = control))
  # Do the quadrature manually, so we can save intermediate results
  utils::capture.output(outeropt <- aghq::optimize_theta(ffouter,startingvalue$theta,control,...))

  # Create the grid
  S <- length(outeropt$mode) # Dimension
  thegrid <- mvQuad::createNIGrid(dim = S,type = "GHe",level = k)
  m <- outeropt$mode
  H <- outeropt$hessian
  mvQuad::rescale(thegrid,m = m,C = Matrix::forceSymmetric(solve(H)),dec.type=2)

  thetaorder <- paste0('theta',1:S)

  nodesandweights <- cbind(mvQuad::getNodes(thegrid),mvQuad::getWeights(thegrid))
  colnames(nodesandweights) <- c(thetaorder,"weights")
  nodesandweights <- as.data.frame(nodesandweights)

  # Compute the log-posterior at the integration points, saving the summaries
  if (length(thetaorder) == 1) {
    distinctthetas <- data.frame(theta1 = nodesandweights[ ,thetaorder])
  } else {
    distinctthetas <- nodesandweights[ ,thetaorder]
  }
  lp <- numeric(nrow(distinctthetas))
  # modesandhessians <- dplyr::as_tibble(distinctthetas) %>%
  #   tibble::add_column(
  #     mode = vector(mode = 'list',length = nrow(distinctthetas)),
  #     H = vector(mode = 'list',length = nrow(distinctthetas)),
  #     logpost = numeric(nrow(distinctthetas))
  #   )

  modesandhessians <- distinctthetas
  modesandhessians$mode <- vector(mode = 'list',length = nrow(distinctthetas))
  modesandhessians$H <- vector(mode = 'list',length = nrow(distinctthetas))
  modesandhessians$logpost <- numeric(nrow(distinctthetas))

  for (i in 1:length(lp)) {
    theta <- as.numeric(distinctthetas[i,thetaorder])
    # Re-use the starting values
    Wstart <- get_Wstart(theta,envtouse)
    # Do the Laplace approx
    ffinner <- list(
      fn = function(W) ffa$fn(W,theta),
      gr = function(W) ffa$gr(W,theta),
      he = function(W) ffa$he(W,theta)
    )

    utils::capture.output(lap <- laplace_approximation(ffinner,Wstart,control = default_control(method = control$inner_method,negate=FALSE)))
    mtmp <- lap$optresults$mode
    names(mtmp) <- paste0("W",1:length(mtmp))
    modesandhessians[i,'mode'] <- list(list(mtmp))
    modesandhessians[i,'H'] <- list(list(lap$optresults$hessian))
    modesandhessians[i,'logpost'] <- as.numeric(lap$lognormconst)

    lp[i] <- as.numeric(lap$lognormconst)
  }

  # Get the normalization constant
  ww <- nodesandweights$weights

  lognormconst <- logsumexp(log(ww) + lp)
  if (control$onlynormconst) return(lognormconst)

  nodesandweights$logpost <- lp
  nodesandweights$logpost_normalized <- lp - lognormconst

  normalized_posterior <- list(
    nodesandweights = nodesandweights,
    grid = thegrid,
    lognormconst = lognormconst
  )
  # Do the rest of the aghq
  d <- length(startingvalue$theta)
  marginals <- vector(mode = "list",length = d)
  for (j in 1:d) marginals[[j]] <- marginal_posterior(outeropt,k,j)

  out <- list(
    normalized_posterior = normalized_posterior,
    marginals = marginals,
    optresults = outeropt,
    modesandhessians = modesandhessians,
    control = control,
    transformation = transformation
  )
  class(out) <- c("marginallaplace","aghq")
  out
}

#' AGHQ-normalized marginal Laplace approximation from a TMB function template
#'
#' Implement the algorithm from \code{aghq::marginal_laplace()}, but making use of
#' \code{TMB}'s automatic Laplace approximation. This function takes a function
#' list from \code{TMB::MakeADFun()} with a non-empty set of \code{random} parameters,
#' in which the \code{fn} and \code{gr} are the unnormalized marginal Laplace
#' approximation and its gradient. It then calls \code{aghq::aghq()} and formats
#' the resulting object so that its contents and class match the output of
#' \code{aghq::marginal_laplace()} and are hence suitable for post-processing
#' with \code{summary}, \code{aghq::sample_marginal()}, and so on.
#'
#' @details Because \code{TMB} does not yet have the Hessian of the log marginal Laplace
#' approximation implemented, a numerically-differentiated jacobian of the gradient
#' is used via \code{numDeriv::jacobian()}. You can turn this off (using \code{ff$he()} instead,
#' which you'll have to modify yourself) using \code{default_control_tmb(numhessian = FALSE)}.
#'
#' @param ff The output of calling \code{TMB::MakeADFun()} with \code{random} set
#' to a non-empty subset of the parameters. **VERY IMPORTANT**: \code{TMB}'s
#' automatic Laplace approximation requires you to write your template implementing
#' the **negated** log-posterior. Therefore, this list that you input here
#' will contain components \code{fn}, \code{gr} and \code{he} that implement the
#' **negated** log-posterior and its derivatives. This is **opposite**
#' to every other comparable function in the \code{aghq} package, and is done
#' here to emphasize compatibility with \code{TMB}.
#' @param transformation Optional. Do the quadrature for parameter \code{theta}, but
#' return summaries and plots for parameter \code{g(theta)}. This applies to the \code{theta}
#' parameters only, not the \code{W} parameters.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#' @param control A list of control parameters. See \code{?default_control} for details. Valid options are:
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
#' \item \code{negate}: default \code{TRUE}. See \code{?default_control_tmb}. Assumes that your \code{TMB} function
#' template computes the **negated** log-posterior, which it must if you're using \code{TMB}'s automatic
#' Laplace approximation, which you must be if you're using this function!}.
#'
#' @inheritParams aghq
#'
#' @return If \code{k > 1}, an object of class \code{marginallaplace}
#' (and inheriting from class \code{aghq}) of the same
#' structure as that returned by \code{aghq::marginal_laplace()}, with \code{plot}
#' and \code{summary} methods, and suitable for input into \code{aghq::sample_marginal()}
#' for drawing posterior samples.
#'
#' @family quadrature
#'
#' @export
#'
marginal_laplace_tmb <- function(ff,k,startingvalue,transformation = default_transformation(),optresults = NULL,basegrid = NULL,control = default_control_tmb(),...) {

  validate_control(control,type='tmb')
  validate_transformation(transformation)
  transformation <- make_transformation(transformation)

  # Get names from TMB function template
  thetanames <- NULL
  if (exists('par',ff)) thetanames <- make.unique(names(ff$par),sep='')

  # Hessian
  if (control$numhessian) {
    ff$he <- function(theta) numDeriv::jacobian(ff$gr,theta,method = 'Richardson')
  }
  ## Do aghq ##
  # The aghq
  quad <- aghq(ff = ff,k = k,transformation = transformation,startingvalue = startingvalue,optresults = optresults,basegrid = basegrid,control = control,...)
  if (control$onlynormconst) return(quad) # NOTE: control was passed to aghq here so quad should itself just be a number

  ## Add on the info needed for marginal Laplace ##
  # Add on the quad point modes and curvatures
  distinctthetas <- quad$normalized_posterior$nodesandweights[ ,grep('theta',colnames(quad$normalized_posterior$nodesandweights))]
  if (!is.data.frame(distinctthetas)) distinctthetas <- data.frame(theta1 = distinctthetas)

  modesandhessians <- distinctthetas
  if (is.null(thetanames)) {
    thetanames <- colnames(distinctthetas)
  } else {
    colnames(modesandhessians) <- thetanames
    colnames(quad$normalized_posterior$nodesandweights)[grep('theta',colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
  }
  modesandhessians$mode <- vector(mode = 'list',length = nrow(distinctthetas))
  modesandhessians$H <- vector(mode = 'list',length = nrow(distinctthetas))

  # if (is.null(thetanames)) {
  #   thetanames <- colnames(distinctthetas)
  # } else {
  #   colnames(modesandhessians)[colnames(modesandhessians) == colnames(distinctthetas)] <- thetanames
  #   colnames(quad$normalized_posterior$nodesandweights)[grep('theta',colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
  # }

    for (i in 1:nrow(distinctthetas)) {
      # Get the theta
      theta <- as.numeric(modesandhessians[i,thetanames])
      # Set up the mode and hessian of the random effects. This happens when you run
      # the TMB objective with a particular theta
      ff$fn(theta)
      # Now pull the mode and hessian. Have to be careful about scoping
      mm <- ff$env$last.par
      modesandhessians[i,'mode'] <- list(list(mm[ff$env$random]))
      H <- ff$env$spHess(mm,random = TRUE)
      H <- rlang::duplicate(H) # Somehow, TMB puts all evals of spHess in the same memory location.
      modesandhessians[i,'H'] <- list(list(H))
    }

  quad$modesandhessians <- modesandhessians

  class(quad) <- c("marginallaplace","aghq")
  quad
}

#' Summary statistics for models using marginal Laplace approximations
#'
#' The \code{summary.marginallaplace} calls \code{summary.aghq}, but also computes
#' summary statistics of the random effects, by drawing from their approximate
#' posterior using \code{aghq::sample_marginal} with the specified number
#' of samples.
#'
#' @param object Object inheriting from **both** classes \code{aghq} and \code{marginallaplace},
#' for example as returned by \code{aghq::marginal_laplace} or \code{aghq::marginal_laplace_tmb}.
#' @param M Number of samples to use to compute summary statistics of the random effects.
#' Default \code{1000}. Lower runs faster, higher is more accurate.
#' @param max_print Sometimes there are a lot of random effects. If there are more random
#' effects than \code{max_print}, the random effects aren't summarized, and a note is printed
#' to this effect. Default \code{30}.
#' @param ... not used.
#'
#' @return A list containing an object of class \code{aghqsummary} (see \code{summary.aghq}).
#'
#' @family quadrature
#'
#' @examples
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' objfunc2dmarg <- function(W,theta) objfunc2d(c(W,theta))
#' objfunc2dmarggr <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::grad(fn,W)
#' }
#' objfunc2dmarghe <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::hessian(fn,W)
#' }
#'
#' funlist2dmarg <- list(
#'   fn = objfunc2dmarg,
#'   gr = objfunc2dmarggr,
#'   he = objfunc2dmarghe
#' )
#'
#' themarginallaplace <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0))
#' summary(themarginallaplace)
#' @export
#'
summary.marginallaplace <- function(object,M=1e03,max_print=30,...) {

  p <- length(object$modesandhessians$mode[[1]])
  summ <- summary.aghq(object,...)
  if (p > max_print) {
    cat(paste0("There are ",p," random effects, but max_print = ",max_print,", so not computing their summary information.\nSet max_print higher than ",p," if you would like to summarize the random effects.\n"))
    return(summ)
  }

  samps <- aghq::sample_marginal(object,M,...)
  means <- apply(samps$samps,1,mean)
  medians <- apply(samps$samps,1,stats::median)
  sds <- apply(samps$samps,1,stats::sd)
  quants <- t(apply(samps$samps,1,stats::quantile,probs = c(.025,.975)))

  modes <- with(object,mapply(modesandhessians$mode,exp(normalized_posterior$nodesandweights$logpost_normalized)*normalized_posterior$nodesandweights$weights,FUN = '*'))
  if (is.array(modes)) {
    modes <- apply(modes,1,sum)
  } else {
    modes <- sum(modes)
  }

  randomeffectsummary <- data.frame(
    mean = means,
    median = medians,
    mode = modes,
    sd = sds,
    `2.5%` = quants[ ,1],
    `97.5%` = quants[ ,2]
  )
  colnames(randomeffectsummary) <- colnames(summ$summarytable)
  rownames(randomeffectsummary) <- NULL

  randomeffectsummary$variable <- names(object$modesandhessians$mode[[1]])

  out <- list(
    aghqsummary = summ,
    randomeffectsummary = randomeffectsummary,
    info = c("M" = M,"max_print" = max_print)
  )
  class(out) <- "marginallaplacesummary"
  out
}

#' Summary statistics for models using marginal Laplace approximations
#'
#' The \code{summary.marginallaplace} calls \code{summary.aghq}, but also computes
#' summary statistics of the random effects, by drawing from their approximate
#' posterior using \code{aghq::sample_marginal} with the specified number
#' of samples.
#'
#' @param x Object of class \code{marginallaplacesummary} returned by calling
#' \code{summary} on an object of class \code{marginallaplace}.
#' @param ... not used.
#'
#' @return Nothing. Prints contents.
#'
#' @family quadrature
#'
#' @examples
#' logfteta2d <- function(eta,y) {
#'   # eta is now (eta1,eta2)
#'   # y is now (y1,y2)
#'   n <- length(y)
#'   n1 <- ceiling(n/2)
#'   n2 <- floor(n/2)
#'   y1 <- y[1:n1]
#'   y2 <- y[(n1+1):(n1+n2)]
#'   eta1 <- eta[1]
#'   eta2 <- eta[2]
#'   sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
#'     sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
#' }
#' set.seed(84343124)
#' n1 <- 5
#' n2 <- 5
#' n <- n1+n2
#' y1 <- rpois(n1,5)
#' y2 <- rpois(n2,5)
#
#' objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
#' objfunc2dmarg <- function(W,theta) objfunc2d(c(W,theta))
#' objfunc2dmarggr <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::grad(fn,W)
#' }
#' objfunc2dmarghe <- function(W,theta) {
#'   fn <- function(W) objfunc2dmarg(W,theta)
#'   numDeriv::hessian(fn,W)
#' }
#'
#' funlist2dmarg <- list(
#'   fn = objfunc2dmarg,
#'   gr = objfunc2dmarggr,
#'   he = objfunc2dmarghe
#' )
#'
#' themarginallaplace <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0))
#' summary(themarginallaplace)
#' @export
#'
print.marginallaplacesummary <- function(x,...) {
  cat("\n==========================================================\n\n")
  cat(paste0("Fixed effects:\n"))
  print(x$aghqsummary)
  cat("==========================================================\n\n")
  cat(paste0("Random effects, based on ",x$info['M']," approximate posterior samples:\n"))
  print(x$randomeffectsummary)
}