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
#' @inheritParams optimize_theta
#' @inheritParams normalize_logpost
#' @param optresults Optional. A list of the results of the optimization of the log
#' posterior, formatted according to the output of \code{aghq::optimize_theta}. The
#' \code{aghq::aghq} function handles the optimization for you; passing this list
#' overrides this, and is useful for when you know your optimization is too difficult to be
#' handled by general-purpose software. See the software paper for several examples of this.
#' If you're unsure whether this option is needed for your problem then it probably is not.
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
aghq <- function(ff,k,startingvalue,optresults = NULL,control = default_control()) {
  # Optimization
  if (is.null(optresults)) optresults <- optimize_theta(ff,startingvalue,control)

  # Normalization
  normalized_posterior <- normalize_logpost(optresults,k)

  # Marginals
  d <- length(startingvalue)
  marginals <- vector(mode = "list",length = d)
  for (j in 1:d) marginals[[j]] <- marginal_posterior(optresults,k,j)

  out <- list(
    normalized_posterior = normalized_posterior,
    marginals = marginals,
    optresults = optresults
  )
  class(out) <- "aghq"
  out
}

#' Summary statistics computed using AGHQ
#'
#' The \code{summary.aghq} method computes means, standard deviations, and
#' quantiles and the associated print method
#' prints these along with diagnostic and other information about
#' the quadrature.
#'
#'
#' @export
#'
summary.aghq <- function(aghqobj) {
  d <- length(aghqobj$optresults$mode)

  # Moments
  themeans <- compute_moment(aghqobj$normalized_posterior,function(x) x)
  thesds <- numeric(d)
  for (j in 1:d) thesds[j] <- sqrt(compute_moment(aghqobj$normalized_posterior,function(x) (x - themeans[j])^2)[j])
  names(thesds) <- names(themeans)

  themoments <- cbind(themeans,thesds)
  colnames(themoments) <- c('mean','sd')
  themoments <- as.data.frame(themoments)

  # Quantiles
  thequants <- vector(mode = 'list',length = d)
  for (j in 1:d) thequants[[j]] <- compute_quantiles(aghqobj$marginals[[j]],c(.025,.5,.975))
  names(thequants) <- paste0("theta",1:d)
  thequants <- t(as.data.frame(thequants))
  colnames(thequants)[2] <- 'median'

  thesummary <- cbind(themoments,thequants)

  out <- list()
  class(out) <- "aghqsummary"
  out$mode <- aghqobj$optresults$mode
  out$hessian <- aghqobj$optresults$hessian
  out$covariance <- solve(out$hessian)
  out$cholesky <- chol(out$covariance)
  out$quadpoints <- as.numeric(aghqobj$normalized_posterior$grid$level)
  out$dim <- length(out$quadpoints)
  out$summarytable <- thesummary

  out
}

print.aghqsummary <- function(summ) {
  cat("AGHQ on a",summ$dim,"dimensional posterior with ",summ$quadpoints,"quadrature points\n\n")
  cat("The posterior mode is:",summ$mode,"\n\n")
  cat("The posterior Hessian at the mode is:\n")
  print(as.matrix(summ$hessian))
  cat("\n")
  cat("The covariance matrix used for the quadrature is...\n")
  print(as.matrix(summ$covariance))
  cat("\n")
  cat("...and its Cholesky is:\n")
  print(as.matrix(summ$cholesky))
  cat("\n")
  cat("Here are some moments and quantiles for theta:\n\n")
  print(summ$summarytable)
  cat("\n")
}


plot.aghq <- function(aghqobj) {
  d <- length(aghqobj$marginals)
  # Compute pdf and cdf
  pdfandcdf <- vector(mode = 'list',length = d)
  for (j in 1:d) pdfandcdf[[j]] <- compute_pdf_and_cdf(aghqobj$marginals[[j]])

  par(mfrow = c(d,2))
  for (j in 1:d) {
    plot(pdfandcdf[[j]]$pdf ~ pdfandcdf[[j]]$theta,
         type = 'l',
         main = paste0("Marg. Post., theta",j),
         xlab = paste0("theta",j),
         ylab = "Density")
    plot(pdfandcdf[[j]]$cdf ~ pdfandcdf[[j]]$theta,
         type = 'l',
         main = paste0("Marg. CDF, theta",j),
         xlab = paste0("theta",j),
         ylab = "CDF")
  }
}