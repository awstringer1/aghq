### Normalization ###
# This script contains a function for normalizing
# the joint posterior using AGHQ

#' Normalize the joint posterior using AGHQ
#'
#' This function takes in the optimization results from \code{aghq::optimize_theta()}
#' and returns a list with the quadrature points, weights, and normalization
#' information. Like \code{aghq::optimize_theta()}, this is designed for use only within
#' \code{aghq::aghq}, but is exported for debugging and documented in case you want to
#' modify it somehow, or something.
#'
#' @param optresults The results of calling \code{aghq::optimize_theta()}: see return value of that function.
#' @param k Integer, the number of quadrature points to use. I suggest at least 3. k = 1 corresponds to a Laplace
#' approximation.
#' @param whichfirst Integer between 1 and the dimension of the parameter space, default 1.
#' The user shouldn't have to worry about this: it's used internally to re-order the parameter vector
#' before doing the quadrature, which is useful when calculating marginal posteriors.
#'
#' @return A list with elements:
#' \itemize{
#' \item{\code{nodesandweights}: }{a dataframe containing the nodes and weights for the adaptive quadrature rule, with the un-normalized and normalized log posterior evaluated at the nodes.}
#' \item{\code{thegrid}: }{a \code{NIGrid} object from the \code{mvQuad} package, see \code{?mvQuad::createNIGrid}.}
#' \item{\code{lognormconst}: }{the actual result of the quadrature: the log of the normalizing constant of the posterior.}
#' }
#'
#' @family quadrature
#'
#' @examples
#' # Same setup as optimize_theta
#' logfteta <- function(eta,y) {
#'   sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
#' }
#' set.seed(84343124)
#' y <- rpois(10,5) # Mode should be sum(y) / (10 + 1)
#' truemode <- log((sum(y) + 1)/(length(y) + 1))
#
#' objfunc <- function(x) logfteta(x,y)
#' funlist <- list(
#'   fn = objfunc,
#'   gr = function(x) numDeriv::grad(objfunc,x),
#'   he = function(x) numDeriv::hessian(objfunc,x)
#' )
#' opt_sparsetrust <- optimize_theta(funlist,1.5)
#' opt_trust <- optimize_theta(funlist,1.5,control = list(method = "trust"))
#' opt_bfgs <- optimize_theta(funlist,1.5,control = list(method = "BFGS"))
#'
#' # Quadrature with 3, 5, and 7 points using sparse trust region optimization:
#' norm_sparse_3 <- normalize_logpost(opt_sparsetrust,3,1)
#' norm_sparse_5 <- normalize_logpost(opt_sparsetrust,5,1)
#' norm_sparse_7 <- normalize_logpost(opt_sparsetrust,7,1)
#'
#' # Quadrature with 3, 5, and 7 points using dense trust region optimization:
#' norm_trust_3 <- normalize_logpost(opt_trust,3,1)
#' norm_trust_5 <- normalize_logpost(opt_trust,5,1)
#' norm_trust_7 <- normalize_logpost(opt_trust,7,1)
#'
#' # Quadrature with 3, 5, and 7 points using BFGS optimization:
#' norm_bfgs_3 <- normalize_logpost(opt_bfgs,3,1)
#' norm_bfgs_5 <- normalize_logpost(opt_bfgs,5,1)
#' norm_bfgs_7 <- normalize_logpost(opt_bfgs,7,1)
#'
#' @export
#'
normalize_logpost <- function(optresults,k,whichfirst = 1) {
  if (as.integer(k) != k) stop(paste0("Please provide an integer k, the number of quadrature points. You provided ",k,"which does not satisfy as.integer(k) == k"))
  # Create the grid
  S <- length(optresults$mode) # Dimension
  thegrid <- mvQuad::createNIGrid(dim = S,type = "GHe",level = k)
  # Reorder the mode and Hessian so that "whichfirst" is first
  # This does not change the normalizing constant of the joint,
  # but is necessary to compute marginals later.
  idxorder <- c(whichfirst,(1:S)[-whichfirst])
  m <- optresults$mode[idxorder]
  H <- optresults$hessian[idxorder,idxorder]
  mvQuad::rescale(thegrid,m = m,C = Matrix::forceSymmetric(solve(H)),dec.type=2) # forceSymmetric for numerical asymmetries

  nodesandweights <- cbind(mvQuad::getNodes(thegrid),mvQuad::getWeights(thegrid))
  colnames(nodesandweights) <- c(stringr::str_c("theta",idxorder),"weights")
  nodesandweights <- as.data.frame(nodesandweights)

  # Compute the log-posterior at the integration points
  thetaorder <- paste0('theta',1:S)
  if (length(idxorder) == 1) {
    nodesandweights$logpost <- sapply(nodesandweights[ ,thetaorder],optresults$ff$fn)
  } else{
    nodesandweights$logpost <- apply(nodesandweights[ ,thetaorder],1,optresults$ff$fn)
  }

  # Get the normalization constant
  ww <- nodesandweights$weights
  pp <- nodesandweights$logpost

  lognormconst <- matrixStats::logSumExp(log(ww) + pp)
  nodesandweights$logpost_normalized <- nodesandweights$logpost - lognormconst

  list(
    nodesandweights = nodesandweights,
    grid = thegrid,
    lognormconst = lognormconst
  )
}


#' Normalize the joint posterior using a Laplace Approximation
#'
#' A Laplace approximation is AGHQ with k = 1 quadrature point. This function
#' is simply a wrapper for \code{aghq::normalize_logpost(...,k=1)}, for users
#' who are more comfortable with Laplace approximations than with the AGHQ methodology.
#' But it's the same thing :)
#'
#' @inherit normalize_logpost params return examples
#'
#' @export
#'
laplace_approximation <- function(optresults) normalize_logpost(optresults,1,1)