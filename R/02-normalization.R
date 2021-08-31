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
#' @param basegrid Optional. Provide an object of class \code{NIGrid} from the \code{mvQuad}
#' package, representing the base quadrature rule that will be adapted. This is only
#' for users who want more complete control over the quadrature, and is not necessary
#' if you are fine with the default option which basically corresponds to
#' \code{mvQuad::createNIGrid(length(theta),'GHe',k,'product')}.
#' @param ndConstruction Create a multivariate grid using a product or sparse construction?
#' Passed directly to \code{mvQuad::createNIGrid()}, see that function for further details. Note
#' that the use of sparse grids within \code{aghq} is currently **experimental** and not supported
#' by tests. In particular, calculation of marginal posteriors is known to fail currently.
#' @param ... Additional arguments to be passed to \code{optresults$ff}, see \code{?optimize_theta}.
#'
#' @return If k > 1, a list with elements:
#' \itemize{
#' \item{\code{nodesandweights}: }{a dataframe containing the nodes and weights for the adaptive quadrature rule, with the un-normalized and normalized log posterior evaluated at the nodes.}
#' \item{\code{thegrid}: }{a \code{NIGrid} object from the \code{mvQuad} package, see \code{?mvQuad::createNIGrid}.}
#' \item{\code{lognormconst}: }{the actual result of the quadrature: the log of the normalizing constant of the posterior.}
#' }
#'
#' If k = 1, then the method returns
#' a numeric value representing the log of the normalizing constant computed using
#' a Laplace approximation.
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
#' @importFrom Matrix determinant
#'
#' @export
#'
normalize_logpost <- function(optresults,k,whichfirst = 1,basegrid = NULL,ndConstruction = "product",...) {
  if (as.integer(k) != k) stop(paste0("Please provide an integer k, the number of quadrature points. You provided ",k,"which does not satisfy as.integer(k) == k"))
  if (k == 1) {
    # Laplace approx: just return the normalizing constant
    return(optresults$ff$fn(optresults$mode,...) - as.numeric(.5 * determinant(optresults$hessian,logarithm = TRUE)$modulus) + .5*dim(optresults$hessian)[1]*log(2*pi))
  }
  # Create the grid
  S <- length(optresults$mode) # Dimension
  if (!is.null(basegrid)) {
    thegrid <- basegrid
    # Check
    if (thegrid$dim != S) stop(paste0("Your startingvalue has dimension ",S,", but the grid you supplied has dimension ",thegrid$dim))
    if (!all(thegrid$features$initial.domain %in% c(-Inf,Inf))) stop("When supplying your own basegrid, you still have to choose a rule corresponding to a Gaussian kernel for the method to make sense. You chose a rule with initial domain not equal to (-Inf,Inf). Check the mvQuad::createNIGrid documentation for a list of available rules which have domain of integration (-Inf,Inf).")
  } else {
    thegrid <- mvQuad::createNIGrid(dim = S,type = "GHe",level = k,ndConstruction = ndConstruction,...)
  }
  # Reorder the mode and Hessian so that "whichfirst" is first
  # This does not change the normalizing constant of the joint,
  # but is necessary to compute marginals later.
  idxorder <- c(whichfirst,(1:S)[-whichfirst])
  m <- optresults$mode[idxorder]
  H <- optresults$hessian[idxorder,idxorder]
  mvQuad::rescale(thegrid,m = m,C = Matrix::forceSymmetric(solve(H)),dec.type=2) # forceSymmetric for numerical asymmetries

  nodesandweights <- cbind(mvQuad::getNodes(thegrid),mvQuad::getWeights(thegrid))
  colnames(nodesandweights) <- c(paste0("theta",idxorder),"weights")
  nodesandweights <- as.data.frame(nodesandweights)

  # Compute the log-posterior at the integration points
  thetaorder <- paste0('theta',1:S)
  if (length(idxorder) == 1) {
    nodesandweights$logpost <- sapply(nodesandweights[ ,thetaorder],optresults$ff$fn,...)
  } else{
    nodesandweights$logpost <- apply(nodesandweights[ ,thetaorder],1,optresults$ff$fn,...)
  }

  # Get the normalization constant
  ww <- nodesandweights$weights
  pp <- nodesandweights$logpost

  # lognormconst <- logsumexp(log(ww) + pp)
  # Account for negative weights (doesn't happen with GHQ but happens for e.g. sparse rules)
  lognormconst <- logdiffexp(
    logsumexp(log(ww[ww>0]) + pp[ww>0]),
    logsumexp(log(-ww[ww<0]) + pp[ww<0])
  )

  nodesandweights$logpost_normalized <- nodesandweights$logpost - lognormconst

  list(
    nodesandweights = nodesandweights,
    grid = thegrid,
    lognormconst = lognormconst
  )
}

#' Obtain the log-normalizing constant from a fitted quadrature object
#'
#' Quick helper S3 method to retrieve the log normalizing constant from an object
#' created using the aghq package. Methods for a list (returned by \code{aghq::normalize_posterior})
#' and for objects of class \code{aghq}, \code{laplace}, and \code{marginallaplace}.
#'
#' @param obj A list returned by \code{aghq::normalize_posterior} or an object of class \code{aghq}, \code{laplace}, or \code{marginallaplace}.
#' @param ... Not used
#'
#' @return A number representing the natural logarithm of the approximated normalizing constant.
#'
#' @family quadrature
#'
#' @export
#'
get_log_normconst <- function(obj,...) UseMethod("get_log_normconst")
#' @rdname get_log_normconst
#' @export
get_log_normconst.default <- function(obj,...) obj$lognormconst
#' @rdname get_log_normconst
#' @export
get_log_normconst.aghq <- function(obj,...) get_log_normconst(obj$normalized_posterior)
#' @rdname get_log_normconst
#' @export
get_log_normconst.laplace <- function(obj,...) obj$lognormconst
#' @rdname get_log_normconst
#' @export
get_log_normconst.marginallaplace <- function(obj,...) get_log_normconst(obj$normalized_posterior)

