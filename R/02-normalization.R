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
#' opt_trust <- optimize_theta(funlist,1.5,control = default_control(method = "trust"))
#' opt_bfgs <- optimize_theta(funlist,1.5,control = default_control(method = "BFGS"))
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
    # thegrid <- basegrid
    # This seems to be the only way to use basegrid without modifying it outside the function, due to how the mvQuad package implements this.
    thegrid <- with(basegrid,mvQuad::createNIGrid(dim = dim,type = type,level = level,ndConstruction = ndConstruction))
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

#' Nested, sparse Gaussian quadrature in AGHQ
#'
#' Compute a whole sequence of normalizing constants
#' for \code{1,3,5,...,k} points,
#' using only the function evaluations from the \code{k}-point rule.
#'
#' @param optresults The results of calling \code{aghq::optimize_theta()}: see return value of that function.
#' The dimension of the parameter \code{p} will be taken from \code{optresults$mode}.
#' @param p Dimension of the variable of integration.
#' @param k Integer, the number of quadrature points to use.
#' @param ndConstruction Create a multivariate grid using a product or sparse construction?
#' Passed directly to \code{mvQuad::createNIGrid()}, see that function for further details.
#' @param ... Additional arguments to be passed to \code{optresults$ff}, see \code{?optimize_theta}.
#'
#' @return For \code{get_quadtable}, a pre-computed table of nodes for the \code{k}-point rule,
#' with weights for the points from each of the \code{1,3,...,k}-point rules, for passing to
#' \code{nested_quadrature}. For \code{nested_quadrature} and \code{adaptive_nested_quadrature}, a named numeric vector of \code{optresults$fn}
#' values for each \code{k}.
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
#'
#' @rdname nested
#' @export
#'
nested_quadrature <- function(optresults,k,ndConstruction = 'product',...) {
  # Compute the quadrature rule at 1,3,5,...,k points using only
  # the function evaluations required to compute it at k points.
  p <- length(optresults$mode)
  quadtable <- get_quadtable(p,k,ndConstruction)
  # TODO: speed and memory comparisons/optimization. Or just store the tables...
  # quadtable <- as(as.matrix(quadtable),'sparseMatrix')
  quadtable[is.na(quadtable)] <- 0
  colSums(quadtable[grep('w',colnames(quadtable))]*apply(quadtable[grep('V',colnames(quadtable))],1,optresults$fn,...))
}
#' @rdname nested
#' @export
adaptive_nested_quadrature <- function(optresults,k,ndConstruction = 'product',...) {
  # Compute the quadrature rule at 1,3,5,...,k points using only
  # the function evaluations required to compute it at k points.
  p <- length(optresults$mode)
  quadtable <- get_quadtable(p,k,ndConstruction)
  # TODO: speed and memory comparisons/optimization. Or just store the tables...
  # quadtable <- as(as.matrix(quadtable),'sparseMatrix')
  quadtable[is.na(quadtable)] <- 0
  # TODO: adaptation
  0
}
#' @rdname nested
#' @export
get_quadtable <- function(p,k,ndConstruction = 'product',...) {
  if (k%%2==0) {
    warning("Nested quadrature usually produces the same grid for adjacent even and odd k. You provided even k. Setting k = k+1.")
    k <- k+1
  }
  # k is now odd
  gg <- mvQuad::createNIGrid(p,'nHe',k,ndConstruction = ndConstruction)
  nn <- as.data.frame(mvQuad::getNodes(gg))
  nn$w <- mvQuad::getWeights(gg)[ ,1]
  colnames(nn)[colnames(nn)=='w'] <- paste0('w',k)

  # recursion
  if (k > 1) {
    nnl <- get_quadtable(p,k-2)
    mergevars <- paste0('V',1:p)
    nnm <- merge(nn,nnl,all=TRUE,by = mergevars)
    return(nnm)
  }
  nn
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
get_log_normconst.numeric <- function(obj,...) obj
#' @rdname get_log_normconst
#' @export
get_log_normconst.aghq <- function(obj,...) get_log_normconst(obj$normalized_posterior)
#' @rdname get_log_normconst
#' @export
get_log_normconst.laplace <- function(obj,...) obj$lognormconst
#' @rdname get_log_normconst
#' @export
get_log_normconst.marginallaplace <- function(obj,...) get_log_normconst(obj$normalized_posterior)

#' Obtain the nodes and weights table from a fitted quadrature object
#'
#' Quick helper S3 method to retrieve the quadrature nodes and weights from an object
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
get_nodesandweights <- function(obj,...) UseMethod("get_nodesandweights")
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.default <- function(obj,...) obj$nodesandweights
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.list <- function(obj,...) obj$nodesandweights
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.data.frame <- function(obj,...) obj
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.aghq <- function(obj,...) get_nodesandweights(obj$normalized_posterior)
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.laplace <- function(obj,...) {
  # This is never actually created, so create it here
  out <- as.data.frame(rbind(c(obj$optresults$mode,
                               exp(-as.numeric(.5 * determinant(obj$optresults$hessian,logarithm = TRUE)$modulus) + .5*dim(obj$optresults$hessian)[1]*log(2*pi)),
                               obj$optresults$ff$fn(obj$optresults$mode),
                               obj$optresults$ff$fn(obj$optresults$mode) - obj$lognormconst)))
  names(out) <- c(paste0('theta',1:length(obj$optresults$mode)),'weights','logpost','logpost_normalized')
  out
}
#' @rdname get_nodesandweights
#' @export
get_nodesandweights.marginallaplace <- function(obj,...) get_nodesandweights(obj$normalized_posterior)

#' Obtain the number of quadrature nodes used from an aghq object
#'
#' Quick helper S3 method to retrieve the number of quadrature points used when creating an aghq object.
#'
#' @param obj Object of class \code{aghq} returned by \code{aghq::aghq}.
#' @param ... Not used
#'
#' @return A numeric vector of length 1 containing \code{k}, the number of quadrature points used.
#'
#' @family quadrature
#'
#' @export
#'
get_numquadpoints <- function(obj,...) as.numeric(obj$normalized_posterior$grid$level)[1]

#' Obtain the parameter dimension from an aghq object
#'
#' Quick helper S3 method to retrieve the parameter dimension from an aghq object.
#'
#' @param obj Object of class \code{aghq} returned by \code{aghq::aghq}.
#' @param ... Not used
#'
#' @return A numeric vector of length 1 containing \code{p}, the parameter dimension.
#'
#' @family quadrature
#'
#' @export
#'
get_param_dim <- function(obj,...) UseMethod('get_param_dim')
#' @rdname get_param_dim
#' @export
get_param_dim.aghq <- function(obj,...) length(obj$optresults$mode)

#' Obtain the optimization results from an aghq object
#'
#' Quick helper S3 method to retrieve the mode and Hessian from an aghq object. The
#' full results of calling \code{aghq::optimize_theta} are stored in \code{obj$optresults}.
#'
#' @param obj Object of class \code{aghq} returned by \code{aghq::aghq}.
#' @param ... Not used
#'
#' @return A named list with elements:
#' \itemize{
#' \item{\code{mode}: a numeric vector of length \code{dim(theta)} containing the posterior mode.}
#' \item{\code{hessian}: a numeric matrix of dimension \code{dim(theta) x dim(theta)} containing the negative Hessian of the log-posterior evaluated at the mode.}
#' }
#' For objects of class \code{marginallaplace}, a third list item \code{modesandhessians} is
#' a \code{data.frame} containing
#' the mode and Hessian of the \code{W} parameters evaluated at each adapted quadrature point.
#'
#' @family quadrature
#'
#' @export
#'
get_opt_results <- function(obj,...) UseMethod('get_opt_results')
#' @rdname get_opt_results
#' @export
get_opt_results.aghq <- function(obj,...) list(mode = obj$optresults$mode,hessian = obj$optresults$hessian)
#' @rdname get_opt_results
#' @export
get_opt_results.marginallaplace <- function(obj,...) list(mode = obj$optresults$mode,hessian = obj$optresults$hessian)

#' Obtain the mode from an aghq object
#'
#' Quick helper method to retrieve the mode from an aghq object. Just
#' calls \code{aghq::get_opt_results}.
#'
#' @param obj Object of class \code{aghq} returned by \code{aghq::aghq}.
#' @param ... Not used
#'
#' @return A numeric vector of length \code{dim(theta)} containing the posterior mode.
#'
#' @family quadrature
#'
#' @export
#'
get_mode <- function(obj,...) get_opt_results(obj,...)$mode

#' Obtain the Hessian from an aghq object
#'
#' Quick helper method to retrieve the Hessian from an aghq object. Just
#' calls \code{aghq::get_opt_results}.
#'
#' @param obj Object of class \code{aghq} returned by \code{aghq::aghq}.
#' @param ... Not used
#'
#' @return A numeric matrix of dimension \code{dim(theta) x dim(theta)} containing the negative Hessian of the log-posterior evaluated at the mode.
#'
#' @family quadrature
#'
#' @export
#'
get_hessian <- function(obj,...) get_opt_results(obj,...)$hessian

