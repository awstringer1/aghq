### Summary Statistics ###
# This script contains functions for computing summary
# statistics from aghq output. Marginal posteriors,
# moments, and quantiles.



#' Marginal Posteriors
#'
#' Compute the marginal posterior for a given parameter using AGHQ.
#'
#' @inheritParams normalize_logpost
#' @param j Integer between 1 and the dimension of the parameter space. Which
#' index of the parameter vector to compute the marginal posterior for.
#'
#' @return a tbl_df/tbl/data.frame containing the normalized log marginal posterior
#' for theta_j evaluated at the original quadrature points. Has columns
#' \code{"thetaj","logpost_normalized","weights"}, where \code{j} is the \code{j} you specified.
#'
#' @family summaries
#'
#' @examples
#' ## A 2d example ##
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
#
#' opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
#'
#' # Now actually do the marginal posteriors
#' marginal_posterior(opt_sparsetrust_2d,3,1)
#' marginal_posterior(opt_sparsetrust_2d,3,2)
#' marginal_posterior(opt_sparsetrust_2d,7,2)
#'
#' @export
#'
marginal_posterior <- function(optresults,k,j) {
  normresults <- normalize_logpost(optresults,k,whichfirst = j)
  nodesandweights <- normresults$nodesandweights
  thetagridfull <- normresults$grid
  lognormconstorig <- normresults$lognormconst

  S <- length(optresults$mode)
  # If S = 1 just return the nodesandweights, properly formatted
  if (S == 1) {
    out <- nodesandweights[ ,c("theta1","logpost_normalized","weights")]
    colnames(out) <- c("theta1","logmargpost","w")
    return(out)
  }
  idxorder <- c(j,(1:S)[-j])

  # In the right order...
  m <- thetagridfull$features$m
  C <- thetagridfull$features$C

  gg <- mvQuad::createNIGrid(1,"GHe",k) # Get the nodes and weights from a 1-d GH rule
  ww <- mvQuad::getWeights(gg)
  wwE <- do.call(expand.grid,rep(list(ww),S)) # Expand these into a grid to get a S-d GH rule
  colnames(wwE) <- paste0("theta",idxorder,"W")
  cholinvH <- chol(C)
  diagcholinvH <- diag(cholinvH) # Get the elements of the diagonal of the inverse cholesky, for multiplying the weights
  nn <- mvQuad::getNodes(thetagridfull)
  colnames(nn) <- paste0("theta",idxorder)

  thetaj <- paste0("theta",j)
  thetaminusj <- paste0("theta",(1:S)[-j])

  nodesandweightsfactored <- merge(cbind(nn,wwE),nodesandweights,by = c(thetaj,thetaminusj),all = FALSE)

  out <- nodesandweightsfactored
  out['id'] <- seq(1,nrow(out),by = 1)
  # Use apply or sapply depending on dimension, because base R data frame indexing
  # does not return a consistent result
  if (S == 2) {
    out['logw'] <- sapply(out[ ,paste0(thetaminusj,"W")],function(x) sum(log(x)) + sum(log(diagcholinvH[-1])))
  } else {
    out['logw'] <- apply(out[ ,paste0(thetaminusj,"W")],1,function(x) sum(log(x)) + sum(log(diagcholinvH[-1])))
  }

  # tapply doesn't preserve the order of its arguments. This caused a bug
  tapres <- tapply(out$logw + out$logpost_normalized,out[ ,thetaj],logsumexp)
  tt <- unique(out[[thetaj]])
  ttord <- match(tt,names(tapres))
  out <- data.frame(theta = tt,logmargpost = as.numeric(tapres[ttord]))
  colnames(out)[colnames(out) == 'theta'] <- thetaj

  out$w <- as.numeric(ww * diagcholinvH[1])
  out
}

#' Interpolate the Marginal Posterior
#'
#' Build a Lagrange polynomial interpolant of the marginal posterior, for plotting
#' and for computing quantiles.
#'
#' @param margpost The output of \code{aghq::marginal_posterior}. See the documentation for that function.
#'
#' @return A function of \code{theta} which computes the log interpolated normalized marginal posterior.
#'
#' @family summaries
#'
#' @export
#'
interpolate_marginal_posterior <- function(margpost) {
  # Unname the theta
  colnames(margpost)[grep("theta",colnames(margpost))] <- "theta"

  as.function(polynom::poly.calc(x = margpost$theta,y = margpost$logmargpost))
}

#' Compute moments
#'
#' Compute the moment of any function ff using AGHQ.
#'
#' @param obj Object of class \code{aghq} output by \code{aghq::aghq()}, or its \code{normalized_posterior} element. See \code{?aghq}.
#' @param ff Any R function which takes in a numeric vector and returns a numeric vector.
#' @param ... Used to pass additional argument \code{ff}.
#'
#' @return A numeric vector containing the moment(s) of ff with respect to the joint
#' distribution being approximated using AGHQ.
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
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#
#' opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
#' norm_sparse_2d_7 <- normalize_logpost(opt_sparsetrust_2d,7,1)
#'
#' # ff = function(x) 1 should return 1,
#' # the normalizing constant of the (already normalized) posterior:
#' compute_moment(norm_sparse_2d_7)
#' # Compute the mean of theta1 and theta2
#' compute_moment(norm_sparse_2d_7,ff = function(x) x)
#' # Compute the mean of lambda1 = exp(theta1) and lambda2 = exp(theta2)
#' lambdameans <- compute_moment(norm_sparse_2d_7,ff = function(x) exp(x))
#' lambdameans
#' # Compare them to the truth:
#' (sum(y1) + 1)/(length(y1) + 1)
#' (sum(y2) + 1)/(length(y2) + 1)
#' # Compute the standard deviation of lambda1
#' lambda1sd <- sqrt(compute_moment(norm_sparse_2d_7,ff = function(x) (exp(x) - lambdameans[1])^2))[1]
#' # ...and so on.
#' @family summaries
#' @export
#'
compute_moment <- function(obj,...) {
  UseMethod("compute_moment")
}
#' @rdname compute_moment
#' @method compute_moment default
#' @export
compute_moment.default <- function(obj,ff = function(x) 1,...) {
  nodesandweights <- obj$nodesandweights

  whereistheta <- grep('theta',colnames(nodesandweights))

  lengthof_f <- length(ff(nodesandweights[1,whereistheta]))

  if (lengthof_f == 1) {
    out <- sum(ff(nodesandweights[ ,whereistheta])* exp(nodesandweights$logpost_normalized) * nodesandweights$weights)
  } else {
    out <- apply(nodesandweights[ ,whereistheta],1,ff)
    out <- apply(out,1,function(x) sum(x * exp(nodesandweights$logpost_normalized) * nodesandweights$weights))
  }

  unname(out)
}
#' @rdname compute_moment
#' @method compute_moment aghq
#' @export
compute_moment.aghq <- function(obj,ff = function(x) 1,...) compute_moment(obj$normalized_posterior,ff)

#' Density and Cumulative Distribution Function
#'
#' Compute the density and cumulative distribution function of the approximate posterior.
#' The density is approximated on a find grid using a polynomial interpolant.
#' The CDF can't be computed exactly (if it could, you wouldn't be using quadrature!),
#' so a fine grid is laid down and the CDF is approximated at each grid point
#' using a simpler integration rule and a polynomial interpolant. This method tends
#' to work well, but won't always.
#'
#' @param obj Either the output of \code{aghq::aghq()}, its list of marginal distributions
#' (element \code{marginals}), or an individual \code{data.frame} containing one of
#' these marginal distributions as output by \code{aghq::marginal_posterior()}.
#' @param finegrid Optional, a grid of values on which to compute the CDF. The default makes
#' use of the values in \code{margpost} but if the results are unsuitable, you may wish to
#' modify this manually.
#' @param transformation Optional. A list containing two functions, \code{fromtheta}
#' and \code{totheta}, which accept and return numeric vectors, defining a parameter transformation for which you would
#' also like the pdf calculated for. See examples. May also have an element \code{jacobian},
#' a function which takes a numeric vector and computes the jacobian of the transformation; if
#' not provided, this is done using \code{numDeriv::jacobian}.
#' @param ... Used to pass additional arguments.
#'
#' @return A tbl_df/tbl/data.frame with columns \code{theta}, \code{pdf} and \code{cdf} corresponding
#' to the value of the parameter and its estimated PDF and CDF at that value.
#'
#' @family summaries
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
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#
#' opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
#' margpost <- marginal_posterior(opt_sparsetrust_2d,3,1) # margpost for theta1
#' thepdfandcdf <- compute_pdf_and_cdf(margpost)
#' with(thepdfandcdf,{
#'   plot(pdf~theta,type='l')
#'   plot(cdf~theta,type='l')
#' })
#'
#' @export
#'
compute_pdf_and_cdf <- function(obj,...) UseMethod("compute_pdf_and_cdf")
#' @rdname compute_pdf_and_cdf
#' @export
compute_pdf_and_cdf.default <- function(obj,transformation = NULL,finegrid = NULL,...) {

  margpostinterp <- interpolate_marginal_posterior(obj)

  thetacol <- colnames(obj)[grep("theta",colnames(obj))]

  if (is.null(finegrid)) {
    rn <- range(obj[[thetacol]])
    rnl <- diff(rn)
    thetarange <- c(min(rn) - rnl/2,max(rn) + rnl/2)
    finegrid <- c(seq(thetarange[1],thetarange[2],length.out=1000))
  }

  out <- data.frame(
    theta = finegrid,
    pdf = exp(margpostinterp(finegrid)),
    cdf = cumsum(exp(margpostinterp(finegrid)) * c(0,diff(finegrid)))
  )

  if (!is.null(transformation)) {
    if (is.null(transformation$jacobian)) {
      transformation$jacobian <- function(theta) {
        out <- numeric(length(theta))
        for (i in 1:length(theta)) {
          out[i] <- det(abs(numDeriv::jacobian(transformation$totheta,transformation$fromtheta(theta[i]))))
        }
        out
      }
    }
    out$transparam <- transformation$fromtheta(out$theta)
    out$pdf_transparam <- out$pdf * transformation$jacobian(out$theta)
  }
  out
}
#' @rdname compute_pdf_and_cdf
#' @export
compute_pdf_and_cdf.list <- function(obj,...) {
  out <- list()
  for (i in 1:length(obj)) out[[i]] <- compute_pdf_and_cdf(obj[[i]],...)
  out
}
#' @rdname compute_pdf_and_cdf
#' @export
compute_pdf_and_cdf.aghq <- function(obj,...) compute_pdf_and_cdf(obj$marginals,...)

#' Quantiles
#'
#' Compute marginal quantiles using AGHQ. This function works by first approximating
#' the CDF using \code{aghq::compute_pdf_and_cdf} and then inverting the approximation numerically.
#'
#' @inheritParams compute_pdf_and_cdf
#' @param q Numeric vector of values in (0,1). The quantiles to compute.
#'
#' @return A named numeric vector containing the quantiles you asked for, for the
#' variable whose marginal posterior you provided.
#'
#' @family summaries
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
#' funlist2d <- list(
#'   fn = objfunc2d,
#'   gr = function(x) numDeriv::grad(objfunc2d,x),
#'   he = function(x) numDeriv::hessian(objfunc2d,x)
#' )
#
#' opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
#' margpost <- marginal_posterior(opt_sparsetrust_2d,3,1) # margpost for theta1
#' etaquant <- compute_quantiles(margpost)
#' etaquant
#' # lambda = exp(eta)
#' exp(etaquant)
#' # Compare to truth
#' qgamma(.025,1+sum(y1),1+n1)
#' qgamma(.975,1+sum(y1),1+n1)
#'
#'
#'
#' @export
#'
compute_quantiles <- function(obj,...) UseMethod("compute_quantiles")
#' @rdname compute_quantiles
#' @export
compute_quantiles.default <- function(obj,q = c(.025,.975),...) {
  pdfandcdf <- compute_pdf_and_cdf(obj)
  out <- numeric(length(q))
  names(out) <- paste0(as.character(100 * q),"%")

  for (i in 1:length(q)) {
    out[i] <- pdfandcdf$theta[max(which(pdfandcdf$cdf < q[i]))]
  }
  out
}
#' @rdname compute_quantiles
#' @export
compute_quantiles.list <- function(obj,...) {
  out <- list()
  for (i in 1:length(obj)) out[[i]] <- compute_quantiles(obj[[i]],...)
  out
}
#' @rdname compute_quantiles
#' @export
compute_quantiles.aghq <- function(obj,...) compute_quantiles(obj$marginals,...)

#' Exact independent samples from an approximate posterior distribution
#'
#' Draws samples from an approximate marginal distribution for general posteriors
#' approximated using \code{aghq}, or from the mixture-of-Gaussians approximation to the variables that were
#' marginalized over in a marginal Laplace approximation fit using \code{aghq::marginal_laplace}
#' or \code{aghq::marginal_laplace_tmb}.
#'
#' @param quad Object from which to draw samples.
#' An object inheriting from class \code{marginallaplace}
#' (the result of running \code{aghq::marginal_laplace} or \code{aghq::marginal_laplace_tmb}),
#' or an object inheriting from class \code{aghq} (the result of running \code{aghq::aghq()}).
#' Can also provide a \code{data.frame} returned by \code{aghq::compute_pdf_and_cdf} in which
#' case samples are returned for \code{transparam} if \code{transformation} is provided,
#' and for \code{param} if \code{transformation = NULL}.
#' @param M Numeric, integer saying how many samples to draw
#' @param transformation Optional.
#' A list containing function \code{fromtheta()} which accepts and returns numeric vectors,
#' defining a parameter transformation for which you would like samples to be taken.
#' See \code{?compute_pdf_and_cdf}. Note that unlike there, where this operation is
#' a bit more complicated, here all is done is samples are taken on the original
#' scale and then \code{transformation$fromtheta()} is called on them before returning.
#' @param ... Used to pass additional arguments.
#'
#' @family sampling
#'
#' @return If run on a \code{marginallaplace} object, a list containing elements:
#' \itemize{
#' \item{\code{samps}: }{ \code{d x M} matrix where \code{d = dim(W)} and each column is a sample
#' from \code{pi(W|Y,theta)}}
#' \item{\code{theta}: }{\code{M x S} tibble where \code{S = dim(theta)} containing the value of \code{theta} for
#' each sample}
#' \item{\code{thetasamples}: }{A list of \code{S} numeric vectors each of length
#' \code{M} where the \code{j}th element is a sample from \code{pi(theta_{j}|Y)}}. These are samples
#' from the **marginals**, NOT the **joint**. Sampling from the joint is a much more difficult
#' problem and how to do so in this context is an active area of research.
#' }
#' If run on an \code{aghq} object, then a list with just the \code{thetasamples} element. It still
#' returns a list to maintain output consistency across inputs.
#'
#' If, for some reason, you don't want to do the sampling from \code{pi(theta|Y)}, you can manually
#' set \code{quad$marginals = NULL}. Note that this sampling is typically *very* fast
#' and so I don't know why you would need to not do it but the option is there if you like.
#'
#' If, again for some reason, you just want samples from one marginal distribution using inverse CDF,
#' you can just do \code{compute_quantiles(quad$marginals[[1]],runif(M))}.
#'
#' @details For objects of class \code{aghq} or their marginal distribution components,
#' sampling is done using the inverse CDF method, which is just \code{compute_quantiles(quad$marginals[[1]],runif(M))}.
#'
#' For marginal Laplace approximations (\code{aghq::marginal_laplace()}): this method samples from the posterior and returns a vector that is ordered
#' the same as the "\code{W}" variables in your marginal Laplace approximation. See Algorithm 1 in
#' Stringer et. al. (2021, https://arxiv.org/abs/2103.07425) for the algorithm; the details of sampling
#' from a Gaussian are described in the reference(s) therein, which makes use of the (sparse)
#' Cholesky factors. These are computed once for each quadrature point and stored.
#'
#' For the marginal Laplace approximations where the "inner" model is handled entirely by \code{TMB}
#' (\code{aghq::marginal_laplace_tmb}), the interface here is identical to above,
#' with the order of the "\code{W}" vector being determined by \code{TMB}. See the
#' \code{names} of \code{ff$env$last.par}, for example (where \code{ff} is your
#' template obtained from a call to \code{TMB::MakeADFun}.
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
# themargsamps <- aghq::sample_marginal(themarginallaplace,10)
#'
#' @export
#'
sample_marginal <- function(quad,...) UseMethod("sample_marginal")
#' @rdname sample_marginal
#' @export
sample_marginal.aghq <- function(quad,M,transformation = NULL,...) {
  out <- list()
  if (is.null(quad$marginals)) return(out)
  for (i in 1:length(quad$marginals)) out[[i]] <- unname(compute_quantiles(quad$marginals[[i]],stats::runif(M)))

  if (!is.null(transformation)) {
    if (is.null(transformation$fromtheta)) warning("transformation provided but transformation$fromtheta appears NULL.\n")
    for (i in 1:length(out)) out[[i]] <- transformation$fromtheta(out[[i]])
  }
  out
}
#' @rdname sample_marginal
#' @export
sample_marginal.marginallaplace <- function(quad,M,transformation = NULL,...) {
  K <- as.numeric(quad$normalized_posterior$grid$level)[1]
  d <- dim(quad$modesandhessians$H[[1]])[1]
  simlist <- quad$modesandhessians
  simlist$L <- lapply(simlist$H,function(h) chol(Matrix::forceSymmetric(h),perm = FALSE))
  simlist$lambda <- exp(quad$normalized_posterior$nodesandweights$logpost_normalized) * quad$normalized_posterior$nodesandweights$weights

  # Sample from the multinomial
  if (M ==1) {
    k <- which(stats::rmultinom(M,1,simlist$lambda) == 1)
  } else {
    k <- apply(stats::rmultinom(M,1,simlist$lambda),2,function(x) which(x == 1))
  }
  tt <- table(k) # Values are number of samples to draw for each k

  # Big Gaussian mixture matrix
  Z <- lapply(split(matrix(stats::rnorm(M*d),nrow = M),k),matrix,nrow = d)

  samps <- mapply(
    function(.x,.y) as.numeric(solve(simlist$L[[as.numeric(.y)]],.x)) + do.call(cbind,rep(list(simlist$mode[[as.numeric(.y)]]),ncol(.x))),
    Z,
    names(Z)
  )

  # Order them properly
  ord <- numeric(length(k))
  cumtab <- cumsum(c(0,table(k)))
  cumtab <- cumtab[-length(cumtab)]
  cnt <- numeric(length(unique(k)))
  names(cnt) <- sort(unique(k))
  for (i in 1:length(k)) {
    wc <- which(names(cnt) == k[i])
    cnt[wc] <- cnt[wc] + 1
    ord[i] <- cumtab[wc] + cnt[wc]
  }

  samps <- Reduce(cbind,samps)
  samps <- samps[ ,ord]

  theta <- simlist[k,paste0('theta',seq(1,length(grep('theta',colnames(simlist)))))]

  # In one dimension, R's indexing is not type consistent
  if (!is.matrix(samps)) {
    samps <- rbind(samps)
    rownames(samps) <- NULL
  }

  if (!inherits(theta,"data.frame")) theta <- data.frame(theta1 = theta)

  out <- list(
    samps = samps,
    theta = theta
  )
  # Add the marginals from theta|Y, with possible transformation
  class(quad) <- "aghq"
  out$thetasamples <- sample_marginal(quad,M,transformation,...)
  out
}