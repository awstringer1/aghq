### Summary Statistics ###
# This script contains functions for computing summary
# statistics from aghq output. Marginal posteriors,
# moments, and quantiles.



#' Marginal Posteriors
#'
#' Compute the marginal posterior for a given parameter using AGHQ. This function is
#' mostly called within \code{aghq()}.
#'
#' @inheritParams normalize_logpost
#' @param j Integer between 1 and the dimension of the parameter space. Which
#' index of the parameter vector to compute the marginal posterior for.
#' @param qq Numeric vector of length \code{>=1} giving the points at which to evaluate the marginal posterior.
#' The default, \code{NULL}, chooses these points in a 'clever' way, see Details.
#' @param quad Object returned by \code{aghq::aghq}.
#' @param method Method for computing the quadrature points used to approximate moment.
#' One of 'reuse' (default) or 'correct'. See details.
#' The default SHOULD be 'correct'; it is currently set to 'reuse' to maintain compatibility of
#' results with previous versions. This will be switched in a future major release.
#'
#' @return a data.frame containing the normalized log marginal posterior
#' for theta_j evaluated at the original quadrature points. Has columns
#' \code{"thetaj","logpost_normalized","weights"}, where \code{j} is the \code{j} you specified.
#'
#' @details If \code{qq=NULL}, then it is set to the unique values in an adapted GHQ grid computed
#' assuming that \code{j=1} (there is nothing special about this procedure, it's just a way to provide
#' an apparently sensible default).
#'
#' If \code{method='reuse'}, then the parameter vector is reordered so \code{j=1}, and the
#' approximate marginal is computed by first computing the whole AGHQ grid, then summing over the other
#' indices. This is an outdated method that does not have any theory pertaining to it, and is included for
#' backwards compatibility. It does not use \code{qq} if supplied.
#'
#' If \code{method='correct'} then the theoretically-justified approximation from Section 2.4 of the 'Stochastic Convergence Rates...'
#' paper is returned.
#'
#' \code{method='auto'} currently chooses \code{'reuse'} for backwards compatibility, but this will be
#' changed in a future release.
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
marginal_posterior <- function(...) UseMethod("marginal_posterior")
#' @rdname marginal_posterior
#' @export
marginal_posterior.aghq <- function(quad,j,qq=NULL,method = c('auto','reuse','correct'),...) {
  method <- method[1]
  if (method=='auto') method <- 'reuse'
  if (method=='reuse') return(marginal_posterior.list(quad$optresults,get_numquadpoints(quad),j,...))

  # marginal_posterior.list(quad$optresults,get_numquadpoints(quad),j,...)
  S <- get_param_dim(quad)
  idxorder <- c(j,(1:S)[-j])
  thetaminusj <- (1:S)[-j]
  if (is.null(qq)) {
    # If q not provided, choose some clever q's and then re-call the function at those q's
    # Choose evaluation points
    mm <- quad$optresults$mode[idxorder]
    HH <- quad$optresults$hessian[idxorder,idxorder]
    LL <- t(solve(chol(HH)))
    gg <- mvQuad::createNIGrid(S,'GHe',get_numquadpoints(quad))
    nn <- as.matrix(mvQuad::getNodes(gg))
    qqq <- unique(sweep(LL%*%t(nn),1,mm,'+')[1, ])

    out <- vector(mode='list',length=length(qqq))
    for (i in 1:length(qqq)) {
      # Call marginal_posterior at each qq[j]
      out[[i]] <- marginal_posterior.aghq(quad = quad,j = j,qq = qqq[i],method = 'correct',...)
    }
    out <- Reduce(rbind,out)
  } else {
    # Get the name
    cname <- colnames(get_nodesandweights(quad))[j]
    if (S==1) {
      out <- data.frame(qq,quad$optresults$ff$fn(qq) - get_log_normconst(quad))
      colnames(out) <- c(cname,'logmargpost')
      return(out)
    }
    # Get the value
    fn <- function(theta) quad$optresults$ff$fn(splice(theta,qq,j))
    gr <- function(theta) quad$optresults$ff$gr(splice(theta,qq,j))[-j]
    he <- function(theta) quad$optresults$ff$he(splice(theta,qq,j))[-j,-j]
    ffm <- list(fn=fn,gr=gr,he=he)
    newcontrol <- quad$control
    newcontrol$onlynormconst <- TRUE
    lognumerator <- get_log_normconst(aghq(ffm,get_numquadpoints(quad),get_mode(quad)[-j],control = newcontrol))
    out <- data.frame(qq,lognumerator-get_log_normconst(quad))
    colnames(out) <- c(cname,'logmargpost')
  }
  out
}
#' @rdname marginal_posterior
#' @export
marginal_posterior.list <- function(optresults,k,j,basegrid = NULL,ndConstruction = 'product',...) {

  # If using sparse grids, marginals are currently not supported
  dummyout <- data.frame(theta1 = 0,logmargpost = 0,w = 0)
  if (!is.null(basegrid)) {
    if (basegrid$ndConstruction != 'product') return(dummyout)
  }
  if (ndConstruction != 'product') return(dummyout)

  normresults <- normalize_logpost(optresults,k,whichfirst = j,basegrid=basegrid,ndConstruction = ndConstruction)
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
#' @param method The method to use. Default is a \code{k} point polynomial interpolant using \code{polynom::poly.calc()}.
#' This has been observed to result in unstable behaviour for larger numbers of quadrature points \code{k},
#' which is of course undesirable. If \code{k > 3}, you can set \code{method = 'spline'} to use \code{splines::interpSpline()} instead,
#' which uses cubic B-Splines. These should always be better than a straight polynomial, except don't work
#' when \code{k < 4} which is why they aren't the default. If you try and set \code{method = 'spline'} with
#' \code{k < 4} it will be changed back to polynomial, with a warning.
#'
#' @return A function of \code{theta} which computes the log interpolated normalized marginal posterior.
#'
#' @family summaries
#'
#' @export
#'
interpolate_marginal_posterior <- function(margpost,method = c('auto','polynomial','spline')) {
  # Unname the theta
  colnames(margpost)[grep("theta",colnames(margpost))] <- "theta"

  goodmethods <- c('auto','polynomial','spline')
  method <- method[1]
  if (!is.character(method)) {
    stop(paste0("'method' should be one of",goodmethods,", instead you supplied a ",class(method)," object"))
  } else {
    if (!(method %in% goodmethods)) {
      stop(paste0("'method' should be one of",goodmethods,", instead you provided ",method))
    }
  }

  # If method is auto, assign spline if k >= 4 and polynomial if k < 4
  k <- nrow(margpost)
  if (method == 'auto') {
    if (k >=4) {
      method <- 'spline'
    } else {
      method <- 'polynomial'
    }
  }

  if (method == 'spline' & k < 4) {
    warning("You asked to use a cubic B-Spline interpolant for the marginal posteriors, however you have k < 4 so it doesn't work. Using a polynomial interpolant instead.\n")
    method <- "polynomial"
  }

  if (method == 'polynomial') {
    # if (k > 3 & verbose) warning("Polynomial interpolation not recommended with more than 3 quadrature points. Try default_control(method = 'auto') or 'spline'")
    out <- as.function(polynom::poly.calc(x = margpost$theta,y = margpost$logmargpost))
  } else if (method == 'spline') {
    ss <- with(margpost,splines::interpSpline(theta,logmargpost,bSpline = TRUE,sparse = TRUE))
    out <- function(x) as.numeric(stats::predict(ss,x)$y)
  } else {
    stop(paste0("Unrecognized interpolation method ",method,", should be one of 'spline' or 'polynomial'.\n"))
  }
  out
}

#' Compute moments
#'
#' Compute the moment of any function ff using AGHQ.
#'
#' @param obj Object of class \code{aghq} output by \code{aghq::aghq()}. See \code{?aghq}.
#' @param ff Any R function which takes in a numeric vector and returns a numeric vector. Exactly one of \code{ff} or \code{gg} must be provided. If both are provided, \code{aghq::compute_moment()} will use \code{gg},
#' without warning.
#' @param gg The output of, or an object which may be input to \code{aghq::make_moment_function()}. See documentation of that function. Exactly one of \code{ff} or \code{gg} must be provided. If both are provided, \code{aghq::compute_moment()} will use \code{gg},
#' without warning.
#' @param nn A numeric scalar. Compute the approximate moment of this order, \code{E(theta^nn|Y)}. See details.
#' If \code{nn} is provided, \code{compute_moment} will use it over \code{ff} or \code{gg},
#' without warning.
#' @param type Either \code{'raw'} (default) or \code{'central'}, see details.
#' @param method Method for computing the quadrature points used to approximate moment. One of \code{'reuse'} (default) or \code{'correct'}. See details. The default SHOULD be \code{'correct'}; it is currently
#' set to \code{'reuse'} to maintain compatibility of results with previous versions. This will be switched in a future major release.
#' @param ... Used to pass additional argument \code{ff}.
#'
#' @return A numeric vector containing the moment(s) of ff with respect to the joint
#' distribution being approximated using AGHQ.
#'
#' @details If mutliple of \code{nn}, \code{gg}, and \code{ff} are provided, then \code{compute_moment}
#' will use \code{nn}, \code{gg}, or \code{ff}, in that order, without warning.
#'
#' There are several approximations available. The "best" one is obtained by specifying \code{gg}
#' and using \code{method = 'correct'}. This recomputes the mode and curvature for the
#' function \code{g(theta)posterior(theta)}, and takes the ratio of the AGHQ approximation
#' to this function to the AGHQ approximation to the marginal likelihood. This obtains the
#' same relative rate of convergence as the AGHQ approximation to the marginal likelihood. It
#' may take a little extra time, and only works for **positive, scalar-valued** functions \code{g}.
#'
#' \code{method = 'reuse'} re-uses the AGHQ adapted points and weights. It's faster than the
#' correct method, because it does not involve any new optimization, it's just a weighted sum.
#' No convergence theory. Seems to work ok in "practice". "Works" for arbitrary \code{g}.
#'
#' Specifying \code{ff} instead of \code{gg} automatically uses \code{method = 'reuse'}. This
#' interface is provided for backwards compatibility mostly. However, one advantage is that
#' it allows for **vector-valued** functions, in which case it just returns the corresponding
#' vector of approximate moments. Also, it only requires the adapted nodes and weights, not
#' the ability to evaluate the log-posterior and its derivatives, although this is unlikely
#' to be a practical concern.
#'
#' Specifying a numeric value \code{nn} will return the moment \code{E(theta^nn|Y)}.
#' This automatically does some internal shifting to get the evaluations away from zero,
#' to avoid the inherent problem of multi-modal "posteriors" that occurs when the posterior
#' mode is near zero, and account for the fact that some of the new adapted quadrature points
#' may be negative. So, the actual return value is \code{E(theta^nn + a|Y) - a} for a cleverly-chosen
#' value \code{a}.
#'
#' Finally, \code{type='raw'} computes raw moments \code{E(g(theta)|Y)}, where \code{type='central'}
#' computes central moments, \code{E(g(theta - E(g(theta)|Y))|Y)}. See examples.
#'
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
#' quad <- aghq(funlist2d,7,c(0,0))
#'
#' @family summaries moments
#' @export
#'
compute_moment <- function(obj,...) {
  UseMethod("compute_moment")
}
#' @rdname compute_moment
# #' @method compute_moment default
#' @export
compute_moment.list <- function(obj,ff = function(x) 1,gg = NULL,nn = NULL,type = c('raw','central'),method = c("auto","reuse","correct"),...) {
  valid_types <- c('raw','central')
  type <- type[1]
  if (!(type %in% valid_types)) stop(paste0("type should be one of: ",valid_types,", but you provided: ",type,".\n"))
  # If the user tries to pass a list but method is not "reuse", print a warning
  method <- method[1]
  valid_methods <- c("auto","reuse","correct")
  if (!(method %in% valid_methods)) stop(paste0("method should be one of: ",valid_methods,", but you provided: ",method,".\n"))
  if (method == 'auto') method <- "reuse"
  if (method != "reuse") warning(paste0("Method = 'reuse' is the only possible option for computing moments when you have passed an object of class 'list'. You selected method = ",method,". You could try passing an object of class 'aghq' to 'compute_moment' if you want to use the more advanced options.\n" ))

  nodesandweights <- get_nodesandweights(obj)
  # Prepare the numeric moment. Nothing fancy to do here, since method = 'reuse'
  if (!is.null(nn)) {
    if (type == 'central') {
      if (nn == 1) return(0) # Central moment of first order is zero by definition
      mm <- compute_moment(obj,nn=1,type='raw')
    } else {
      mm <- 0
    }
    ff <- function(theta) (theta-mm)^nn
  }
  if (is.null(gg) & is.null(ff)) stop("You provided NULL for ff and gg in compute_moment. You have to provide one of these.\n")
  if (!is.null(gg)) {
    # Use gg over ff if both provided.
    if (!validate_moment(gg)) stop("validate_moment(gg) failed inside compute_moment.\n")
    ff <- function(theta) exp(gg$fn(theta))
  }
  # whereistheta <- grep('theta',colnames(nodesandweights))
  whereistheta <- 1:(ncol(nodesandweights)-3)

  lengthof_f <- length(ff(nodesandweights[1,whereistheta]))

  if (lengthof_f == 1) {
    out <- sum(ff(nodesandweights[ ,whereistheta])* exp(nodesandweights$logpost_normalized) * nodesandweights$weights)
  } else {
    out <- apply(nodesandweights[ ,whereistheta],1,ff)
    out <- apply(out,1,function(x) sum(x * exp(nodesandweights$logpost_normalized) * nodesandweights$weights))
  }

  # If shifting was applied, unshift
  if (!is.null(get_shift(gg))) out <- out - get_shift(gg)

  unname(out)
}
#' @rdname compute_moment
# #' @method compute_moment aghq
#' @export
compute_moment.aghq <- function(obj,ff = function(x) 1,gg = NULL,nn = NULL,type = c('raw','central'),method = c("auto","reuse","correct"),...) {
  valid_types <- c('raw','central')
  type <- type[1]
  if (!(type %in% valid_types)) stop(paste0("type should be one of: ",valid_types,", but you provided: ",type,".\n"))
  valid_methods <- c("auto","reuse","correct")
  method <- method[1]
  if (!(method %in% valid_methods)) stop(paste0("method should be one of: ",valid_methods,", but you provided: ",method,".\n"))
  if (method == 'auto') method <- "reuse"
  if (method == 'reuse') return(compute_moment.list(obj=obj$normalized_posterior,ff=ff,gg=gg,nn=nn,type=type,...))

  # Proceeding as if method = 'correct'
  if (is.null(nn) & is.null(gg) & is.null(ff)) stop("You provided NULL for ff, gg, and nn in compute_moment. You have to provide one of these.\n")
  if (!is.null(nn)) {
    if (type == 'central') {
      if (nn == 1) return(0) # Central moment of first order is zero by definition
      mm <- compute_moment(obj,nn=1,type='raw',method = 'correct')
    } else {
      mm <- 0
    }
    p <- get_param_dim(obj)
    out <- numeric(p)
    for (j in 1:p) {
      gg <- make_numeric_moment_function(nn,j,quad = obj,centre = mm,shift = NULL)
      out[j] <- compute_moment(obj,gg = gg,method = 'correct',type = 'raw') # Type=raw because already centred
      out[j] <- out[j] - get_shift(gg)
    }
  } else if (is.null(gg)) {
    if (!validate_moment(ff)) stop("Function ff provided to compute_moment fails validation. Run validate_moment(ff) to see why.\n")
    p <- get_param_dim(obj)
    q <- length(ff(get_mode(obj)))

    if (q>1) {
      # if q>1, run compute_moment on each coordinate
      out <- numeric(q)
      for (j in 1:q) out[j] <- compute_moment(obj=obj,ff=function(theta) ff(theta)[j],gg=NULL,method='correct')
    } else {
      # if q=1, construct gg and then call compute_moment on that
      out <- compute_moment(obj=obj,ff=NULL,gg=make_moment_function(function(theta) log(ff(theta))),method='correct')
    }
  } else {
    if(!validate_moment(gg)) stop("Function gg provided to compute_moment fails validation. Run validate_moment(gg) to see why.\n")
    gg <- make_moment_function(gg)
    # out <- compute_moment.list(obj=obj$normalized_posterior,ff=function(x) exp(gg$fn(x)),gg=NULL,method='reuse')
    newff <- list(
      fn = function(theta) obj$optresults$ff$fn(theta) + gg$fn(theta),
      gr = function(theta) obj$optresults$ff$gr(theta) + gg$gr(theta),
      he = function(theta) obj$optresults$ff$he(theta) + gg$he(theta)
    )
    newcontrol <- obj$control
    newcontrol$onlynormconst <- TRUE
    out <- exp(get_log_normconst(aghq(ff = newff,k = get_numquadpoints(obj),startingvalue = get_mode(obj),control = newcontrol)) - get_log_normconst(obj))
  }
  out
}
# compute_moment.aghq <- function(obj,ff = function(x) 1,...) compute_moment(obj,ff,...)
#' @rdname compute_moment
#' @export
compute_moment.default <- function(obj,ff = function(x) 1,gg = NULL,method = c("auto","reuse","correct"),...) {
  stop(paste0("Unrecognized object of class: ",class(obj)," passed to comupute_moment.\n"))
}

#' Correct the marginals of a fitted aghq object
#'
#' The default method of computation for aghq objects computes approximate marginals using an outdated algorithm
#' with no known theoretical properties. The updated algorithm computes pointwise approximate marginals that
#' satisfy the same rate of convergence as the original approximate marginal likelihood. This function takes
#' a fitted aghq object and recomputes its marginals using the updated algorithm
#'
#' @param quad An object of class \code{aghq} returned by \code{aghq::aghq()}.
#' @param ... Not used
#'
#' @return An object of class \code{aghq} equal to the provided object, but with its
#' \code{marginals} component replaced with one calculated using \code{method='correct'}.
#' See \code{marginal_posterior}.
correct_marginals <- function(quad,...) {
  p <- length(quad$marginals)
  margout <- vector(mode='list',length = p)
  for (j in 1:p) margout[[j]] <- marginal_posterior(quad,j,method = 'correct')
  quad$marginals <- margout
  quad
}

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
#' @param transformation Optional. Calculate pdf/cdf for a transformation of the parameter
#' whose posterior was normalized using adaptive quadrature.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#' @param interpolation Which method to use for interpolating the marginal posterior, \code{'polynomial'} (default)
#' or \code{'spline'}? If \code{k > 3} then the polynomial may be unstable and you should use the spline, but the spline
#' doesn't work *unless* \code{k > 3} so it's not the default. See \code{interpolate_marginal_posterior()}.
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
compute_pdf_and_cdf.default <- function(obj,transformation = default_transformation(),finegrid = NULL,interpolation = 'auto',...) {

  # if (!is.null(transformation)) transformation <- Map(match.fun,transformation)
  validate_transformation(transformation)
  transformation <- make_transformation(transformation)

  margpostinterp <- interpolate_marginal_posterior(obj,interpolation)

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

  # if (!is.null(transformation)) {
  #   if (is.null(transformation$jacobian)) {
  #     transformation$jacobian <- function(theta) {
  #       out <- numeric(length(theta))
  #       for (i in 1:length(theta)) {
  #         out[i] <- det(abs(numDeriv::jacobian(transformation$totheta,transformation$fromtheta(theta[i]))))
  #       }
  #       out
  #     }
  #   }
  out$transparam <- transformation$fromtheta(out$theta)
  out$pdf_transparam <- out$pdf * transformation$jacobian(out$theta)
  # }
  out
}
#' @rdname compute_pdf_and_cdf
#' @export
compute_pdf_and_cdf.list <- function(obj,transformation = default_transformation(),...) {
  out <- list()
  for (i in 1:length(obj)) out[[i]] <- compute_pdf_and_cdf(obj[[i]],transformation = transformation,...)
  out
}
#' @rdname compute_pdf_and_cdf
#' @export
compute_pdf_and_cdf.aghq <- function(obj,transformation = obj$transformation,...) compute_pdf_and_cdf(obj$marginals,transformation = transformation,...)

#' Quantiles
#'
#' Compute marginal quantiles using AGHQ. This function works by first approximating
#' the CDF using \code{aghq::compute_pdf_and_cdf} and then inverting the approximation numerically.
#'
#' @inheritParams compute_pdf_and_cdf
#' @param q Numeric vector of values in (0,1). The quantiles to compute.
#' @param transformation Optional. Calculate marginal quantiles for a transformation of the parameter
#' whose posterior was normalized using adaptive quadrature.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#' Note that since \code{g} has to be monotone anyways, this just returns \code{sort(g(q))} instead of \code{q}.
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
compute_quantiles <- function(obj,q = c(.025,.975),transformation = default_transformation(),...) UseMethod("compute_quantiles")
#' @rdname compute_quantiles
#' @export
compute_quantiles.default <- function(obj,q = c(.025,.975),transformation = default_transformation(),interpolation = 'auto',...) {
  # Compute the pdf and cdf WITHOUT the transformation
  # It's better to just transform the quantiles.
  pdfandcdf <- compute_pdf_and_cdf(obj,interpolation = interpolation)
  out <- numeric(length(q))
  increasing <- TRUE

  validate_transformation(transformation)
  transformation <- make_transformation(transformation)

  for (i in 1:length(q)) out[i] <- pdfandcdf$theta[max(which(pdfandcdf$cdf < q[i]))]

  # if (!is.null(transformation)) {
  #   if (is.null(transformation$fromtheta)) warning("transformation provided but transformation$fromtheta appears NULL.\n")
  #   transformation <- Map(match.fun,transformation)

    # Have to check if it's increasing or decreasing so can reverse order if necessary
  increasing <- transformation$fromtheta(min(out)) <= transformation$fromtheta(max(out))
  for (i in 1:length(out)) out[i] <- transformation$fromtheta(out[i])
  # }

  if (!increasing) out <- rev(out)
  names(out) <- paste0(as.character(100 * q),"%")
  out
}
#' @rdname compute_quantiles
#' @export
compute_quantiles.list <- function(obj,q = c(.025,.975),transformation = default_transformation(),...) {
  out <- list()
  for (i in 1:length(obj)) out[[i]] <- compute_quantiles(obj[[i]],q,transformation,...)
  out
}
#' @rdname compute_quantiles
#' @export
compute_quantiles.aghq <- function(obj,q = c(.025,.975),transformation = obj$transformation,...) compute_quantiles(obj$marginals,q,transformation,...)

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
#' @param transformation Optional. Draw samples for a transformation of the parameter
#' whose posterior was normalized using adaptive quadrature.
#' \code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
#' or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
#' @param interpolation Which method to use for interpolating the marginal posteriors
#' (and hence to draw samples using the inverse CDF method), \code{'auto'} (choose for you), \code{'polynomial'}
#' or \code{'spline'}? If \code{k > 3} then the polynomial may be unstable and you should use the spline, but the spline
#' doesn't work *unless* \code{k > 3} so it's not the default. The default of \code{'auto'} figures this out for you.
#' See \code{interpolate_marginal_posterior()}.
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
#' If \code{getOption('mc.cores',1L) > 1}, the Cholesky decompositions of the Hessians are computed
#' in parallel using \code{parallel::mcapply}, for the Gaussian approximation involved for objects of class \code{marginallaplace}. This step is slow
#' so may be sped up by parallelization, if the matrices are sparse (and hence the operation is just slow, but not memory-intensive).
#' Uses the \code{parallel} package so is not available on Windows.
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
sample_marginal <- function(quad,M,transformation = default_transformation(),interpolation = 'auto',...) UseMethod("sample_marginal")
#' @rdname sample_marginal
#' @export
sample_marginal.aghq <- function(quad,M,transformation = quad$transformation,interpolation = 'auto',...) {
  out <- list()
  if (is.null(quad$marginals)) return(out)
  for (i in 1:length(quad$marginals)) out[[i]] <- unname(compute_quantiles(quad$marginals[[i]],stats::runif(M),interpolation = interpolation,transformation = transformation))

  # if (!is.null(transformation)) {
  #   if (is.null(transformation$fromtheta)) warning("transformation provided but transformation$fromtheta appears NULL.\n")
  #   transformation <- Map(match.fun,transformation)
  #   for (i in 1:length(out)) out[[i]] <- transformation$fromtheta(out[[i]])
  # }
  out
}
#' @rdname sample_marginal
#' @export
sample_marginal.marginallaplace <- function(quad,M,transformation = quad$transformation,interpolation = 'auto',...) {
  numcores = getOption('mc.cores',1L)
  K <- as.numeric(quad$normalized_posterior$grid$level)[1]
  d <- dim(quad$modesandhessians$H[[1]])[1]
  simlist <- quad$modesandhessians
  # Avoid use of parallel computing on Windows
  if (.Platform$OS.type == 'windows') numcores <- 1
  if (numcores > 1) {
    # mclapply does not preserve the order of its arguments
    # simlist$L <- parallel::mclapply(simlist$H,function(h) chol(Matrix::forceSymmetric(h),perm = FALSE),mc.cores = numcores)
    simlist$L <- parallel::mclapply(simlist$H,function(h) Matrix::Cholesky(as(Matrix::forceSymmetric(h),'sparseMatrix'),perm = TRUE,LDL=FALSE),mc.cores = numcores)

  } else {
    # simlist$L <- lapply(simlist$H,function(h) chol(Matrix::forceSymmetric(h),perm = FALSE))
    simlist$L <- lapply(simlist$H,function(h) Matrix::Cholesky(as(Matrix::forceSymmetric(h),'sparseMatrix'),perm = TRUE,LDL=FALSE))

  }
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
    # function(.x,.y) as.numeric(solve(simlist$L[[as.numeric(.y)]],.x)) + do.call(cbind,rep(list(simlist$mode[[as.numeric(.y)]]),ncol(.x))),
    function(.x,.y) as.numeric(Matrix::solve(simlist$L[[as.numeric(.y)]],Matrix::solve(simlist$L[[as.numeric(.y)]],.x,system="Lt"),system='Pt')) + do.call(cbind,rep(list(simlist$mode[[as.numeric(.y)]]),ncol(.x))),
    Z,
    names(Z),
    SIMPLIFY = FALSE
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

  # BUG FIX: these lines were ad-hoc, and didn't return the correct columns
  # theta <- simlist[k, paste0("theta", seq(1, length(grep("theta",colnames(simlist)))))]
  # theta <- simlist[k, 1:(ncol(simlist)-5)]
  d <- length(quad$optresults$mode)
  thetanames <- colnames(quad$normalized_posterior$nodesandweights)[1:d]
  theta <- simlist[k,thetanames]
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
  out$thetasamples <- sample_marginal(quad,M,transformation,interpolation,...)
  out
}

#' Marginal Parameter Transformations
#'
#' These functions make it easier for the user to represent marginal parameter transformations
#' for which inferences are to be made. Suppose quadrature is done on the posterior for parameter \code{theta},
#' but interest lies in parameter \code{lambda = g(theta)} for smooth, monotone, univariate
#' \code{g}. This interface lets the user provide \code{g}, \code{g^-1}, and (optionally)
#' the jacobian \code{dtheta/dlambda}, and \code{aghq} will do quadrature on the \code{theta} scale
#' but report summaries on the \code{lambda} scale. See a note in the Details below about
#' multidimensional parameters.
#'
#' @param ... Used to pass arguments to methods.
#' @param fromtheta Function \code{g: R^p -> R^p}, where \code{p = dim(theta)}.
#' Must take vector \code{theta_1...theta_p} and return vector \code{g_1(theta_1)...g_p(theta_p)}, i.e.
#' only independent/marginal transformations are allowed (but these are the only ones
#' of interest, see below). For \code{j=1...p}, the parameter of
#' inferential interest is \code{lambda_j = g_j(theta_j)} and the parameter whose posterior is being
#' normalized via \code{aghq} is \code{theta_j}.  Passed to \code{match.fun}.
#' @param totheta Inverse function \code{g^-1(theta)}. Specifically, takes vector
#' \code{g_1(theta_1)...g_p(theta_p)} and returns vector \code{theta_1...theta_p}.
#' @param jacobian (optional) Function taking \code{theta} and returning the absolute value of the determinant of
#' the Jacobian \code{dtheta/dg(theta)}. If not provided, a numerically differentiated Jacobian is used as
#' follows: \code{numDeriv::jacobian(totheta,fromtheta(theta))}. Passed to \code{match.fun}.
#' @param translist A list with elements \code{totheta}, \code{fromtheta}, and, optionally, \code{jacobian}.
#' @param transobj An object of class \code{aghqtrans}. Just returns this object. This is for internal
#' compatibility.
#'
#' @return Object of class \code{aghqtrans}, which is simply a list with elements \code{totheta},
#' \code{fromtheta}, and \code{jacobian}. Object is suitable for checking with \code{aghq::validate_transformation}
#' and for inputting into any function in \code{aghq} which takes a \code{transformation} argument.
#'
#' @details Often, the scale on which quadrature is done is not the scale on which the user
#' wishes to make inferences. For example, when a parameter \code{lambda>0} is
#' of interest, the posterior for \code{theta = log(lambda)} may be better approximated
#' by a log-quadratic than that for \code{lambda}, so running \code{aghq} on the
#' likelihood and prior for \code{theta} may lead to faster and more stable optimization
#' as well as more accurate estimates. But, interest is still in the original parameter
#' \code{lambda = exp(theta)}.
#'
#' These considerations are by no means unique to the use of quadrature-based approximate
#' Bayesian inferences. However, when using (say) \code{MCMC}, inferences for summaries
#' of transformations of the parameter are just as easy as for the un-transformed parameter.
#' When using quadrature, a little bit more work is needed.
#'
#' The \code{aghq} package provides an interface for computing
#' posterior summaries of smooth, monotonic parameter transformations. If quadrature
#' is done on parameter \code{theta} and \code{g(theta)} is a univariate, smooth, monotone function,
#' then inferences are made for \code{lambda = g(theta)}. In the case that \code{theta} is
#' \code{p}-dimensional, \code{p > 1}, the supplied function \code{g} is understood to
#' take in \code{theta_1...theta_p} and return \code{g_1(theta_1)...g_p(theta_p)}. The
#' jacobian is diagonal.
#'
#' To reiterate, all of this discussion applies only to *marginal* parameter transformations.
#' For the full joint parameter, the only summary statistics you can even calculate at all
#' (at present?) are moments, and you can already calculate the moment of any function \code{h(theta)}
#' using \code{aghq::compute_moment}, so no additional interface is needed here.
#'
#' @family transformations
#'
#' @examples
#'
#' make_transformation('log','exp')
#' make_transformation(log,exp)
#'
#' @export
make_transformation <- function(...) UseMethod('make_transformation')
#' @rdname make_transformation
#' @export
make_transformation.aghqtrans <- function(transobj,...) transobj
#' @rdname make_transformation
#' @export
make_transformation.list <- function(translist,...) make_transformation.default(translist$totheta,translist$fromtheta,translist$jacobian,...)
#' @rdname make_transformation
#' @export
make_transformation.default <- function(totheta,fromtheta,jacobian = NULL,...) {
  totheta <- match.fun(totheta)
  fromtheta <- match.fun(fromtheta)

  trans <- list(totheta=totheta,fromtheta=fromtheta)
  if (is.null(jacobian)) {
    trans$jacobian <- function(theta) {
      out <- numeric(length(theta))
      for (i in 1:length(theta)) {
        out[i] <- det(abs(numDeriv::jacobian(totheta,fromtheta(theta[i]))))
      }
      out
    }
  } else {
    trans$jacobian <- match.fun(jacobian)
  }
  class(trans) <- 'aghqtrans'
  trans
}

#' Validate a transformation object
#'
#' Routine for checking whether a given transformation is valid.
#'
#' @param ... Used to pass arguments to methods.
#' @param trans A transformation object of class \code{aghqtrans} returned by \code{make_transformation}.
#' @param translist A list. Will be checked, passed to \code{aghqtrans}, and then checked again.
#' @param checkinverse Default \code{FALSE}, do not check that \code{totheta(fromtheta(theta)) = theta}. Otherwise,
#' a vector of values for which to perform that check. No default values are provided, since \code{validate_transformation}
#' has no way of determining the domain and range of \code{totheta} and \code{fromtheta}. This argument is
#' used internally in \code{aghq} package functions, with cleverly chosen check values.
#'
#' @details This function checks that:
#' \itemize{
#' \item{The supplied object contains elements \code{totheta}, \code{fromtheta}, and \code{jacobian}, and that they are all functions,}
#' \item{If \code{checkinverse} is a vector of numbers, then it checks that \code{totheta(fromtheta(checkinverse)) == checkinverse}.}
#' }
#' In addition, if a \code{list} is provided, the function first checks that it contains the right elements,
#' then passes it to \code{make_transformation}, then checks that.
#'
#' This function throws an informative error messages when checks don't pass or themselves throw errors.
#'
#' @return \code{TRUE} if the function runs to completion without throwing an error.
#'
#' @family transformations
#'
#' @examples
#'
#' t <- make_transformation(log,exp)
#' validate_transformation(t)
#' t2 <- list(totheta = log,fromtheta = exp)
#' validate_transformation(t2)
#' \dontrun{
#' t3 <- make_transformation(log,log)
#' checkvals <- exp(exp(rnorm(10)))
#' # Should throw an informative error because log isn't the inverse of log.
#' validate_transformation(t3,checkinverse = checkvals)
#' }
#' @export
validate_transformation <- function(...) UseMethod('validate_transformation')
#' @rdname validate_transformation
#' @export
validate_transformation.aghqtrans <- function(trans,checkinverse = FALSE,...) {
  # Check names
  valid_names <- c("totheta","fromtheta","jacobian")
  if (!all(names(trans) == valid_names)) {
    stop(paste0("Transformation object should have names: ",valid_names,". The provided object has names: ",names(trans),"."))
  }
  # Check functions
  for (fn in names(trans)) {
    if (!is.function(trans[[fn]])) stop("Transformation object should be a list of functions, but element ",fn," inherits from class ",class(fn),".")
  }
  # Check inverse
  if (checkinverse[1] | (!is.logical(checkinverse[1]) & checkinverse[1]==0)) {
    # Compute the values and check they are all numeric
    fromthetavals <- trans$fromtheta(checkinverse)
    if (any(is.infinite(fromthetavals))) stop("trans$fromtheta(checkinverse) produced infinite values. Check the domain and range of totheta and fromtheta.")
    if (any(is.na(fromthetavals))) stop("trans$fromtheta(checkinverse) produced NA values. Check the domain and range of totheta and fromtheta.")
    if (any(is.nan(fromthetavals))) stop("trans$fromtheta(checkinverse) produced NaN values. Check the domain and range of totheta and fromtheta.")
    tothetavals <- trans$totheta(fromthetavals)
    if (any(is.infinite(tothetavals))) stop("trans$totheta(trans$fromtheta(checkinverse)) produced infinite values. Check the domain and range of totheta and fromtheta.")
    if (any(is.na(tothetavals))) stop("trans$totheta(trans$fromtheta(checkinverse)) produced NA values. Check the domain and range of totheta and fromtheta.")
    if (any(is.nan(tothetavals))) stop("trans$totheta(trans$fromtheta(checkinverse)) produced NaN values. Check the domain and range of totheta and fromtheta.")
    if (!all(tothetavals == checkinverse)) stop("Elements totheta and fromtheta do not appear to be each others' inverse functions. Please check this.")
  }
  TRUE
}
#' @rdname validate_transformation
#' @export
validate_transformation.list <- function(translist,checkinverse = FALSE,...) {
  # First check it has the required elements, then create a transformation argument
  # and check that
  valid_names <- c("totheta","fromtheta") # Don't need it to have a jacobian
  if (all(valid_names %in% names(translist))) {
    return(validate_transformation(trans = make_transformation(translist)))
  }
  stop(paste0("Transformation object should have names: ",valid_names,". The provided object has names: ",names(translist),"."))
}
#' @rdname validate_transformation
#' @export
validate_transformation.default <- function(...) FALSE

#' Default transformation
#'
#' Default (identity) transformation object. Default argument in package functions
#' which accept transformations, and useful for user inspection.
#'
#' @family transformations
#'
#' @examples
#'
#' default_transformation()
#'
#' @export
default_transformation <- function() make_transformation(totheta = force,fromtheta = force)


#' Moments of Positive Functions
#'
#' Given an object \code{quad} of class \code{aghq} returned by \code{aghq::aghq()}, \code{aghq::compute_moment()}
#' will compute the moment of a positive function \code{g(theta)} of parameter \code{theta}. The present function,
#' \code{aghq::make_moment_function()}, assists the user in constructing the appropriate input to \code{aghq::compute_moment()}.
#'
#' @param ... Used to pass arguments to methods.
#' @param gg LOGARITHM of function `R^p -> R^+` of which the moment is to be computed along with its two derivatives. So for example providing gg = function(x) x will compute the moment of exp(x).
#' Provided either as a \code{function}, a \code{list}, an \code{aghqtrans} object, or an \code{aghqmoment} object. See details.
#'
#' @return Object of class \code{aghqmoment}, which is a list with elements \code{fn},
#' \code{gr}, and \code{he}, exactly like the input to \code{aghq::aghq()} and related functions. Here \code{gg$fn} is
#' \code{log(gg(theta))}, \code{gg$gr} is its gradient, and \code{gg$he} its Hessian.
#' Object is suitable for checking with \code{aghq::validate_moment()}
#' and for inputting into \code{aghq::compute_moment()}.
#'
#' @details The approximation of moments of positive functions implemented in \code{aghq::compute_moment()}
#' achieves the same asymptotic rate of convergence as the marginal likelihood. This involves computing a new mode and
#' Hessian depending on the original posterior mode and Hessian, and `g`. These computations are handled by \code{aghq::compute_moment()},
#' re-using information from the original quadrature when feasible.
#'
#' Computation of moments is defined only for scalar-valued functions, with a vector moment just defined as a vector of moments. Consequently,
#' the user may input to \code{aghq:compute_moment()} a function \code{g: R^p -> R^q+} for any \code{q}, and that function will return the corresponding
#' vector of moments. This is handled within \code{aghq::compute_moment()}. The \code{aghq::make_moment_function()} interface accepts the logarithm of \code{gg: R^p -> R^+}, i.e.
#' a multivariable, scalar-valued positive function. This is mostly to keep first and second derivatives as 1d and 2d arrays (i.e. the gradient and the Hessian);
#' I deemed it too confusing for the user and the code-base to work with Jacobians and 2nd derivative tensors (if you're confused just reading this, there you go!).
#' But, see \code{aghq::compute_moment()} for how to handle the very common case where the *same* trasnformation is desired of all parameter coordinates; for example
#' when all parameters are on the log-scale and you want to compute \code{E(exp(theta))} for *vector* \code{theta}.
#'
#' If \code{gg} is a \code{function} or a \code{character} (like \code{'exp'}) it is first passed to \code{match.fun}, and then the output
#' object is constructed using \code{numDeriv::grad()} and \code{numDeriv::hessian()}. If \code{gg} is a \code{list} then it is assumed to
#' have elements \code{fn}, \code{gr}, and \code{he} of the correct form, and these elements are extracted and then passed back to \code{make_moment_function()}.
#' If \code{gg} is an object of class \code{aghqtrans} returned by \code{aghq::make_transformation()}, then \code{gg$fromtheta}
#' is passed back to \code{make_moment_function()}. If \code{gg} is an object of class \code{aghqtrans} then it is just returned.
#'
#' @family moments
#'
#' @examples
#' # E(exp(x))
#' mom1 <- make_moment_function(force) # force = function(x) x
#' mom2 <- make_moment_function('force')
#' mom3 <- make_moment_function(list(fn=function(x) x,gr=function(x) 1,he = function(x) 0))
#'
#' @export
make_moment_function <- function(...) UseMethod('make_moment_function')
#' @rdname make_moment_function
#' @export
make_moment_function.aghqmoment <- function(gg,...) gg
#' @rdname make_moment_function
#' @export
make_moment_function.aghqtrans <- function(gg,...) make_moment_function.default(function(theta) log(gg$from_theta(theta)),...)
#' @rdname make_moment_function
#' @export
make_moment_function.function <- function(gg,...) {
  gg <- match.fun(gg)
  fn <- gg
  gr <- function(theta) numDeriv::grad(fn,theta)
  he <- function(theta) numDeriv::hessian(fn,theta)
  make_moment_function.list(list(fn=fn,gr=gr,he=he),...)
}
#' @rdname make_moment_function
#' @export
make_moment_function.character <- function(gg,...) {
  gg <- match.fun(gg)
  make_moment_function.function(gg,...)
}
#' @rdname make_moment_function
#' @export
make_moment_function.list <- function(gg,...) {
  # Check it has the correct form before passing to default method
  nms <- c("fn","gr","he")
  if (!all(nms%in%names(gg)))
    stop(paste0("make_moment_function requires a list including elements: ",nms,".\nThe provided list has names: ",names(gg),".\nIt does not have elements: ",setdiff(names(gg),nms),".\n"))
  if (!all(Reduce(c,lapply(gg[nms],inherits,what='function')))) {
    classes <- lapply(gg[nms],class)
    names(classes) <- nms
    stop(paste0("make_moment_function requires a list of functions. The elements of your list are: ",classes,".\n"))
  }
  make_moment_function.default(gg,...)
}
#' @rdname make_moment_function
#' @export
make_moment_function.default <- function(gg,...) {
  # Nothing left to do!
  class(gg) <- "aghqmoment"
  gg
}

#' Validate a moment function object
#'
#' Routine for checking whether a given moment function object is valid.
#'
#' @param ... Used to pass arguments to methods.
#' @param moment An object to check if it is a valid moment function or not. Can be an object of class \code{aghqmoment} returned by \code{aghq::make_moment_function()},
#' or any object that can be passed to \code{aghq::make_moment_function()}.
#' @param checkpositive Default \code{FALSE}, do not check that \code{gg$fn(theta) > 0}. Otherwise,
#' a vector of values for which to perform that check. No default values are provided, since \code{validate_moment}
#' has no way of determining the domain and range of \code{gg$fn}. This argument is
#' used internally in \code{aghq} package functions, with cleverly chosen check values.
#'
#' @details This function checks that:
#' \itemize{
#' \item{The supplied object contains elements \code{fn}, \code{gr}, and \code{he}, and that they are all functions,}
#' \item{If \code{checkpositive} is a vector of numbers, then it checks that \code{gg$fn(checkpositive)} is not \code{-Inf}, \code{NA}, or \code{NaN}. (It actually uses \code{is.infinite} for the first.)}
#' }
#' In addition, if a \code{list} is provided, the function first checks that it contains the right elements,
#' then passes it to \code{make_moment_function}, then checks that. If a \code{function} or a \code{character} is provided,
#' it checks that \code{match.fun} works, and returns any errors or warnings from doing so in a clear way.
#'
#' This function throws an informative error messages when checks don't pass or themselves throw errors.
#'
#' @return \code{TRUE} if the function runs to completion without throwing an error.
#'
#' @family moments
#'
#' @examples
#'
#' mom1 <- make_moment_function(exp)
#' mom2 <- make_moment_function('exp')
#' mom3 <- make_moment_function(list(fn=function(x) x,gr=function(x) 1,he = function(x) 0))
#' validate_moment(mom1)
#' validate_moment(mom2)
#' validate_moment(mom3)
#' \dontrun{
#' mombad1 <- list(exp,exp,exp) # No names
#' mombad2 <- list('exp','exp','exp') # List of not functions
#' mombad3 <- make_moment_function(function(x) -exp(x)) # Not positive
#' validate_moment(mombad1)
#' validate_moment(mombad2)
#' validate_moment(mombad3)
#' }
#' @export
validate_moment <- function(...) UseMethod('validate_moment')
#' @rdname validate_moment
#' @export
validate_moment.aghqmoment <- function(moment,checkpositive = FALSE,...) {
  # Check names
  valid_names <- c("fn","gr","he","shift")
  if (!all(names(moment) %in% valid_names)) {
    stop(paste0("Moment object should have names: ",paste0(valid_names,sep=","),
                ". The provided object has names: ",paste0(names(moment),sep=","),"."))  }
  # Check functions
  function_names <- c("fn","gr","he")
  for (fn in names(moment)) {
    if (fn %in% function_names & !is.function(moment[[fn]])) stop("Moment object should be a list of functions, but element ",fn," inherits from class ",class(fn),".")
  }
  # Check positive
  if (checkpositive[1] | (!is.logical(checkpositive[1]) & checkpositive[1]==0)) {
    # Compute the values and check they are all numeric
    suppressWarnings(momentvals <- moment$fn(checkpositive))
    if (any(is.infinite(momentvals))) stop("moment$fn(checkpositive) produced infinite values. Check that the function used to create your moment object returns only positive values.")
    if (any(is.na(momentvals))) stop("moment$fn(checkpositive) produced NA values. Check that the function used to create your moment object returns only positive values.")
    if (any(is.nan(momentvals))) stop("moment$fn(checkpositive) produced NaN values. Check that the function used to create your moment object returns only positive values.")
  }
  TRUE
}
#' @rdname validate_moment
#' @export
validate_moment.list <- function(moment,checkpositive = FALSE,...) {
  # First check it has the required elements, then create a moment object
  # and check that
  valid_names <- c("fn","gr","he")
  if (all(valid_names %in% names(moment))) {
    return(validate_moment(make_moment_function(moment)))
  }
  stop(paste0("Moment object should have names: ",paste0(valid_names,sep=","),
              ". The provided object has names: ",paste0(names(moment),sep=","),"."))
}
#' @rdname validate_moment
#' @export
validate_moment.function <- function(moment,checkpositive = FALSE,...) {
  e <- tryCatch(match.fun(moment),warning = function(w) w,error = function(e) e)
  if (inherits(e,'condition')) stop(paste0("Tried calling match.fun on your moment function, and returned the following condition: ",e,".\n"))
  validate_moment(make_moment_function(moment))
}
#' @rdname validate_moment
#' @export
validate_moment.character <- function(moment,checkpositive = FALSE,...) {
  e <- tryCatch(match.fun(moment),warning = function(w) w,error = function(e) e)
  if (inherits(e,'condition')) stop(paste0("Tried calling match.fun on your moment function, and returned the following condition: ",e,".\n"))
  validate_moment(make_moment_function(moment))
}
#' @rdname validate_moment
#' @export
validate_moment.default <- function(moment,...) {
  if (is.null(moment)) stop("NULL function provided to validate_moment.\n")
  FALSE
}

#' Compute numeric moments
#'
#' Create a function suitable for computation of numeric moments. This function is
#' used internally by \code{compute_moment} when the user chooses \code{nn}, and is
#' unlikely to need to be called by a user directly.
#'
#' @param nn Order of moment to be computed, see \code{nn} argument of \code{compute_moment}.
#' @param quad Optional, object of class \code{aghq}, only used if \code{shift} is not \code{NULL}.
#' @param j Numeric, positive integer, index of parameter vector to compute the numeric
#' moment for.
#' @param centre Numeric scalar, added to \code{shift} to ensure that central moments
#' remain far from zero.
#' @param shift Numeric scalar, amount by which to shift \code{theta}. The function that this
#' outputs is \code{g(theta) = (theta)^nn + shift}, and \code{shift} is returned with the
#' object so that it may later be subtracted. Default of \code{NULL} chooses this value
#' internally.
#' @param gg Object of class \code{aghqmoment}. Returns the \code{shift} applied to
#' the moment function. Returns \code{0} if no shift applied.
#' @param ... Not used.
#'
#' @return Object of class \code{aghqmoment}, see \code{make_moment_function}
#'
#' @export
#'
make_numeric_moment_function <- function(nn,j,quad = NULL,centre = 0,shift = NULL,...) {
  if (nn==0) return(make_moment_function(function(x) 1))
  if (nn < 0) stop("Moments of order < 0 currently not supported by the numeric moment interface. You could construct this yourself using make_moment_function.")
  if (!is.null(quad)) {
    if (!inherits(quad,'aghq')) stop(paste0("You must provide an aghq quadrature object. Your quad object has class: ",class(quad),".\n"))
  }
  if (as.integer(j) != j) stop(paste0("You must provide an integer index, j. You provided: ",j,"\n."))
  p <- get_param_dim(quad)
  if (j < 1 | j > p) stop(paste0("You must provide an index j such that 1 <= j <= dim(theta). You provided: ",j,", but dim(theta) = ",p,"\n."))
  if (centre==0 & 'center' %in% names(list(...))) centre <- list(...)$center # For the Americans
  if (is.null(shift)) {
    if (is.null(quad)) {
      shift <- 0 # Can't do much if no information provided
    } else {
      # Get the shift
      nodes <- as.numeric(get_nodesandweights(quad)[ ,j]) # Take the jth column
      buffer <- 10*diff(range(nodes)) # Pretty arbitrary, until someone comes up with something better
      shift <- -1*min(nodes) + buffer + centre # This will be ADDED
    }
  }
  # Create the function
  fn <- function(theta) log((theta-centre)^nn + shift) # Supposed to be positive...
  gr <- function(theta) ( nn * (theta - centre)^(nn-1) ) / ( (theta-centre)^nn + shift )
  if (nn == 1) {
    # Avoid division by 0 when theta = centre
    he <- function(theta) -1 / ( (theta-centre)^nn + shift )^2
  } else {
    he <- function(theta) ( -nn*(theta-centre)^(nn-1) ) / ( (theta-centre)^nn + shift )^2 + ( nn*(nn-1)*(theta-centre)^(nn-2) ) / ( (theta-centre)^nn + shift )
  }
  gg <- make_moment_function(list(fn=fn,gr=gr,he=he))
  gg$shift <- shift
  gg
}
#' @rdname make_numeric_moment_function
#' @export
get_shift <- function(gg) {
  if (is.null(gg$shift)) return(0)
  gg$shift
}


