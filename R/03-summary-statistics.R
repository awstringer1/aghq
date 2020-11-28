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
    out <- dplyr::as_tibble(out)
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

  nodesandweightsfactored <- cbind(nn,wwE) %>% dplyr::inner_join(nodesandweights,by = c(thetaj,thetaminusj))

  out <- nodesandweightsfactored %>%
    dplyr::mutate(id = 1:(dplyr::n())) %>%
    tidyr::pivot_longer(tidyselect::one_of(paste0(thetaminusj,"W")),names_to = "v",values_to = "w") %>%
    dplyr::group_by(.data[[thetaj]],.data$id) %>%
    dplyr::summarize(logw = sum(log(.data$w)) + sum(log(.env$diagcholinvH[-1])),
                     logpost_normalized = mean(.data$logpost_normalized)) %>%
    dplyr::summarize(logmargpost = matrixStats::logSumExp(.data$logw + .data$logpost_normalized))


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
  colnames(margpost)[stringr::str_detect(colnames(margpost),"theta")] <- "theta"
  as.function(polynom::poly.calc(x = margpost$theta,y = margpost$logmargpost))
}

#' Compute moments
#'
#' Compute the moment of any function ff using AGHQ.
#'
#' @param normalized_posterior The output of \code{aghq::normalize_logpost}. See the documentation for that function.
#' @param ff Any R function which takes in a numeric vector and returns a numeric vector.
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
#'
#' @importFrom rlang .data .env
#' @family summaries
#' @export
#'
compute_moment <- function(normalized_posterior,ff = function(x) 1) {
  nodesandweights <- normalized_posterior$nodesandweights

  whereistheta <- stringr::str_detect(colnames(nodesandweights),'theta')

  lengthof_f <- length(ff(nodesandweights[1,whereistheta]))

  if (lengthof_f == 1) {
    out <- sum(ff(nodesandweights[ ,whereistheta])* exp(nodesandweights$logpost_normalized) * nodesandweights$weights)
  } else {
    out <- apply(nodesandweights[ ,whereistheta],1,ff) %>%
      apply(1,function(x) sum(x * exp(nodesandweights$logpost_normalized) * nodesandweights$weights))
  }

  # thetastodo <- colnames(nodesandweights)[stringr::str_detect(colnames(nodesandweights),"theta")]
  #
  # out <- numeric(length(thetastodo))
  # names(out) <- thetastodo
  # for (theta in thetastodo) {
  #   out[theta] <- sum(ff(nodesandweights[[theta]]) * nodesandweights$weights * exp(nodesandweights$logpost_normalized))
  # }

  unname(out)
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
#' @param margpost The output of \code{aghq::marginal_posterior}. See the documentation for that function.
#' @param finegrid Optional, a grid of values on which to compute the CDF. The default makes
#' use of the values in \code{margpost} but if the results are unsuitable, you may wish to
#' modify this manually.
#' @param transformation Optional. A list containing two functions, \code{fromtheta}
#' and \code{totheta}, which accept and return numeric vectors, defining a parameter transformation for which you would
#' also like the pdf calculated for. See examples. May also have an element \code{jacobian},
#' a function which takes a numeric vector and computes the jacobian of the transformation; if
#' not provided, this is done using \code{numDeriv::jacobian}.
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
#' par(mfrow = c(1,2))
#' with(thepdfandcdf,{
#'   plot(pdf~theta,type='l')
#'   plot(cdf~theta,type='l')
#' })
#'
#' @export
#'
compute_pdf_and_cdf <- function(margpost,transformation = NULL,finegrid = NULL) {

  margpostinterp <- interpolate_marginal_posterior(margpost)

  thetacol <- colnames(margpost)[stringr::str_detect(colnames(margpost),"theta")]
  if (is.null(finegrid)) {
    rn <- range(margpost[[thetacol]])
    rnl <- diff(rn)
    thetarange <- c(min(rn) - rnl/2,max(rn) + rnl/2)
    finegrid <- c(seq(thetarange[1],thetarange[2],length.out=1000))
  }

  out <- dplyr::tibble(theta = finegrid,
                 pdf = exp(margpostinterp(.data$theta)),
                 cdf = cumsum(.data$pdf * c(0,diff(.data$theta))))

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
    out <- out %>%
      dplyr::mutate(transparam = transformation$fromtheta(.data[["theta"]]),
             pdf_transparam = .data[["pdf"]] * transformation$jacobian(.data[["theta"]]))
  }
  out
}

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
compute_quantiles <- function(margpost,q = c(.025,.975)) {
  pdfandcdf <- compute_pdf_and_cdf(margpost)
  out <- numeric(length(q))
  names(out) <- paste0(as.character(100 * q),"%")

  for (i in 1:length(q)) {
    out[i] <- pdfandcdf$theta[max(which(pdfandcdf$cdf < q[i]))]
  }
  out
}

# TODO: examples
#