% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-summary-statistics.R
\name{compute_quantiles}
\alias{compute_quantiles}
\alias{compute_quantiles.default}
\alias{compute_quantiles.list}
\alias{compute_quantiles.aghq}
\title{Quantiles}
\usage{
compute_quantiles(
  obj,
  q = c(0.025, 0.975),
  transformation = default_transformation(),
  ...
)

\method{compute_quantiles}{default}(
  obj,
  q = c(0.025, 0.975),
  transformation = default_transformation(),
  interpolation = "auto",
  ...
)

\method{compute_quantiles}{list}(
  obj,
  q = c(0.025, 0.975),
  transformation = default_transformation(),
  ...
)

\method{compute_quantiles}{aghq}(
  obj,
  q = c(0.025, 0.975),
  transformation = obj$transformation,
  ...
)
}
\arguments{
\item{obj}{Either the output of \code{aghq::aghq()}, its list of marginal distributions
(element \code{marginals}), or an individual \code{data.frame} containing one of
these marginal distributions as output by \code{aghq::marginal_posterior()}.}

\item{q}{Numeric vector of values in (0,1). The quantiles to compute.}

\item{transformation}{Optional. Calculate marginal quantiles for a transformation of the parameter
whose posterior was normalized using adaptive quadrature.
\code{transformation} is either: a) an \code{aghqtrans} object returned by \code{aghq::make_transformation},
or b) a list that will be passed to that function internally. See \code{?aghq::make_transformation} for details.
Note that since \code{g} has to be monotone anyways, this just returns \code{sort(g(q))} instead of \code{q}.}

\item{...}{Used to pass additional arguments.}

\item{interpolation}{Which method to use for interpolating the marginal posterior, \code{'polynomial'} (default)
or \code{'spline'}? If \code{k > 3} then the polynomial may be unstable and you should use the spline, but the spline
doesn't work \emph{unless} \code{k > 3} so it's not the default. See \code{interpolate_marginal_posterior()}.}
}
\value{
A named numeric vector containing the quantiles you asked for, for the
variable whose marginal posterior you provided.
}
\description{
Compute marginal quantiles using AGHQ. This function works by first approximating
the CDF using \code{aghq::compute_pdf_and_cdf} and then inverting the approximation numerically.
}
\examples{
logfteta2d <- function(eta,y) {
  # eta is now (eta1,eta2)
  # y is now (y1,y2)
  n <- length(y)
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  y1 <- y[1:n1]
  y2 <- y[(n1+1):(n1+n2)]
  eta1 <- eta[1]
  eta2 <- eta[2]
  sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
    sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
}
set.seed(84343124)
n1 <- 5
n2 <- 5
n <- n1+n2
y1 <- rpois(n1,5)
y2 <- rpois(n2,5)
objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
funlist2d <- list(
  fn = objfunc2d,
  gr = function(x) numDeriv::grad(objfunc2d,x),
  he = function(x) numDeriv::hessian(objfunc2d,x)
)
opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
margpost <- marginal_posterior(opt_sparsetrust_2d,3,1) # margpost for theta1
etaquant <- compute_quantiles(margpost)
etaquant
# lambda = exp(eta)
exp(etaquant)
# Compare to truth
qgamma(.025,1+sum(y1),1+n1)
qgamma(.975,1+sum(y1),1+n1)



}
\seealso{
Other summaries: 
\code{\link{compute_pdf_and_cdf}()},
\code{\link{interpolate_marginal_posterior}()},
\code{\link{marginal_posterior}()}
}
\concept{summaries}
