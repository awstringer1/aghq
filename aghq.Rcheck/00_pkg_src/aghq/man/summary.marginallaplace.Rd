% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04-aghq.R
\name{summary.marginallaplace}
\alias{summary.marginallaplace}
\title{Summary statistics for models using marginal Laplace approximations}
\usage{
\method{summary}{marginallaplace}(object, M = 1000, max_print = 30, ...)
}
\arguments{
\item{object}{Object inheriting from \strong{both} classes \code{aghq} and \code{marginallaplace},
for example as returned by \code{aghq::marginal_laplace} or \code{aghq::marginal_laplace_tmb}.}

\item{M}{Number of samples to use to compute summary statistics of the random effects.
Default \code{1000}. Lower runs faster, higher is more accurate.}

\item{max_print}{Sometimes there are a lot of random effects. If there are more random
effects than \code{max_print}, the random effects aren't summarized, and a note is printed
to this effect. Default \code{30}.}

\item{...}{not used.}
}
\value{
A list containing an object of class \code{aghqsummary} (see \code{summary.aghq})
plus the following objects:
}
\description{
The \code{summary.marginallaplace} calls \code{summary.aghq}, but also computes
summary statistics of the random effects, by drawing from their approximate
posterior using \code{aghq::sample_marginal} with the specified number
of samples.
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
objfunc2dmarg <- function(W,theta) objfunc2d(c(W,theta))
objfunc2dmarggr <- function(W,theta) {
  fn <- function(W) objfunc2dmarg(W,theta)
  numDeriv::grad(fn,W)
}
objfunc2dmarghe <- function(W,theta) {
  fn <- function(W) objfunc2dmarg(W,theta)
  numDeriv::hessian(fn,W)
}

funlist2dmarg <- list(
  fn = objfunc2dmarg,
  gr = objfunc2dmarggr,
  he = objfunc2dmarghe
)

themarginallaplace <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0))
summary(themarginallaplace)
}
\seealso{
Other quadrature: 
\code{\link{aghq}()},
\code{\link{get_log_normconst}()},
\code{\link{laplace_approximation}()},
\code{\link{marginal_laplace_tmb}()},
\code{\link{marginal_laplace}()},
\code{\link{normalize_logpost}()},
\code{\link{optimize_theta}()},
\code{\link{plot.aghq}()},
\code{\link{print.aghqsummary}()},
\code{\link{print.aghq}()},
\code{\link{print.laplacesummary}()},
\code{\link{print.laplace}()},
\code{\link{print.marginallaplacesummary}()},
\code{\link{summary.aghq}()},
\code{\link{summary.laplace}()}
}
\concept{quadrature}
