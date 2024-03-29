% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02-normalization.R
\name{get_param_dim}
\alias{get_param_dim}
\alias{get_param_dim.aghq}
\title{Obtain the parameter dimension from an aghq object}
\usage{
get_param_dim(obj, ...)

\method{get_param_dim}{aghq}(obj, ...)
}
\arguments{
\item{obj}{Object of class \code{aghq} returned by \code{aghq::aghq}.}

\item{...}{Not used}
}
\value{
A numeric vector of length 1 containing \code{p}, the parameter dimension.
}
\description{
Quick helper S3 method to retrieve the parameter dimension from an aghq object.
}
\seealso{
Other quadrature: 
\code{\link{aghq}()},
\code{\link{get_hessian}()},
\code{\link{get_log_normconst}()},
\code{\link{get_mode}()},
\code{\link{get_nodesandweights}()},
\code{\link{get_numquadpoints}()},
\code{\link{get_opt_results}()},
\code{\link{laplace_approximation}()},
\code{\link{marginal_laplace_tmb}()},
\code{\link{marginal_laplace}()},
\code{\link{nested_quadrature}()},
\code{\link{normalize_logpost}()},
\code{\link{optimize_theta}()},
\code{\link{plot.aghq}()},
\code{\link{print.aghqsummary}()},
\code{\link{print.aghq}()},
\code{\link{print.laplacesummary}()},
\code{\link{print.laplace}()},
\code{\link{print.marginallaplacesummary}()},
\code{\link{summary.aghq}()},
\code{\link{summary.laplace}()},
\code{\link{summary.marginallaplace}()}
}
\concept{quadrature}
