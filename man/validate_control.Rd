% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00-misc-support-functions.R
\name{validate_control}
\alias{validate_control}
\title{Validate a control list}
\usage{
validate_control(control, type = c("aghq", "marglaplace", "tmb"), ...)
}
\arguments{
\item{control}{A list.}

\item{type}{One of \code{c('aghq','marglapace','tmb')}. The type of control object to validate. Will
basically validate against the arguments required by \code{aghq}, \code{marginal_laplace}, and \code{marginal_laplace_tmb},
respectively.}

\item{...}{Not used.}
}
\value{
Logical, \code{TRUE} if the list of control arguments is valid, else \code{FALSE}.
}
\description{
This function checks that the names and value types for a supplied \code{control} list
are valid and are unlikely to cause further errors within \code{aghq} and related functions.
Users should not have to worry about this and should just use \code{default_control()} and related
constructors.
}
\details{
To users reading this: just use \code{default_control()}, \code{default_control_marglaplace()}, or \code{default_control_tmb()}
as appropriate, to ensure that your control arguments are correct. This function just exists to provide more
descriptive error messages in the event that an incompatible list is provided.
}
\examples{
validate_control(default_control())
validate_control(default_control_marglaplace(),type = "marglaplace")
validate_control(default_control_tmb(),type = "tmb")

}
