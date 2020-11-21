### Misc functions ###
# This script contains miscillaneous, non-exported functions for supporting the functions
# in the AGHQ package.

default_control <- function() {
  list(
    method = c("sparse_trust","trust","BFGS")
  )
}

default_control_marglaplace <- function() {
  list(
    method = c("BFGS","sparse_trust","trust"),
    inner_method = c("sparse_trust","trust","BFGS")
  )
}