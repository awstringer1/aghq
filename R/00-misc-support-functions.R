### Misc functions ###
# This script contains miscillaneous, non-exported functions for supporting the functions
# in the AGHQ package.

default_control <- function() {
  list(
    method = c("sparse_trust","trust","BFGS")
  )
}