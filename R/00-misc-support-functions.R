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

# Logsumexp
# From: https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp
# Accessed on: 2021/03/27 10:49AM
logsumexp <- function(l) {
  n <- length(l)
  L <- sort(l, decreasing = TRUE)
  S <- rep(L[1], n)
  for (k in 1:(n-1)) {
    S[k+1] <- max(L[k+1], S[k]) + log1p(exp(-abs(L[k+1] - S[k])))
  }
  S[n]
}