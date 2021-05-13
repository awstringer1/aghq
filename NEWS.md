# aghq 0.2.0

## New features

- Added new function `marginal_laplace_tmb()` which improves compatibility with the `TMB` package.

- Turned all the main functions for computing posterior summaries into `S3` methods
with methods for objects of class `aghq`.

## Bug fixes

- Fixed `default_control_tmb()` which was missing the `ndConstruction = 'product'` default argument
which caused the function to throw a cryptic error.

## Other

- Removed external dependencies:
  
  - `matrixStats`: the only function I had used was `logSumExp`, which I replaced
  with a slower but completely base `R` version from https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp (accessed on 2021/03/27) in order to remove this dependency. This does not have a meaningful effect on the performance of `aghq`.
  
  - All 'tidyverse' subsidiary packages.
