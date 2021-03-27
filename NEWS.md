# aghq 0.2.0

## New features

## Bug fixes

## Other

- Removed external dependencies:
  
  - `matrixStats`: the only function I had used was `logSumExp`, which I replaced
  with a slower but completely base `R` version from https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp (accessed on 2021/03/27) in order to remove this dependency. This does not have a meaningful effect on the performance of `aghq`.