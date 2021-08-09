# aghq 0.3.0 

## New features

- Added support for doing the multiple required Cholesky decompositions in parallel in `sample_marginal.marginallaplace`.

- Switched from using `chol` to using `Matrix::Cholesky` with `perm = TRUE` inside `sample_marginal.marginallaplace`. This
uses fill-reducing permutations (described in the `Matrix` package documentation) and users have reported immense speedup
in certain applications, notably where the hessian of the Gaussian variables has some parts which are highly dense and other
parts which are highly sparse.

- Added spline-based interpolation to `interpolate_marginal_posterior` and all downstream functions. Now the calculation
of marginals doesn't get less stable as more quadrature points are added. Added package `splines` to `Imports`
to support this; since polynomial interpolation almost always gives unstable answers when the number of quadrature
points is even moderate, I consider this a necessary `Import`.

- Added the option `interpolation` to `default_control()` (see documentation).
Default option of `auto` designed to always give stable marginal posterior interpolation
regardless of the number of quadrature points.


## Bug fixes

- `aghq::laplace_approximation()` had a typo and was returning the wrong value. This has been fixed and tests added for its accuracy based on an example with a known answer.

- Fixed an issue in the optimization where the `trustOptim` package was not being checked
for, and this was throwing a cryptic error. Now, it throws a less cryptic error.

- Changed default options for optimization in all functions to `optim(...method = 'BFGS')`, in case
users do not have the `trust` or `trustOptim` packages installed.

- Added `expm1` to `logdiffexp()` to improve numerical stability.

## Other



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
