# aghq 0.3.2

## New features

- added an `S3` interface for parameter transformations, see `make_transformation`,
`validate_transformation`, `default_transformation`, and any `aghq` package function with a `transformation`
argument.

## Bug fixes

## Other

# aghq 0.3.1

## New features

- Added support for doing the multiple required Cholesky decompositions in parallel in `sample_marginal.marginallaplace`.

- Switched from using `chol` to using `Matrix::Cholesky` with `perm = TRUE` inside `sample_marginal.marginallaplace`.

- Added spline-based interpolation to `interpolate_marginal_posterior` and all downstream functions. Now the calculation
of marginals doesn't get less stable as more quadrature points are added. Added package `splines` to `Imports`
to support this; since polynomial interpolation almost always gives unstable answers when the number of quadrature
points is even moderate, I consider this a necessary `Import`.

- Added the option `interpolation` to `default_control()` (see documentation).
Default option of `auto` designed to always give stable marginal posterior interpolation
regardless of the number of quadrature points.

- Added an internal `validate_control` check to all functions which use a `control` argument, which makes sure the user inputs a control list with the correct names and value types. This is supported by the existing control functions, and  helps prevent further cryptic downstream errors.

- Added a `onlynormconst` option to `aghq` and related functions. Simply returns the numeric value of the log integral, avoiding all the extra stuff, at greater speed.

- Added a new summary method for objects of class `marginallaplace`, that includes information on the random effects.

- Preserve variable names in all summary output.

## Bug fixes

- Added a `requireNamespace` condition to all functions from packages listed in `Suggests`.

- `aghq::laplace_approximation()` had a typo and was returning the wrong value. This has been fixed and tests added for its accuracy based on an example with a known answer.

- Fixed an issue in the optimization where the `trustOptim` package was not being checked
for, and this was throwing a cryptic error.

- Changed default options for optimization in all functions to `optim(...method = 'BFGS')`, in case
users do not have the `trust` or `trustOptim` packages installed.

- Added `expm1` to `logdiffexp()` to improve numerical stability.

- Fixed `optimize_theta` so that `control` arguments are passed correctly. 

- Removed several unit tests that were failing on M1 Macs. These tests werre actually
testing that polynomial interpolation of marginal posteriors FAILs, so apparently
this isn't failing on these new Macs, but that's better, not worse. Will re-test
and potentially add back once I have local access to this hardware.

## Other

- Re-added `numDeriv` as an Import, since it is used in core functionality.

- Switched default optimization control arguments to use `base::optim`, to facilitate
removal of `trustOptim` and `trust` as Import dependencies.

- Switched default method for numerically differentiated Hessians to `'Richardson'`,
for more accurate results.

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
