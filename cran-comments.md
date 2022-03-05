## Test environments
* local R installation on Mac OS Catalina 12.2.1, arm64 M1 processor, R 4.1.1
* ubuntu 16.04 (on travis-ci), (oldrelease, release, devel)
* win-builder (oldrelease, release, devel)

## R CMD check results

There were no ERRORs or WARNINGs or NOTEs on any platform when all suggested packages were installed.

There was one NOTE when 'trustOptim' and 'trust' were not installed:

> checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'trustOptim', 'trust'
  
I have added checks for these packages' installation inside `optimize_theta`, where
they are used, and automatically switch the user-requested method to use `stats::optim`,
with a warning, if they are not installed. I have added code to skip unit tests involving
these packages when they are not installed. Both of these packages are available from
`CRAN` and should not be a problem for users to obtain.

## Notes
