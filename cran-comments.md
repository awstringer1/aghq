## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), (oldrelease, release, devel)
* win-builder (oldrelease, release, devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

There were no ERRORs or WARNINGs on any platform.

1 NOTE when checking on Windows:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alex Stringer <alex.stringer@mail.utoronto.ca>'

Re-submission, addressing helpful comments from CRAN

> Possibly mis-spelled words in DESCRIPTION:
  AGHQ (8:33)

This word is not misspelled (it's an acronym) and I don't think there is anything
I can do on my end about the other parts of the NOTE.

> Please do not start the description with "This package", package name,
title or similar.

> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Thank you for providing guidance on how to properly format the DESCRIPTION
file. I have implemented these changes.

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      pipe.Rd:  \value
      
Thank you for catching this omission. This came from my including of the "pipe"
operator from the `magrittr` package, and I forgot to include the `\value` tag
in the documentation. I fixed it.

> Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.: 04-aghq.R

Thank you for this advice, I have removed the call to `par()` from the plot 
method.

> Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.: inst/doc/aghq.R, compute_pdf_and_cdf.Rd
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

> old <- options(digits = 3)
...
options(old)

Thank you for this advice, I have removed the call to `par()` from the example
in the `compute_pdf_and_cdf` documentation.

Thank you!
