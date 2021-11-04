## Test environments
* local R installation on Mac OS Catalina 10.15.7, Intel processor, R 4.1.1
* ubuntu 16.04 (on travis-ci), (oldrelease, release, devel)
* win-builder (oldrelease, release, devel)

## R CMD check results

There were no ERRORs or WARNINGs on any platform.

  There was one NOTE:
  
  Maintainer: ‘Alex Stringer <alex.stringer@uwaterloo.ca>’
  
  New submission
  
  Package was archived on CRAN
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2021-10-22 as check problems were not
      corrected in time.

## Notes

I have removed the offending URLs from the previous submission, and changed the license to the standard GPL-3 thereby removing that NOTE.

The package was previously removed from CRAN because of a failed check pertaining to the non-conditional
use of packages listed in `Suggests`. This has been fixed. I apologize for not correcting this
before the stated deadline and thank the maintainers for considering this submission.
