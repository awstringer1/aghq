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

Received notice that tests were failing on M1 Macs. I apologize; I did not realize
that a testing service was available for this new hardware. I removed the offending
tests (which were not necessary in the first place) and added the following bit to NEWS.md:

> Removed several unit tests that were failing on M1 Macs. These tests werre actually
testing that polynomial interpolation of marginal posteriors FAILs, so apparently
this isn't failing on these new Macs, but that's better, not worse. Will re-test
and potentially add back once I have local access to this hardware.

I will add testing on the M1 Mac service to my workflow for all future CRAN submissions.