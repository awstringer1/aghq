## Test environments
* local R installation on Mac OS Catalina 10.15.7, Intel processor, R 4.1.1
* ubuntu 16.04 (on travis-ci), (oldrelease, release, devel)
* win-builder (oldrelease, release, devel)

## R CMD check results

There were no ERRORs or WARNINGs on any platform.

There were no NOTEs on the Mac or Ubuntu platforms.

There was one NOTE on the Windows platform:

  Maintainer: 'Alex Stringer <alex.stringer@mail.utoronto.ca>'
  
  New submission
  
  Package was archived on CRAN
  
  Non-FOSS package license (file LICENSE)
  
  Possibly misspelled words in DESCRIPTION:
    AGHQ (8:9)
    aghq (13:9)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2021-10-22 as check problems were not
      corrected in time.
  
  Found the following (possibly) invalid URLs:
    URL: https://cran.rstudio.com/web/packages/grand-total/aghq/index.html
      From: README.md
      Status: 404
      Message: Not Found
      CRAN URL not in canonical form
    URL: https://travis-ci.com/awstringer1/aghq (moved to https://www.travis-ci.com/awstringer1/aghq)
      From: README.md
      Status: 404
      Message: Not Found
    Canonical CRAN.R-project.org URLs use https.

The LICENSE file corresponds to the GNU General Public License.

The misspelled words are acronyms pertaining to the method and package name.

The invalid URLs corrspond to badges on the Github page for this package, and are working as intended.

## Notes

The package was previously removed from CRAN because of a failed check pertaining to the non-conditional
use of packages listed in `Suggests`. This has been fixed. I apologize for not correcting this
before the stated deadline and thank the maintainers for considering this submission.
