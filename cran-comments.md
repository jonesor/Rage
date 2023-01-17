## R CMD check results

## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel), fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

❯ On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Owen Jones <jones@biology.sdu.dk>’
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1111/2041-210X.13289
      From: inst/doc/a03_LifeHistoryTraits.html
      Status: 503
      Message: Service Unavailable
    URL: https://doi.org/10.1111/j.2041-210X.2010.00087.x
      From: inst/doc/a03_LifeHistoryTraits.html
      Status: 503
      Message: Service Unavailable

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

* This is a new release.
