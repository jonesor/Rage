## R CMD check results

* This is a new release.
* Includes improvements to functionality and some method corrections.

## Test environments
- R-hub Windows Server 2022, R-devel, 64 bit
- R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
- R-hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results

── Rage 1.4.0: NOTE

  Build ID:   Rage_1.4.0.tar.gz-55334a159687416cbb72b3b4ede9df0c
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  6h 16m 30.2s ago
  Build time: 4m 45.1s

❯ checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── Rage 1.4.0: NOTE

  Build ID:   Rage_1.4.0.tar.gz-0575e5a1c59f4ad2ac57eac4f69449cc
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  6h 16m 30.5s ago
  Build time: 1h 48m 44.6s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Owen Jones <jones@biology.sdu.dk>’
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1111/2041-210X.13289
      From: inst/doc/a03_LifeHistoryTraits.html
      Status: 403
      Message: Forbidden
    URL: https://doi.org/10.1111/j.2041-210X.2010.00087.x
      From: inst/doc/a03_LifeHistoryTraits.html
      Status: 403
      Message: Forbidden
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1111/2041-210X.13792
      From: inst/CITATION
      Status: Forbidden
      Message: 403

0 errors ✔ | 0 warnings ✔ | 1 note ✖

── Rage 1.4.0: NOTE

  Build ID:   Rage_1.4.0.tar.gz-ccb7f7da53a64a8ca7ed52eec418f69e
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  6h 16m 30.6s ago
  Build time: 1h 41m 59.6s

❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable

0 errors ✔ | 0 warnings ✔ | 1 note ✖

* This is a new release.
