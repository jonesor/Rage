## R CMD check results

* This is a minor update. Adding two new functions and deprecating an old one. Minor bug fixes.

* NOTES include: 
1) Maintainer Owen Jones <jones@biology.sdu.dk>
2) Possibly invalid URLs - these have been verified as correct.
3) package 'V8' unavailable
4) Found the following files/directories:''NULL'', 'lastMiKTeXException'

I believe these are not problematic notes.

## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [28s] NOTE
  Maintainer: 'Owen Jones <jones@biology.sdu.dk>'
  
  Found the following (possibly) invalid URLs:
    URL: http://rich-iannone.github.io/DiagrammeR/ (moved to https://rich-iannone.github.io/DiagrammeR/)
      From: man/plot_life_cycle.Rd
      Status: 200
      Message: OK
    URL: https://jonesor.github.io/Rage/articles/a03_AgeFromStage.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a04_TernaryPlots.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a05_LifeHistoryTraits.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found

❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

❯ On windows-x86_64-devel (r-devel)
  checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

❯ On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... [10s/63s] NOTE
  Maintainer: ‘Owen Jones <jones@biology.sdu.dk>’
  
  Found the following (possibly) invalid URLs:
    URL: http://rich-iannone.github.io/DiagrammeR/ (moved to https://rich-iannone.github.io/DiagrammeR/)
      From: man/plot_life_cycle.Rd
      Status: 200
      Message: OK
    URL: https://jonesor.github.io/Rage/articles/a03_AgeFromStage.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a04_TernaryPlots.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a05_LifeHistoryTraits.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found

❯ On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [12s/80s] NOTE
  Maintainer: ‘Owen Jones <jones@biology.sdu.dk>’
  
  Found the following (possibly) invalid URLs:
    URL: http://rich-iannone.github.io/DiagrammeR/ (moved to https://rich-iannone.github.io/DiagrammeR/)
      From: man/plot_life_cycle.Rd
      Status: 200
      Message: OK
    URL: https://jonesor.github.io/Rage/articles/a03_AgeFromStage.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a04_TernaryPlots.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found
    URL: https://jonesor.github.io/Rage/articles/a05_LifeHistoryTraits.html
      From: inst/doc/a01_GettingStarted.html
      Status: 404
      Message: Not Found

0 errors ✔ | 0 warnings ✔ | 6 notes ✖
