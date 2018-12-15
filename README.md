
<!-- README.md is generated from README.Rmd. Please edit that file -->
Rage
====

[![Travis-CI Build Status](https://travis-ci.org/jonesor/Rage.svg?branch=master)](https://travis-ci.org/jonesor/Rage) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jonesor/Rage?branch=master&svg=true)](https://ci.appveyor.com/project/jonesor/Rage)

An R package for matrix population model analysis. Note this package is at an early stage of development, and may contain bugs.

Installation
------------

Install from GitHub with:

``` r
devtools::install_github("jonesor/Rage")
#
# or
#
install.packages("remotes") # smaller and quicker to install than devtools
remotes::install_github("jonesor/Rage")
```

Usage
-----

``` r
library(Rage)
```

First we need a matrix population model separated into the growth/survival component and sexual reproduction component (possibly also a clonal component, but we'll skip that here).

``` r
matU <- rbind(c(0.1,   0,   0,   0), # growth/survival component
              c(0.6, 0.2, 0.1,   0),
              c(  0, 0.5, 0.5, 0.1),
              c(  0,   0, 0.3, 0.8))

matF <- rbind(c(  0,   0, 0.2, 0.6), # sexual reproduction component
              c(  0,   0, 0.1, 0.2),
              c(  0,   0,   0,   0),
              c(  0,   0,   0,   0))
```

Now we can get to work. Here are just a few example analyses. <br>

Calculate life expectancy:

``` r
lifeExpectancy(matU)
#> [1] 5.749
```

Calculate net reproductive rate:

``` r
R0(matU, matF)
#> [1] 2.464
```

Produce a life table:

``` r
makeLifeTable(matU, matF, nSteps = 15)
#>     x     lx      dx     qx     Lx     Tx     ex     mx    lxmx
#> 1   0 1.0000 0.30000 0.3000 0.8500 4.2910 4.2910 0.0000 0.00000
#> 2   1 0.7000 0.21000 0.3000 0.5950 3.4410 4.9158 0.0000 0.00000
#> 3   2 0.4900 0.08700 0.1776 0.4465 2.8460 5.8082 0.1837 0.09000
#> 4   3 0.4030 0.05490 0.1362 0.3755 2.3995 5.9542 0.3573 0.14400
#> 5   4 0.3481 0.04263 0.1225 0.3268 2.0240 5.8144 0.4731 0.16470
#> 6   5 0.3055 0.03542 0.1160 0.2878 1.6972 5.5560 0.5457 0.16668
#> 7   6 0.2700 0.03031 0.1122 0.2549 1.4094 5.2192 0.5903 0.15940
#> 8   7 0.2397 0.02637 0.1100 0.2266 1.1545 4.8158 0.6174 0.14802
#> 9   8 0.2134 0.02318 0.1087 0.2018 0.9280 4.3492 0.6339 0.13525
#> 10  9 0.1902 0.02051 0.1078 0.1799 0.7262 3.8184 0.6438 0.12244
#> 11 10 0.1697 0.01822 0.1074 0.1606 0.5463 3.2195 0.6498 0.11025
#> 12 11 0.1515 0.01622 0.1071 0.1434 0.3857 2.5466 0.6534 0.09896
#> 13 12 0.1352 0.01446 0.1069 0.1280 0.2424 1.7920 0.6555 0.08866
#> 14 13 0.1208 0.01290 0.1068 0.1143 0.1143 0.9466 0.6568 0.07934
#> 15 14 0.1079      NA     NA     NA     NA     NA 0.6576 0.07095
```

Collapse the matrix population model to two stages, pre-reproductive and reproductive, and conduct perturbation analyses on the vital rates of the collapsed matrix:

``` r
# pre-reproductive stages are 1 and 2, reproductive stages are 3 and 4
mpm_collapse <- collapseMatrix(matU, matF, collapse = list(1:2, 3:4))

# perturbation analysis (S- represent sensitivities, and E- elasticities)
vitalRatePerturbation(mpm_collapse$matU, mpm_collapse$matF)
#>   SSurvival SGrowth SShrinkage SReproduction SClonality ESurvival EGrowth
#> 1     1.312  0.2707    -0.3888        0.2698          0         1 0.09701
#>   EShrinkage EReproduction EClonality
#> 1    -0.2368        0.1643          0
```
