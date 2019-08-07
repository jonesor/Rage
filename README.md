
<!-- README.md is generated from README.Rmd. Please edit that file -->
Rage
====

[![Travis-CI Build Status](https://travis-ci.org/jonesor/Rage.svg?branch=devel)](https://travis-ci.org/jonesor/Rage) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jonesor/Rage?branch=devel&svg=true)](https://ci.appveyor.com/project/jonesor/Rage) [![Coverage status](https://codecov.io/gh/jonesor/Rage/branch/devel/graph/badge.svg)](https://codecov.io/github/jonesor/Rage?branch=devel)

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
life_expect(matU)
#> [1] 5.749
```

Calculate net reproductive rate:

``` r
net_repro_rate(matU, matF)
#> [1] 2.464
```

Produce a life table:

``` r
mpm_to_table(matU, matF, xmax = 15)
#>     x      lx      dx     qx     Lx     Tx     ex     mx    lxmx
#> 1   0 1.00000 0.30000 0.3000 0.8500 4.3932 4.3932 0.0000 0.00000
#> 2   1 0.70000 0.21000 0.3000 0.5950 3.5432 5.0617 0.1837 0.12857
#> 3   2 0.49000 0.08700 0.1776 0.4465 2.9482 6.0167 0.3573 0.17509
#> 4   3 0.40300 0.05490 0.1362 0.3755 2.5017 6.2076 0.4731 0.19068
#> 5   4 0.34810 0.04263 0.1225 0.3268 2.1261 6.1078 0.5457 0.18994
#> 6   5 0.30547 0.03542 0.1160 0.2878 1.7993 5.8903 0.5903 0.18031
#> 7   6 0.27005 0.03031 0.1122 0.2549 1.5116 5.5974 0.6174 0.16674
#> 8   7 0.23974 0.02637 0.1100 0.2266 1.2567 5.2418 0.6339 0.15197
#> 9   8 0.21337 0.02318 0.1087 0.2018 1.0301 4.8279 0.6438 0.13737
#> 10  9 0.19018 0.02051 0.1078 0.1799 0.8283 4.3554 0.6498 0.12358
#> 11 10 0.16967 0.01822 0.1074 0.1606 0.6484 3.8215 0.6534 0.11086
#> 12 11 0.15146 0.01622 0.1071 0.1434 0.4878 3.2209 0.6555 0.09929
#> 13 12 0.13524 0.01446 0.1069 0.1280 0.3445 2.5472 0.6568 0.08883
#> 14 13 0.12079 0.01290 0.1068 0.1143 0.2165 1.7922 0.6576 0.07943
#> 15 14 0.10789 0.01151 0.1067 0.1021 0.1021 0.9466 0.6581 0.07100
#> 16 15 0.09638      NA     NA     NA     NA     NA 0.6583 0.06345
```

Collapse the matrix population model to two stages, pre-reproductive and reproductive, and conduct perturbation analyses on the vital rates of the collapsed matrix:

``` r
# pre-reproductive stages are 1 and 2, reproductive stages are 3 and 4
collapsed_mpm <- mpm_collapse(matU, matF, collapse = list(1:2, 3:4))

# perturbation analysis (summed elasticities for various vital-rates)
perturb_vr(collapsed_mpm$matU, collapsed_mpm$matF, type = "elasticity")
#> $survival
#> [1] 1
#> 
#> $growth
#> [1] 0.09701
#> 
#> $shrinkage
#> [1] -0.2368
#> 
#> $fecundity
#> [1] 0.1643
#> 
#> $clonality
#> [1] 0
```
