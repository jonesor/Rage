## Resubmission

This is a resubmission. 
The previous version was Archived on 2021-07-04 as CRAN check issues were not corrected in time. The issue was that a package from `Suggests` was not used conditionally in a vignette.

In this version I have made the necessary corrections. Specifically I have:

* Removed the vignette that required `ggtern`. I have also removed `ggtern` from `Suggests` in the `DESCRIPTION`

## Test environments

### local

- R version 4.0.5 (2021-03-31)
- Platform: x86_64-apple-darwin17.0 (64-bit)
- Running under: macOS Catalina 10.15.7

### GitHub Actions

#### Microsoft Windows

__R latest release__

- R version 4.1.0 (2021-05-18)
- Platform: x86_64-w64-mingw32/x64 (64-bit)
- Running under: Windows Server x64 (build 17763)

#### Apple OSX

__R latest release__

- R version 4.1.0 (2021-05-18)
- Platform: x86_64-apple-darwin17.0 (64-bit)
- Running under: macOS Catalina 10.15.7

__R development version__

- R Under development (unstable) (2021-07-07 r80605)
- Platform: x86_64-apple-darwin17.0 (64-bit)
- Running under: macOS Catalina 10.15.7

#### Linux (Ubuntu)

__R latest release__

- R version 4.1.0 (2021-05-18)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Ubuntu 20.04.2 LTS

__R development version__

- R Under development (unstable) (2021-07-07 r80605)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Ubuntu 20.04.2 LTS


## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

Maintainer: ‘Owen Jones <jones@biology.sdu.dk>’
New submission

## Downstream dependencies

There are no current downstream dependencies.

