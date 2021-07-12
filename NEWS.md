# Rage (development version)

No significant changes yet.

# Rage 1.1.0

* Removed the vignette that made heavy use of `ggtern` package. This package was not available for some builds of R and thus caused problems for CRAN.

# Rage 1.0.0

* Released on CRAN

# Rage 0.2.0

_Released on Github on 25th April 2021_

* Name changes of functions to consistent snake_case.
* Improved documentation with fuller descriptions and executable examples.
* Grouped functions by type in the documentation.
* Improved vignettes by greatly expanding them.
* Added support using stage names in addition to stage number, plus helper function `name_stages` for (re)naming MPM stages.
* Replaced old function `lifeTimeRepEvents` with specific functions for events (`gen_time`, `life_expect`, `longevity`).
* Added functions for manipulating matrices: `mpm_collapse`, `mpm_rearrange`, `mpm_split`, `mpm_standardize`. 
* Added functions to calculate various vital rates `vr`, `vr_mat`, `vr_vec`.
* Added utility functions `repro_stages` and `standard_stages` to identify reproductive stages either logically (TRUE/FALSE) or by the standardised set of reproductive stages (propagule, pre-reproductive, reproductive and post-reproductive).
* Renamed `vitalRates` to `vital_rates`.
* Renamed (and modified) old function `makeLifeTable` to create new function, `mpm_to_table`. 
* Added `lifetable_convert` function to convert between types of life table (hazard, survivorship and survival probability).
* Replaced `matrixElementPerturbation` and `vitalRatePerturbation` functions with enhanced perturbation functions: `perturb_matrix`, `perturb_stochastic`, `perturb_trans`, `perturb_vr`.
* Added utility functions (`utils.R`) to do various tasks like check validity of matrices, calculate mean matrices, calculate matrix inverse.
* Expanded use of unit tests for all functions.
* Updated DESCRIPTION with contributors
* Added build checks via continuous integration on Travis, Appveyor and Github actions (including weekly checks).
* Added machine-readable codemeta-data information (`codemeta.json`)



# Rage 0.1.0

_Released on Github on 14th December 2018_

First (pre) release package. Functions include: `R0`, `dEntropy`, `kEntropy`, `lifeTimeRepEvents`, `longevity`, `makeLifeTable`, `matrixElementPerturbation`, `plotLifeCycle`, `qsdConverge`, `reprodStages`, `standardizedVitalrates`, `vitalRatePerturbation`, `vitalRates`.

