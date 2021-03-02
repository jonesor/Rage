# Rage (development version)

* Name changes of functions to consistent snake_case.
* Improved documentation with fuller descriptions and executable examples.
* Grouped functions by type in the documentation.
* Improved vignettes by greatly expanding them.
* Replaced old function `lifeTimeRepEvents` with specific functions for events (`gen_time`, `life_expect`, `longevity`).
* Added functions for manipulating matrices: `mpm_collapse`, `mpm_rearrange`, `mpm_split`, `mpm_standardize`. 
* Added functions to calculate various vital rates `vr`, `vr_mat`, `vr_vec`.
* Added utility functions `repro_stages` and `standard_stages` to identify reproductive stages either logically (TRUE/FALSE) or by the standardised set of reproductive stages (propagule, pre-reproductive, reproductive and post-reproductive).
* Renamed (and modified) old function `makeLifeTable` to create new function, `mpm_to_table`. 
* Added `lifetable_convert` function to convert between types of life table (hazard, survivorship and survival probability).
* Replaced `matrixElementPerturbation` and `vitalRatePerturbation` functions with enhanced perturbation functions: `perturb_matrix`, `perturb_stochastic`, `perturb_trans`, `perturb_vr`.
* Added utility functions (`utils.R`) to do various tasks like check validity of matrices, calculate mean matrices, calculate matrix inverse.
* Expanded use of unit tests for all functions.
* Updated DESCRIPTION with contributors
* Added build checks via continuous integration on Travis, Appveyor and Github actions (including weekly checks).


# Rage 0.1.0

First (pre) release package. Functions include: `R0`, `dEntropy`, `kEntropy`, `lifeTimeRepEvents`, `longevity`, `makeLifeTable`, `matrixElementPerturbation`, `plotLifeCycle`, `qsdConverge`, `reprodStages`, `standardizedVitalrates`, `vitalRatePerturbation`, `vitalRates`.

