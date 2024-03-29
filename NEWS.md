# Rage (development version)

# Rage 1.7.0

* `plot_life_cycle` now includes a logical argument `show_gv` which allows users to print the `graphViz` code constructed by the function to the screen. This is useful for use elsewhere, or for debugging plots.
* The function `entropy_k` is being deprecated. We change the function name to `life_elas` and recommend use of the new functions `entropy_k_age` and `entropy_k_stage` for Leslie and Lefkovitch matrices respectively. See the article by de Vries et al (2023) for the rationale for this change.

# Rage 1.6.0

* Bug fix in `plot_life_cycle` to ensure that nodes are placed in the correct order. Additional argument `node_order` also allows nodes to be placed in any arbitrary order.
* Changes to `life_expect_var` and `life_expect_mean` to ensure that outputs are consistent (as vectors).
* Improved test coverage (still >95%)


# Rage 1.5.1

* Patch to fix test errors on older systems (r-oldrel-macos-arm64). No changes to main package code. 

# Rage 1.5.0

* added `is_leslie_matrix` function to detect a Leslie matrix.
* added `leslie_collapse` function to aggregate a Leslie matrix to a desired dimension. This function is similar to `mpm_collapse` but it will ensure that the returned aggregated matrix is also in the form of a Leslie matrix. Code contributed by Richard A. Hinrichsen (@hinrich62).
* Error in `life_expect_var` (#193) corrected (contribution from Chrissy Hernández, @chrissy3815). Note that `life_expect_mean` and `life_expect_var` now function slightly differently. Be aware!
* performance enhancements to speed up code and make it more robust.
* increased unit test coverage to >95%


# Rage 1.4.0

* `entropy_d` `entropy_k` `shape_rep` and `shape_surv` can now accept matrices directly. In previous versions, these functions required input of `lx` and/or `mx` trajectories, but now the functions can optionally use `mpm_to_...` functions to create these vectors internally (#174).
* added `remove_final` argument to `mpm_to_table` function. This allows users to optionally remove the final row of the life table to prevent the artificial inflation of mortality/hazard caused by the fact that the final age-class is assumed to be closed (and hence all individuals must die).

# Rage 1.3.0

* Fixed an error in the `entropy_d` function, which calculates Demetrius's entropy.

* Added two methods to calculate generation time in `gen_time`: average parent-offspring age difference & expected age at first reproduction (#183).

# Rage 1.2.0

* Improved documentation of `mpm_to_table` and related functions.
* Standardised the `lx_crit`, `conv` and `xmax` argument values across age-from-stage functions.
* Minor edits to documentation for stylistic consistency.
* Removed dependency on the package `popbio`.
* Added vignette with suggested quality control.


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

