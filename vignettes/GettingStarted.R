## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(digits = 4)
#ragg_png = function(..., res = 192) {
#  ragg::agg_png(..., res = res, units = "in")
#}
#knitr::opts_chunk$set(dev = "ragg_png", fig.ext = "png")
# knitr::opts_chunk$set(dev = "png", dev.args = list(type = "cairo-png"))

## -----------------------------------------------------------------------------
library(Rage)  # load Rage
data(mpm1)     # load data object 'mpm1'
mpm1           # display the contents

## -----------------------------------------------------------------------------
life_expect(matU = mpm1$matU, start = 1)  # life expectancy from "seed" stage
life_expect(matU = mpm1$matU, start = 2)  # life expectancy from "small" stage
longevity(matU = mpm1$matU, start = 2, lx_crit = 0.05)  # post-germination years until survivorship falls below 5%

## -----------------------------------------------------------------------------
mature_age(matU = mpm1$matU, matR = mpm1$matF, start = 2)     # post-germination years to first reproduction
mature_prob(matU = mpm1$matU, matR = mpm1$matF, start = 2)    # post-germination Pr(survival to first repro)
net_repro_rate(matU = mpm1$matU, matR = mpm1$matF)            # net reproductive rate (aggregate)
gen_time(matU = mpm1$matU, matR = mpm1$matF)                  # generation time (aggregate)

## -----------------------------------------------------------------------------
lx <- mpm_to_lx(matU = mpm1$matU, start = 2)
px <- mpm_to_px(matU = mpm1$matU, start = 2)
hx <- mpm_to_hx(matU = mpm1$matU, start = 2)
mx <- mpm_to_mx(matU = mpm1$matU, matR = mpm1$matF, start = 2)

# then calculate life history traits
entropy_k(lx)       # Keyfitz' entropy
entropy_d(lx, mx)   # Demetrius' entropy
shape_surv(lx)      # shape of survival/mortality trajectory
shape_rep(mx)       # shape of fecundity trajectory

## ---- warning=FALSE, message=FALSE--------------------------------------------
# derive life table from MPM
lt <- mpm_to_table(mpm1$matU, start = 2)

# calculate time to QSD
(q <- qsd_converge(mpm1$matU, start = 2))

# plot mortality trajectory w/ vertical line at time to QSD
par(mar = c(4.5, 4.5, 1, 1))
plot(qx ~ x, data = lt, type = "l", ylim = c(0, 0.65))
abline(v = q, lty = 2)

## -----------------------------------------------------------------------------
# calculate the shape of the survival/mortality trajectory
shape_surv(lt$lx)       # based on full lx trajectory
shape_surv(lt$lx[1:q])  # based on lx trajectory prior to the QSD

## -----------------------------------------------------------------------------
vr_vec_survival(mpm1$matU)
vr_vec_growth(mpm1$matU, exclude = c(1, 5))
vr_vec_shrinkage(mpm1$matU, exclude = 5)
vr_vec_stasis(mpm1$matU)
vr_vec_dorm_enter(mpm1$matU, dorm_stages = 5)
vr_vec_dorm_exit(mpm1$matU, dorm_stages = 5)
vr_vec_reproduction(mpm1$matU, mpm1$matF)

## -----------------------------------------------------------------------------
# derive full MPM (matA)
mpm1$matA <- mpm1$matU + mpm1$matF

# calculate stable stage distribution at equilibrium using popbio::stable.stage
library(popbio)
w <- popbio::stable.stage(mpm1$matA)

# calculate MPM-specific vital rates
vr_survival(mpm1$matU, exclude_col = c(1, 5), weights_col = w)
vr_growth(mpm1$matU, exclude = c(1, 5), weights_col = w)
vr_shrinkage(mpm1$matU, exclude = c(1, 5), weights_col = w)
vr_stasis(mpm1$matU, exclude = c(1, 5), weights_col = w)
vr_dorm_enter(mpm1$matU, dorm_stages = 5, weights_col = w)
vr_dorm_exit(mpm1$matU, dorm_stages = 5, weights_col = w)
vr_fecundity(mpm1$matU, mpm1$matF, weights_col = w)

## -----------------------------------------------------------------------------
# matrix element perturbation
perturb_matrix(mpm1$matA, type = "sensitivity")

# vital rate perturbation
# (we use as.data.frame here for prettier printing)
as.data.frame(perturb_vr(mpm1$matU, mpm1$matF, type = "sensitivity"))

# transition type perturbation
as.data.frame(perturb_trans(mpm1$matU, mpm1$matF, type = "sensitivity"))

## -----------------------------------------------------------------------------
# collapse 'small', 'medium', and 'large' stages into single stage class
col1 <- mpm_collapse(mpm1$matU, mpm1$matF, collapse = list(1, 2:4, 5))
col1$matA

## -----------------------------------------------------------------------------
# compare population growth rate of original and collapsed MPM (preserved)
popbio::lambda(mpm1$matA)
popbio::lambda(col1$matA)

# compare net reproductive rate of original and collapsed MPM (not preserved)
net_repro_rate(mpm1$matU, mpm1$matF)
net_repro_rate(col1$matU, col1$matF)

