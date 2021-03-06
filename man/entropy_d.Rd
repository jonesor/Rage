% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_d.R
\name{entropy_d}
\alias{entropy_d}
\title{Calculate Demetrius' entropy from trajectories of age-specific survivorship
and fecundity}
\usage{
entropy_d(lx, mx)
}
\arguments{
\item{lx}{Age-specific survivorship trajectory (a vector of
monotonically-declining values in the interval [0,1]).}

\item{mx}{Age-specific fecundity trajectory (a vector of non-negative values).}
}
\value{
Demetrius' entropy.
}
\description{
This function calculates Demetrius' entropy from vectors of age-specific
survivorship (lx) and fecundity (mx).
}
\details{
#' @section Warning:
Note that this function may produce unexpected results if used on partial
survivorship and fecundity trajectories. In addition, it is sensitive to the
length of the these vectors. We direct users to the functions
`\code{\link{shape_surv}}` and `\code{\link{shape_rep}}` which are relatively
robust to these issues.
}
\examples{
data(mpm1)

# derive trajectories of lx and mx, starting from stage 2
lx <- mpm_to_lx(mpm1$matU, start = 2)
mx <- mpm_to_mx(mpm1$matU, mpm1$matF, start = 2)

# calculate Demetrius' entropy
entropy_d(lx, mx)

}
\references{
Demetrius, L., & Gundlach, V. M. 2014. Directionality theory and
  the entropic principle of natural selection. Entropy 16: 5428-5522.
}
\seealso{
Other life history traits: 
\code{\link{entropy_k}()},
\code{\link{gen_time}()},
\code{\link{life_expect_mean}()},
\code{\link{longevity}()},
\code{\link{net_repro_rate}()},
\code{\link{repro_maturity}},
\code{\link{shape_rep}()},
\code{\link{shape_surv}()}
}
\author{
Roberto Salguero-Gomez <rob.salguero@zoo.ox.ac.uk>

Patrick Barks <patrick.barks@gmail.com>
}
\concept{life history traits}
