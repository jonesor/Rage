% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mature_.R
\name{repro_maturity}
\alias{repro_maturity}
\alias{mature_prob}
\alias{mature_age}
\alias{mature_distrib}
\title{Age of reproductive maturity}
\usage{
mature_prob(matU, matR, start = 1L)

mature_age(matU, matR, start = 1L)

mature_distrib(matU, start = 1L, repro_stages)
}
\arguments{
\item{matU}{The survival component of a matrix population model (i.e., a
square projection matrix reflecting survival-related transitions; e.g.
progression, stasis, and retrogression). Optionally with named rows and
columns indicating the corresponding life stage names.}

\item{matR}{The reproductive component of a matrix population model (i.e., a
square projection matrix reflecting transitions due to reproduction; either
sexual, clonal, or both). Optionally with named rows and columns indicating
the corresponding life stage names.}

\item{start}{The index (or stage name) of the first stage at which the author
considers the beginning of life. Defaults to \code{1}.}

\item{repro_stages}{A vector of stage names or indices indicating which
stages are reproductive. Alternatively, a logical vector of length
\code{ncol(matU)} indicating whether each stage is reproductive
(\code{TRUE}) or not (\code{FALSE}).}
}
\value{
For \code{mature_distrib}, a vector giving the proportion of
  individuals that first reproduce within each stage class. For all others, a
  scalar trait value.
}
\description{
Apply Markov chain approaches to compute age-specific
trajectory of reproduction for individuals in a matrix population model.
Includes functions to calculate the probability of achieving reproductive
maturity (\code{mature_prob}), mean age at first reproduction
(\code{mature_age}), and distribution of individuals first achieving
reproductive maturity among stage class (\code{mature_distrib}).
}
\note{
Note that the units of time in returned values are the same as the
  \code{ProjectionInterval} of the MPM.
}
\examples{
data(mpm1)

mature_prob(mpm1$matU, mpm1$matF, start = 2)
mature_age(mpm1$matU, mpm1$matF, start = 2)

### distribution of first reproductive maturity among stage classes
repstage <- repro_stages(mpm1$matF)
mature_distrib(mpm1$matU, start = 2, repro_stages = repstage)

}
\references{
Caswell, H. 2001. Matrix Population Models: Construction,
  Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
  978-0878930968
}
\seealso{
Other life history traits: 
\code{\link{entropy_d}()},
\code{\link{entropy_k_age}()},
\code{\link{entropy_k_stage}()},
\code{\link{entropy_k}()},
\code{\link{gen_time}()},
\code{\link{life_elas}()},
\code{\link{life_expect_mean}()},
\code{\link{longevity}()},
\code{\link{net_repro_rate}()},
\code{\link{shape_rep}()},
\code{\link{shape_surv}()}
}
\author{
Roberto Salguero-Gomez <rob.salguero@zoo.ox.ac.uk>

Hal Caswell <hcaswell@whoi.edu>

Owen R. Jones <jones@biology.sdu.dk>

Patrick Barks <patrick.barks@gmail.com>
}
\concept{life history traits}
