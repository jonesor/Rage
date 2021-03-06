% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/life_expect.R
\name{life_expect_mean}
\alias{life_expect_mean}
\alias{life_expect_var}
\alias{life_expect}
\title{Calculate mean and variance of life expectancy from a matrix population model}
\usage{
life_expect_mean(matU, start = 1L)

life_expect_var(matU, start = 1L)

life_expect(matU, start = 1L)
}
\arguments{
\item{matU}{The survival component of a MPM (i.e. a square projection matrix
reflecting survival-related transitions; e.g. progression, stasis, and
retrogression). Optionally with named rows and columns indicating the
corresponding life stage names.}

\item{start}{The index (or stage name) of the first stage of the life cycle
which the user considers to be the beginning of life. Defaults to 1.
Alternately, a numeric vector giving the starting population vector (in which
 case \code{length(start)} must match \code{ncol(matU))}. See section
\emph{Starting from multiple stages}.}
}
\value{
Returns life expectancy. If \code{matU} is singular (often indicating
  infinite life expectancy), returns \code{NA}.
}
\description{
Applies Markov chain approaches to obtain mean and variance of life
expectancy from a matrix population model (MPM).
}
\section{Starting from multiple stages}{

Rather than specifying argument \code{start} as a single stage class from
which all individuals start life, it may sometimes be desirable to allow for
multiple starting stage classes. For example, if the user wants to start their
calculation of life expectancy from reproductive maturity (i.e. first
reproduction), they should account for the possibility that there may be
multiple stage classes in which an individual could first reproduce.

To specify multiple starting stage classes, specify argument \code{start} as
the desired starting population vector (\strong{n1}), giving the proportion
of individuals starting in each stage class (the length of \code{start}
should match the number of columns in the relevant MPM).

See function \code{\link{mature_distrib}} for calculating the proportion of
individuals achieving reproductive maturity in each stage class.
}

\examples{
data(mpm1)

# mean life expectancy starting from stage class 2 
life_expect_mean(mpm1$matU, start = 2)
life_expect_mean(mpm1$matU, start = "small")  # equivalent using named life stages

# mean life expectancy starting from first reproduction
rep_stages <- repro_stages(mpm1$matF)
n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = rep_stages)
life_expect_mean(mpm1$matU, start = n1)

# variance of life expectancy from stage class 2
life_expect_var(mpm1$matU, start = 1)


}
\references{
Caswell, H. 2001. Matrix Population Models: Construction,
  Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
  978-0878930968
}
\seealso{
Other life history traits: 
\code{\link{entropy_d}()},
\code{\link{entropy_k}()},
\code{\link{gen_time}()},
\code{\link{longevity}()},
\code{\link{net_repro_rate}()},
\code{\link{repro_maturity}},
\code{\link{shape_rep}()},
\code{\link{shape_surv}()}
}
\author{
Roberto Salguero-Gomez <rob.salguero@zoo.ox.ac.uk>

Hal Caswell <hcaswell@whoi.edu>
}
\concept{life history traits}
