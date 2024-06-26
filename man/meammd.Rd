% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meammd.R
\name{meammd}
\alias{meammd}
\title{MEA-MMD: Multivariate Efficient Approximate Maximum Mean Discrepancy}
\usage{
meammd(
  X,
  Y,
  beta = -0.1,
  pval = TRUE,
  type = c("proj", "dist"),
  numproj = 20,
  nmethod = c(2, 1),
  distpval = c("Hommel", "Fisher"),
  numperm = 200,
  seednum = 0,
  alternative = c("greater", "two.sided"),
  allowzeropval = FALSE
)
}
\arguments{
\item{X}{Matrix (or vector) of observations in first sample.}

\item{Y}{Matrix (or vector) of observations in second sample.}

\item{beta}{kernel parameter. Must be positive; if not, computes
median heuristic in quadratic time for each projection.
Default value
is \code{-0.1}, which will force median heuristic to be used.}

\item{pval}{Boolean for whether to compute p-value or not.}

\item{type}{The type of projection used. Either \code{"proj"} for
random projections (default) or \code{"dist"} for interpoint
distances.}

\item{numproj}{Number of projections (only used if \code{type="proj"}).
Default is \code{20}.}

\item{nmethod}{Norm used for interpoint distances, if \code{type="dist"}.
Needs to be either \code{2} (for two-norm, default) or
\code{1} (for one-norm).}

\item{distpval}{The p-value combination procedure if \code{type="dist"}.
Options are \code{"Hommel"} (default) or \code{"Fisher"}.
The Hommel method is preferred since the Type I error does
not seem to be controlled if the Fisher method is used.}

\item{numperm}{Number of permutations. Default is \code{200}.}

\item{seednum}{Seed number for generating permutations. Default is \code{0},
which means seed is set randomly. For values larger than
\code{0}, results will be reproducible.}

\item{alternative}{A character string specifying the alternative hypothesis,
which must be either \code{"greater"} (default) or
\code{"two.sided"}. In Gretton et al., the
MMD test statistic is specified so that if it is
significantly larger than zero, then the null hypothesis
that the two samples come from the same distribution
should be rejected. For this reason, \code{"greater"}
is recommended. The test will still work
in many cases with \code{"two.sided"} specified, but this
could lead to problems in certain cases.}

\item{allowzeropval}{A boolean, specifying whether we will allow zero
p-values or not. Default is \code{FALSE}; then
a threshold of \code{0.5 / (numperm+1)} is used,
and if the computed p-value is less than this
threshold, it is then set to be this value.
this avoids the possibility of zero p-values.}
}
\value{
A list with the following elements:
\describe{
\item{\code{pval}}{The p-value of the test, if it is
computed (\code{pval=TRUE}). Otherwise,
it is set to \code{NA}.}
\item{\code{stat}}{The statistic of the test, which
is only returned when \code{type="proj"},
otherwise it is set to \code{NA}.}
}
}
\description{
Computes maximum mean discrepancy statistics with Laplacian
or Gaussian kernel.
Suitable for multivariate data. Naive approach, quadratic in number
of observations.
}
\examples{
X <- matrix(c(1:12), ncol=2, byrow=TRUE)
Y <- matrix(c(13:20), ncol=2, byrow=TRUE)
# using the random projections method
mmdList <- meammd(X=X, Y=Y, pval=TRUE, type="proj", numproj=50)

# using the method were distances are computed to the various points 
mmdList <- meammd(X=X, Y=Y, pval=TRUE, type="dist")


}
\references{
Bodenham, D. A., and Kawahara, Y. (2023)
"euMMD: efficiently computing the MMD two-sample test statistic for
univariate data." Statistics and Computing 33.5 (2023): 110.
}
