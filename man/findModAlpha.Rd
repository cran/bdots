\name{find.mod.alpha}
\alias{find.mod.alpha}

\title{ Calculate a corrected alpha}

\description{ Calculates an adjusted alpha to test against for timeseries datasets}

\usage{
find.mod.alpha(rho, n, df, alpha = .05, error.acc = .001,
	grad.diff = ifelse(cores > 3, .5, .1),
	verbose = FALSE, cores = 2, method = "t")
}

\arguments{
  \item{rho}{ Output from doubleGauss.fit or logistic.fit}
  \item{n}{ Total number of tests being performed (i.e. length of the time series)}
  \item{df}{ Degrees of freedom used for the individual t-tests.}
  \item{alpha}{ Familywise error rate.}
  \item{error.acc}{ Acceptable distance between estimated error rate and intended familywise error rate}
  \item{grad.diff}{ Distance between estimates for calculation of Newton's and Hailey's methods}
  \item{verbose}{ If TRUE, prints adjusted alpha estimate for each iteration.}
  \item{cores}{ Number of cores to use. More than 2 cores is not currently supported and may provide non-sensical results.}
  \item{method}{ Multivariate calculation approach. Set to "NormApprox" for normal approximation instead of t.}
}
\details{ Calculates an adjusted alpha to test against for timeseries datasets. This is intended for use with t-statistics. An AR1 autocorrelation structure is assumed.}

\value{ Numeric.}

\note{ There are no further notes}

\examples{
\dontrun{
# Sample Data
x <- numeric(500)
x[1] <- rnorm(1)
for(i in 2:500) x[i] <- rnorm(1, x[i - 1] * .9, 1 - .9 ^ 2)

# Estimate autocorrelation
# Assume the x values are t-statistics based on 50 subjects
rho.est <- ar(x, FALSE, order.max = 1)$ar
alphastar <- find.mod.alpha(rho.est, alpha = .05, n = 500, df = 49)
alphastar
}}

\keyword{ ts}