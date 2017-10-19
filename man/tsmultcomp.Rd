\name{tsmultcomp}
\alias{tsmultcomp}

\title{ Calculate a corrected alpha}

\description{ Calculates an adjusted alpha to test against for timeseries datasets}

\usage{
tsmultcomp(rho, alpha, N, df, maxItr = 1000, tol = 1e-8, verbose = FALSE)
}

\arguments{
  \item{rho}{ Output from doubleGauss.fit or logistic.fit}
  \item{alpha}{ Familywise error rate.}
  \item{N}{ Total number of tests being performed (i.e. length of the time series)}
  \item{df}{ Degrees of freedom used for the individual t-tests.}
  \item{maxItr}{ Maximum number of iterations to reach convergence for adjusted alpha estimate.}
  \item{tol}{ Tolerance level for assessing convergence of the adjusted alpha.}
  \item{verbose}{ If TRUE, prints adjusted alpha estimate for each iteration.}
}
\details{ Calculates an adjusted alpha to test against for timeseries datasets. This is intended for use with t-statistics. An AR1 autocorrelation structure is assumed.}

\value{ Numeric.}

\note{ There are no further notes}

\examples{
# Sample Data
x <- numeric(500)
x[1] <- rnorm(1)
for(i in 2:500) x[i] <- rnorm(1, x[i - 1] * .9, 1 - .9 ^ 2)

# Estimate autocorrelation
# Assume the x values are t-statistics based on 50 subjects
rho.est <- ar(x, FALSE, order.max = 1)$ar
alphastar <- tsmultcomp(rho.est, alpha = .05, N = 500, df = 49)
alphastar
}

\keyword{ ts}