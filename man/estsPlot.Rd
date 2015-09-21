\name{ests.plot}
\alias{ests.plot}

\title{ Plot Parameter Estimates}

\description{ Plots of a histogram of the parameter estimates}

\usage{
ests.plot(part1.list)
}

\arguments{
  \item{part1.list}{ Output from doubleGauss.fit or logistic.fit}

}
\details{ Plots of a histogram of the parameter estimates}

\value{ NULL}

\note{ There are no further notes}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)ests.plot(out.1)
ests.plot(out.1)
}}

\keyword{ htest}