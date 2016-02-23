\name{bdots.write.csv}
\alias{bdots.write.csv}

\title{ Write to CSV}

\description{ Write bootstrapped estimates and confidence intervals to csv}

\usage{
bdots.write.csv(part2.list, file, ...)
}

\arguments{
  \item{part2.list}{list. Output from doubleGauss.boot}
  \item{file}{Name of file to write to}
  \item{...}{Further arguments to write.csv}
}
\details{ rite bootstrapped estimates and confidence intervals to csv}

\value{ NULL }

\note{ There are no further notes}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)
out.2 <- logistic.boot(out.1)
bdots.write.csv(out.2, "CIOutput.csv", row.names = FALSE)
}}

\keyword{ htest}