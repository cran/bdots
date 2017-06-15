\name{Print Fits}
\alias{printFits}

\title{ Print summary of quality of curve fits}

\description{ Print summary of quality of curve fits}

\usage{
printFits(part1.list)
}

\arguments{
  \item{part1.list}{list. Output from doubleGauss.fit}
}
\details{ Print summary of quality of curve fits}

\note{ There are no further notes}

\examples{
\dontrun{
data(ci)
ci.2 <- subset(ci, ci$LookType == "Cohort" | ci$LookType == "Unrelated")
ci.2$Group <- ci.2$protocol
ci.2$Curve <- ifelse(ci.2$LookType == "Cohort", 1, 2)
out.1 <- doubleGauss.fit(ci.2, 4, diffs = TRUE)
out.1 <- doubleGauss.refit(out.1, subj = c(13, 23), group = c(2, 2),
	curves = c(2, 2), cor = c(FALSE, FALSE))
printFits(out.1)
out.2 <- doubleGauss.boot(out.1)
replot(out.2, ylim = c(-0.01, 0.1), bucket.lim = c(0, 0.08))
}}

\keyword{ htest}