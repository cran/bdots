\name{Refitting Step}
\alias{doubleGauss.refit}
\alias{logistic.refit}

\title{ Refit Subjects Individual Curves}

\description{ Refit Subjects from 2 groups with the 6-parameter Double Gaussian or the 4-parameter Logistic.
	Can specify starting parameters}

\usage{
doubleGauss.refit(part1.list, subj, group, curves = NULL,
	params = NULL, cor=NULL, rho.0 = NULL)
logistic.refit(part1.list, subj, group, curves = NULL,
	params = NULL, cor = NULL, rho.0 = NULL)
}

\arguments{
  \item{part1.list}{Output from fitting step}
  \item{subj}{numeric vector. Subject numbers (within their group) that you want to refit}
  \item{group}{numeric vector. Group numbers corresponding to subject vector}
  \item{curves}{numeric vector. Curve numbers if it's a differences of logistics fit}
  \item{params}{list of numeric vectors (length 4 or 6). Parameter estimates for the 6 parameter DoubleGauss in the
	order mu, height, variance 1, variance 2, base 1, base 2. Parameter estimate for the 4 parameter Logistic in the
	order min, peak, slope, crossover}
  \item{cor}{logical vector. If TRUE assumes an correlation structure of AR(rho) and if FALSE assumes no
	correlation structure}
  \item{rho.0}{numeric. assumed autocorrelation of errors for subject's curve}

}
\details{ Refit Subjects from 2 groups with the 6-parameter Double Gaussian or the 4-parameter Logistic.
	Can specify starting parameters}

\value{ List}

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
out.2 <- doubleGauss.boot(out.1)
replot(out.2, ylim = c(-0.01, 0.1), bucket.lim = c(0, 0.08))
}}

\keyword{ htest}