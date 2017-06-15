\name{Boostrap Step}
\alias{doubleGauss.boot}
\alias{logistic.boot}

\title{ Bootstrap on the Fitted Parameters}

\description{ Bootstrap on the fitted parameters, plot the estimates, and highlight the significant regions}

\usage{
doubleGauss.boot(part1.list, seed = new.seed(), alpha = 0.05, paired = FALSE,
	N.iter = 1000, cores = 1, p.adj = "oleson", test.spots = NULL,
	time.test = NULL, test.params = FALSE)
logistic.boot(part1.list, seed = new.seed(), alpha = 0.05, paired = FALSE,
	N.iter = 1000, cores = 1, p.adj = "oleson", test.spots = NULL,
	time.test = NULL, test.params = FALSE)
}

\arguments{
  \item{part1.list}{list. Output from doubleGauss.fit}
  \item{seed}{integer. What to set seed at}
  \item{alpha}{numeric (Between 0 and 1). Probability of familywise Type I Error}
  \item{paired}{boolean. Whether the same subjects are in both data sets}
  \item{N.iter}{numeric (positive integer). Number of bootstrap iterations to run}
  \item{cores}{integer. Number of cores on the localhost to use}
  \item{p.adj}{Options: oleson, fdr, none}
  \item{test.spots}{numeric. Specify specific x-values for testing at}
  \item{time.test}{numeric. Specify individual time points to conduct t-tests at without any p-value correction}
  \item{test.params}{boolean. Whether to test for significant differences in group mean parameter estimates. Performs a 2-sample t-test with equal variance assumption if paired = FALSE. If paired = TRUE, performs a paired t-test.}
}
\details{ Bootstrap on the fitted parameters, plot the estimates, and highlight the significant regions}

\value{ List: for input in replot}

\note{ There are no further notes}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)
out.2 <- logistic.boot(out.1)
replot(out.2, bucket.lim = c(0, 1))

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