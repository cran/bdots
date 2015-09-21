\name{Fitting Step}
\alias{doubleGauss.fit}
\alias{logistic.fit}

\title{ Fit Subjects Individual Curves}

\description{ Fit Subjects from 2 groups with the 6-parameter Double Gaussian or the 4-parameter Logistic}

\usage{
doubleGauss.fit(data, col, concave = TRUE, diffs = FALSE,
	rho.0 = 0.9, cor = TRUE, cores = 1)
logistic.fit(data, col, diffs = FALSE, rho.0 = 0.9, cor = TRUE, cores = 1)
}

\arguments{
  \item{data}{ data.frame. A data.frame with the columns 'Subject', 'Time', and 'Group'. 'Subject' designates
	the subject number (numeric), 'Time' designates the time (numeric, should be ordered from low to high),
	and 'Group' designates the group (should be 2 unique groups)}
  \item{col}{ numeric. The column in the data.frame that corresponds to the eyetracking}
  \item{concave}{boolean. TRUE indicates concave UP, FALSE indicates concave DOWN. Only for Double Gaussian.}
  \item{diffs}{ boolean. If the each group is calculating the difference of 2 logistic curves, set to TRUE.
	In this case, there needs to be a numeric 'Curve' column (with only 1s and 2s) designating the
	secondary curve (2) to subtract from the primary curve (1).}
  \item{rho.0}{ numeric (Between 0 and 1). Assumed autocorrelation of errors for individual subject's curve}
  \item{cor}{ boolean. If TRUE assumes an AR1 autocorrelation structure among the residuals for fitting }
  \item{cores}{integer. Number of cores on the localhost to use}
}
\details{ Fit Subjects from 2 groups with the 6-parameter Double Gaussian or the 4-parameter Logistic}

\value{ List: for input into refit, boot, plot.ests, and plot.subjs functions}

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