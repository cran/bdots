\name{replot}
\alias{replot}

\title{ Replot Boostrapped Output}

\description{ Plot the boostrapped output with different parameters than the default ones.}

\usage{
replot(part2.list, xlim = NULL, ylim = c(0, 1), main = "Curve",
	legend.location = "topleft", bucket.lim = c(0, .9))
}

\arguments{
  \item{part2.list}{ list. Output from doubleGauss.boot or logistic.boot}
  \item{xlim}{ numeric vector (length = 2). Start and end point of x-axis. If NULL, takes the full time course}
  \item{ylim}{ numeric vector (length = 2). Start and end point of y-axis}
  \item{main}{ string. Title}
  \item{legend.location}{ string. Location of the legend}
  \item{bucket.lim}{ numeric vector (length = 2). How far the yellow significant region goes
	on the y axis}
}
\details{ Plot the boostrapped output with different parameters than the default ones.}

\value{ NULL}

\note{ Options for legend.location include "topleft", "top", "topright", "right", "bottomright",
	"bottom", "bottomleft", "left"
}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)
out.2 <- logistic.boot(out.1)
replot(out.2, bucket.lim = c(0, 1))
}}

\keyword{ htest}