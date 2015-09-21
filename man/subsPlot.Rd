\name{subs.plot}
\alias{subs.plot}

\title{ Plot subjects raw data along with function fits}

\description{ Plot subjects raw data along with function fits}

\usage{
subs.plot(part1.list, legend.spot = "topright")
}

\arguments{
  \item{part1.list}{ Output from doubleGauss.fit or logistic.fit}
  \item{legend.spot}{ string. Location of the legend}

}
\details{ Plot subjects raw data along with function fits}

\value{ NULL}

\note{ Options for legend.location include "topleft", "top", "topright", "right", "bottomright",
	"bottom", "bottomleft", "left"
}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)ests.plot(out.1)
subs.plot(out.1)
}}

\keyword{ htest}