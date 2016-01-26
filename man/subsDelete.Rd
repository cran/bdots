\name{subs.delete}
\alias{subs.delete}

\title{ Delete individual subjects from both data and fits}

\description{ Delete individual subjects from both data and fits}

\usage{
subs.delete(part1.list, subj, group)
}

\arguments{
  \item{part1.list}{ Output from doubleGauss.fit or logistic.fit}
  \item{subj}{ numeric. Subject numbers to remove.}
  \item{group}{ numeric. Corresponding group numbers to subjects.}

}
\details{ Delete individual subjects from both data and fits}

\value{ part1.list}

\examples{
\dontrun{
data(ci)
ci.1 <- subset(ci, ci$LookType == "Target")
ci.1$Group <- ci.1$protocol
out.1 <- logistic.fit(ci.1, 4)ests.plot(out.1)
#Remove subject 1 from group 1 and subject 10 from group 2
out.1 <- subjs.delete(out.1, subj = c(1, 10), group = c(1, 2))
}}

\keyword{ htest}