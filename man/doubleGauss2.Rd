% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doubleGauss2.R
\name{doubleGauss2}
\alias{doubleGauss2}
\title{DoubleGauss2 curve function for nlme}
\usage{
doubleGauss2(dat, y, time, params = NULL, concave = TRUE, ...)
}
\arguments{
\item{dat}{subject data to be used}

\item{y}{outcome variable, character vector}

\item{time}{time variable, character vector}

\item{params}{\code{NULL} unless user wants to specify starting parameters for gnls}

\item{concave}{Boolean}

\item{...}{just in case}
}
\description{
DoubleGauss2 function used in fitting nlme curve for observations
}
\details{
User should only have to worry about setting concavity
of this function. Presently only work for time series scaled out to 2000ms

\code{y ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
/ (2 * sig1 ^ 2)) * (ht - base1) + base1)
+ (mu <= time) * (exp(-1 * (time - mu) ^ 2
                         / (2 * sig2 ^ 2)) * (ht - base2) + base2)}
}
