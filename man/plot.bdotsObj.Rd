% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plot.bdotsObj}
\alias{plot.bdotsObj}
\title{Plot a bdotsFit object}
\usage{
\method{plot}{bdotsObj}(x, fitCode = NULL, gridSize = NULL, plotfun = "fits", ...)
}
\arguments{
\item{x}{An object of class 'bdotsObj' returned from \code{bdotsFit}}

\item{fitCode}{Currently not used}

\item{gridSize}{Length one numeric indicating size of plot grid. Default is
2x2. For right now, they are square}

\item{plotfun}{Plot either subject fits or model parameters with "fits" or "pars"}

\item{...}{ignore for now (other args to plot.generic)}
}
\value{
This will return a list of all of the plots rendered.
}
\description{
Plot individual fits or model fit parameters from an object of class
'bdotsObj'. These functions are not very stable
}
\details{
Right now, these functions are a bit unstable and expected to change.
The largest current issue is with the placement of the legend, which cannot
be adjusted. If you are running into issues with seeing things correctly, try
making the "Plots" window in RStudio larger before running this function
}
