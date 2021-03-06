% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMQR.R
\name{plot.MultiQR}
\alias{plot.MultiQR}
\title{Fan plot for \code{MultiQR} objects}
\usage{
\method{plot}{MultiQR}(
  plotdata,
  targetTimes = NULL,
  quantiles = colnames(plotdata),
  ylim = "auto",
  q50_line = F,
  Legend = NULL,
  cols = colorRampPalette(c("dodgerblue1", "cyan")),
  ...
)
}
\arguments{
\item{plotdata}{A \code{MultiQR} object to be plotted.}

\item{targetTimes}{A vector of forecast target times corresponding to \code{plotdata}.}

\item{quantiles}{A charater vector. By default \code{colnames(plotdata)},
i.e. all available quantiles. Alternatively, a subset may be specified here.}

\item{ylim}{Passed to \code{plot}. Defaults to \code{"auto"} which set
ylim equal to the range of \code{plotdata} expanded by 20\%.}

\item{q50_line}{Should \code{q50} be plotted as a distinct line?
Not advised. This can mis-lead users who may interpret it as the
most likely temporal trajectory.}

\item{Legend}{Location of legend to be produced. E.g. "topleft"}

\item{cols}{Color scheme. May be a list of colors with length
\code{>=floor(length(qunatiles)/2)} or function that returns a vector of colors, e.g.
\code{colorRampPalette()}.}

\item{...}{Additional arguments passed to \code{plot()}.}
}
\value{
A plot of a \code{MQR}.
}
\description{
This S3 method for \code{plot} produces a fan plot of a MultiQR object.
}
\details{
If \code{is.null(targetTimes)} row indexes of \code{plotdata} are used.
}
\author{
Jethro Browell, \email{jethro.browell@strath.ac.uk}
}
\keyword{Quantile}
\keyword{Regression,}
\keyword{plot}
