\name{plot.mona}
\alias{plot.mona}
\title{Banner of Monothetic Divisive Hierarchical Clusterings}
\description{
  Creates the banner of a \code{mona} object.
}
\usage{
plot.mona(x, main = paste("Banner of ", deparse(attr(x, "Call"))),
          col = 2, axes = TRUE, adj = 0,
          nmax.lab = 35, max.strlen = 5, \dots)
}
\arguments{
  \item{x}{an object of class \code{"mona"}, typically created by
    \code{\link{mona}(.)}.}
  \item{main}{main title for the plot, with a convenient
    default.  See documentation for \code{main} in \code{\link{plot.default}}.}
  \item{nmax.lab}{integer indicating the number of labels which is
    considered too large for labeling.}
  \item{max.strlen}{positive integer giving the length to which
    strings are truncated in labeling.}
  \item{axes}{logical, indicating if (labeled) axes should be drawn.}
  \item{col,adj,\dots}{graphical parameters (see \code{\link{par}}) may also
    be supplied as arguments to this function.}
}
\section{Side Effects}{
  A banner is plotted on the current graphics device.
}
\details{
  Plots the separation step at which clusters are splitted.  The
  observations are given in the order found by the \code{mona}
  algorithm, the numbers in the \code{step} vector are represented as
  bars between the observations.

  When a long bar is drawn between two observations,
  those observations have the same value for each variable.
  See chapter 7 of Kaufman and Rousseeuw (1990).
}
\note{
  In the banner plot,
  observation labels are only printed when the number of observations is
  limited less than \code{nmax.lab} (35, by default), for readability.
  Moreover, observation labels are truncated to maximally
  \code{max.strlen} (5) characters.
}
\references{see those in \code{\link{plot.agnes}}.}
\seealso{
  \code{\link{mona}}, \code{\link{mona.object}}, \code{\link{par}}.
}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
