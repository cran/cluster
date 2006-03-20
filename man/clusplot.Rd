\name{clusplot}
\alias{clusplot}
\title{Cluster Plot - Generic Function}
\description{
  Draws a 2-dimensional \dQuote{clusplot} on the current graphics device.
  This is a generic function with a default and \code{partition} method.
}
\usage{
clusplot(x, \dots)
}
\arguments{
  \item{x}{an \R object.}
  \item{\dots}{additional arguments for \code{\link{methods}}.
    Graphical parameters (see \code{\link{par}}) may also
    be supplied as arguments to this function.}
}
\section{Side Effects}{
  a 2-dimensional clusplot is created on the current graphics device.
}
\seealso{(for references and examples)
  \code{\link{clusplot.default}}, \code{\link{clusplot.partition}},
  \code{\link{partition.object}}, \code{\link{pam}}, \code{\link{fanny}},
  \code{\link{clara}}.
}
\keyword{cluster}
\keyword{hplot}

