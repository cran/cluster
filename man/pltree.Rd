\name{pltree}
\alias{pltree}
\title{Clustering Trees - Generic Function}
\description{
  Generic function drawing a clustering tree (``dendrogram'') on the
  current graphics device.
  There is Functions with names beginning in \code{pltree.} will be
  methods for this function.
}
\usage{
pltree(x, \dots)
}
\arguments{
  \item{x}{an \R object (for which a \code{pltree} method is defined).}
  \item{...}{methods may have additional arguments.  Graphical
    parameters (see \code{\link{par}}) may also be supplied as arguments
    to this function.}
}
\section{Side Effects}{
  a clustering tree is created on the current graphics device.
}
\seealso{
  \code{\link{twins.object}}, \code{\link{agnes}}, \code{\link{diana}},
  \code{\link{par}}.
}
\keyword{cluster}
\keyword{hplot}

