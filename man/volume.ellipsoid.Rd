\name{volume.ellipsoid}
\alias{volume}
\alias{volume.ellipsoid}
\title{Compute the Volume of Planar Object}
\description{
  Compute the volume of a planar object.
  This is a generic function and a method for \code{ellipsoid} objects.
}
\usage{
\method{volume}{ellipsoid}(object)
}
\arguments{
  \item{object}{an \R object the volume of which is wanted.}
}
\value{
  a number, the volume of the given \code{object}.
}
\seealso{\code{\link{ellipsoidhull}} for spanning ellipsoid computation.}
\examples{
example(ellipsoidhull) # which defines `ellipsoid' object <namefoo>
%% FIXME
}
\keyword{utilities}
