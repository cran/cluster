\name{predict.ellipsoid}
\alias{predict.ellipsoid}
\alias{ellipsoidPoints}
\title{Predict Method for Ellipsoid Objects}
\description{
  Compute points on the ellipsoid boundary, mostly for drawing.
}
\usage{
\method{predict}{ellipsoid}(object, n.out=201, \dots)
ellipsoidPoints(A, d2, loc, n=201)
}
\arguments{
  \item{object}{an object of class \code{ellipsoid}, typically from
    \code{\link{ellipsoidhull}()}.}
  \item{n.out,n}{the number of points to create.}
  \item{A, d2, loc}{arguments of the auxilary \code{ellipsoidPoints}, see below.}
  \item{\dots}{passed to and from methods.}
}
\details{
   Note \code{ellipsoidPoints} is the workhorse function of the
   \code{predict.ellipsoid} method; the class of \code{object} is not
   checked; it must solely have valid components \code{loc} (length \eqn{p}),
   the \eqn{p \times p}{p x p} matrix \code{cov} (corresponding to
   \code{A}) and \code{d2} for the center, the shape (``covariance'')
   matrix and the squared average radius or distance.
}
\value{
  a numeric matrix of dimension \code{n.out} times \eqn{p}.
}
\seealso{\code{\link{ellipsoidhull}}, also for examples.}
\keyword{dplot}
\keyword{utilities}
