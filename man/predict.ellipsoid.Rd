\name{predict.ellipsoid}
\alias{predict.ellipsoid}
\alias{ellipsoidPoints}
\title{Predict Method for Ellipsoid Objects}
\description{
  Compute points on the ellipsoid boundary, mostly for drawing.
}
\usage{
%  method *and* stand alone:
predict.ellipsoid(object, n.out=201, \dots)
\method{predict}{ellipsoid}(object, n.out=201, \dots)
ellipsoidPoints(A, d2, loc, n=201)
}
\arguments{
  \item{object}{an object of class \code{ellipsoid}, typically from
    \code{\link{ellipsoidhull}()}; alternatively any list-like object
    with proper components, see details below.}
  \item{n.out,n}{the number of points to create.}
  \item{A, d2, loc}{arguments of the auxilary \code{ellipsoidPoints},
    see below.}
  \item{\dots}{passed to and from methods.}
}
\details{
   Note \code{ellipsoidPoints} is the workhorse function of
   \code{predict.ellipsoid} a standalone function and method for
   \code{ellipsoid} objects, see \code{\link{ellipsoidhull}}.
   The class of \code{object} is not checked; it must solely have valid
   components \code{loc} (length \eqn{p}), the \eqn{p \times p}{p x p}
   matrix \code{cov} (corresponding to \code{A}) and \code{d2} for the
   center, the shape (``covariance'') matrix and the squared average
   radius or distance.
}
\value{
  a numeric matrix of dimension \code{n.out} times \eqn{p}.
}
\seealso{\code{\link{ellipsoidhull}}.}
\examples{
 ## see also  example(ellipsoidhull)

## Robust vs. L.S. covariance matrix
set.seed(143)
x <- rt(200, df=3)
y <- 3*x + rt(200, df=2)
plot(x,y, main="non-normal data (N=200)")
X <- cbind(x,y)
C.ls <- cov(X) ; m.ls <- colMeans(X)
lines(ellipsoidPoints(C.ls, d2 = 2, loc=m.ls), col="green")
if(require(MASS)) {
  Cxy <- cov.rob(cbind(x,y))
  lines(ellipsoidPoints(Cxy$cov, d2 = 2, loc=Cxy$center), col="red")
}# MASS
}
\keyword{dplot}
\keyword{utilities}
