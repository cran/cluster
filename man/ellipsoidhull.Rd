\name{ellipsoidhull}
\alias{ellipsoidhull}
\alias{print.ellipsoid}
\title{Compute the Ellipsoid Hull or Spanning Ellipsoid of a Point Set}
\description{
  Compute the ``ellipsoid hull'' or ``spanning ellipsoid'', i.e. the
  ellipsoid of minimal volume (`area' in 2D) such that all given points
  lie just inside or on the boundary of the ellipsoid.
}
\usage{
ellipsoidhull(x, tol=0.01, maxit=5000,
              ret.wt = FALSE, ret.sqdist = FALSE, ret.pr = FALSE)
\method{print}{ellipsoid}(x, digits = max(1, getOption("digits") - 2), \dots)
}
\arguments{
  \item{x}{the \eqn{n} \eqn{p}-dimensional points  asnumeric
    \eqn{n\times p}{n x p} matrix.}
  \item{tol}{convergence tolerance for Titterington's algorithm.
    Setting this to much smaller values may drastically increase the number of
    iterations needed, and you may want to increas \code{maxit} as well.}
  \item{maxit}{integer giving the maximal number of iteration steps for
    the algorithm.}
  \item{ret.wt, ret.sqdist, ret.pr}{logicals indicating if additional
    information should be returned, \code{ret.wt} specifying the
    \emph{weights}, \code{ret.sqdist} the \emph{\bold{sq}uared
      \bold{dist}ances} and \code{ret.pr} the final \bold{pr}obabilities
    in the algorithms.}
  \item{digits,\dots}{the usual arguments to \code{\link{print}} methods.}
}
\details{
  The ``spanning ellipsoid'' algorithm is said to stem from
  Titterington(1976), in Pison et al(1999) who use it for
  \code{\link{clusplot.default}}.\cr
  The problem can be seen as a special case of the ``Min.Vol.''
  ellipsoid of which a more more flexible and general implementation is
  \code{\link[lqs]{cov.mve}} in the \code{lqs} package.
}
\value{
  an object of class \code{"ellipsoid"}, basically a \code{\link{list}}
  with several components, comprising at least
  \item{cov}{\eqn{p\times p}{p x p} \emph{covariance} matrix description
    the ellipsoid.}
  \item{loc}{\eqn{p}-dimensional location of the ellipsoid center.}
  \item{d2}{average squared radius.}
  \item{wt}{the vector of weights iff \code{ret.wt} was true.}
  \item{sqdist}{the vector of squared distances iff \code{ret.sqdist} was true.}
  \item{prob}{the vector of algorithm probabilities iff \code{ret.pr} was true.}
  \item{it}{number of iterations used.}
  \item{tol, maxit}{just the input argument, see above.}
  \item{eps}{the achieved tolerance which is the maximal squared radius
    minus \eqn{p}.}
  \item{ierr}{error code as from the algorithm; \code{0} means \emph{ok}.}
  \item{conv}{logical indicating if the converged.  This is defined as
    \code{it < maxit && ierr == 0}.}
}
\references{
  Pison, G., Struyf, A. and Rousseeuw, P.J. (1999)
  Displaying a Clustering with CLUSPLOT,
  \emph{Computational Statistics and Data Analysis}, \bold{30}, 381--392.\cr
  A version of this is available as technical report from
  \url{http://win-www.uia.ac.be/u/statis/abstract/Disclu99.htm}

  D.N. Titterington. (1976)
  Algorithms for computing {D}-optimal design on finite design spaces.  In
  \emph{Proc.\ of the 1976 Conf.\ on Information Science and Systems},
  213--216; John Hopkins University.
}

\author{Martin Maechler did the present class implementation; Rousseeuw
  et al did the underlying code.}
\seealso{\code{\link[base]{chull}} for the convex hull,
  \code{\link{clusplot}} which makes use of this; \code{\link[lqs]{cov.mve}}.}
\examples{
x <- rnorm(100)
xy <- unname(cbind(x, rnorm(100) + 2*x + 10))
exy <- ellipsoidhull(xy)
exy # >> calling print.ellipsoid()

plot(xy)
lines(predict(exy))
points(rbind(exy$loc), col = "red", cex = 3, pch = 13)

exy <- ellipsoidhull(xy, tol = 1e-7, ret.wt = TRUE, ret.sq = TRUE)
str(exy) # had small `tol', hence many iterations
(ii <- which(zapsmall(exy $ wt) > 1e-6)) # only about 4 to 6 points
round(exy$wt[ii],3); sum(exy$wt[ii]) # sum to 1
}
\keyword{dplot}
\keyword{hplot}% << ? chull has "hplot" as well.
