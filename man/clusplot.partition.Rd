\name{clusplot.partition}
\alias{clusplot.partition}
\title{Bivariate Clusplot of a Partitioning Object}
\description{
  Clusplot (Clustering Plot) method for an object of class \code{partition}.
}
\usage{
\method{clusplot}{partition}(x, main = NULL, \dots)
}
\arguments{
  \item{x}{
    an object of class \code{"partition"}, e.g. created by the functions
    \code{\link{pam}}, \code{\link{clara}}, or \code{\link{fanny}}.}
  \item{main}{title for the plot; when \code{NULL} (by default), a title
    is constructed, using \code{x$call}.}
  \item{\dots}{all optional arguments available for the
    \code{\link{clusplot.default}} function (except for the \code{diss}
    one) may also be supplied to this function.  Graphical parameters
    (see \code{\link{par}}) may also be supplied as arguments to this
    function.}
}
\value{
  An invisible list with components
  \item{Distances}{
    When option lines is 1 or 2 we optain a k by k matrix (k is the number of
    clusters). The element at row j and column s is the distance between
    ellipse j and ellipse s.  If lines=0, then the value of this
    component is NA.
  }
  \item{Shading}{
    A vector of length k (where k is the number of clusters), containing
    the amount of shading per cluster.  Let y be a vector where element
    i is the ratio between the number of objects in cluster i and the
    area of ellipse i.  When the cluster i is a line segment, y[i] and
    the density of the cluster are set to NA.  Let z be the sum of all
    the elements of y without the NAs.  Then we put shading = y/z *37 + 3.
  }
}
\details{
  This \code{clusplot.partition()} method relies on
  \code{\link{clusplot.default}}.

  If the clustering algorithms \code{pam}, \code{fanny} and \code{clara}
  are applied to a data matrix of observations-by-variables then a
  clusplot of the resulting clustering can always be drawn.  When the
  data matrix contains missing values and the clustering is performed
  with \code{\link{pam}} or \code{\link{fanny}}, the dissimilarity
  matrix will be given as input to \code{clusplot}.  When the clustering
  algorithm \code{\link{clara}} was applied to a data matrix with NAs
  then clusplot will replace the missing values as described in
  \code{\link{clusplot.default}}, because a dissimilarity matrix is not
  available.
}
\seealso{\code{\link{clusplot.default}} for references;
  \code{\link{partition.object}}, \code{\link{pam}},
  \code{\link{pam.object}}, \code{\link{clara}},
  \code{\link{clara.object}}, \code{\link{fanny}},
  \code{\link{fanny.object}}, \code{\link{par}}.
}
\examples{
## generate 25 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
clusplot(pam(x, 2))
## add noise, and try again :
x4 <- cbind(x, rnorm(25), rnorm(25))
clusplot(pam(x4, 2))
}
\keyword{cluster}
\keyword{hplot}
