\name{pltree.twins}
\alias{pltree.twins}
\title{Clustering Tree of a Hierarchical Clustering}
\description{
  Draws a clustering tree (dendrogramm) of a \code{twins} object, i.e.,
  hierarchical clustering, typically resulting from \code{agnes} or
  \code{diana}.
}
\usage{
pltree.twins(x, main = paste("Dendrogram of ", deparse(call)),
             labels = NULL, ylab = "Height", \dots)
}
\arguments{
  \item{x}{
    an object of class \code{"twins"}, typically created by either
    \code{\link{agnes}()} or \code{\link{diana}()}.}
  \item{main}{main title.}
  \item{labels}{labels to use; the default is constructed from \code{x}.}
  \item{ylab}{label for y-axis.}
  \item{\dots}{
    Graphical parameters (see \code{\link{par}}) may also
    be supplied as arguments to this function.}
}
\value{
  a NULL value is returned.
}
\details{
  Creates a plot of a clustering tree given a \code{twins} object.  The
  leaves of the tree are the original observations.  In case of an
  agglomerative clustering, two branches come together at the distance
  between the two clusters being merged.  For a divisive clustering, a
  branch splits up at the diameter of the cluster being splitted.
}
\seealso{
  \code{\link{agnes}}, \code{\link{agnes.object}}, \code{\link{diana}},
  \code{\link{diana.object}}, \code{\link{hclust}}, \code{\link{par}},
  \code{\link{plot.agnes}}, \code{\link{plot.diana}}.
}
\examples{
}
\keyword{cluster}
\keyword{hplot}



