% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{pltree.twins}
\alias{pltree.twins}
\title{
Clustering Tree of an Agglomerative or a Divisive Hierarchical Clustering
}
\description{
Creates a clustering tree of a \code{twins} object.
}
\usage{
pltree.twins(x, \dots)


}
\arguments{
\item{x}{
an object of class \code{"twins"}, created by the functions \code{agnes} or \code{diana}.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.


}}
\value{
a NULL value is returned.


}
\section{Side Effects}{
A clustering tree is plotted on the current graphics device.


}
\details{
Creates a plot of a clustering tree given a \code{twins} object.
The leaves of the tree are the original observations. In case of an agglomerative
clustering, two branches come together
at the distance between the two clusters being merged. For a divisive
clustering, a branch splits up at the diameter of the cluster being splitted.


}
\seealso{
\code{\link{agnes}}, \code{\link{agnes.object}}, \code{\link{diana}}, \code{\link{diana.object}}, \code{\link{hclust}}, \code{\link{par}}, 
\code{\link{plot.agnes}}, \code{\link{plot.diana}}, \code{\link{plclust}}.
}
\keyword{libcluster}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
