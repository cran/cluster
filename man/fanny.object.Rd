% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{fanny.object}
\alias{fanny.object}
\title{
Fuzzy Analysis Object
}
\description{
These are objects of class
"fanny"
They represent a fuzzy clustering of a dataset.


}
\section{GENERATION}{
This class of objects is returned from \code{fanny}.


}
\section{METHODS}{
The \code{"fanny"} class has methods for the following generic functions:
\code{print}, \code{summary}.


}
\section{INHERITANCE}{
The class \code{"fanny"} inherits from \code{"partition"}.
Therefore, the generic functions \code{plot} and \code{clusplot} can be used on a 
\code{fanny} object.


}
\section{STRUCTURE}{
A legitimate \code{fanny} object is a list with the following components:


}
\arguments{
\item{objective}{
the objective function and the number of iterations the \code{fanny} algorithm
needed to reach this minimal value.


}
\item{membership}{
matrix containing the memberships for each pair consisting of an
observation and a cluster.


}
\item{coeff}{
Dunn's partition coefficient F(k) of the clustering, where k is the number
of clusters.
F(k) is the sum of all squared membership coefficients,
divided by the number of observations. Its value is always between 1/k and 1.
The normalized form of the coefficient is also given. It is defined as
(F(k) - 1/k) / (1 - 1/k), and ranges between 0 and 1.
A low value of Dunn's coefficient indicates a very fuzzy clustering,
whereas a value close to 1 indicates a near-crisp clustering.


}
\item{clustering}{
the clustering vector of the nearest crisp clustering. A vector with length
equal to the number of observations, giving for each observation the number 
of the cluster to which it has the largest membership.


}
\item{silinfo}{
list with all information necessary to construct a silhouette plot of the
nearest crisp clustering. This list is only available when 1 < k < n.
The first component is a matrix, with for each observation i the cluster to 
which i belongs, as well as the neighbor cluster of i (the cluster, 
not containing i, for which the average dissimilarity between its observations 
and i is minimal), and the silhouette width of i.
The other two components give the average silhouette width per cluster and
the average silhouette width for the dataset.
See \code{plot.partition} for more information.


}
\item{diss}{
an object of class \code{"dissimilarity"}, representing the total dissimilarity
matrix of the dataset.


}
\item{data}{
a matrix containing the original or standardized measurements, depending 
on the \code{stand} option of the function \code{fanny}. If a dissimilarity matrix was 
given as input structure, then this component is not available.


}}
\seealso{
\code{\link{fanny}}, \code{\link{dissimilarity.object}}, \code{\link{partition.object}}, \code{\link{plot.partition}}.
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
