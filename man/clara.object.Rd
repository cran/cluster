% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{clara.object}
\alias{clara.object}
\title{Clustering Large Applications (CLARA) Object}
\description{
  The objects of class \code{"clara"} represent a partitioning of a large
  dataset into clusters.
}
\section{GENERATION}{
  This class of objects is returned from \code{\link{clara}}.
}
\section{METHODS}{
  The \code{"clara"} class has methods for the following generic functions:
  \code{print}, \code{summary}.
}
\section{INHERITANCE}{
  The class \code{"clara"} inherits from \code{"partition"}.
  Therefore, the generic functions \code{plot} and \code{clusplot} can
  be used on a \code{clara} object.
}
\value{
A legitimate \code{clara} object is a list with the following components:

\item{sample}{
labels or case numbers of the observations in the best sample, that is,
the sample used by the \code{clara} algorithm for the final partition.
}

\item{medoids}{
the medoids or representative objects of the clusters.
It is a matrix with in each row the coordinates of one medoid.
}

\item{clustering}{
the clustering vector. A vector with length equal to the number of 
observations, giving the number of the cluster to which each observation 
belongs.
}

\item{objective}{
  the objective function for the final clustering of the entire dataset.}

\item{clusinfo}{
matrix, each row gives numeric information for one cluster. These are
the cardinality of the cluster (number of observations), and the maximal and
average dissimilarity between the observations in the cluster and the cluster's
medoid. The last column is the maximal dissimilarity
between the observations in the cluster and the cluster's medoid, divided by
the minimal dissimilarity between the cluster's medoid and the medoid of
any other cluster. If this ratio is small, the cluster is well-separated
from the other clusters.
}

\item{silinfo}{
list with all information necessary to construct a silhouette plot of the
clustering of the best sample. This list is only available when 1 < k < n.
The first component is a matrix, with for each observation i in the best 
sample, the cluster to which i belongs,
as well as the neighbor cluster of i (the cluster, not containing
i, for which the average dissimilarity between its observations and i is
minimal), and the silhouette width of i.
The other two components give the average silhouette width per cluster and
the average silhouette width for the best sample.
See \code{plot.partition} for more information.
}

\item{diss}{
an object of class \code{"dissimilarity"}, representing the total dissimilarity
matrix of the dataset.
}
\item{data}{
a matrix containing the original or standardized measurements, depending 
on the \code{stand} option of the function \code{clara}. 
}
}
\seealso{
  \code{\link{clara}}, \code{\link{dissimilarity.object}},
  \code{\link{partition.object}}, \code{\link{plot.partition}}.
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
