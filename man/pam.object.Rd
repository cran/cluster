% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{pam.object}
\alias{pam.object}
\title{
Partitioning Around Medoids Object
}
\description{
These are objects of class
"pam"
They represent a partitioning of a dataset into clusters.


}
\section{GENERATION}{
This class of objects is returned from \code{pam}.


}
\section{METHODS}{
The \code{"pam"} class has methods for the following generic functions:
\code{print}, \code{summary}.


}
\section{INHERITANCE}{
The class \code{"pam"} inherits from \code{"partition"}.
Therefore, the generic functions \code{plot} and \code{clusplot} can be used on a 
\code{pam} object.


}
\section{STRUCTURE}{
A legitimate \code{pam} object is a list with the following components:


}
\arguments{
\item{medoids}{
the medoids or representative objects of the clusters. If a dissimilarity
matrix was given as input to \code{pam}, then a vector of numbers or labels of
observations is given, else \code{medoids} is a matrix with in each row the
coordinates of one medoid.


}
\item{clustering}{
the clustering vector. A vector with length equal to the number of 
observations, giving the number of the cluster to which each observation 
belongs.


}
\item{objective}{
the objective function after the first and second step of the \code{pam}
algorithm.


}
\item{clusinfo}{
matrix, each row gives numerical information for one cluster. These are
the cardinality of the cluster (number of observations), the maximal and 
average dissimilarity between the observations in the cluster and the 
cluster's medoid, the diameter of the cluster (maximal dissimilarity between 
two observations of the cluster), and the separation of the cluster (minimal 
dissimilarity between an observation of the cluster and an observation of 
another cluster).


}
\item{isolation}{
vector with length equal to the number of clusters, specifying which
clusters are isolated clusters (L- or L*-clusters) and which clusters are
not isolated.
A cluster is an L*-cluster iff its diameter is smaller than its separation.
A cluster is an L-cluster iff for each observation i the maximal dissimilarity
between i and any other observation of the cluster is smaller than the minimal
dissimilarity between i and any observation of another cluster.
Clearly each L*-cluster is also an L-cluster.


}
\item{silinfo}{
list with all information necessary to construct a silhouette plot of the
clustering. This list is only available when 1 < k < n.
The first component is a matrix, with for each observation i the cluster to 
which i belongs, as well as the neighbor cluster of i (the cluster, not 
containing i, for which the average dissimilarity between its observations 
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
on the \code{stand} option of the function \code{pam}. If a dissimilarity matrix was 
given as input structure, then this component is not available.


}}
\seealso{
\code{\link{pam}}, \code{\link{dissimilarity.object}}, \code{\link{partition.object}}, \code{\link{plot.partition}}.
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
