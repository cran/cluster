% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{diana.object}
\alias{diana.object}
\title{
Divisive Analysis Object
}
\description{
These are objects of class
"diana"
They represent a divisive hierarchical clustering of a dataset.


}
\section{GENERATION}{
This class of objects is returned from \code{diana}.


}
\section{METHODS}{
The \code{"diana"} class has methods for the following generic functions:
\code{print}, \code{summary}, \code{plot}.


}
\section{INHERITANCE}{
The class \code{"diana"} inherits from \code{"twins"}.
Therefore, the generic function \code{pltree} can be used on a \code{diana} object.


}
\section{STRUCTURE}{
A legitimate \code{diana} object is a list with the following components:


}
\arguments{
\item{order}{
a vector giving a permutation of the original observations to allow for 
plotting, in the sense that the branches of a clustering tree will not cross.


}
\item{order.lab}{
a vector similar to \code{order}, but containing observation labels instead of 
observation numbers. This component is only available if the original
observations were labelled.


}
\item{height}{
a vector with the diameters of the clusters prior to splitting.


}
\item{dc}{
the divisive coefficient, measuring the clustering structure of the dataset.
For each observation i, denote by d(i) the diameter of the last cluster to 
which it belongs (before being split off as a single observation), 
divided by the diameter of the whole dataset.
The \code{dc} is the average of all 1 - d(i). It can also be seen as the average
width (or the percentage filled) of the banner plot.
Because \code{dc} grows with the number of observations, this measure should not
be used to compare datasets of very different sizes.


}
\item{merge}{
an (n-1) by 2 matrix, where n is the number of observations. Row i of \code{merge}
describes the split at step n-i of the clustering. If a number
j in row r is negative, then the single observation |j| is split off at
stage n-r. If j is positive, then the cluster that will be splitted at
stage n-j (described by row j), is split off at stage n-r.


}
\item{diss}{
an object of class \code{"dissimilarity"}, representing the total dissimilarity
matrix of the dataset.


}
\item{data}{
a matrix containing the original or standardized measurements, depending 
on the \code{stand} option of the function \code{agnes}. If a dissimilarity matrix was 
given as input structure, then this component is not available.


}}
\seealso{
\code{\link{agnes}}, \code{\link{diana}}, \code{\link{plot.diana}}, \code{\link{twins.object}}.
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
