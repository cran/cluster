% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{dissimilarity.object}
\alias{dissimilarity.object}
\title{
Dissimilarity Matrix Object
}
\description{
These are objects of class
"dissimilarity"
They represent the dissimilarity matrix of a dataset.


}
\section{GENERATION}{
\code{daisy} returns this class of objects.
Also the functions \code{pam}, \code{clara}, \code{fanny}, \code{agnes}, and \code{diana} return 
a \code{dissimilarity} object, as one component of their return objects.


}
\section{METHODS}{
The \code{"dissimilarity"} class has methods for the following generic functions:
\code{print}.


}
\section{STRUCTURE}{
The dissimilarity matrix is symmetric, and hence is represented as a vector
to save storage space.
For i less than j, the dissimilarity between row i and row j is element
nrow(x)*(i-1) - i*(i-1)/2 + j-i of that vector.
The length of the vector is nrow(x)*(nrow(x)-1)/2,
that is, it is of order nrow(x) squared.
The object has the following attributes:


}
\arguments{
\item{Size}{
the number of observations in the dataset.


}
\item{Metric}{
the metric used for calculating the dissimilarities. Possible values are
"euclidean", "manhattan", "mixed" (if variables of different types were
present in the dataset), and "unspecified".


}
\item{Labels}{
optionally, contains the labels, if any, of the observations of the dataset.


}
\item{NA.message}{
optionally, if a dissimilarity could not be computed, because of too many
missing values for some observations of the dataset.


}}
\seealso{
\code{\link{daisy}}, \code{\link{dist}}, \code{\link{pam}}, \code{\link{clara}}, \code{\link{fanny}}, \code{\link{agnes}}, \code{\link{diana}}.
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
