% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{mona.object}
\alias{mona.object}
\title{
Monothetic Analysis Object
}
\description{
These are objects of class
"mona"
They represent the divisive hierarchical clustering of a dataset with
binary measurements only.


}
\section{GENERATION}{
This class of objects is returned from \code{mona}.


}
\section{METHODS}{
The \code{"mona"} class has methods for the following generic functions:
\code{print}, \code{summary}, \code{plot}.


}
\section{STRUCTURE}{
A legitimate \code{mona} object is a list with the following components:


}
\arguments{
\item{data}{
matrix with the same dimensions as the original data matrix,
but with factors coded as 0 and 1, and all missing values replaced.


}
\item{order}{
a vector giving a permutation of the original observations to allow for 
plotting, in the sense that the branches of a clustering tree will not cross.


}
\item{order.lab}{
a vector similar to \code{order}, but containing observation labels instead of 
observation numbers. This component is only available if the original
observations were labelled.


}
\item{variable}{
vector with length n-1, where n is the number of observations,
specifying the variables used to separate the observations of \code{order}.


}
\item{step}{
vector with length n-1, where n is the number of observations,
specifying the separation steps at which the observations of \code{order} are
separated.
}}
\seealso{
\code{\link{mona}}, \code{\link{plot.mona}}.
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
