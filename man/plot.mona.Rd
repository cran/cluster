% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{plot.mona}
\alias{plot.mona}
\title{
Banner of Monothetic Divisive Hierarchical Clusterings
}
\description{
Creates the banner of a \code{mona} object.
}
\usage{
plot.mona(x, \dots)


}
\arguments{
\item{x}{
an object of class \code{"mona"}, created by the function \code{mona}.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.


}}
\value{
a NULL value is returned.


}
\section{Side Effects}{
A banner is plotted on the current graphics device.


}
\details{
Plots the separation step at which clusters are splitted.
The observations are given in the order found by the \code{mona} algorithm, the
numbers in the \code{step} vector are represented as bars between the observations.
When a long bar is drawn between two observations,
those observations have the same value for each variable.
See chapter 7 of Kaufman and Rousseeuw (1990).
}
\section{NOTE}{
In the banner plot,
observation labels are only printed when the number of observations is 
limited to less than 35, for readability.
Moreover, observation labels are truncated to at most 5 characters.


}
\references{
Kaufman, L. and Rousseeuw, P.J. (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis. }
Wiley, New York.


Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust 
Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis, \bold{26}, 17-37.}
}
\seealso{
\code{\link{mona}}, \code{\link{mona.object}}, \code{\link{par}}.
}
\keyword{libcluster}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
