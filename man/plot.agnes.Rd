% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{plot.agnes}
\alias{plot.agnes}
\title{
Plots of an Agglomerative Hierarchical Clustering
}
\description{
Creates plots for visualizing an \code{agnes} object.
}
\usage{
plot.agnes(x, ask=F, \dots)
}
\arguments{
\item{x}{
an object of class \code{"agnes"}, created by the function \code{agnes}.
}
\item{ask}{
if TRUE, \code{plot.agnes} operates in interactive mode.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.


}}
\value{
a NULL value is returned.


}
\section{Side Effects}{
An appropriate plot is produced on the current graphics device. This can
be one or both of the following choices:
Banner
Clustering tree
}
\details{
When \code{ask=T}, rather than producing each plot sequentially,
\code{plot.agnes} displays a menu listing all the plots that can be produced.
If the menu is not desired but a pause between plots is still wanted
one must set \code{par(ask=T)} before invoking the plot command.


The banner displays the hierarchy of clusters, and is equivalent to a tree.
See Rousseeuw (1986) or chapter 5 of Kaufman and Rousseeuw (1990).
The banner plots distances at which observations and clusters are merged.
The observations are listed in the order found by the \code{agnes} algorithm, 
and the numbers in the \code{height} vector are represented as bars between the 
observations.


The leaves of the clustering tree are the original observations. 
Two branches come together at the distance between the two clusters being 
merged.
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


Rousseeuw, P.J. (1986). A visual display for hierarchical classification,
in 
\emph{Data Analysis and Informatics 4.}
Edited by E. Diday, Y. Escoufier,
L. Lebart, J. Pages, Y. Schektman, and R. Tomassone. North-Holland,
Amsterdam. pp. 743-748.


Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust 
Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis, \bold{26}, 17-37.}
}
\seealso{
\code{\link{agnes}}, \code{\link{agnes.object}}, \code{\link{twins.object}}, \code{\link{par}}.
}
\keyword{libcluster}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
