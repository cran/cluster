\name{plot.diana}
\alias{plot.diana}
\title{Plots of a Divisive Hierarchical Clustering}
\description{
  Creates plots for visualizing a \code{diana} object.
}
\usage{
plot.diana(x, ask = FALSE, which.plots = NULL,
         main = paste("Banner of ", deparse(attr(x, "Call"))),
         sub  = paste("Divisive Coefficient = ", round(x$dc, digits = 2)),
         adj = 0, nmax.lab = 35, max.strlen = 5, \dots)
}
\arguments{
  \item{x}{an object of class \code{"diana"}, created by the function
    \code{diana}.}
  \item{ask}{if TRUE, \code{plot.diana} operates in interactive mode.}

  \item{\dots}{Graphical parameters (see \code{\link{par}}) may also
    be supplied as arguments to this function.}
}
\section{Side Effects}{
An appropriate plot is produced on the current graphics device. This can
be one or both of the following choices:
Banner
Clustering tree
}
\details{
When \code{ask= TRUE}, rather than producing each plot sequentially,
\code{plot.diana} displays a menu listing all the plots that can be produced.
If the menu is not desired but a pause between plots is still wanted
one must set \code{par(ask= TRUE)} before invoking the plot command.

The banner displays the hierarchy of clusters, and is equivalent to a tree.
See Rousseeuw (1986) or chapter 6 of Kaufman and Rousseeuw (1990).
The banner plots the diameter of each cluster being splitted.
The observations are listed in the order found by the \code{diana} algorithm, and 
the numbers in the \code{height} vector are represented as bars between the 
observations.


The leaves of the clustering tree are the original observations.
A branch splits up at the diameter of the cluster being splitted.
}
\note{
In the banner plot,
observation labels are only printed when the number of observations is 
limited to less than 35, for readability.

Moreover, observation labels are truncated to at most 5 characters.
}
\references{
Kaufman, L. and Rousseeuw, P.J. (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis.}
Wiley, New York.

Rousseeuw, P.J. (1986). A visual display for hierarchical classification,
in
\emph{Data Analysis and Informatics 4.}
Edited by E. Diday, Y. Escoufier,
L. Lebart, J. Pages, Y. Schektman, and R. Tomassone. North-Holland,
Amsterdam. pp. 743-748.

Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust 
Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis}, \bold{26}, 17-37.
}
\seealso{
  \code{\link{diana}}, \code{\link{diana.object}},
  \code{\link{twins.object}}, \code{\link{par}}. 
}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
