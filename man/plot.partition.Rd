% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{plot.partition}
\alias{plot.partition}
\title{
Plot of a Partition of the Data Set
}
\description{
Creates plots for visualizing a \code{partition} object.
}
\usage{
plot.partition(x, ask=F, \dots)
}
\arguments{
\item{x}{
an object of class \code{"partition"}, e.g. created by the functions \code{pam},
\code{clara}, and \code{fanny}.
}
\item{ask}{
if TRUE, \code{plot.partition} operates in interactive mode.


All optional arguments available to the function \code{clusplot.default}
(except for the \code{diss} option) may also be supplied to this function.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.
}}
\value{
a NULL value is returned.
}
\section{Side Effects}{
An appropriate plot is produced on the current graphics device. This can
be one or both of the following choices:
Clusplot
Silhouette plot
}
\details{
When \code{ask=T}, rather than producing each plot sequentially,
\code{plot.partition} displays a menu listing all the plots that can be produced.
If the menu is not desired but a pause between plots is still wanted
one must set \code{par(ask=T)} before invoking the plot command.


The clusplot of a cluster partition consists of a two-dimensional 
representation of the observations, in which the clusters are 
indicated by ellipses. (See clusplot.partition for more details.)


The silhouette plot of a nonhierarchical clustering is fully described in 
Rousseeuw (1987) and in chapter 2 of Kaufman and Rousseeuw (1990).
For each observation i, a bar is drawn, representing the silhouette width s(i)
of the observation. Observations are grouped per cluster, starting with 
cluster 1 at the top. Observations with a large s(i) (almost 1) are very well 
clustered, a small s(i) (around 0) means that the observation lies between 
two clusters, and observations with a negative s(i) are probably placed in 
the wrong cluster.
A clustering can be performed for several values of \code{k} (the number of
clusters). Finally, choose the value of \code{k} with the largest overall
average silhouette width.


The silhouette width is computed as follows:
Put a(i) = average dissimilarity between i and all other points of the
cluster to which i belongs. For all clusters C, put d(i,C) = average
dissimilarity of i to all observations of C. The smallest of these d(i,C) is
denoted as b(i), and can be seen as the dissimilarity between i and its
neighbor cluster. Finally, put s(i) = ( b(i) - a(i) ) / max( a(i), b(i) ).
The overall average silhouette width is then simply the average of s(i) over
all observations i.
}
\section{NOTE}{
In the silhouette plot,
observation labels are only printed when the number of observations is 
limited to less than 40, for readability.
Moreover, observation labels are truncated to at most 5 characters.


}
\references{
Kaufman, L. and Rousseeuw, P.J. (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis. }
Wiley, New York.


Rousseeuw, P.J. (1987). Silhouettes: A graphical aid to the interpretation
and validation of cluster analysis. 
\emph{J. Comput. Appl. Math., }
\bold{20}, 53-65.


Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust 
Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis, \bold{26}, 17-37.}
}
\seealso{
\code{\link{partition.object}}, \code{\link{clusplot.partition}}, \code{\link{clusplot.default}},
\code{\link{pam}}, \code{\link{pam.object}}, \code{\link{clara}}, \code{\link{clara.object}}, \code{\link{fanny}},
\code{\link{fanny.object}}, \code{\link{par}}.
}
\examples{
## generate 25 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
plot(pam(x, 2))
}
\keyword{libcluster}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
