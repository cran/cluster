% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{clusplot.partition}
\alias{clusplot.partition}
\title{
Bivariate clusplot of a partitioning object
}
\description{
Creates the clusplot of a \code{partition} object.
}
\usage{
clusplot.partition(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{"partition"}, e.g. created by the functions \code{pam},
\code{clara}, or \code{fanny}.


All optional arguments available to the function clusplot.default (except for the \code{diss} option) may also be supplied to this function.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.
}
\value{
An invisible list with components:
}
\item{Distances}{
When option lines is 1 or 2 we optain a k by k matrix (k is the number of 
clusters). The element at row j and column s is the distance between ellipse 
j and ellipse s.
If lines=0, then the value of this component is NA.


}
\item{Shading}{
A vector of length k (where k is the number of clusters), containing the 
amount of shading per cluster. Let y be a vector where element i is the ratio 
between the number of objects in cluster i and the area of ellipse i. When 
the cluster i is a line segment, y[i] and the density of the cluster are set 
to NA. Let z be the sum of all the elements of y without the NAs. Then we 
put shading = y/z *37 + 3 
}}
\section{Side Effects}{
a visual display of  the  clustering  is  plotted  on  the current graphics
device.
}
\details{
clusplot uses the functions princomp and  cmdscale.  These
functions   are   data  reduction  techniques.  They  will
represent the data in a bivariate plot.  Ellipses  are
then   drawn  to indicate the clusters. The further layout
of the plot is determined by the optional arguments.


If the clustering algorithms \code{pam}, \code{fanny} and \code{clara} are applied to a data 
matrix of observations-by-variables then a clusplot of the resulting 
clustering can always be drawn.
When the data matrix contains missing values and the clustering is performed 
with \code{pam} or \code{fanny}, the dissimilarity matrix will be given as input to 
\code{clusplot}. When the clustering algorithm \code{clara} was applied to a
data matrix with NAs then clusplot will replace the missing values as 
described in \code{clusplot.default}, because a dissimilarity matrix is not 
available. 
}
\references{
Kaufman, L. and Rousseeuw, P.J. (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis. }
Wiley, New York.


Pison, G., Struyf, A. and Rousseeuw, P.J. (1997). 
Displaying a Clustering with CLUSPLOT,
\emph{Technical Report, University of Antwerp, submitted.}


Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). 
Integrating Robust Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis, \bold{26}, 17-37.}
}
\seealso{
\code{\link{partition.object}}, \code{\link{pam}}, \code{\link{pam.object}}, \code{\link{clara}}, \code{\link{clara.object}}, \code{\link{fanny}},
\code{\link{fanny.object}}, \code{\link{par}}, \code{\link{clusplot.default}}.
}
\examples{
## generate 25 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
clusplot(pam(x, 2))
}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
