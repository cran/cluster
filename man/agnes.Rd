% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{agnes}
\alias{agnes}
\title{
Agglomerative Nesting
}
\description{
Returns a list representing an agglomerative hierarchical clustering of
the dataset.
}
\usage{
agnes(x, diss = F, metric = "euclidean", stand = F, method = "average")


}
\arguments{
\item{x}{
data matrix or dataframe, or dissimilarity matrix, depending on the
value of the \code{diss} argument.


In case of a matrix or dataframe, each row corresponds to an observation,
and each column corresponds to a variable. All variables must be numeric.
Missing values (NAs) are allowed.


In case of a dissimilarity matrix, \code{x} is typically the output of \code{daisy} or
\code{dist}. Also a vector with length n*(n-1)/2 is allowed (where n is the
number of observations), and will be interpreted in the same way as the output
of the above-mentioned functions. Missing values (NAs) are not allowed.


}
\item{diss}{
logical flag: if TRUE, then \code{x} will be considered as a dissimilarity
matrix. If FALSE, then \code{x} will be considered as a matrix of observations
by variables.


}
\item{metric}{
character string specifying the metric to be used for calculating
dissimilarities between observations.
The currently available options are "euclidean" and "manhattan".
Euclidean distances are root sum-of-squares of differences, and
manhattan distances are the sum of absolute differences.
If \code{x} is already a dissimilarity matrix, then this argument will
be ignored.


}
\item{stand}{
logical flag: if TRUE, then the measurements in \code{x} are standardized before
calculating the dissimilarities. Measurements are standardized for each
variable (column), by subtracting the variable's mean value and dividing by 
the variable's mean absolute deviation.
If \code{x} is already a dissimilarity matrix, then this argument
will be ignored.


}
\item{method}{
character string defining the clustering method. The five methods
implemented are "average" (group average method), 
"single" (single linkage), "complete" (complete linkage), 
"ward" (Ward's method), and "weighted" (weighted average linkage).
Default is "average".


}}
\value{
an object of class \code{"agnes"} representing the clustering.
See agnes.object for details.


}
\details{
\code{agnes} is fully described in chapter 5 of Kaufman and Rousseeuw (1990).
Compared to other agglomerative clustering methods such as \code{hclust}, 
\code{agnes} has the following features: (a) it yields the
agglomerative coefficient (see \code{agnes.object})
which measures the amount of clustering structure found; and (b)
apart from the usual tree it also provides the banner, a novel
graphical display (see \code{plot.agnes}).


The \code{agnes}-algorithm constructs a hierarchy of clusterings.
At first, each observation
is a small cluster by itself. Clusters are merged until only one large
cluster remains which contains all the observations.
At each stage the two "nearest" clusters are combined to form one larger
cluster. For \code{method}="average", the distance between two clusters is the
average of the dissimilarities between the points in one cluster and the
points in the other cluster. In \code{method}="single", we use
the smallest dissimilarity between a point in the first cluster
and a point in the second cluster (nearest neighbor method). 
When \code{method}="complete", we use the
largest dissimilarity between a point in the first cluster and a point
in the second cluster (furthest neighbor method).


}
\section{BACKGROUND}{
Cluster analysis divides a dataset into groups (clusters) of observations that
are similar to each other. Hierarchical methods like \code{agnes}, \code{diana}, and
\code{mona} construct a hierarchy of clusterings, with the number of clusters
ranging from one to the number of observations. Partitioning methods like 
\code{pam}, \code{clara}, and \code{fanny} require that the number of clusters be given by
the user.


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
\code{\link{agnes.object}}, \code{\link{daisy}}, \code{\link{diana}}, \code{\link{dist}}, \code{\link{hclust}}, \code{\link{plot.agnes}}, 
\code{\link{twins.object}}.
}
\examples{
data(votes.repub)
agn1 <- agnes(votes.repub, metric = "manhattan", stand = TRUE)
print(agn1)
plot(agn1)
agn2 <- agnes(daisy(votes.repub), diss = TRUE, method = "complete")
plot(agn2)

data(agriculture)
## Plot similar to Figure 7 in ref
\dontrun{plot(agnes(agriculture), ask = TRUE)}
\testonly{plot(agnes(agriculture))}
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
