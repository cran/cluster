% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{diana}
\alias{diana}
\title{
Divisive Analysis
}
\description{
Returns a list representing a divisive hierarchical clustering of
the dataset.
}
\usage{
diana(x, diss = F, metric = "euclidean", stand = F)


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
matrix. If FALSE, then \code{x} will be considered as a matrix of
observations by variables.


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


}}
\value{
an object of class \code{"diana"} representing the clustering.
See diana.object for details.


}
\details{
\code{diana} is fully described in chapter 6 of Kaufman and Rousseeuw (1990).
It is probably unique in computing a divisive hierarchy, whereas most
other software for hierarchical clustering is agglomerative.
Moreover, \code{diana} provides (a) the divisive coefficient
(see \code{diana.object}) which measures the amount of clustering structure
found; and (b) the banner, a novel graphical display
(see \code{plot.diana}).


The \code{diana}-algorithm constructs a hierarchy of clusterings,
starting with one large
cluster containing all n observations. Clusters are divided until each cluster
contains only a single observation.
At each stage, the cluster with the largest diameter is selected.
(The diameter of a cluster is the largest dissimilarity between any
two of its observations.)
To divide the selected cluster, the algorithm first looks for its most
disparate observation (i.e., which has the largest average dissimilarity to the
other observations of the selected cluster). This observation initiates the
"splinter group". In subsequent steps, the algorithm reassigns observations
that are closer to the "splinter group" than to the "old party". The result
is a division of the selected cluster into two new clusters.


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
\code{\link{agnes}}, \code{\link{diana.object}}, \code{\link{daisy}}, \code{\link{dist}}, \code{\link{plot.diana}}, \code{\link{twins.object}}.
}
\examples{
data(votes.repub)
dv <- diana(votes.repub, metric = "manhattan", stand = TRUE)
print(dv)
plot(dv)

data(agriculture)
## Plot similar to Figure 8 in ref
\dontrun{plot(diana(agriculture), ask = TRUE)}
\testonly{plot(diana(agriculture))}
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
