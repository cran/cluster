% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{fanny}
\alias{fanny}
\title{
Fuzzy Analysis
}
\description{
Returns a list representing a fuzzy clustering of the data into \code{k}
clusters.
}
\usage{
fanny(x, k, diss = F, metric = "euclidean", stand = F)
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
\item{k}{
integer, the number of clusters.
It is required that 0 < k < n/2 where n is the number of observations. 


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


}}
\value{
an object of class \code{"fanny"} representing the clustering.
See \code{fanny.object} for details.
}
\details{
In a fuzzy clustering, each observation is "spread out" over the various
clusters. Denote by u(i,v) the membership of observation i to cluster v.
The memberships are nonnegative, and for a fixed observation i they sum to 1.
The particular method \code{fanny} stems from chapter 4 of
Kaufman and Rousseeuw (1990).
Compared to other fuzzy clustering methods, \code{fanny} has the following
features: (a) it also accepts a dissimilarity matrix; (b) it is
more robust to the \code{spherical cluster} assumption; (c) it provides
a novel graphical display, the silhouette plot (see \code{plot.partition}).


Fanny aims to minimize the objective function
  \deqn{\sum_{v=1}^k
    \frac{\sum_{i=1}^n\sum_{j=1}^n u_{iv}^2 u_{jv}^2 d(i,j)}{
      2 \sum_{j=1}^n u_{jv}^2}}{
    SUM_v (SUM_(i,j) u(i,v)^2 u(j,v)^2 d(i,j)) / (2 SUM_j u(j,v)^2)}
where n is the number of observations, k is the number of clusters and
d(i,j) is the dissimilarity between observations i and j.


}
\section{BACKGROUND}{
Cluster analysis divides a dataset into groups (clusters) of observations that
are similar to each other. Partitioning methods like \code{pam}, \code{clara}, and
\code{fanny} require that the number of clusters be given by the user.
Hierarchical methods like \code{agnes}, \code{diana}, and \code{mona} construct a
hierarchy of clusterings, with the number of clusters ranging from one to
the number of observations.


}
\references{
    Kaufman, L. and Rousseeuw, P.J. (1990). 
    \emph{Finding Groups in Data: An Introduction to Cluster Analysis}.
    Wiley, New York.

    Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996):
    Clustering in an Object-Oriented Environment.
    \emph{Journal of Statistical Software}, \bold{1}.
    \url{http://www.stat.ucla.edu/journals/jss/}
    
    Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating
    Robust Clustering Techniques in S-PLUS,
    \emph{Computational Statistics and Data Analysis}, \bold{26}, 17-37.
}
\seealso{
\code{\link{fanny.object}}, \code{\link{daisy}}, \code{\link{partition.object}}, \code{\link{plot.partition}}, \code{\link{dist}}.
}
\examples{
## generate 25 objects, divided into two clusters, and 3 objects lying
## between those clusters. 
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)),
           cbind(rnorm(3,3.5,0.5), rnorm(3,3.5,0.5)))
fannyx <- fanny(x, 2)
fannyx
summary(fannyx)
plot(fannyx)

data(ruspini)
## Plot similar to Figure 6 in Stryuf et al (1996)
plot(fanny(ruspini, 5))
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
