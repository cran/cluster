% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{clara}
\alias{clara}
\title{
Clustering Large Applications
}
\description{
Returns a list representing a clustering of the data into \code{k} clusters.
}
\usage{
clara(x, k, metric = "euclidean", stand = F, samples = 5, 
      sampsize = 40 + 2 * k)


}
\arguments{
\item{x}{
data matrix or dataframe, each row corresponds to an observation,
and each column corresponds to a variable. All variables must be numeric.
Missing values (NAs) are allowed.


}
\item{k}{
integer, the number of clusters.
It is required that 0 < k < n where n is the number of observations. 


}
\item{metric}{
character string specifying the metric to be used for calculating
dissimilarities between observations.
The currently available options are "euclidean" and "manhattan".
Euclidean distances are root sum-of-squares of differences, and
manhattan distances are the sum of absolute differences.


}
\item{stand}{
logical flag: if TRUE, then the measurements in \code{x} are standardized before
calculating the dissimilarities. Measurements are standardized for each
variable (column), by subtracting the variable's mean value and dividing by 
the variable's mean absolute deviation.


}
\item{samples}{
integer, number of samples to be drawn from the dataset.


}
\item{sampsize}{
integer, number of observations in each sample. \code{sampsize} should be higher 
than the number of clusters (\code{k}) and at most the number of observations 
(nrow(\code{x})).


}}
\value{
an object of class \code{"clara"} representing the clustering.
See clara.object for details.


}
\details{
\code{clara} is fully described in chapter 3 of Kaufman and Rousseeuw (1990).
Compared to other partitioning methods such as \code{pam}, it can deal with
much larger datasets. Internally, this is achieved by considering
sub-datasets of fixed size, so that the time and storage requirements
become linear in nrow(\code{x}) rather than quadratic.


Each sub-dataset is partitioned into \code{k} clusters using the same 
algorithm as in the \code{pam} function.
Once \code{k} representative objects have been selected from the
sub-dataset, each observation of the entire dataset is assigned
to the nearest medoid.
The sum of the dissimilarities of the observations to their closest medoid, is
used as a measure of the quality of the clustering. The sub-dataset
for which the sum is minimal, is retained.
A further analysis is carried out on the final partition.
Each sub-dataset is forced to contain the medoids obtained from the best 
sub-dataset until then.
Randomly drawn observations are added to this set until \code{sampsize}
has been reached.


}
\section{BACKGROUND}{
Cluster analysis divides a dataset into groups (clusters) of observations that
are similar to each other. Partitioning methods like \code{pam}, \code{clara}, and
\code{fanny} require that the number of clusters be given by the user.
Hierarchical methods like \code{agnes}, \code{diana}, and \code{mona} construct a
hierarchy of clusterings, with the number of clusters ranging from one to
the number of observations.


}
\section{NOTE}{
For small datasets (say with fewer than 200 observations),
the function \code{pam} can be used directly.


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
\code{\link{clara.object}}, \code{\link{pam}}, \code{\link{partition.object}}, \code{\link{plot.partition}}.
}
\examples{
## generate 500 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(200,0,8), rnorm(200,0,8)),
           cbind(rnorm(300,50,8), rnorm(300,50,8)))
clarax <- clara(x, 2)
clarax
clarax$clusinfo
plot(clarax)

## `xclara' is an artificial data set with 3 clusters of 1000 bivariate
## objects each.
data(xclara)
## Plot similar to Figure 5 in Struyf et al (1996)
\dontrun{plot(clara(xclara, 3), ask = TRUE)}
\testonly{plot(clara(xclara, 3))}
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
