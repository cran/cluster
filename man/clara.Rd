\name{clara}
\alias{clara}
\title{Clustering Large Applications}
\description{
  Computes a \code{"clara"} object, a list representing a clustering of
  the data into \code{k} clusters.
}
\usage{
clara(x, k, metric = "euclidean", stand = FALSE, samples = 5,
      sampsize = 40 + 2 * k)
}
\arguments{
  \item{x}{
    data matrix or data frame, each row corresponds to an observation,
    and each column corresponds to a variable.  All variables must be numeric.
    Missing values (NAs) are allowed. }
  \item{k}{integer, the number of clusters.
    It is required that \eqn{0 < k < n} where \eqn{n} is the number of
    observations (i.e., n = \code{nrow(x)}).}
  \item{metric}{
    character string specifying the metric to be used for calculating
    dissimilarities between observations.
    The currently available options are "euclidean" and "manhattan".
    Euclidean distances are root sum-of-squares of differences, and
    manhattan distances are the sum of absolute differences.
  }
  \item{stand}{logical, indicating if the measurements in \code{x} are
    standardized before calculating the dissimilarities.  Measurements
    are standardized for each variable (column), by subtracting the
    variable's mean value and dividing by the variable's mean absolute
    deviation.
  }
  \item{samples}{integer, number of samples to be drawn from the dataset.}
  \item{sampsize}{integer, number of observations in each
    sample. \code{sampsize} should be higher than the number of clusters
    (\code{k}) and at most the number of observations (n = \code{nrow(x)}).}
}
\value{
  an object of class \code{"clara"} representing the clustering.  See
  \code{\link{clara.object}} for details.
}
\details{
  \code{clara} is fully described in chapter 3 of Kaufman and Rousseeuw (1990).
  Compared to other partitioning methods such as \code{pam}, it can deal with
  much larger datasets.  Internally, this is achieved by considering
  sub-datasets of fixed size (\code{sampsize}) such that the time and
  storage requirements become linear in \eqn{n} rather than quadratic.

  Each sub-dataset is partitioned into \code{k} clusters using the same
  algorithm as in \code{\link{pam}}.\cr
  Once \code{k} representative objects have been selected from the
  sub-dataset, each observation of the entire dataset is assigned
  to the nearest medoid.

  The sum of the dissimilarities of the observations to their closest
  medoid is used as a measure of the quality of the clustering.  The
  sub-dataset for which the sum is minimal, is retained.  A further
  analysis is carried out on the final partition.

  Each sub-dataset is forced to contain the medoids obtained from the
  best sub-dataset until then.  Randomly drawn observations are added to
  this set until \code{sampsize} has been reached.
}
\note{
%% mostly by Martin Mächler :
  The random sampling is implemented with a \emph{very} simple scheme
  (with period \eqn{2^{16} = 65536}) inside the Fortran code,
  independently of \R's random number generation, and as a matter of
  fact, deterministically.

  The storage requirement of \code{clara} computation (for small
  \code{k}) is about
  \eqn{O(n \times p) + O(j^2)}{O(n * p) + O(j^2)} where
  \eqn{j = \code{sampsize}}, and \eqn{(n,p) = \code{dim(x)}}.
  The CPU computing time (again neglecting small \code{k}) is about
  \eqn{O(n \times p \times j^2 \times N)}{O(n * p * j^2 * N)}, where
  \eqn{N = \code{samples}}.

  For ``small'' datasets, the function \code{\link{pam}} can be used
  directly.  What can be considered \emph{small}, is really a function
  of available computing power, both memory (RAM) and speed.
  Originally (1990), ``small'' meant less than 100 observations;
  later, the authors said \emph{``small (say with fewer than 200
  observations)''}..
}
\seealso{
  \code{\link{agnes}} for background and references;
  \code{\link{clara.object}}, \code{\link{pam}},
  \code{\link{partition.object}}, \code{\link{plot.partition}}.
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
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
