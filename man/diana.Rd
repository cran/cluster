\name{diana}
\alias{diana}
\title{DIvisive ANAlysis Clustering}
\description{
  Computes a divisive hierarchical clustering of the dataset
  returning an object of class \code{diana}.
}
\usage{
diana(x, diss = inherits(x, "dist"), metric = "euclidean", stand = FALSE)
}
\arguments{
  \item{x}{
    data matrix or data frame, or dissimilarity matrix or object,
    depending on the value of the \code{diss} argument.

    In case of a matrix or data frame, each row corresponds to an observation,
    and each column corresponds to a variable.  All variables must be numeric.
    Missing values (\code{\link{NA}}s) \emph{are} allowed.

    In case of a dissimilarity matrix, \code{x} is typically the output
    of \code{\link{daisy}} or \code{\link{dist}}.  Also a vector of
    length n*(n-1)/2 is allowed (where n is the number of observations),
    and will be interpreted in the same way as the output of the
    above-mentioned functions. Missing values (NAs) are \emph{not} allowed.
  }
  \item{diss}{
    logical flag: if TRUE (default for \code{dist} or
    \code{dissimilarity} objects), then \code{x} will be considered as a
    dissimilarity matrix.  If FALSE, then \code{x} will be considered as
    a matrix of observations by variables.
  }
  \item{metric}{
    character string specifying the metric to be used for calculating
    dissimilarities between observations.\cr
    The currently available options are "euclidean" and
    "manhattan".  Euclidean distances are root sum-of-squares of
    differences, and manhattan distances are the sum of absolute
    differences.  If \code{x} is already a dissimilarity matrix, then
    this argument will be ignored.
  }
  \item{stand}{logical; if true, the measurements in \code{x} are
    standardized before calculating the dissimilarities.  Measurements
    are standardized for each variable (column), by subtracting the
    variable's mean value and dividing by the variable's mean absolute
    deviation.  If \code{x} is already a dissimilarity matrix, then this
    argument will be ignored.}
}
\value{
  an object of class \code{"diana"} representing the clustering.  See
  \code{?\link{diana.object}} for details.
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
contains only a single observation.\cr
At each stage, the cluster with the largest diameter is selected.
(The diameter of a cluster is the largest dissimilarity between any
two of its observations.)\cr
To divide the selected cluster, the algorithm first looks for its most
disparate observation (i.e., which has the largest average dissimilarity to the
other observations of the selected cluster). This observation initiates the
"splinter group". In subsequent steps, the algorithm reassigns observations
that are closer to the "splinter group" than to the "old party". The result
is a division of the selected cluster into two new clusters.
}
\seealso{
  \code{\link{agnes}} also for background and references;
  \code{\link{diana.object}}, \code{\link{daisy}}, \code{\link{dist}},
  \code{\link{plot.diana}}, \code{\link{twins.object}}.
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
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
