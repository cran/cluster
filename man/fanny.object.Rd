\name{fanny.object}
\alias{fanny.object}
\title{Fuzzy Analysis (FANNY) Object}
\description{
  The objects of class \code{"fanny"} represent a fuzzy clustering of a
  dataset.
}
\section{GENERATION}{
  These objects are returned from \code{\link{fanny}}.
}
\section{METHODS}{
  The \code{"fanny"} class has methods for the following generic functions:
  \code{print}, \code{summary}.
}
\section{INHERITANCE}{
  The class \code{"fanny"} inherits from \code{"partition"}.
  Therefore, the generic functions \code{plot} and \code{clusplot} can
  be used on a \code{fanny} object.
}
\value{
  A legitimate \code{fanny} object is a list with the following components:
  \item{membership}{
    matrix containing the memberships for each pair consisting of an
    observation and a cluster.
  }
  \item{coeff}{
    Dunn's partition coefficient \eqn{F(k)} of the clustering, where
    \eqn{k} is the number of clusters. \eqn{F(k)} is the sum of all
    squared membership coefficients, divided by the number of
    observations.  Its value is always between \eqn{1/k} and 1.

    The normalized form of the coefficient is also given. It is defined as
    \eqn{(F(k) - 1/k) / (1 - 1/k)}, and ranges between 0 and 1.
    A low value of Dunn's coefficient indicates a very fuzzy clustering,
    whereas a value close to 1 indicates a near-crisp clustering.
  }
  \item{clustering}{
    the clustering vector of the nearest crisp clustering, see
    \code{\link{partition.object}}.}
  \item{objective}{
    the objective function and the number of iterations the \code{fanny}
    algorithm needed to reach this minimal value.
  }
  \item{diss}{
    an object of class \code{"dissimilarity"}, see
    \code{\link{partition.object}}.}
  \item{call}{generating call, see \code{\link{partition.object}}.}
  \item{silinfo}{
    list with silhouette information of the nearest crisp clustering, see
    \code{\link{partition.object}}.}
  \item{data}{matrix, possibibly standardized, or NULL, see
    \code{\link{partition.object}}.}
}
\seealso{
  \code{\link{fanny}}, \code{\link{dissimilarity.object}},
  \code{\link{partition.object}}, \code{\link{plot.partition}}.
}
\keyword{cluster}

