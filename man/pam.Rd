%-*- mode: Rd; kept-old-versions: 12;  kept-new-versions: 20; -*-
\name{pam}
\alias{pam}
\title{Partitioning Around Medoids}
\description{
  Partitioning (clustering) of the data into \code{k} clusters ``around
  medoids'', a more robust version of K-means.
}
\usage{
pam(x, k, diss = inherits(x, "dist"), metric = "euclidean", stand = FALSE)
}
\arguments{
  \item{x}{
    data matrix or data frame, or dissimilarity matrix or object,
    depending on the value of the \code{diss} argument.

    In case of a matrix or data frame, each row corresponds to an
    observation, and each column corresponds to a variable.  All
    variables must be numeric.  Missing values (\code{\link{NA}}s)
    \emph{are} allowed---as long as every pair of observations has at
    least one case not missing.

    In case of a dissimilarity matrix, \code{x} is typically the output
    of \code{\link{daisy}} or \code{\link{dist}}.  Also a vector of
    length n*(n-1)/2 is allowed (where n is the number of observations),
    and will be interpreted in the same way as the output of the
    above-mentioned functions. Missing values (NAs) are \emph{not}
    allowed.
  }
  \item{k}{positive integer specifying the number of clusters, less than
    the number of observations.}
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
  an object of class \code{"pam"} representing the clustering.  See
  \code{?\link{pam.object}} for details.
}
\details{
  \code{pam} is fully described in chapter 2 of Kaufman and Rousseeuw
  (1990).  Compared to the k-means approach in \code{kmeans}, the
  function \code{pam} has the following features: (a) it also accepts a
  dissimilarity matrix; (b) it is more robust because it minimizes a sum
  of dissimilarities instead of a sum of squared euclidean distances;
  (c) it provides a novel graphical display, the silhouette plot (see
  \code{plot.partition}) (d) it allows to select the number of clusters
  using \code{\link{silhouette}()$avg.width} on the result (or its
  component \code{silinfo}, see also \code{\link{pam.object}}).

  The \code{pam}-algorithm is based on the search for \code{k}
  representative objects or medoids among the observations of the
  dataset.  These observations should represent the structure of the
  data.  After finding a set of \code{k} medoids, \code{k} clusters are
  constructed by assigning each observation to the nearest medoid.  The
  goal is to find \code{k} representative objects which minimize the sum
  of the dissimilarities of the observations to their closest
  representative object.
  \cr
  The algorithm first looks for a good initial set of medoids (this is
  called the \bold{build} phase).  Then it finds a local minimum for the
  objective function, that is, a solution such that there is no single
  switch of an observation with a medoid that will decrease the
  objective (this is called the \bold{swap} phase).
}
\note{
  For datasets larger than (say) 200 observations, \code{pam} will take a lot of
  computation time.  Then the function \code{\link{clara}} is preferable.
}
\seealso{
  \code{\link{agnes}} for background and references;
  \code{\link{pam.object}}, \code{\link{clara}}, \code{\link{daisy}},
  \code{\link{partition.object}}, \code{\link{plot.partition}},
  \code{\link{dist}}.
}
\examples{
## generate 25 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
pamx <- pam(x, 2)
pamx
summary(pamx)
plot(pamx)

pam(daisy(x, metric = "manhattan"), 2, diss = TRUE)

data(ruspini)
## Plot similar to Figure 4 in Stryuf et al (1996)
\dontrun{plot(pam(ruspini, 4), ask = TRUE)}
\testonly{plot(pam(ruspini, 4))}
}
\keyword{cluster}
