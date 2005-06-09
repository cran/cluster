\name{fanny}
\alias{fanny}
\title{Fuzzy Analysis Clustering}
\description{
  Computes a fuzzy clustering of the data into \code{k} clusters.
}
\usage{
fanny(x, k, diss = inherits(x, "dist"),
      memb.exp = 2, metric = "euclidean", stand = FALSE,
      maxit = 500, tol = 1e-15)
}
\arguments{
  \item{x}{
    data matrix or data frame, or dissimilarity matrix, depending on the
    value of the \code{diss} argument.

    In case of a matrix or data frame, each row corresponds to an observation,
    and each column corresponds to a variable. All variables must be numeric.
    Missing values (NAs) are allowed.

    In case of a dissimilarity matrix, \code{x} is typically the output
    of \code{\link{daisy}} or \code{\link{dist}}.  Also a vector of
    length n*(n-1)/2 is allowed (where n is the number of observations),
    and will be interpreted in the same way as the output of the
    above-mentioned functions.  Missing values (NAs) are not allowed.
  }
  \item{k}{integer giving the desired number of clusters.  It is
    required that \eqn{0 < k < n/2} where \eqn{n} is the number of
    observations.}
  \item{diss}{
    logical flag: if TRUE (default for \code{dist} or
    \code{dissimilarity} objects), then \code{x} is assumed to be a
    dissimilarity matrix.  If FALSE, then \code{x} is treated as
    a matrix of observations by variables.
  }
  \item{memb.exp}{number \eqn{r} strictly larger than 1 specifying the
    \emph{membership exponent} used in the fit criterion; see the
    \sQuote{Details} below. Default: \code{2} which used to be hardwired
    inside FANNY.}
  \item{metric}{
    character string specifying the metric to be used for calculating
    dissimilarities between observations.
    The currently available options are "euclidean" and "manhattan".
    Euclidean distances are root sum-of-squares of differences, and
    manhattan distances are the sum of absolute differences.
    If \code{x} is already a dissimilarity matrix, then this argument will
    be ignored.
  }
  \item{stand}{logical; if true, the measurements in \code{x} are
    standardized before calculating the dissimilarities.  Measurements
    are standardized for each variable (column), by subtracting the
    variable's mean value and dividing by the variable's mean absolute
    deviation.  If \code{x} is already a dissimilarity matrix, then this
    argument will be ignored.}
  \item{maxit, tol}{maximal number of iterations and default tolerance
    for convergence (relative convergence of the fit criterion) for the
    FANNY algorithm.  The defaults \code{maxit = 500} and \code{tol =
      1e-15} used to be hardwired inside the algorithm.}
}
\value{
  an object of class \code{"fanny"} representing the clustering.
  See \code{\link{fanny.object}} for details.
}
\details{
  In a fuzzy clustering, each observation is ``spread out'' over the various
  clusters. Denote by u(i,v) the membership of observation i to cluster v.
  The memberships are nonnegative, and for a fixed observation i they sum to 1.
  The particular method \code{fanny} stems from chapter 4 of
  Kaufman and Rousseeuw (1990) (see the references in
  \code{\link{daisy}}) and has been extended to allow user specified
  \code{memb.exp}.

  Fanny aims to minimize the objective function
  \deqn{\sum_{v=1}^k
    \frac{\sum_{i=1}^n\sum_{j=1}^n u_{iv}^r u_{jv}^r d(i,j)}{
      2 \sum_{j=1}^n u_{jv}^r}}{%
    SUM_[v=1..k] (SUM_(i,j) u(i,v)^r u(j,v)^r d(i,j)) / (2 SUM_j u(j,v)^r)}
  where \eqn{n} is the number of observations, \eqn{k} is the number of
  clusters, \eqn{r} is the membership exponent \code{memb.exp} and
  \eqn{d(i,j)} is the dissimilarity between observations \eqn{i} and \eqn{j}.
  \cr Note that \eqn{r \to 1}{r -> 1} gives increasingly crisper
  clusterings whereas \eqn{r \to \infty}{r -> Inf} leads to complete
  fuzzyness.  K\&R(1990), p.191 note that values too close to 1 can lead
  to slow convergence.

  Compared to other fuzzy clustering methods, \code{fanny} has the following
  features: (a) it also accepts a dissimilarity matrix; (b) it is
  more robust to the \code{spherical cluster} assumption; (c) it provides
  a novel graphical display, the silhouette plot (see
  \code{\link{plot.partition}}).
}
\seealso{
  \code{\link{agnes}} for background and references;
  \code{\link{fanny.object}}, \code{\link{partition.object}},
  \code{\link{plot.partition}}, \code{\link{daisy}}, \code{\link{dist}}.
}
\examples{
## generate 10+15 objects in two clusters, plus 3 objects lying
## between those clusters.
x <- rbind(cbind(rnorm(10, 0, 0.5), rnorm(10, 0, 0.5)),
           cbind(rnorm(15, 5, 0.5), rnorm(15, 5, 0.5)),
           cbind(rnorm( 3,3.2,0.5), rnorm( 3,3.2,0.5)))
fannyx <- fanny(x, 2)
## Note that observations 26:28 are "fuzzy" (closer to # 2):
fannyx
summary(fannyx)
plot(fannyx)

(fan.x.15 <- fanny(x, 2, memb.exp = 1.5)) # 'crispier' for obs. 26:28
(fanny(x, 2, memb.exp = 3))               # more fuzzy in general

data(ruspini)
## Plot similar to Figure 6 in Stryuf et al (1996)
plot(fanny(ruspini, 5))
}
\keyword{cluster}
