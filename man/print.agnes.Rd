\name{print.agnes}
\alias{print.agnes}
\title{Print Method for AGNES Objects}
\description{
  Prints the ordering of objects, distances between merging clusters,
  agglomerative coefficient, and merge structure of an \code{agnes} object.

  This is a method for the function \code{\link{print}()} for objects
  inheriting from class \code{\link{agnes}}.
}
\usage{
\method{print}{agnes}(x, \dots)
}
\arguments{
  \item{x}{an agnes object.}
  \item{\dots}{potential further arguments (require by generic).}
}
\seealso{
  \code{\link{agnes}}, \code{\link{agnes.object}}, \code{\link{print}},
  \code{\link{print.default}}.
}
\keyword{cluster}
\keyword{print}
