\title{Print and Summary Methods for Dissimilarity Objects}
\name{print.dissimilarity}
\alias{print.dissimilarity}
\alias{summary.dissimilarity}
\alias{print.summary.dissimilarity}
\description{
  Print or summarize the distances and the attributes of a
  \code{dissimilarity} object.

  These are methods for the functions \code{print()} and \code{summary()} for
  \code{dissimilarity} objects.  See \code{print}, \code{print.default},
  or \code{summary} for the general behavior of these.
}
\usage{
  \method{print}{dissimilarity}(x, \dots)
\method{summary}{dissimilarity}(object, \dots)
\method{print}{summary.dissimilarity}(x, \dots)
}
\arguments{
  \item{x, object}{a \code{dissimilarity} object or a
    \code{summary.dissimilarity} one for \code{print.summary.dissimilarity()}.}
  \item{\dots}{potential further arguments (require by generic).}
}
\seealso{
  \code{\link{daisy}}, \code{\link{dissimilarity.object}},
  \code{\link{print}}, \code{\link{print.default}}.
}
\examples{
 ## See  example(daisy)

 sd <- summary(daisy(matrix(rnorm(100), 20,5)))
 sd # -> print.summary.dissimilarity(.)
 str(sd)
}
\keyword{cluster}
\keyword{print}


