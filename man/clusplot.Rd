% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{clusplot}
\alias{clusplot}
\title{
Clusplot - Generic Function
}
\description{
Creates a  2-dimensional clusplot on the current graphics device.
This is a generic function.
Functions with names beginning in \code{clusplot.} will be
methods for this function.
\code{partition}.
}
\usage{
clusplot(x, \dots)
}
\arguments{
\item{x}{
an S-PLUS object.
}
\item{...}{
methods may have additional arguments.
Graphical parameters (see \code{\link{par}}) may also
be supplied as arguments to this function.
}}
\section{Side Effects}{
a 2-dimensional clusplot is created on the current graphics device.
}
\seealso{
\code{\link{clusplot.default}}, \code{\link{clusplot.partition}}, \code{\link{partition.object}}, \code{\link{pam}}, \code{\link{fanny}}, \code{\link{clara}}.


}
\keyword{libcluster}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
