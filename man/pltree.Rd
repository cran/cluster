% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{pltree}
\alias{pltree}
\title{
Clustering Trees - Generic Function
}
\description{
Creates a clustering tree on the current graphics device.
This is a generic function.
Functions with names beginning in \code{pltree.} will be
methods for this function.
\code{twins}.
}
\usage{
pltree(x, \dots)
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
a clustering tree is created on the current graphics device.
}
\seealso{
\code{\link{twins.object}}, \code{\link{agnes}}, \code{\link{diana}}, \code{\link{par}}.
}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
