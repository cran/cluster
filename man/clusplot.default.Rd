\name{clusplot.default}
\alias{clusplot.default}
\title{Bivariate Cluster Plot (clusplot)}
\description{
  Creates a bivariate plot visualizing a partition (clustering) of the data. All
  observation are represented by points in the plot, using principal
  components or multidimensional scaling. Around each cluster an ellipse
  is drawn. 
}
\usage{
clusplot.default(x, clus, diss = FALSE, cor = TRUE, stand = FALSE,
                 lines = 2, shade = FALSE, color = FALSE,
                 labels= 0, plotchar = TRUE,
                 col.p = "dark green", col.txt = col.p,
                 span = TRUE, xlim = NULL, ylim = NULL, \dots)
}
\arguments{
  \item{x}{matrix or dataframe, or dissimilarity matrix, depending on
    the value of the \code{diss} argument.

    In case of a matrix (alike), each row corresponds to an observation,
    and each column corresponds to a variable.  All variables must be
    numeric. Missing values (\code{\link{NA}}s) are allowed.  They are
    replaced by the median of the corresponding variable.  When some
    variables or some observations contain only missing values, the
    function stops with a warning message.  

    In case of a dissimilarity matrix, \code{x} is the output of
    \code{\link{daisy}} or \code{\link{dist}} or a symmetric matrix.  Also,
    a vector of length \eqn{n*(n-1)/2} is allowed (where \eqn{n} is the
    number of observations), and will be interpreted in the same way as
    the output of the above-mentioned functions.  Missing values (NAs)
    are not allowed.
  }
  \item{clus}{
    a vector of length n representing a clustering of \code{x}.  For
    each observation the vector lists the number or name of the cluster
    to which it has been assigned. \code{clus} is often the clustering
    component of the output of \code{\link{pam}}, \code{\link{fanny}} or
    \code{\link{clara}}.}
  \item{diss}{
    logical indicating if \code{x} will be considered as a dissimilarity
    matrix or a matrix of observations by variables (see \code{x}
    arugment above).}
  \item{cor}{
    logical flag, only used when working with a data matrix (\code{diss
      = FALSE}). If TRUE, then the variables are scaled to unit variance.}
  \item{stand}{
    logical flag: if true, then the representations of the n observations in the
    2-dimensional plot are standardized.
  }
  \item{lines}{
    integer out of \code{0, 1, 2}, used to obtain an idea of the
    distances between ellipses.  The distance between two ellipses E1
    and E2 is measured along the line connecting the centers \eqn{m1}
    and \eqn{m2} of the two ellipses. 

    In case E1 and E2 overlap on the line through \eqn{m1} and \eqn{m2},
    no line is drawn.  Otherwise, the result depends on the value of 
    \code{lines}: If
    \describe{
      \item{lines = 0,}{no distance lines will appear on the plot;}
      \item{lines = 1,}{the line segment between \eqn{m1} and \eqn{m2} is drawn;}
      \item{lines = 2,}{a line segment between the boundaries of E1 and
	E2 is drawn (along the line connecting \eqn{m1} and \eqn{m2}).}
    }
  }
  \item{shade}{
    logical flag: if TRUE, then the ellipses are shaded in relation to their
    density. The density is the number of points in the cluster divided by the 
    area of the ellipse.
  }
  \item{color}{
    logical flag: if TRUE, then the ellipses are colored with respect to their
    density. With increasing density, the colors are light blue, light
    green, red and purple.  To see these colors on the graphics device, an
    appropriate color scheme should be selected (we recommend a white
    background).}
  \item{labels}{
    integer code, currently one of 0,1,2,3 and 4.  If
    \describe{
      \item{labels= 0,}{no labels are placed in the plot;}
      \item{labels= 1,}{points and ellipses can be identified in the plot (see
	\code{\link{identify}});}
      \item{labels= 2,}{all points and ellipses are labelled in the plot;}
      \item{labels= 3,}{only the points are labelled in the plot;}
      \item{labels= 4,}{only the ellipses are labelled in the plot.}
    }
    The levels of the vector \code{clus} are taken as labels for the
    clusters.  The labels 
    of the points are the rownames of \code{x} if \code{x} is matrix like.
    Otherwise (\code{diss = TRUE}), \code{x} is a vector, point labels
    can be attached to \code{x} as a "Labels" attribute
    (\code{attr(x,"Labels")}), as is done for the output of
    \code{\link{daisy}}.  
    
    A possible \code{\link{names}} attribute of \code{clus} will not
    be taken into account. 
  }
  \item{plotchar}{
    logical flag: if TRUE, then the plotting symbols differ for points belonging
    to different clusters.
  }
  \item{span}{
    logical flag: if TRUE, then each cluster is represented by the ellipse with
    smallest area containing all its points. (This is a special case of the
    minimum volume ellipsoid.)\cr
    If FALSE, the ellipse is based on the mean and covariance matrix of the 
    same points, often yielding a much larger ellipse.
    
    There are also some special cases:  When a cluster consists of only
    one point, a tiny circle is drawn around it.  When the points of a
    cluster fall on a straight line, \code{span=FALSE} draws a narrow
    ellipse around it and \code{span=TRUE} gives the exact line segment.
  }
  \item{col.p}{color code used for the observation points.}
  \item{col.txt}{color code for used for the labels.}

  \item{xlim, ylim}{length 2 vectors giving the x- and y- ranges as in
    \code{\link{plot.default}}.}
  \item{\dots}{Further graphical parameters may also be supplied, see
    \code{\link{par}}.} 
}% End Arguments

\value{
  An invisible list with components:
  \item{Distances}{
    When \code{lines} is 1 or 2 we optain a k by k matrix (k is the number of
    clusters).  The element in \code{[i,j]} is the distance between ellipse
    i and ellipse j.\cr
    If \code{lines = 0}, then the value of this component is NA.
  }
  \item{Shading}{
    A vector of length k (where k is the number of clusters), containing the
    amount of shading per cluster. Let y be a vector where element i is the
    ratio between the number of points in cluster i and the area of ellipse i. 
    When the cluster i is a line segment, y[i] and the density of the cluster are 
    set to  NA. Let z be the sum of all the elements of y without the NAs.
    Then we put shading = y/z *37 + 3 .
  }
}

\section{Side Effects}{
  a visual display of the clustering is plotted on the current graphics device.
}
\details{
  \code{clusplot} uses the functions \code{\link[mva]{princomp}} and
  \code{\link[mva]{cmdscale}}.  These functions are 
  data reduction techniques. They will represent the data in a bivariate plot.
  Ellipses are then drawn to indicate the clusters.  The further layout of the
  plot is determined by the optional arguments.
}
\note{
  When we have 4 or fewer clusters, then the \code{color=TRUE} gives
  every cluster a different color.  When there are more than 4 clusters,
  clusplot uses the function \code{\link[cluster]{pam}} to cluster the
  densities into 4 groups such that ellipses with nearly the same
  density get the same color.
}
\references{
  Kaufman, L. and Rousseeuw, P.J. (1990).
  \emph{Finding Groups in Data: An Introduction to Cluster Analysis.}
  Wiley, New York.

  Pison, G., Struyf, A. and Rousseeuw, P.J. (1997). 
  Displaying a Clustering with CLUSPLOT,
  \emph{Technical Report, University of Antwerp}, submitted.

  Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). 
  Integrating Robust Clustering Techniques in S-PLUS,
  \emph{Computational Statistics and Data Analysis}, \bold{26}, 17-37.
}
\seealso{
  \code{\link{princomp}}, \code{\link{cmdscale}}, \code{\link{pam}},
  \code{\link{clara}}, \code{\link{daisy}}, \code{\link{par}},
  \code{\link{identify}}, \code{\link{cov.mve}},
  \code{\link{clusplot.partition}}.
}
\examples{
## plotting votes.diss(dissimilarity) in a bivariate plot and
## partitioning into 2 clusters
data(votes.repub)
votes.diss <- daisy(votes.repub)
votes.clus <- pam(votes.diss, 2, diss = TRUE)$clustering
clusplot(votes.diss, votes.clus, diss = TRUE, shade = TRUE)

if(interactive()) #  uses identify() *interactively* :
  clusplot(votes.diss, votes.clus, diss = TRUE, shade = TRUE,
           labels = 1)

## plotting iris (dataframe) in a 2-dimensional plot and partitioning
## into 3 clusters.
data(iris)
iris.x <- iris[, 1:4]
clusplot(iris.x, pam(iris.x, 3)$clustering, diss = FALSE,
         plotchar = TRUE, color = TRUE)
}
\keyword{cluster}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
