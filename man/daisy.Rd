% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{daisy}
\alias{daisy}
\title{
Dissimilarity Matrix Calculation
}
\description{
Returns a matrix containing all the
pairwise dissimilarities (distances) between observations in the dataset.
The original variables may be of mixed types.


}
\usage{
daisy(x, metric = "euclidean", stand = F, type = list())


}
\arguments{
\item{x}{
data matrix or dataframe. Dissimilarities will be computed between the
rows of \code{x}.
Columns of class \code{numeric} will be recognized as interval scaled variables,
columns of class \code{factor} will be recognized as nominal variables,
and columns of class \code{ordered} will be recognized as ordinal variables.
Other variable types should be specified with the \code{type} argument.
Missing values (NAs) are allowed.


}
\item{metric}{
character string specifying the metric to be used.
The currently available options are "euclidean" and "manhattan".
Euclidean distances are root sum-of-squares of differences, and
manhattan distances are the sum of absolute differences.
If not all columns of \code{x} are numeric, then this argument
will be ignored.


}
\item{stand}{
logical flag: if TRUE, then the measurements in \code{x} are standardized before
calculating the dissimilarities. Measurements are standardized for each
variable (column), by subtracting the variable's mean value and dividing by 
the variable's mean absolute deviation.
If not all columns of \code{x} are numeric, then this argument
will be ignored.


}
\item{type}{
list containing some (or all) of the types of the variables
(columns) in \code{x}. The list may contain the following components: \code{ordratio}
(ratio scaled variables to be treated as ordinal variables), \code{logratio}
(ratio scaled variables that must be logarithmically transformed), \code{asymm}
(asymmetric binary variables). Each component's value is a vector,
containing the names or the numbers of the corresponding columns of \code{x}.
Variables not mentioned in the \code{type} list are interpreted as usual
(see argument \code{x}).


}}
\value{
an object of class \code{"dissimilarity"} containing the dissimilarities among
the rows of x. This is typically the input for the functions \code{pam},
\code{fanny}, \code{agnes} or \code{diana}. See dissimilarity.object for details.


}
\details{
\code{daisy} is fully described in chapter 1 of Kaufman and Rousseeuw (1990).
Compared to \code{dist} whose input must be numeric variables, the main
feature of \code{daisy} is its ability to handle other variable types as well
(e.g. nominal, ordinal, asymmetric binary) even when different types occur
in the same dataset.


In the \code{daisy} algorithm,
missing values in a row of x are not included in the
dissimilarities involving that row. If all variables are interval scaled,
the metric is "euclidean", and ng is the number of columns in which
neither row i and j have NAs, then the dissimilarity d(i,j) returned is
sqrt(ncol(x)/ng) times the Euclidean distance between the two vectors
of length ng shortened to exclude NAs. The rule is similar for the
"manhattan" metric, except that the coefficient is ncol(x)/ng.
If ng is zero, the dissimilarity is NA.


When some variables have a type other than interval scaled, the
dissimilarity between two rows is the weighted sum of the contribution of
each variable.
The weight becomes zero when that variable is missing in either or both
rows, or when the variable is asymmetric binary and both values are
zero. In all other situations, the weight of the variable is 1.
The contribution of nominal or binary variable a to the total dissimilarity
is zero if both values are different, else it
is equal to 1. The contribution of other variables is the absolute
difference of both values, divided by the total range of that variable.
Ordinal variables are first converted to ranks.
If nok is the number of nonzero weights, the dissimilarity is
multiplied by the factor 1/nok and thus ranges between 0 and 1.
If nok is zero, the dissimilarity is NA.


}
\section{BACKGROUND}{
Dissimilarities are used as inputs to cluster analysis and 
multidimensional scaling. The choice of metric may have a
large impact.


}
\references{
Kaufman, L. and Rousseeuw, P.J. (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis. }
Wiley, New York.


Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust 
Clustering Techniques in S-PLUS,
\emph{Computational Statistics and Data Analysis, \bold{26}, 17-37.}


}
\seealso{
\code{\link{dissimilarity.object}}, \code{\link{dist}}, \code{\link{pam}}, \code{\link{fanny}}, \code{\link{clara}}, \code{\link{agnes}}, \code{\link{diana}}.
}
\examples{
data(agriculture)
## Example 1 in ref
## Compute the dissimilarities using Euclidean metric and without
## standardization 
daisy(agriculture, metric = "euclidean", stand = FALSE)

data(flower)
## Example 2 in ref
daisy(flower, type = list(asymm = 3))
daisy(flower, type = list(asymm = c(1, 3), ordratio = 7))
}
\keyword{libcluster}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
