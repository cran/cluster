% %Z%%M% version %I% created %G%
% %Z%%Q%
\name{mona}
\alias{mona}
\title{
Monothetic Analysis
}
\description{
Returns a list representing a divisive hierarchical clustering of
a dataset with binary variables only.
}
\usage{
mona(x)
}
\arguments{
\item{x}{
data matrix or data frame in which each row corresponds to an observation,
and each column corresponds to a variable. All variables must be binary.
A limited number of missing values (NAs) is allowed. Every observation
must have at least one value different from NA. No variable should have
half of its values missing. There must be at least one variable which has
no missing values. A variable with all its non-missing values identical,
is not allowed.


}}
\value{
an object of class \code{"mona"} representing the clustering.
See \code{mona.object} for details.
}
\details{
\code{mona} is fully described in chapter 7 of Kaufman and Rousseeuw (1990).
It is "monothetic" in the sense that each division is based on a
single (well-chosen) variable, whereas most other hierarchical methods
(including \code{agnes} and \code{diana}) are "polythetic", i.e. they use
all variables together.


The \code{mona}-algorithm constructs a hierarchy of clusterings,
starting with one large
cluster. Clusters are divided until all observations in the same cluster have
identical values for all variables.
At each stage, all clusters are divided according to the values of one
variable. A cluster is divided into one cluster with all observations having
value 1 for that variable, and another cluster with all observations having
value 0 for that variable.


The variable used for splitting a cluster is the variable with the maximal
total association to the other variables, according to the observations in the
cluster to be splitted. The association between variables f and g
is given by a(f,g)*d(f,g) - b(f,g)*c(f,g), where a(f,g), b(f,g), c(f,g),
and d(f,g) are the numbers in the contingency table of f and g.
[That is, a(f,g) (resp. d(f,g)) is the number of observations for which f and g
both have value 0 (resp. value 1); b(f,g) (resp. c(f,g)) is the number of
observations for which f has value 0 (resp. 1) and g has value 1 (resp. 0).]
The total association of a variable f is the sum of its associations to all
variables.


This algorithm does not work with missing values, therefore the data are
revised, e.g. all missing values are filled in. To do this, the same measure
of association between variables is used as in the algorithm. When variable
f has missing values, the variable g with the largest absolute association
to f is looked up. When the association between f and g is positive,
any missing value of f is replaced by the value of g for the same
observation. If the association between f and g is negative, then any missing
value of f is replaced by the value of 1-g for the same
observation.
}
\section{BACKGROUND}{
Cluster analysis divides a dataset into groups (clusters) of observations that
are similar to each other. Hierarchical methods like \code{agnes}, \code{diana}, and
\code{mona} construct a hierarchy of clusterings, with the number of clusters
ranging from one to the number of observations. Partitioning methods like 
\code{pam}, \code{clara}, and \code{fanny} require that the number of clusters be given by
the user.
}
\references{
    Kaufman, L. and Rousseeuw, P.J. (1990). 
    \emph{Finding Groups in Data: An Introduction to Cluster Analysis}.
    Wiley, New York.

    Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996):
    Clustering in an Object-Oriented Environment.
    \emph{Journal of Statistical Software}, \bold{1}.
    \url{http://www.stat.ucla.edu/journals/jss/}

    Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating
    Robust Clustering Techniques in S-PLUS,
    \emph{Computational	Statistics and Data Analysis}, \bold{26},
    17-37.
}
\seealso{
\code{\link{mona.object}}, \code{\link{plot.mona}}.
}
\examples{
data(animals)
ma <- mona(animals)
ma
## Plot similar to Figure 10 in Struyf et al (1996)
plot(ma)
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.
