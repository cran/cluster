\name{flower}
\alias{flower}
\title{Flower Characteristics}
\usage{data(flower)}
\description{8 characteristics for 18 popular flowers.}
\format{
  A data frame with 18 observations on 8 variables:
  \tabular{rll}{
    [ , 1] \tab factor  \tab winters \cr
    [ , 2] \tab factor  \tab shadow \cr
    [ , 3] \tab factor  \tab tubers \cr
    [ , 4] \tab factor  \tab color \cr
    [ , 5] \tab ordered \tab soil \cr
    [ , 6] \tab ordered \tab preference \cr
    [ , 7] \tab numeric \tab height \cr
    [ , 8] \tab numeric \tab distance
  }

  \describe{
    \item{winters}{is binary and indicates whether the plant may be left
      in the garden when it freezes.}

    \item{shadow}{is binary and shows whether the plant needs to stand
      in the shadow.}

    \item{tubers}{is asymmetric binary and distinguishes between plants
      with tubers and plants that grow in any other way.}

    \item{color}{is nominal and specifies the flower's color (1 = white,
      2 = yellow, 3 = pink, 4 = red, 5 = blue).}

    \item{soil}{is ordinal and indicates whether the plant grows in dry
      (1), normal (2), or wet (3) soil.}

    \item{preference}{is ordinal and gives someone's preference ranking
      going from 1 to 18.}

    \item{height}{is interval scaled, the plant's height in centimeters.}

    \item{distance}{is interval scaled, the distance in centimeters that
      should be left between the plants.}
  }
}
\source{The reference below.}
\references{
  Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996):
  Clustering in an Object-Oriented Environment.
  \emph{Journal of Statistical Software}, \bold{1}.
  \url{http://www.stat.ucla.edu/journals/jss/}
}
\examples{
data(flower)
## Example 2 in ref
daisy(flower, type = list(asymm = 3))
daisy(flower, type = list(asymm = c(1, 3), ordratio = 7))
}
\keyword{datasets}
