library(cluster)

data(animals)
(mani <- mona(animals))

str(mani)

set.seed(253)
require(MASS)
n <- 512; p <- 3
Sig <- diag(p); Sig[] <- 0.8 ^ abs(col(Sig) - row(Sig))
x3 <- mvrnorm(n, rep(0,p), Sig) >= 0
x <- cbind(x3, rbinom(n, size=1, prob = 1/2))

(mx <- mona(x))
str(mx)

