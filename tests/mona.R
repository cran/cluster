library(cluster)

data(animals)
(mani <- mona(animals))

str(mani)

if(require(MASS)) {

    if(paste(R.version$major, R.version$minor, sep=".") >= 1.7)  RNGversion(1.6)
    set.seed(253)
    n <- 512; p <- 3
    Sig <- diag(p); Sig[] <- 0.8 ^ abs(col(Sig) - row(Sig))
    x3 <- mvrnorm(n, rep(0,p), Sig) >= 0
    x <- cbind(x3, rbinom(n, size=1, prob = 1/2))

    (mx <- mona(x))
    str(mx)
}
