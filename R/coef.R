#### R-interface to  Agglomerative / Divisive coefficient
####
coef.twins <- function(object, ...)
{
    if(inherits(object, "agnes"))
	object$ac
    else if(inherits(object, "diana"))
	object$dc
    else
	stop("invalid 'twins' object")
}

coef.hclust <- function(object, ...)
{
    ## Author: Martin Maechler, Date: 27 Nov 2004
    ## Now "really" using $merge _and_ $height -- assuming they match!
    ht  <- object$height
    mrg <- object$merge
    nh <- length(ht)
    stopifnot(nh > 0, is.matrix(mrg), dim(mrg) == c(nh,2),
              is.numeric(ht), is.numeric(mrg),
              !is.unsorted(ht))# then they match with merge
    ## stopifnot(all.equal(1:n, sort(-mrg[mrg < 0])))

    1 - sum(rowSums(mrg < 0) * ht) / max(ht) / (nh+1)
}

if(FALSE){ ##-- experiments

coef.hclust <- function(object, ...)
{
    ## Purpose: Compute agglomerative coefficient from hclust
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Nov 2004
    nh <- length(ht <- object$height)
    stopifnot(nh > 0, is.numeric(ht))
    .Fortran("bncoef",
	     n =  as.integer(nh + 1),
	     ban= as.double(c(0., ht)),
	     cf = double(1),
	     PACKAGE = "cluster")$cf
}

### The same with R-only code -- vectorized
coef2.hclust <- function(object, ...)
{
    ## Author: Martin Maechler, Date: 27 Nov 2004
    nh <- length(ht <- object$height)
    stopifnot(nh > 0, is.numeric(ht))
    1- mean(c(ht[1], pmin(ht[1:(nh-1)], ht[2:nh]), ht[nh]))/max(ht)
}

m.hclust <- function(object, ...)
{
    ## Author: Martin Maechler, Date: 27 Nov 2004
    ## Now "really" using $merge _and_ $height -- assuming they match!
    ## and slowly implementing the algorithm
    nh <- length(ht  <- object$height)
    nm <- nrow(mrg <- object$merge)
    n <- nh+1 # number of objects
    stopifnot(nh > 0, is.matrix(mrg), nm == nh,
              is.numeric(ht), is.numeric(mrg),
              all.equal(1:n, sort(-mrg[mrg < 0])))

    max.h <- max(ht)
    m <- numeric(n)
    for(i in 1:n) { ## compute m[i]: the height at which i-th object is merged
        j <- which(mrg == -i, arr.ind = TRUE)[,1]
        m[i] <- ht[j]
    }
    m <- m / max(ht) ## == (m[i]) as in ../man/agnes.object.Rd
    m
    ## 1 - mean(m.hclust())  should give *the* coefficient
}


}
