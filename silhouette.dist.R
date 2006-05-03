## silhouette.dist.R
##
## Author: Romain Francois <Romain.Francois@inria.fr>
##
###############################################################################

## MM: What I don't like is that  *order* of the arguments
## --  the first argument should be the important one, i.e. the clustering;
##     this is the case for all other silhouette.*() methods.
##    --> rather make this as one of the two methods of  silhouette.default()

##--> first ask Romain if that's fine with him...

silhouette.dist <- function(x, clustering, ...)
{
    cll <- match.call()

    n <- length(clustering)
    if(!all(clustering == round(clustering)))
        stop("'clustering' must only have integer codes")
    k <- length(unique(clustering))
    if(k <= 1 || k >= n) # silhouette undefined for trivial clusterings
        return(NA)

    wds <- matrix(NA, n,3, dimnames =
                  list(names(clustering), c("cluster","neighbor","sil_width")))

    out <- .C('sildist',
              d = as.numeric(x),
              as.integer(n),
              as.integer(clustering),
              as.integer(k),
              diC =    numeric(n*k),
              counts = integer(n*k),
              bi = numeric(n),
              ai = numeric(n),
              si = numeric(n),
              neighbor = integer(n),
              DUP = FALSE, PACKAGE = "cluster")

    wds[,1] <- clustering
    wds[,2] <- out$neighbor
    wds[,3] <- out$si

    attr(wds, "Ordered") <- FALSE
    attr(wds, "call") <- cll
    class(wds) <- "silhouette"
    wds
}
