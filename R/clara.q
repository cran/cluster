#### CLARA := Clustering LARge Applications
####
#### Note that the algorithm is O(n), but O(ns^2) where ns == sampsize

### FIXME :
##  should not necessarily keep data in result, because "large" !
##  OTOH, data is used for clusplot.partition() 

## Note:  ./plotpart.q	is also working with clara() objects 


clara <- function(x, k, metric = "euclidean", stand = FALSE,
		  samples = 5, sampsize = 40 + 2 * k)
{
    meanabsdev <- function(y)
	mean(abs(y - mean(y, na.rm = TRUE)), na.rm = TRUE)
    upper.to.lower.tri.inds <- function(n)
    {
	return(unlist(lapply(0:(n - 2), function(x, n)
			     cumsum(x:(n - 2)), n = n)) +
	       rep(1 + cumsum(0:(n - 2)), (n - 1):1))
    }
    ## check type of input matrix and values of input numbers
    if((!is.data.frame(x) && !is.numeric(x)) ||
       (!all(sapply(x, data.class) == "numeric")))
	stop(message = "x is not a numeric dataframe or matrix.")
    x <- data.matrix(x)
    n <- nrow(x)
    if((k <- as.integer(k)) < 1 || k > n - 1)
	stop("The number of cluster should be at least 1 and at most n-1." )
    if((sampsize <- as.integer(sampsize)) < k)
	stop(paste(c("'sampsize' should be at least", k,
		     "(number of clusters)"), collapse = " "))
    if(n < sampsize)
	stop(paste(c("Number of objects is", n,
		     ", should be at least", sampsize, "(sampsize)"),
		   collapse = " "))
    namx <- dimnames(x)[[1]]
    ## standardize, if necessary
    x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
    ## put info about metric, size and NAs in arguments for the Fortran call
    jp <- ncol(x2)
    jtmd  <- ifelse(is.na(rep(1, n) %*% x2), -1, 1)
    mdata <- is.na(min(x2))

    ## FIXME: The following will go wrong as soon as  min(x2) < -5e15
    valmisdat <- min(x2, na.rm = TRUE) - 0.5
    x2[is.na(x2)] <- valmisdat

    x3 <- as.double(as.vector(t(x2)))# transposing LARGE x ..not efficient ....
    ## call Fortran routine
    res <- .Fortran("clara",
		    as.integer(n),
		    as.integer(jp),
		    as.integer(k),
		    clu = x3,# transpose (x [n * jp] )
		    nran  = as.integer(samples),
		    nsam  = sampsize,
		    dis	  = double(1 + (sampsize * (sampsize - 1))/2),
		    mdata = as.integer(mdata),
		    valmd = rep(as.double(valmisdat), jp),
		    jtmd  = as.integer(jtmd),
		    ndyst = as.integer(if(metric == "manhattan") 2 else 1),
		    integer(sampsize),
		    integer(sampsize),
		    sample = integer(sampsize),# = nbest
		    integer(k),
		    med = integer(k),# = nrx
		    double(k),
		    double(k),
		    double(k),
		    avdis  = double(k),# = ttbes
		    maxdis = double(k),# = rdbes
		    ratdis = double(k),# = rabes
		    size  = integer(k),# = mtt
		    obj	  = double(1),
		    avsil = double(k),
		    ttsil = double(1),
		    silinf = matrix(0, sampsize, 4),
		    jstop = as.integer(0),
		    double(sampsize),
		    double(sampsize),
		    double(sampsize),
		    integer(sampsize),
		    integer(sampsize),
		    integer(sampsize),
		    integer(sampsize),
		    integer(sampsize),
		    integer(sampsize),
		    PACKAGE = "cluster")	
    ## give a warning when errors occured
    if(res$jstop == 1)
	stop("For each sample at least one object was found which could not be assigned to a cluster (because of missing values).")
    if(res$jstop == 2)
	stop("Each of the random samples contains objects between which no distance can be computed.")
    sildim <- res$silinf[, 4]	
    ## adapt Fortran output to S:
    ## convert lower matrix, read by rows, to upper matrix, read by rows.
    disv <- res$dis[-1]
    disv[disv == -1] <- NA
    disv <- disv[upper.to.lower.tri.inds(sampsize)]
    class(disv) <- "dissimilarity"
    attr(disv, "Size") <- sampsize
    attr(disv, "Metric") <- metric	
    attr(disv, "Labels") <- namx[res$sample]
    ## add labels to Fortran output
    res$med <- x[res$med, ]
    res$clu <- matrix(res$clu, nrow = n, ncol = jp, byrow = TRUE)[, 1]
    if(length(namx) != 0) {
	sildim <- namx[sildim]
	res$sample <- namx[res$sample]
	names(res$clu) <- namx
    }
    ## add dimnames to Fortran output
    clusinf <- cbind(size = res$size, "max_diss" = res$maxdis,
		     "av_diss" = res$avdis, isolation = res$ratdis)
    if(k != 1) {
	dimnames(res$silinf) <- list(sildim,
				     c("cluster", "neighbor", "sil_width", ""))
	clustering <- list(sample = res$sample, medoids = res$med, 
			   clustering = res$clu, objective = res$obj,
			   clusinfo = clusinf,
			   silinfo = list(width = res$silinf[, -4], 
			   clus.avg.widths = res$avsil[1:k],
			   avg.width = res$ttsil),
			   diss = disv)
    }
    else {
	clustering <- list(sample = res$sample, medoids = res$med, 
			   clustering = res$clu, objective = res$obj,
			   clusinfo = clusinf, diss = disv)
    }
    x2[x2 == valmisdat] <- NA
    clustering$data <- x2
    class(clustering) <- c("clara", "partition")
    attr(clustering, "Call") <- sys.call()
    clustering
}

print.clara <- function(x, ...)
{
    cat("Best sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("Medoids:\n");			print(x$medoids, ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    cat("Objective function:\n");	print(x$objective, ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

summary.clara <- function(x, ...)
{
    class(x) <- "summary.clara"
    x
}

print.summary.clara <- function(x, ...)
{
    cat("Best sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("Medoids:\n");			print(x$medoids, ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    cat("Objective function:\n");	print(x$objective, ...)
    cat("\nNumerical information per cluster:\n")
    print(x$clusinfo, ...)
    if(length(x$silinfo) != 0) {
	cat("\nSilhouette plot information for best sample:\n")
	print(x$silinfo[[1]], ...)
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of best sample:\n")
	print(x$silinfo[[3]], ...)
    }
    cat("\n")
    print(x$diss, ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

