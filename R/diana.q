diana <- function(x, diss = FALSE, metric = "euclidean", stand = FALSE)
{
    meanabsdev <- function(y)
    {
	mean(abs(y - mean(y, na.rm = TRUE)), na.rm = TRUE)
    }
    size <- function(d)
    {
	discr <- 1 + 8 * length(d)
	sqrtdiscr <- round(sqrt(discr))
	if(round(sqrtdiscr)^2 != discr)
	    return(0)
	(1 + sqrtdiscr)/2
    }
    lower.to.upper.tri.inds <- function(n)
    {
	return(c(0, unlist(lapply(2:(n - 1), function(x, n)
				  cumsum(c(0, (n - 2):(n - x))), n = n))) +
	       rep(1:(n - 1), 1:(n - 1)))
    }
    upper.to.lower.tri.inds <- function(n)
    {
	return(unlist(lapply(0:(n - 2), function(x, n)
			     cumsum(x:(n - 2)), n = n)) +
	       rep(1 + cumsum(0:(n - 2)), (n - 1):1))
    }
    if(diss) {
	##check type of input vector
	if(is.na(min(x)))
	    stop("NA-values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") {
	    if(!is.numeric(x) || size(x) == 0)
		stop("x is not of class dissimilarity and can not be converted to this class."
		     )
	    ##convert input vector to class "dissimilarity"
	    class(x) <- "dissimilarity"
	    attr(x, "Size") <- size(x)
	    attr(x, "Metric") <- "unspecified"
	}
	n <- attr(x, "Size")
	dv <- x[lower.to.upper.tri.inds(n)]
	##prepare arguments for the Fortran call
	dv <- c(0, dv)
	jp <- 1
	valmd <- double(1)
	jtmd <- integer(1)
	ndyst <- 0
	x2 <- double(n)
	jdyss <- 1
	dv2 <- double(1 + (n * (n - 1))/2)
    }
    else {
	##check type of input matrix
	if((!is.data.frame(x) && !is.numeric(x)) ||
	   (!all(sapply(x, data.class) == "numeric")))
	    stop("x is not a numeric dataframe or matrix.")
	x <- data.matrix(x) # standardize, if necessary
	x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
	ndyst <- if(metric == "manhattan") 2 else 1
	n <- nrow(x2)
	jp <- ncol(x2)
	jtmd <- ifelse(is.na(rep(1, n) %*% x2), -1, 1)
	valmisdat <- min(x2, na.rm = TRUE) - 0.5
	x2[is.na(x2)] <- valmisdat
	valmd <- rep(valmisdat, jp)
	jdyss <- 0
	dv <- double(1 + (n * (n - 1))/2)
	dv2 <- double(1 + (n * (n - 1))/2)
    }
    jalg <- 2
    ## call Fortran routine
    storage.mode(dv) <- "double"
    storage.mode(dv2) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(valmd) <- "double"
    storage.mode(jtmd) <- "integer"
    merge <- matrix(0, n - 1, 2)
    storage.mode(merge) <- "integer"
    res <- .Fortran("twins",
		    as.integer(n),
		    as.integer(jp),
		    x2,
		    dv,
		    dis = dv2,
		    ok = as.integer(jdyss),
		    valmd,
		    jtmd,
		    as.integer(ndyst),
		    as.integer(jalg),
		    as.integer(0),
		    integer(n),
		    ner = integer(n),
		    ban = double(n),
		    dc = as.double(0),
		    merge = merge,
		    PACKAGE = "cluster")
    if(!diss) {
	## give warning if some dissimilarities are missing.
	if(res$ok == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.\n"
		 )
	## adapt Fortran output to S:
	## convert lower matrix, read by rows, to upper matrix, read by rows.
	disv <- res$dis[-1]
	disv[disv == -1] <- NA
	disv <- disv[upper.to.lower.tri.inds(n)]
	class(disv) <- "dissimilarity"
	attr(disv, "Size") <- nrow(x)
	attr(disv, "Metric") <- metric
	attr(disv, "Labels") <- dimnames(x)[[1]]
	##add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0) {
	    order.lab <- dimnames(x)[[1]][res$ner]
	}
    }
    else {
	disv <- x
	##add labels to Fortran output
	if(length(attr(x, "Labels")) != 0) {
	    order.lab <- attr(x, "Labels")[res$ner]
	}
    }
    clustering <- list(order = res$ner, height = res$ban[-1], dc = res$dc,
		       merge = res$merge, diss = disv)
    if(exists("order.lab"))
	clustering$order.lab <- order.lab
    if(!diss) {
	x2[x2 == valmisdat] <- NA
	clustering$data <- x2
    }
    class(clustering) <- c("diana", "twins")
    attr(clustering, "Call") <- sys.call()
    clustering
}

print.diana <- function(x, ...)
{
    cat("Merge:\n")
    print(x$merge, ...)
    cat("Order of objects:\n")
    print(if (length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Height:\n")
    print(x$height, ...)
    cat("Divisive coefficient:\n")
    print(x$dc, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.diana <- function(object, ...)
{
    class(object) <- "summary.diana"
    object
}

print.summary.diana <- function(x, ...)
{
    cat("Merge:\n")
    print(x$merge, ...)
    cat("Order of objects:\n")
    print(if (length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Height:\n")
    print(x$height, ...)
    cat("Divisive coefficient:\n")
    print(x$dc, ...)
    cat("\n")
    print(x$diss, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}
