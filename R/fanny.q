fanny <- function(x, k, diss = FALSE, metric = "euclidean", stand = FALSE)
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
    if(diss) {
	## check type of input vector
	if(is.na(min(x)))
	    stop("NA-values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") {
	    if(!is.numeric(x) || size(x) == 0)
		stop("x is not of class dissimilarity and can not be converted to this class." )	
	    ## convert input vector to class "dissimilarity"
	    class(x) <- "dissimilarity"
	    attr(x, "Size") <- size(x)
	    attr(x, "Metric") <- "unspecified"
	}
	## prepare arguments for the Fortran call
	n <- attr(x, "Size")
	if((k < 1) || (k > floor(n/2) - 1))
	    stop("The number of cluster should be at least 1 and at most n/2 - 1." )
	dv <- c(x, 0)
	jp <- 1
	valmd <- double(1)
	jtmd <- integer(1)
	ndyst <- 0
	x2 <- double(n)
	jdyss <- 1
    }
    else {
	##check type of input matrix 
	if((!is.data.frame(x) && !is.numeric(x)) ||
	   (!all(sapply(x, data.class) == "numeric")))
	    stop("x is not a numeric dataframe or matrix.")
	x <- data.matrix(x)	
	## standardize, if necessary
	x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
	##put info about metric, size and NAs in arguments for the Fortran call
	ndyst <- if(metric == "manhattan") 2 else 1
	n <- nrow(x2)
	if((k < 1) || (k > floor(n/2) - 1))
	    stop("The number of cluster should be at least 1 and at most n/2 - 1."
		 )
	jp <- ncol(x2)
	jtmd <- ifelse(is.na(rep(1, n) %*% x2), -1, 1)
	valmisdat <- min(x2, na.rm = TRUE) - 0.5
	x2[is.na(x2)] <- valmisdat
	valmd <- rep(valmisdat, jp)
	jdyss <- 0
	dv <- double(1 + (n * (n - 1))/2)
    }
    ##call Fortran routine
    storage.mode(dv) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(valmd) <- "double"
    storage.mode(jtmd) <- "integer"
    res <- .Fortran("fanny",
		    as.integer(n),
		    as.integer(jp),
		    as.integer(k),
		    x2,
		    dis = dv,
		    ok = as.integer(jdyss),
		    valmd,
		    jtmd,
		    as.integer(ndyst),
		    integer(n),
		    integer(n),
		    integer(n),
		    double(n),
		    p = matrix(0, n, k),
		    matrix(0, n, k),
		    avsil = double(k),
		    integer(k),
		    double(k),
		    double(k),
		    double(n),
		    ttsil = as.double(0),
		    eda = as.double(0),
		    edb = as.double(0),
		    obj = double(2),
		    clu = integer(n),
		    silinf = matrix(0, n, 4),
		    as.double(1e-15),
		    PACKAGE = "cluster")
    sildim <- res$silinf[, 4]
    if(diss) {
	disv <- x	
	##add labels to Fortran output
	if(length(attr(x, "Labels")) != 0) {
	    sildim <- attr(x, "Labels")[sildim]
	    dimnames(res$p) <- list(attr(x, "Labels"), NULL)
	    names(res$clu) <- attr(x, "Labels")
	}
    }
    else {
	##give warning if some dissimilarities are missing.
	if(res$ok == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.")	
	disv <- res$dis[ - (1 + (n * (n - 1))/2)]
	disv[disv == -1] <- NA
	class(disv) <- "dissimilarity"
	attr(disv, "Size") <- nrow(x)
	attr(disv, "Metric") <- metric
	attr(disv, "Labels") <- dimnames(x)[[1]]	
	##add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0) {
	    sildim <- dimnames(x)[[1]][sildim]
	    dimnames(res$p) <- list(dimnames(x)[[1]], NULL)
	    names(res$clu) <- dimnames(x)[[1]]
	}
    }
    ##add dimnames to Fortran output
    names(res$obj) <- c("iterations", "objective")
    res$coeff <- c(res$eda, res$edb)
    names(res$coeff) <- c("dunn_coeff", "normalized")
    if(k != 1) {
	dimnames(res$silinf) <- list(sildim,
				     c("cluster", "neighbor", "sil_width", ""))
	clustering <- list(membership = res$p, coeff = res$coeff, 
			   clustering = res$clu, objective = res$obj,
			   silinfo = 
			   list(widths = res$silinf[, -4],
				clus.avg.widths = res$avsil[1:k],
				avg.width = res$ttsil),
			   diss = disv)
    }
    else {
	clustering <- list(membership = res$p, coeff = res$coeff, 
			   clustering = res$clu, objective = res$obj,
			   diss = disv)
    }
    if(!diss) {
	x2[x2 == valmisdat] <- NA
	clustering$data <- x2
    }
    class(clustering) <- c("fanny", "partition")
    attr(clustering, "Call") <- sys.call()
    clustering
}

print.fanny <- function(x, ...)
{
    print(x$objective, ...)
    cat("Membership coefficients:\n")
    print(x$membership, ...)
    cat("Coefficients:\n")
    print(x$coeff, ...)
    cat("Closest hard clustering:\n")
    print(x$clustering, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.fanny <- function(x, ...)
{
    object <- x
    class(object) <- "summary.fanny"
    object
}

print.summary.fanny <- function(x, ...)
{
    print(x$objective, ...)
    cat("Membership coefficients:\n")
    print(x$membership, ...)
    cat("Coefficients:\n")
    print(x$coeff, ...)
    cat("Closest hard clustering:\n")
    print(x$clustering, ...)
    if(length(x$silinfo) != 0) {
	cat("\nSilhouette plot information:\n")
	print(x$silinfo[[1]], ...)
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of total data set:\n")
	print(x$silinfo[[3]], ...)
    }
    cat("\n")
    print(x$diss, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

