#### PAM : Partitioning Around Medoids
#### --- $Id: pam.q,v 1.11 2002/03/04 10:44:45 maechler Exp maechler $
pam <- function(x, k, diss = inherits(x, "dist"),
                metric = "euclidean", stand = FALSE)
{
    if(diss) {
	## check type of input vector
	if(is.na(min(x)))
	    stop("NA-values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") {
	    if(!is.numeric(x) || is.na(sizeDiss(x)))
		stop("x is not of class dissimilarity and can not be converted to this class." )
	    ## convert input vector to class "dissimilarity"
	    class(x) <- "dissimilarity"
	    attr(x, "Size") <- sizeDiss(x)
	    attr(x, "Metric") <- "unspecified"
	}
	## adapt S dissimilarities to Fortran:
	## convert upper matrix, read by rows, to lower matrix, read by rows.
	n <- attr(x, "Size")
	dv <- x[lower.to.upper.tri.inds(n)]
	## prepare arguments for the Fortran call
	dv <- c(0, dv)
	jp <- 1
	valmd <- double(1)
	jtmd <- integer(1)
	ndyst <- 0
	x2 <- double(n)
    }
    else {
        ## check input matrix and standardize, if necessary
	x <- data.matrix(x)
	if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")
        x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
	## put info about metric, size and NAs in arguments for the Fortran call
	ndyst <- if(metric == "manhattan") 2 else 1
	n <- nrow(x2)
	jp <- ncol(x2)
	jtmd <- as.integer(ifelse(is.na(rep(1, n) %*% x2), -1, 1))
	valmisdat <- min(x2, na.rm=TRUE) - 0.5 #(double) VALue for MISsing DATa
	x2[is.na(x2)] <- valmisdat
	valmd <- rep(valmisdat, jp)
	dv <- double(1 + (n * (n - 1))/2)
    }
    if((k <- as.integer(k)) < 1 || k >= n)
        stop("The number of clusters `k' must be in {1,2, .., n-1}.")
    ## call Fortran routine
    storage.mode(dv) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(valmd) <- "double"
    res <- .Fortran("pam",
		    as.integer(n),
		    as.integer(jp),
		    k,
		    x = x2,
		    dys = dv,
                    jdyss = as.integer(diss),# 0/1
		    valmd,
		    jtmd,
		    as.integer(ndyst),
		    integer(n),
		    integer(n),
		    integer(n),
		    double(n),
		    double(n),
		    avsil = double(n),
		    double(n),
		    ttsil = as.double(0),
		    med = integer(k),
		    obj = double(2),
		    clu = integer(n),
		    clusinf = matrix(0., k, 5),
		    silinf = matrix(0., n, 4),
		    isol = integer(k),
		    PACKAGE = "cluster")
    sildim <- res$silinf[, 4]
    if(diss) {
	disv <- x
	## add labels to Fortran output
	if(length(attr(x, "Labels")) != 0) {
	    sildim <- attr(x, "Labels")[sildim]
	    names(res$clu) <- attr(x, "Labels")
	    res$med <- attr(x, "Labels")[res$med]
	}
    }
    else {
	## give warning if some dissimilarities are missing.
	if(res$jdyss == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.")
	## adapt Fortran output to S:
	## convert lower matrix, read by rows, to upper matrix, read by rows.
	disv <- res$dys[-1]
	disv[disv == -1] <- NA
	disv <- disv[upper.to.lower.tri.inds(n)]
	class(disv) <- "dissimilarity"
	attr(disv, "Size") <- nrow(x)
	attr(disv, "Metric") <- metric
	attr(disv, "Labels") <- dimnames(x)[[1]]
	## add labels to Fortran output
	res$med <- x[res$med,  ]
	if(length((dimnames(x)[[1]])) != 0) {
	    sildim <- dimnames(x)[[1]][sildim]
	    names(res$clu) <- dimnames(x)[[1]]
	}
    }
    ## add dimnames to Fortran output
    names(res$obj) <- c("build", "swap")
    res$isol <- factor(res$isol, levels = c(0, 1, 2),
		       labels = c("no", "L", "L*"))
    names(res$isol) <- 1:k
    dimnames(res$clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
					  "diameter", "separation"))
    ## construct S object
    clustering <-
        list(medoids = res$med, clustering = res$clu,
             objective = res$obj, isolation = res$isol,
             clusinfo = res$clusinf,
             silinfo = if(k != 1) {
                 dimnames(res$silinf) <-
                     list(sildim, c("cluster", "neighbor", "sil_width", ""))
                 list(widths = res$silinf[, -4],
                      clus.avg.widths = res$avsil[1:k],
                      avg.width = res$ttsil)
             },
             diss = disv,
             call = match.call())
    if(!diss) {
	x2[x2 == valmisdat] <- NA
	clustering$data <- x2
    }
    class(clustering) <- c("pam", "partition")
    clustering
}

print.pam <- function(x, ...)
{
    cat("Medoids:\n")
    print(x$medoids, ...)
    cat("Clustering vector:\n")
    print(x$clustering, ...)
    cat("Objective function:\n")
    print(x$objective, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.pam <- function(object, ...)
{
    class(object) <- "summary.pam"
    object
}

print.summary.pam <- function(x, ...)
{
    cat("Medoids:\n")
    print(x$medoids, ...)
    cat("Clustering vector:\n")
    print(x$clustering, ...)
    cat("Objective function:\n")
    print(x$objective, ...)
    cat("\nNumerical information per cluster:\n")
    print(x$clusinfo, ...)
    cat("\nIsolated clusters:\n")
    cat("L-clusters: ")
    print(names(x$isolation[x$isolation == "L"]), quote = FALSE, ...)
    cat("L*-clusters: ")
    print(names(x$isolation[x$isolation == "L*"]), quote = FALSE, ...)
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

