#### PAM : Partitioning Around Medoids
#### --- $Id: pam.q,v 1.17 2003/04/30 14:24:00 maechler Exp $
pam <- function(x, k, diss = inherits(x, "dist"),
		metric = "euclidean", stand = FALSE, cluster.only = FALSE,
                keep.diss = !diss && !cluster.only && n < 100,
                keep.data = !diss && !cluster.only)
{
    if(diss <- as.logical(diss)) {
	## check type of input vector
	if(any(is.na(x)))
	    stop("NA-values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") {
	    if(!is.numeric(x) || is.na(sizeDiss(x)))
		stop("x is not of class dissimilarity and can not be converted to this class." )
	    ## convert input vector to class "dissimilarity"
	    class(x) <- ..dClass
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
	mdata <- FALSE
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
	if((mdata <- any(inax <- is.na(x2)))) { # TRUE if x[] has any NAs
	    jtmd <- as.integer(ifelse(apply(inax, 2, any), -1, 1))
	    ## VALue for MISsing DATa
	    valmisdat <- 1.1* max(abs(range(x2, na.rm=TRUE)))
	    x2[inax] <- valmisdat
	    valmd <- rep(valmisdat, jp)
	}
	dv <- double(1 + (n * (n - 1))/2)
    }
    if((k <- as.integer(k)) < 1 || k >= n)
	stop("Number of clusters `k' must be in {1,2, .., n-1}; hence n >= 2")
    ## call Fortran routine
    storage.mode(dv) <- "double"
    storage.mode(x2) <- "double"
    res <- .C("pam",
	      as.integer(n),
	      as.integer(jp),
	      k,
	      x = x2,
	      dys = dv,
	      jdyss = as.integer(diss),
	      if(mdata)valmd else double(1),
	      if(mdata) jtmd else integer(1),
	      as.integer(ndyst),
	      integer(n),		# nsend[]
	      logical(n),		# nrepr[]
	      integer(if(cluster.only) 1 else n), # nelem[]
	      double(n),		# radus[]
	      double(n),		# damer[]
	      avsil = double(n),	# `ttd'
	      double(n),		# separ[]
	      ttsil = as.double(0),
	      obj = c(as.double(cluster.only), 0.),# in & out!
	      med = integer(if(cluster.only) 1 else k),
	      clu = integer(n),
	      clusinf = if(cluster.only) 0. else matrix(0., k, 5),
	      silinf  = if(cluster.only) 0. else matrix(0., n, 4),
	      isol = integer(if(cluster.only) 1 else k),
	      DUP = FALSE, # care!!
	      PACKAGE = "cluster")

    xLab <- if(diss) attr(x, "Labels") else dimnames(x)[[1]]
    if(length(xLab) > 0)
        names(res$clu) <- xLab
    if(cluster.only)
        return(res$clu)

    ## Else, usually
    sildim <- res$silinf[, 4]
    if(diss) {
	disv <- x
	## add labels to Fortran output
	if(length(xLab) > 0) {
	    sildim <- xLab[sildim]
	    res$med <- xLab[res$med]
	}
    }
    else {
	## give warning if some dissimilarities are missing.
	if(res$jdyss == -1)
	    stop("No clustering performed, NAs in the computed dissimilarity matrix.")
        if(keep.diss) {
            ## adapt Fortran output to S:
            ## convert lower matrix, read by rows, to upper matrix, read by rows.
            disv <- res$dys[-1]
            disv[disv == -1] <- NA
            disv <- disv[upper.to.lower.tri.inds(n)]
            class(disv) <- ..dClass
            attr(disv, "Size") <- nrow(x)
            attr(disv, "Metric") <- metric
            attr(disv, "Labels") <- dimnames(x)[[1]]
        }
	## add labels to Fortran output
	res$med <- x[res$med,  , drop =FALSE]
	if(length(xLab) > 0)
	    sildim <- xLab[sildim]
    }
    ## add dimnames to Fortran output
    names(res$obj) <- c("build", "swap")
    res$isol <- factor(res$isol, levels = 0:2, labels = c("no", "L", "L*"))
    names(res$isol) <- 1:k
    dimnames(res$clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
					  "diameter", "separation"))
    ## construct S object
    r <-
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
	     diss = if(keep.diss)disv,
	     call = match.call())
    if(keep.data && !diss) {
	if(mdata) x2[x2 == valmisdat] <- NA
	r$data <- x2
    }
    class(r) <- c("pam", "partition")
    r
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
    cat("Medoids:\n");			print(x$medoids, ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    cat("Objective function:\n");	print(x$objective, ...)
    cat("\nNumerical information per cluster:\n"); print(x$clusinfo, ...)
    cat("\nIsolated clusters:\n L-clusters: ")
    print(names(x$isolation[x$isolation == "L"]), quote = FALSE, ...)
    cat(" L*-clusters: ")
    print(names(x$isolation[x$isolation == "L*"]), quote = FALSE, ...)
    if(length(x$silinfo) != 0) {
	cat("\nSilhouette plot information:\n")
	print(x$silinfo[[1]], ...)
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of total data set:\n")
	print(x$silinfo[[3]], ...)
    }
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

