### $Id: clara.q,v 1.15 2002/09/09 09:36:53 maechler Exp $

#### CLARA := Clustering LARge Applications
####
#### Note that the algorithm is O(n), but O(ns^2) where ns == sampsize

### FIXME :
##  should not necessarily keep data in result, because "large" !
##  OTOH, data is used for clusplot.partition() :
## Note:  ./plotpart.q	is also working with clara() objects

clara <- function(x, k, metric = "euclidean", stand = FALSE,
		  samples = 5, sampsize = 40 + 2 * k,
		  trace = 0, keepdata = TRUE)
{
    ## check type of input matrix and values of input numbers
    x <- data.matrix(x)
    if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")
    n <- nrow(x)
    if((k <- as.integer(k)) < 1 || k > n - 1)
	stop("The number of cluster should be at least 1 and at most n-1." )
    if((sampsize <- as.integer(sampsize)) < max(2,k))
	stop(paste(c("'sampsize' should be at least", max(2,k),
		     " = max(2, number of clusters)"), collapse = " "))
    if(n < sampsize)
	stop(paste(c("Number of objects is", n,
		     ", should be at least", sampsize, "(sampsize)"),
		   collapse = " "))
    if((samples <- as.integer(samples)) < 1)
	stop("'samples' should be at least 1")

    namx <- dimnames(x)[[1]]
    ## standardize, if necessary
    data <- x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
    ## put info about metric, size and NAs in arguments for the Fortran call
    jp <- ncol(x2)

    if((mdata <- any(inax <- is.na(x2)))) { # TRUE if x[] has any NAs
	jtmd <- as.integer(ifelse(apply(inax, 2, any), -1, 1))
	## VALue for MISsing DATa
	valmisdat <- 1.1* max(abs(range(x2, na.rm=TRUE)))
	x2[inax] <- valmisdat
    }

    if((trace <- as.integer(trace)))
	cat("calling .C(\"clara\", *):\n")
    res <- .C("clara",
	      n,
	      jp,
	      k,
	      clu = as.double(t(x2)),# transposing LARGE x[n,jp] ..fixme?..
	      nran  = samples,
	      nsam  = sampsize,
	      dis   = double(1 + (sampsize * (sampsize - 1))/2),
	      mdata = as.integer(mdata),
	      valmd = if(mdata) rep(valmisdat, jp) else -1.,
	      jtmd  = if(mdata) jtmd else integer(1),
	      ndyst = as.integer(if(metric == "manhattan") 2 else 1),
	      integer(sampsize),# = nrepr
	      integer(sampsize),# = nsel
	      sample= integer(sampsize),# = nbest
	      integer(k),		# = nr
	      med = integer(k),		# = nrx
	      double(k),		# = radus
	      double(k),		# = ttd
	      double(k),		# = ratt
	      avdis  = double(k),	# = ttbes
	      maxdis = double(k),	# = rdbes
	      ratdis = double(k),	# = rabes
	      size  = integer(k),	# = mtt
	      obj   = double(1),
	      avsil = double(k),
	      ttsil = double(1),
	      silinf = matrix(0, sampsize, 4),
	      jstop = integer(1),
	      trace = trace,
	      tmp  = double (3 * sampsize),
	      itmp = integer(6 * sampsize),
	      DUP = FALSE,
	      PACKAGE = "cluster")
    ## give a warning when errors occured
    if(res$jstop == 1)
	stop("For each sample at least one object was found which\n",
	     " could not be assigned to a cluster (because of missing values).")
    if(res$jstop == 2)
	stop("Each of the random samples contains objects between which\n",
	     " no distance can be computed.")
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
    res$clu <- as.integer(matrix(res$clu, nrow= n, ncol= jp, byrow= TRUE)[, 1])
    if(!is.null(namx)) {
	sildim <- namx[sildim]
	res$sample <- namx[res$sample]
	names(res$clu) <- namx
    }
    ## add dimnames to Fortran output
    r <- list(sample = res$sample, medoids = res$med,
	      clustering = res$clu, objective = res$obj,
	      clusinfo = cbind(size = res$size, "max_diss" = res$maxdis,
	      "av_diss" = res$avdis, isolation = res$ratdis),
	      diss = disv, call = match.call())
    if(k > 1) {
	dimnames(res$silinf) <- list(sildim,
				     c("cluster", "neighbor", "sil_width", ""))
	r$silinfo <- list(widths = res$silinf[, -4],
			  clus.avg.widths = res$avsil,
			  avg.width = res$ttsil)
    }
    if(keepdata) r$data <- data
    class(r) <- c("clara", "partition")
    r
}

print.clara <- function(x, ...)
{
    cat("Call:	", deparse(x$call),
	"\nMedoids:\n");		print(x$medoids, ...)
    cat("Objective function:\t ", format(x$objective, ...),"\n",
	"Clustering vector: \t", sep=""); str(x$clustering, vec.len = 7)
    cat("Cluster sizes:	    \t", x$clusinfo[,1],
	"\nBest sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

summary.clara <- function(object, ...)
{
    class(object) <- "summary.clara"
    object
}

print.summary.clara <- function(x, ...)
{
    cat("Object of class `clara' from call:\n", deparse(x$call),
	"\nMedoids:\n");		print(x$medoids, ...)
    cat("Objective function:\t ", format(x$objective, ...),
	"\nNumerical information per cluster:\n")
    print(x$clusinfo, ...)
    if(has.sil <- !is.null(x$silinfo)) {
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of best sample:",
	    format(x$silinfo[[3]], ...), "\n")
    }
    cat("\nBest sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    if(has.sil) {
	cat("\nSilhouette plot information for best sample:\n")
	print(x$silinfo[[1]], ...)
    }
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

