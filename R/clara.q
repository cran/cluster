"clara" <-
function(x, k, metric = "euclidean", stand = F, samples = 5, sampsize = 40 + 2 *
	k)
{
	meanabsdev <- function(y)
	{
		mean(abs(y - mean(y, na.rm = T)), na.rm = T)
	}
	upper.to.lower.tri.inds <- function(n)
	{
        	return(unlist(lapply(0:(n - 2), function(x, n)
	        	cumsum(x:(n - 2)), n = n)) +
			rep(1 + cumsum(0:(n - 2)), (n - 1):1))
	}
#check type of input matrix and values of input numbers
	if((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x, data.class) ==
		"numeric")))
		stop(message = "x is not a numeric dataframe or matrix.")
	x <- data.matrix(x)
	n <- nrow(x)
	if((k < 1) || (k > n - 1))
		stop(message = "The number of cluster should be at least 1 and at most n-1."
			)
	if(sampsize < k) {
		warning <- paste(c("'sampsize' should be at least", k,
			"(number of clusters)"), collapse = " ")
		stop(message = warning)
	}
	if(n < sampsize) {
		warning <- paste(c("Number of objects is", n,
			", should be at least", sampsize, "(sampsize)"),
			collapse = " ")
		stop(message = warning)
	}
#standardize, if necessary
	if(stand) {
		x2 <- scale(x, scale = apply(x, 2, meanabsdev))
	}
	else x2 <- x
#put info about metric, size and NAs in arguments for the Fortran call
	if(metric == "manhattan")
		ndyst <- 2
	else ndyst <- 1
	n <- nrow(x2)
	jp <- ncol(x2)
	jtmd <- ifelse(is.na(rep(1, n) %*% x2), -1, 1)
	mdata <- ifelse(is.na(min(x2)), 1, 0)
	valmisdat <- min(x2, na.rm = T) - 0.5
	x2[is.na(x2)] <- valmisdat
	valmd <- rep(valmisdat, jp)
	jdyss <- 0
	x3 <- as.vector(t(x2))	
#call Fortran routine
	storage.mode(x3) <- "double"
	storage.mode(valmd) <- "double"
	storage.mode(jtmd) <- "integer"
	res <- .Fortran("clara",
		as.integer(n),
		as.integer(jp),
		as.integer(k),
		clu = x3,
		as.integer(samples),
		as.integer(sampsize),
		dis = double(1 + (sampsize * (sampsize - 1))/2),
		as.integer(mdata),
		valmd,
		jtmd,
		as.integer(ndyst),
		integer(sampsize),
		integer(sampsize),
		sample = integer(sampsize),
		integer(k),
		med = integer(k),
		double(k),
		double(k),
		double(k),
		avdis = double(k),
		maxdis = double(k),
		ratdis = double(k),
		size = integer(k),
		obj = as.double(0),
		avsil = double(k),
		ttsil = as.double(0),
		silinf = matrix(0, sampsize, 4),
		stop = as.integer(0),
		double(sampsize),
		double(sampsize),
		double(sampsize),
		integer(sampsize),
		integer(sampsize),
		integer(sampsize),
		integer(sampsize),
		integer(sampsize),
		integer(sampsize))	
#give a warning when errors occured
	if(res$stop == 1)
		stop(message = "For each sample at least one object was found which could not be assigned to a cluster (because of missing values)."
			)
	if(res$stop == 2)
		stop(message = "Each of the random samples contains objects between which no distance can be computed."
			)
	sildim <- res$silinf[, 4]	
#adapt Fortran output to S-Plus:
#convert lower matrix, read by rows, to upper matrix, read by rows.
	disv <- res$dis[-1]
	disv[disv == -1] <- NA
	disv <- disv[upper.to.lower.tri.inds(sampsize)]
	class(disv) <- "dissimilarity"
	attr(disv, "Size") <- sampsize
	attr(disv, "Metric") <- metric	
	attr(disv, "Labels") <- dimnames(x)[[1]][res$sample]
#add labels to Fortran output
	res$med <- x[res$med,  ]
	res$clu <- matrix(res$clu, nrow = n, ncol = jp, byrow = T)[, 1]
	if(length(dimnames(x)[[1]]) != 0) {
		sildim <- dimnames(x)[[1]][sildim]
		res$sample <- dimnames(x)[[1]][res$sample]
		names(res$clu) <- dimnames(x)[[1]]
	}
#add dimnames to Fortran output
	clusinf <- cbind(res$size, res$maxdis, res$avdis, res$ratdis)
	dimnames(clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
		"isolation"))
	if(k != 1) {
		dimnames(res$silinf) <- list(sildim, c("cluster", "neighbor", 
			"sil_width", ""))
		clustering <- list(sample = res$sample, medoids = res$med, 
			clustering = res$clu, objective = res$obj, clusinfo = 
			clusinf, silinfo = list(width = res$silinf[, -4], 
			clus.avg.widths = res$avsil[1:k], avg.width = res$ttsil
			), diss = disv)
	}
	else {
		clustering <- list(sample = res$sample, medoids = res$med, 
			clustering = res$clu, objective = res$obj, clusinfo = 
			clusinf, diss = disv)
	}
	x2[x2 == valmisdat] <- NA
	clustering$data <- x2
	class(clustering) <- c("clara", "partition")
	attr(clustering, "Call") <- sys.call()
	clustering
}

"print.clara" <-
function(x, ...)
{
	cat("Best sample:\n")
	print(x$sample, quote = F, ...)
	cat("Medoids:\n")
	print(x$medoids, ...)
	cat("Clustering vector:\n")
	print(x$clustering, ...)
	cat("Objective function:\n")
	print(x$objective, ...)
	cat("\nAvailable arguments:\n")
	print(names(x), ...)
	invisible(x)
}

"summary.clara" <-
function(x)
{
	object <- x
	class(object) <- "summary.clara"
	object
}

"print.summary.clara" <-
function(x, ...)
{
	cat("Best sample:\n")
	print(x$sample, quote = F, ...)
	cat("Medoids:\n")
	print(x$medoids, ...)
	cat("Clustering vector:\n")
	print(x$clustering, ...)
	cat("Objective function:\n")
	print(x$objective, ...)
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
	cat("\nAvailable arguments:\n")
	print(names(x), ...)
	invisible(x)
}

