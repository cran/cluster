"pam" <-
function(x, k, diss = F, metric = "euclidean", stand = F)
{
	meanabsdev <- function(y)
	{
		mean(abs(y - mean(y, na.rm = T)), na.rm = T)
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
#check type of input vector
		if(is.na(min(x))) stop(message =
				"NA-values in the dissimilarity matrix not allowed."
				)
		if(data.class(x) != "dissimilarity") {
			if(!is.numeric(x) || size(x) == 0) stop(message =
				  "x is not of class dissimilarity and can not be converted to this class."
				  )	
#convert input vector to class "dissimilarity"
			class(x) <- "dissimilarity"
			attr(x, "Size") <- size(x)
			attr(x, "Metric") <- "unspecified"
		}
#adapt S-Plus dissimilarities to Fortran:
#convert upper matrix, read by rows, to lower matrix, read by rows.
		n <- attr(x, "Size")
		if((k < 1) || (k > n - 1))
			stop(message = 
				"The number of cluster should be at least 1 and at most n-1."
				)
		dv <- x[lower.to.upper.tri.inds(n)]	
#prepare arguments for the Fortran call
		dv <- c(0, dv)
		jp <- 1
		valmd <- double(1)
		jtmd <- integer(1)
		ndyst <- 0
		x2 <- double(n)
		jdyss <- 1
	}
	else {
#check type of input matrix
		if((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x, 
			data.class) == "numeric"))) stop(message = 
				"x is not a numeric dataframe or matrix.")
		x <- data.matrix(x)	
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
		if((k < 1) || (k > n - 1))
			stop(message = 
				"The number of cluster should be at least 1 and at most n-1."
				)
		jp <- ncol(x2)
		jtmd <- ifelse(is.na(rep(1, n) %*% x2), -1, 1)
		valmisdat <- min(x2, na.rm = T) - 0.5
		x2[is.na(x2)] <- valmisdat
		valmd <- rep(valmisdat, jp)
		jdyss <- 0
		dv <- double(1 + (n * (n - 1))/2)
	}
#call Fortran routine
	storage.mode(dv) <- "double"
	storage.mode(x2) <- "double"
	storage.mode(valmd) <- "double"
	storage.mode(jtmd) <- "integer"
	res <- .Fortran("pam",
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
		double(n),
		avsil = double(n),
		double(n),
		ttsil = as.double(0),
		med = integer(k),
		obj = double(2),
		clu = integer(n),
		clusinf = matrix(0, k, 5),
		silinf = matrix(0, n, 4),
		isol = integer(k))
	sildim <- res$silinf[, 4]
	if(diss) {
		disv <- x	
#add labels to Fortran output
		if(length(attr(x, "Labels")) != 0) {
			sildim <- attr(x, "Labels")[sildim]
			names(res$clu) <- attr(x, "Labels")
			res$med <- attr(x, "Labels")[res$med]
		}
	}
	else {
#give warning if some dissimilarities are missing.
		if(res$ok == -1)
		          stop(message =
				"No clustering performed, NA-values in the dissimilarity matrix.\n"
				)	
#adapt Fortran output to S-Plus:
#convert lower matrix, read by rows, to upper matrix, read by rows.
		disv <- res$dis[-1]
		disv[disv == -1] <- NA
		disv <- disv[upper.to.lower.tri.inds(n)]
		class(disv) <- "dissimilarity"
		attr(disv, "Size") <- nrow(x)
		attr(disv, "Metric") <- metric
		attr(disv, "Labels") <- dimnames(x)[[1]]	
#add labels to Fortran output
		res$med <- x[res$med,  ]
		if(length((dimnames(x)[[1]])) != 0) {
			sildim <- dimnames(x)[[1]][sildim]
			names(res$clu) <- dimnames(x)[[1]]
		}
	}
#add dimnames to Fortran output
	names(res$obj) <- c("build", "swap")
	res$isol <- factor(res$isol, levels = c(0, 1, 2), labels = c("no", "L",
		"L*"))
	names(res$isol) <- 1:k
	dimnames(res$clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
		"diameter", "separation"))
	if(k != 1) {
		dimnames(res$silinf) <- list(sildim, c("cluster", "neighbor", 
			"sil_width", ""))	#construct S-Plus object
		clustering <- list(medoids = res$med, clustering = res$clu, 
			objective = res$obj, isolation = res$isol, clusinfo = 
			res$clusinf, silinfo = list(widths = res$silinf[, -4], 
			clus.avg.widths = res$avsil[1:k], avg.width = res$ttsil
			), diss = disv)
	}
	else {
		clustering <- list(medoids = res$med, clustering = res$clu, 
			objective = res$obj, isolation = res$isol, clusinfo = 
			res$clusinf, diss = disv)
	}
	if(!diss) {
		x2[x2 == valmisdat] <- NA
		clustering$data <- x2
	}
	class(clustering) <- c("pam", "partition")
	attr(clustering, "Call") <- sys.call()
	clustering
}

"print.pam" <-
function(x, ...)
{
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

"summary.pam" <- 
function(x)
{
	object <- x
	class(object) <- "summary.pam"
	object
}

"print.summary.pam" <- 
function(x, ...)
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
	print(names(x$isolation[x$isolation == "L"]), quote = F, ...)
	cat("L*-clusters: ")
	print(names(x$isolation[x$isolation == "L*"]), quote = F, ...)
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
	cat("\nAvailable arguments:\n")
	print(names(x), ...)
	invisible(x)
}

