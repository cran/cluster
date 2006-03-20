#### $Id: fanny.q 2886 2005-06-03 06:21:13Z maechler $
fanny <- function(x, k, diss = inherits(x, "dist"),
                  memb.exp = 2, metric = "euclidean", stand = FALSE,
                  maxit = 500, tol = 1e-15)
{
    if((diss <- as.logical(diss))) {
	## check type of input vector
	if(any(is.na(x))) stop(..msg$error["NAdiss"])
	if(data.class(x) != "dissimilarity") { # try to convert to
	    if(!is.null(dim(x))) {
		x <- as.dist(x)         # or give an error
	    } else {
		## possibly convert input *vector*
		if(!is.numeric(x) || is.na(n <- sizeDiss(x)))
		    stop(..msg$error["non.diss"])
		attr(x, "Size") <- n
	    }
	    class(x) <- dissiCl
	    if(is.null(attr(x,"Metric"))) attr(x, "Metric") <- "unspecified"
	}
	## prepare arguments for the Fortran call
	n <- attr(x, "Size")
	dv <- as.double(c(x, 0))
	jp <- 1
	mdata <- FALSE
	ndyst <- 0
	x2 <- double(n)
	jdyss <- 1
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
	jdyss <- 0
    }
    if((k <- as.integer(k)) < 1 || k > n%/%2 - 1)
	stop("'k' (number of clusters) must be in {1,2, .., n/2 -1}")
    if(length(memb.exp) != 1 || (memb.exp <- as.double(memb.exp)) < 1
       || memb.exp == Inf)
        stop("'memb.exp' must be a finite number > 1")
    if((maxit <- as.integer(maxit)[1]) < 0)
        stop("'maxit' must be non-negative integer")

    ## call Fortran routine
    storage.mode(x2) <- "double"
    res <- .Fortran("fanny",
		    as.integer(n),
		    as.integer(jp),
		    k,
		    x2,
		    dis = dv,
		    ok = as.integer(jdyss),
		    if(mdata)valmd else double(1),
		    if(mdata) jtmd else integer(jp),
		    as.integer(ndyst),
		    integer(n),
		    integer(n),
		    integer(n),
		    double(n),
		    p = matrix(0., n, k),
		    matrix(0., n, k),
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
		    silinf = matrix(0., n, 4),
                    memb.exp = memb.exp,
		    tol = as.double(tol),
                    maxit = maxit,
		    PACKAGE = "cluster")
    sildim <- res$silinf[, 4]
    if(diss) {
	disv <- x
	## add labels to Fortran output
	if(length(attr(x, "Labels")) != 0) {
	    sildim <- attr(x, "Labels")[sildim]
	    dimnames(res$p) <- list(attr(x, "Labels"), NULL)
	    names(res$clu) <- attr(x, "Labels")
	}
    }
    else {
	## give warning if some dissimilarities are missing.
	if(res$ok == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.")
	disv <- res$dis[ - (1 + (n * (n - 1))/2)]
	disv[disv == -1] <- NA
	class(disv) <- dissiCl
	attr(disv, "Size") <- nrow(x)
	attr(disv, "Metric") <- metric
	attr(disv, "Labels") <- dimnames(x)[[1]]
	## add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0) {
	    sildim <- dimnames(x)[[1]][sildim]
	    dimnames(res$p) <- list(dimnames(x)[[1]], NULL)
	    names(res$clu) <- dimnames(x)[[1]]
	}
    }
    if(!(converged <- res$maxit > 0)) {
        warning(sprintf(
            "FANNY algorithm has not converged in 'maxit' = %d iterations",
                        maxit))
    }

    ## add dimnames to Fortran output
    names(res$obj) <- c("iterations", "objective")
    res$coeff <- c(res$eda, res$edb)
    names(res$coeff) <- c("dunn_coeff", "normalized")

    r <- list(membership = res$p, coeff = res$coeff,
              memb.exp = memb.exp,
              clustering = res$clu,
              # 'obj*': also containing iterations for back compatibility:
              objective = c(res$obj, "tolerance" = res$tol),
              convergence = c(res$obj[1], converged = converged),
              diss = disv, call = match.call())
    if(k != 1) {
	dimnames(res$silinf) <- list(sildim,
				     c("cluster", "neighbor", "sil_width", ""))
	r$silinfo <- list(widths = res$silinf[, -4],
                          clus.avg.widths = res$avsil[1:k],
                          avg.width = res$ttsil)
    }
    if(!diss) {
	if(mdata) x2[x2 == valmisdat] <- NA
	r$data <- x2
    }
    class(r) <- c("fanny", "partition")
    r
}

## non-exported:
.print.fanny <- function(x, digits = getOption("digits"), ...) {
    cat("Fuzzy Clustering object of class 'fanny' :")
    print(formatC(cbind(" " = c("m.ship.expon." = x$memb.exp,
                        x$objective[c("objective", "tolerance")],
                        x$convergence)),
                  digits = digits),
	  quote = FALSE, ...)
    cat("Membership coefficients:\n");  print(x$membership, ...)
    cat("Coefficients:\n");             print(x$coeff, ...)
    cat("Closest hard clustering:\n");  print(x$clustering, ...)
}

print.fanny <- function(x, digits = getOption("digits"), ...)
{
    .print.fanny(x, digits = digits, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.fanny <- function(object, ...)
{
    class(object) <- "summary.fanny"
    object
}

print.summary.fanny <- function(x, digits = getOption("digits"), ...)
{
    .print.fanny(x, digits = digits, ...)
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
