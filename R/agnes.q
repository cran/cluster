#### $Id: agnes.q,v 1.17 2004/03/11 16:26:40 maechler Exp maechler $
agnes <- function(x, diss = inherits(x, "dist"), metric = "euclidean",
		  stand = FALSE, method = "average",
                  keep.diss = n < 100, keep.data = !diss)
{
    METHODS <- c("average", "single", "complete", "ward", "weighted")
    ## hclust has more;  1    2         3           4       5
    meth <- pmatch(method, METHODS)
    if(is.na(meth))
	stop("invalid clustering method")
    if(meth == -1)
	stop("ambiguous clustering method")
    if((diss <- as.logical(diss))) {
	## check type of input vector
	if(any(is.na(x))) stop(..msg$error["NAdiss"])
	if(data.class(x) != "dissimilarity") { # try to convert to
	    if(!is.null(dim(x))) {
		x <- as.dist(x) # or give an error
	    } else {
		## possibly convert input *vector*
		if(!is.numeric(x) || is.na(n <- sizeDiss(x)))
		    stop(..msg$error["non.diss"])
		attr(x, "Size") <- n
	    }
	    class(x) <- dissiCl
	    if(is.null(attr(x,"Metric"))) attr(x, "Metric") <- "unspecified"
	}
	n <- attr(x, "Size")
	dv <- x[lower.to.upper.tri.inds(n)]
	## prepare arguments for the Fortran call
	dv <- c(0., dv)# "double", 1st elem. "only for Fortran" (?)
	jp <- 1
	mdata <- FALSE
	ndyst <- 0
	x2 <- double(1)
	jdyss <- 1 # distances on _input_
    }
    else {
	## check input matrix and standardize, if necessary
	x <- data.matrix(x)
	if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")
	x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
        storage.mode(x2) <- "double"
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
	jdyss <- 0 # distances to be computed
    }
    if(keep.diss) jdyss <- jdyss + 10
    ## call Fortran routine
    res <- .Fortran("twins",
		    as.integer(n),
		    as.integer(jp),
		    x2,
		    dv,
		    dis = double(if(keep.diss) length(dv) else 1),
		    ok = as.integer(jdyss),
		    if(mdata) valmd else double(1),
		    if(mdata) jtmd else integer(1),
		    as.integer(ndyst),
		    as.integer(1),# jalg = 1 <==> AGNES
		    meth,# integer
		    integer(n),
		    ner = integer(n),
		    ban = double(n),
		    ac = as.double(0),
		    merge = matrix(0:0, n - 1, 2), # integer
                    DUP = FALSE,
		    PACKAGE = "cluster")
    if(!diss) {
	##give warning if some dissimilarities are missing.
	if(res$ok == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.\n" )
        if(keep.diss) {
            ## adapt Fortran output to S:
            ## convert lower matrix,read by rows, to upper matrix, read by rows.
            disv <- res$dis[-1]
            disv[disv == -1] <- NA
            disv <- disv[upper.to.lower.tri.inds(n)]
            class(disv) <- dissiCl
            attr(disv, "Size") <- nrow(x)
            attr(disv, "Metric") <- metric
            attr(disv, "Labels") <- dimnames(x)[[1]]
        }
	##add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0)
	    order.lab <- dimnames(x)[[1]][res$ner]
    }
    else {
        if(keep.diss) disv <- x
	##add labels to Fortran output
	if(length(attr(x, "Labels")) != 0)
	    order.lab <- attr(x, "Labels")[res$ner]
    }
    clustering <- list(order = res$ner, height = res$ban[-1], ac = res$ac,
		       merge = res$merge, diss = if(keep.diss)disv,
                       call = match.call(), method = METHODS[meth])
    if(exists("order.lab"))
	clustering$order.lab <- order.lab
    if(keep.data && !diss) {
	if(mdata) x2[x2 == valmisdat] <- NA
	clustering$data <- x2
    }
    class(clustering) <- c("agnes", "twins")
    clustering
}

summary.agnes <- function(object, ...)
{
    class(object) <- "summary.agnes"
    object
}

print.agnes <- function(x, ...)
{
    cat("Call:	", deparse(x$call),
	"\nAgglomerative coefficient: ", format(x$ac, ...),
	"\nOrder of objects:\n")
    print(if(length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Height (summary):\n");		print(summary(x$height), ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

print.summary.agnes <- function(x, ...)
{
    ## a bit more than print.agnes() ..
    cat("Object of class `agnes' from call:\n", deparse(x$call),
	"\nAgglomerative coefficient: ", format(x$ac, ...),
	"\nOrder of objects:\n")
    print(if(length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Merge:\n");			print(x$merge, ...)
    cat("Height:\n");			print(x$height, ...)
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

