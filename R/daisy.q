#### $Id: daisy.q,v 1.13 2003/06/03 13:39:28 maechler Exp $
daisy <-
function(x, metric = c("euclidean","manhattan"), stand = FALSE, type = list())
{
    ## check type of input matrix
    if(length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x)))
	stop("x is not a dataframe or a numeric matrix.")
    if(length(type)) {
	tA <- type$asymm
	tS <- type$symm
	if((!is.null(tA) || !is.null(tS))) {
	    d.bin <- as.data.frame(x[,c(tA,tS), drop = FALSE])
	    if(!all(sapply(lapply(d.bin, function(y)
				  levels(as.factor(y))), length) == 2))
		stop("at least one binary variable has more than 2 levels.")
	    ## Convert factors to integer, such that ("0","1") --> (0,1):
	    if(any(is.f <- sapply(d.bin, is.factor)))
		d.bin[is.f] <- lapply(d.bin[is.f],
				      function(f) as.integer(as.character(f)))
	    if(!all(sapply(d.bin, function(y)
			   is.logical(y) ||
			   all(sort(unique(as.numeric(y[!is.na(y)])))==0:1))))
		stop("at least one binary variable has values other than 0,1, and NA")
	}
    }
    ## transform variables and construct `type' vector
    n <- dx[1]# nrow
    p <- dx[2]# ncol
    if(is.data.frame(x)) {
	type2 <- sapply(x, data.class)
	x <- data.matrix(x)
    } else type2 <- rep("numeric", p)
    if(length(type)) {
	if(!is.list(type)) stop("invalid `type'; must be named list")
	tT <- type$ ordratio
	tL <- type$ logratio
	x[, names(type2[tT])] <- unclass(as.ordered(x[, names(type2[tT])]))
	x[, names(type2[tL])] <- log10(		    x[, names(type2[tL])])
	type2[tA] <- "A"
	type2[tS] <- "S"
	type2[tT] <- "T" # was "O" (till 2000-12-14) accidentally !
    }
    type2[tI <- type2 %in% c("numeric", "integer") ] <- "I"
    if(any(tI) && any(iBin <- apply(x[,tI, drop = FALSE],2,
				    function(v) length(table(v)) == 2)))
	warning("binary variable(s) ", paste(which(tI)[iBin], collapse=","),
		" treated as interval scaled")

    type2[type2 == "ordered"] <- "O"
    type2[type2 == "factor"] <- "N"
    if(any(ilog <- type2 == "logical")) {
	warning("setting `logical' variable",if(sum(ilog)>1)"s " else " ",
		which(ilog), " to type `asymm'")
	type2[ilog] <- "A"
    }
    ## standardize, if necessary
    if(all(type2 == "I")) {
	if(stand)
	    x <- scale(x, scale = apply(x, 2,
			  function(y)
			  mean(abs(y - mean(y, na.rm = TRUE)), na.rm = TRUE)))
	jdat <- 2
	metric <- match.arg(metric)
	ndyst <- if(metric == "manhattan") 2 else 1
    }
    else { ## mixed case
	if(!missing(metric))
	    warning("`metric' is not used with mixed variables")
	colmin	 <- apply(x, 2, min, na.rm = TRUE)
	colrange <- apply(x, 2, max, na.rm = TRUE) - colmin
	x <- scale(x, center = colmin, scale = colrange)
	jdat <- 1
	ndyst <- 0
    }
    ##	type2 <- paste(type2, collapse = "")
    typeCodes <- c('A','S','N','O','I','T')
    type3 <- match(type2, typeCodes)# integer
    if(any(ina <- is.na(type3)))
	stop("invalid type", type2[ina],"  for column numbers", which(is.na))
    if((mdata <- any(inax <- is.na(x)))) { # TRUE if x[] has any NAs
	jtmd <- as.integer(ifelse(apply(inax, 2, any), -1, 1))
	## VALue for MISsing DATa
	valmisdat <- 1.1* max(abs(range(x, na.rm=TRUE)))
	x[inax] <- valmisdat
	valmd <- rep(valmisdat, p)
    }
    ## call Fortran routine
    storage.mode(x) <- "double"
    disv <- .Fortran("daisy",
		     n,
		     p,
		     x,
		     if(mdata)valmd else double(1),
		     if(mdata) jtmd else integer(1),
		     as.integer(jdat),
		     type3,		# vtype
		     as.integer(ndyst),
		     dis = double((n * (n - 1))/2),
                     NAOK = TRUE,# only to allow "+- Inf"
		     DUP = FALSE,
		     PACKAGE = "cluster")$dis
    ## adapt Fortran output to S:
    ## convert lower matrix, read by rows, to upper matrix, read by rows.
    disv[disv == -1] <- NA
    full <- matrix(0, n, n)
    full[!lower.tri(full, diag = TRUE)] <- disv
    disv <- t(full)[lower.tri(full)]
    ## give warning if some dissimilarities are missimg
    if(any(is.na(disv))) attr(disv, "NA.message") <-
	"NA-values in the dissimilarity matrix !"
    ## construct S object -- "dist" methods are *there* !
    class(disv) <- ..dClass
    attr(disv, "Labels") <- dimnames(x)[[1]]
    attr(disv, "Size") <- n
    attr(disv, "Metric") <- ifelse(!ndyst, "mixed", metric)
    if(!ndyst) attr(disv, "Types") <- typeCodes[type3]
    disv
}

print.dissimilarity <- function(x, ...)
{
    cat("Dissimilarities :\n")
    print(as.vector(x), ...)
    cat("\n")
    if(!is.null(attr(x, "na.message")))
	cat("Warning : ", attr(x, "NA.message"), "\n")
    cat("Metric : ", attr(x, "Metric"),
	if(!is.null(aT <- attr(x,"Types")))
	paste(";  Types =", paste(aT, collapse=", ")), "\n")
    cat("Number of objects : ", attr(x, "Size"), "\n", sep="")
    invisible(x)
}

summary.dissimilarity <- function(object, ...)
{
    sx <- summary(as.vector(object), ...)
    at <- attributes(object)
    r <- c(list(summ = sx, n = length(object)), at[names(at) != "class"])
    class(r) <- "summary.dissimilarity"
    r
}

print.summary.dissimilarity <- function(x, ...)
{
    cat(x$n, "dissimilarities, summarized :\n")
    print(x$summ, ...)
    cat("Metric : ", x $ Metric,
	if(!is.null(aT <- x $ Types))
	paste(";  Types =", paste(aT, collapse=", ")), "\n")
    cat("Number of objects : ", x $ Size, "\n", sep="")
    if(!is.null(x $ na.message))
	cat("Warning : ", x $ NA.message, "\n")
    invisible(x)
}
