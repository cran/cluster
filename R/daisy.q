"daisy" <-
function(x, metric = "euclidean", stand = F, type = list())
{
	meanabsdev <- function(y)
	{
		mean(abs(y - mean(y, na.rm = T)), na.rm = T)
	}
	levs <- function(y)
	{
		levels(as.factor(y))
	}
#check type of input matrix
	if(!is.data.frame(x) && !is.numeric(x))
		stop(message = "x is not a dataframe or a numeric matrix.")
	if(!is.null(type$asymm)) if(!all(sapply(lapply(as.data.frame(x[, type$
			asymm]), levs), length) == 2))
			stop(message =
				"asymmetric binary variable has more than 2 levels."
				)
#transform variables and construct `type' vector
	type2 <- sapply(x, data.class)
	x <- data.matrix(x)
	x[, names(type2[type$ordratio])] <- codes(as.ordered(x[, names(type2[
		type$ordratio])]))
	x[, names(type2[type$logratio])] <- log10(x[, names(type2[type$logratio
		])])
	type2[type$asymm] <- "A"
	type2[type$ordratio] <- "O"
	type2[type2 == "numeric"] <- "I"
	type2[type2 == "ordered"] <- "O"
	type2[type2 == "factor"] <- "N"
#standardize, if necessary
	if(all(type2 == "I")) {
		if(stand) {
			x <- scale(x, scale = apply(x, 2, meanabsdev))
		}
		jdat <- 2
		if(metric == "manhattan")
			ndyst <- 2
		else ndyst <- 1
	}
	else {
		colmin <- apply(x, 2, min, na.rm = T)
		colextr <- apply(x, 2, max, na.rm = T) - colmin
		x <- scale(x, center = colmin, scale = colextr)
		jdat <- 1
		ndyst <- 0
	}
#	type2 <- paste(type2, collapse = "")
#put info about NAs in arguments for the Fortran call
	jtmd <- ifelse(is.na(rep(1, nrow(x)) %*% x), -1, 1)
	valmisdat <- min(x, na.rm = T) - 0.5
	x[is.na(x)] <- valmisdat
	valmd <- rep(valmisdat, ncol(x))
#call Fortran routine
	storage.mode(x) <- "double"
	storage.mode(valmd) <- "double"
	storage.mode(jtmd) <- "integer"
        type3 <- as.integer(match(type2, c('A','S','N','O','I','T')))
	res <- .Fortran("daisy",
		as.integer(nrow(x)),
		as.integer(ncol(x)),
		x,
		valmd,
		jtmd,
		as.integer(jdat),
		type3,
		as.integer(ndyst),
		dis = double(1 + (nrow(x) * (nrow(x) - 1))/2))
#adapt Fortran output to S-Plus:
#convert lower matrix, read by rows, to upper matrix, read by rows.
	disv <- res$dis[-1]
	disv[disv == -1] <- NA
	full <- matrix(0, nrow(x), nrow(x))
	full[!lower.tri(full, diag = T)] <- disv
	disv <- t(full)[lower.tri(full)]
#give warning if some dissimilarities are missimg
	if(is.na(min(disv))) attr(disv, "NA.message") <-
			"NA-values in the dissimilarity matrix !"
#construct S-Plus object
	class(disv) <- "dissimilarity"
	attr(disv, "Labels") <- dimnames(x)[[1]]
	attr(disv, "Size") <- nrow(x)
	attr(disv, "Metric") <- ifelse(ndyst == 0, "mixed", metric)
	disv
}

"print.dissimilarity" <-
function(x, ...)
{
	cat("Dissimilarities :\n")
	print(as.vector(x), ...)
	cat("\n")
	if(!is.null(attr(x, "na.message")))
		cat("Warning : ", attr(x, "NA.message"), "\n")
	cat("Metric : ", attr(x, "Metric"), "\n")
	cat("Number of objects : ", attr(x, "Size"), "\n")
	invisible(x)
}
