"mona" <-
function(x)
{
	levs <- function(y)
	{
		levels(as.factor(y))
	}
#check type of input matrix
	if(!is.matrix(x) && !is.data.frame(x))
	    stop("x must be a matrix or data frame.")
	if(!all(sapply(lapply(as.data.frame(x), levs), length) == 2))
		stop(message = "All variables must be binary (factor with 2 levels)."
			)
	n <- nrow(x)
	jp <- ncol(x)	
#change levels of input matrix
	x2 <- apply(as.matrix(x), 2, factor)
	x2[x2 == "1"] <- "0"
	x2[x2 == "2"] <- "1"
	x2[x2 == "NA"] <- "2"
##	x2 <- paste(x2, collapse = "")	
#call Fortran routine
##	storage.mode(x2) <- "character"
	storage.mode(x2) <- "integer"        
	res <- .Fortran("mona",
		as.integer(n),
		as.integer(jp),
		x2 = x2,
		error = as.integer(0),
		nban = integer(n),
		ner = integer(n),
		integer(n),
		lava = integer(n),
		integer(jp))	
#give a warning when errors occured
	if(res$error == 1)
		stop(message = "No clustering performed, an object was found with all values missing."
			)
	if(res$error == 2)
		stop(message = "No clustering performed, a variable was found with at least 50% missing values."
			)
	if(res$error == 3)
		stop(message = "No clustering performed, a variable was found with all non missing values identical."
			)
	if(res$error == 4)
		stop(message = "No clustering performed, all variables have at least one missing value."
			)
	res$x2 <- matrix(as.numeric(substring(res$x2, 1:nchar(res$x2), 1:nchar(
		res$x2))), n, jp)
	dimnames(res$x2) <- dimnames(x)	
#add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0)
		order.lab <- dimnames(x)[[1]][res$ner]
	if(length(dimnames(x)[[2]]) != 0) {
		lava <- as.character(res$lava)
		lava[lava != "0"] <- dimnames(x)[[2]][res$lava]
		lava[lava == "0"] <- "NULL"
		res$lava <- lava
	}
#construct S-Plus object
	clustering <- list(data = res$x2, order = res$ner, variable = res$lava[
		-1
		], step = res$nban[-1])
	if(exists("order.lab"))
		clustering$order.lab <- order.lab
	class(clustering) <- "mona"
	attr(clustering, "Call") <- sys.call()
	clustering
}

"print.mona" <-
function(x, ...)
{
	cat("Revised data:\n")
	print(x$data, quote = F, ...)
	cat("Order of objects:\n")
	if (length(x$order.lab) != 0)
		print(x$order.lab, quote = F, ...)
	else
		print(x$order, quote = F, ...)
	cat("Variable used:\n")
	print(x$variable, quote = F, ...)
	cat("Separation step:\n")
	print(x$step, ...)
	cat("\nAvailable arguments:\n")
	print(names(x), ...)
	invisible(x)
}

"summary.mona" <- 
function(x)
{
	object <- x
	class(object) <- "summary.mona"
	object
}

"print.summary.mona" <- 
function(x, ...)
{
	print.mona(x, ...)
	invisible(x)
}

