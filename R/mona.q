mona <- function(x)
{
    levs <- function(y) levels(as.factor(y))

    ## check type of input matrix
    if(!is.matrix(x) && !is.data.frame(x))
        stop("x must be a matrix or data frame.")
    if(!all(sapply(lapply(as.data.frame(x), levs), length) == 2))
        stop(message = "All variables must be binary (factor with 2 levels).")
    n <- nrow(x)
    jp <- ncol(x)
    ##change levels of input matrix
    x2 <- apply(as.matrix(x), 2, factor)
    x2[x2 == "1"] <- "0"
    x2[x2 == "2"] <- "1"
    x2[x2 == "NA"] <- "2"
    ##	x2 <- paste(x2, collapse = "")
    ##	storage.mode(x2) <- "character"
    ## call Fortran routine
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
                    integer(jp),
                    PACKAGE = "cluster")
    ## give a warning when errors occured
    if(res$error != 0) {
        ch <- "No clustering performed,"
        switch(res$error,
               ## 1 :
               stop(paste(ch,"an object was found with all values missing.")),
               ## 2 :
               stop(paste(ch,"a variable was found with at least 50% missing values.")),
               ## 3 :
               stop(paste(ch,"a variable was found with all non missing values identical.")),
               ## 4 :
               stop(paste(ch,"all variables have at least one missing value."))
               )
    }
    res$x2 <- matrix(as.numeric(substring(res$x2,
                                          1:nchar(res$x2), 1:nchar(res$x2))),
                     n, jp)
    dimnames(res$x2) <- dnx <- dimnames(x)
    ## add labels to Fortran output
    if(length(dnx[[2]]) != 0) {
        lava <- as.character(res$lava)
        lava[lava != "0"] <- dnx[[2]][res$lava]
        lava[lava == "0"] <- "NULL"
        res$lava <- lava
    }
    ## construct "mona" object
    clustering <- list(data = res$x2, order = res$ner,
                       variable = res$lava[ -1 ], step = res$nban[-1])
    if(length(dnx[[1]]) != 0)
        clustering$order.lab <- dnx[[1]][res$ner]
    class(clustering) <- "mona"
    attr(clustering, "Call") <- sys.call()
    clustering
}

print.mona <- function(x, ...)
{
    cat("Revised data:\n")
    print(x$data, quote = FALSE, ...)
    cat("Order of objects:\n")
    print(if (length(x$order.lab) != 0) x$order.lab else x$order,
          quote = FALSE, ...)
    cat("Variable used:\n")
    print(x$variable, quote = FALSE, ...)
    cat("Separation step:\n")
    print(x$step, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.mona <- function(object, ...)
{
    class(object) <- "summary.mona"
    object
}

print.summary.mona <- function(x, ...)
{
    print.mona(x, ...)
    invisible(x)
}

