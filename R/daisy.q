daisy <-
function(x, metric = c("euclidean","manhattan"), stand = FALSE, type = list())
{
    ## check type of input matrix
    if(!is.data.frame(x) && !is.numeric(x))
        stop("x is not a dataframe or a numeric matrix.")
    if(!is.null(tA <- type$asymm) &&
       !all(sapply(lapply(as.data.frame(x[,tA]),
                          function(y) levels(as.factor(y))), length) == 2))
        stop("asymmetric binary variable has more than 2 levels.")
    ## transform variables and construct `type' vector
    type2 <- sapply(x, data.class)
    x <- data.matrix(x)
    n <- nrow(x)
    if(length(type) > 0) {
        if(!is.list(type)) stop("invalid `type'; must be named list")
        tT <- type$ ordratio
        tL <- type$ logratio
        x[, names(type2[tT])] <- codes(as.ordered(x[, names(type2[tT])]))
        x[, names(type2[tL])] <- log10(           x[, names(type2[tL])])
        type2[type$asymm] <- "A"
        type2[tT] 	  <- "T" # was "O" (till 2000-12-14) accidentally !
    }
    type2[type2 == "numeric"] <- "I"
    type2[type2 == "ordered"] <- "O"
    type2[type2 == "factor"] <- "N"
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
        colmin   <- apply(x, 2, min, na.rm = TRUE)
        colrange <- apply(x, 2, max, na.rm = TRUE) - colmin
        x <- scale(x, center = colmin, scale = colrange)
        jdat <- 1
        ndyst <- 0
    }
    ## 	type2 <- paste(type2, collapse = "")
    ## put info about NAs in arguments for the Fortran call
    jtmd <- ifelse(is.na(rep(1, n) %*% x), -1, 1)
    valmisdat <- min(x, na.rm = TRUE) - 0.5
    x[is.na(x)] <- valmisdat
    valmd <- rep(valmisdat, ncol(x))
    ## call Fortran routine
    storage.mode(x) <- "double"
    storage.mode(valmd) <- "double"
    storage.mode(jtmd) <- "integer"
    type3 <- as.integer(match(type2, c('A','S','N','O','I','T')))
    res <- .Fortran("daisy",
                    as.integer(n),
                    as.integer(ncol(x)),
                    x,
                    valmd,
                    jtmd,
                    as.integer(jdat),
                    type3,
                    as.integer(ndyst),
                    dis = double(1 + (n * (n - 1))/2),
                    PACKAGE = "cluster")
    ## adapt Fortran output to S:
    ## convert lower matrix, read by rows, to upper matrix, read by rows.
    disv <- res$dis[-1]
    disv[disv == -1] <- NA
    full <- matrix(0, n, n)
    full[!lower.tri(full, diag = TRUE)] <- disv
    disv <- t(full)[lower.tri(full)]
    ## give warning if some dissimilarities are missimg
    if(is.na(min(disv))) attr(disv, "NA.message") <-
        "NA-values in the dissimilarity matrix !"
    ## construct S object -- "dist" methods are *there* !
    class(disv) <- c("dissimilarity", "dist")
    attr(disv, "Labels") <- dimnames(x)[[1]]
    attr(disv, "Size") <- n
    attr(disv, "Metric") <- ifelse(ndyst == 0, "mixed", metric)
    disv
}

print.dissimilarity <- function(x, ...)
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

summary.dissimilarity <- function(x, ...)
{
    cat(length(x), "dissimilarities, summarized :\n")
    print(sx <- summary(as.vector(x), ...))
    cat("\n")
    if(!is.null(attr(x, "na.message")))
        cat("Warning : ", attr(x, "NA.message"), "\n")
    cat("Metric : ", attr(x, "Metric"), "\n",
        "Number of objects : ", attr(x, "Size"), "\n", sep="")
    invisible(sx)
}
