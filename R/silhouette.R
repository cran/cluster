silhouette <- function(x, ...) UseMethod("silhouette")

## Accessor and more:
silhouette.partition <- function(x, ...) {
    r <- x$silinfo$widths
    if(is.null(r))
        stop("invalid partition object")
    attr(r, "Ordered") <- TRUE # (cluster <increasing>, s.i <decreasing>)
    attr(r, "call") <- x$call
    class(r) <- "silhouette"
    r
}

silhouette.default <- function(x, dist, dmatrix, ...) {
    cll <- match.call()
    if(!is.null(cl <- x$clustering)) x <- cl
    n <- length(x)
    if(!all(x == round(x))) stop("`x' must only have integer codes")
    k <- length(clid <- sort(unique(x)))
    if(k <= 1 || k >= n)
        return(NA)
    ## check dist/dmatrix
    if(missing(dist)) {
        if(missing(dmatrix))
            stop("Need either a dissimilarity `dist' or diss.matrix `dmatrix'")
        if(is.null(dm <- dim(dmatrix)) || length(dm) != 2 || !all(n == dm))
            stop("`dmatrix' is not a dissimilarity matrix compatible to `x'")
    } else { # `dist'
        dist <- as.dist(dist) # hopefully
        if(n != attr(dist, "Size"))
            stop("clustering `x' and dissimilarity `dist' are incompatible")
        dmatrix <- as.matrix(dist)# so we can apply(.) below
    }
    wds <- matrix(NA, n,3, dimnames =
                  list(names(x), c("cluster","neighbor","sil_width")))
    for(j in 1:k) { # j-th cluster:
        Nj <- sum(iC <- x == clid[j])
        wds[iC, "cluster"] <- clid[j]
        ## minimal distances to points in all other clusters:
        diC <- rbind(apply(dmatrix[!iC, iC, drop = FALSE], 2,
                           function(r) tapply(r, x[!iC], mean)))# (k-1) x Nj
        ## max.col() breaks ties at random;  rather do not want random
        ##   behavior of silhouette, (but rather "pam" compatible one):
        minC <- apply(diC, 2, which.min)
        ## FIXME minC <- max.col(-t(diC))
        ## FIXME : extend max.col(*, ties.method = "min") {or similar} !
        wds[iC,"neighbor"] <- clid[-j][minC]
        s.i <- if(Nj > 1) {
            a.i <- colSums(dmatrix[iC, iC])/(Nj - 1) # length(a.i)= Nj
            b.i <- diC[cbind(minC, seq(along = minC))]
            ifelse(a.i != b.i, (b.i - a.i) / pmax(b.i, a.i), 0)
        } else 0
        wds[iC,"sil_width"] <- s.i
    }
    attr(wds, "Ordered") <- FALSE
    attr(wds, "call") <- cll
    class(wds) <- "silhouette"
    wds
}


sortSilhouette <- function(object, ...) {
    if(attr(object,"Ordered")) return(object)
    ## Else :
    if(is.null(n <- nrow(object)) || n < 1)
        stop("invalid silhouette structure")
    if(is.null(rownames(object)))
        rownames(object) <- as.character(1:n)
    ## k <- length(clid <- sort(unique(cl <- object[,"cluster"])))# cluster ID s
    cl <- object[,"cluster"]
    r <- object[iOrd <- order(cl, - object[,"sil_width"]) , , drop = FALSE]
    attributes(r) <- attributes(object) # but:
    attr(r,"iOrd") <- iOrd # the ordering
    attr(r,"Ordered") <- TRUE
    r
}

summary.silhouette <- function(object, FUN = mean, ...)
{
    if(ncol(object) != 3) stop("invalid `silhouette' object")
    cl <- object[, "cluster"]
    si <- object[, "sil_width"]
    r <- list(si.summary = summary(si, ...),
              clus.avg.widths = tapply(si, cl, FUN),
              clus.sizes = table(cl),
              avg.width = FUN(si),
              call = attr(object,"call"),
              Ordered = attr(object,"Ordered"))
    class(r) <- "summary.silhouette"
    r
}

print.summary.silhouette <- function(x, ...)
{
    k <- length(csiz <- x$clus.sizes)
    cat("Silhouette of", sum(csiz), "units in", k, "clusters",
        if(!is.null(x$call)) paste("from", deparse(x$call)),
        ":\nCluster sizes and average silhouette widths:\n")
    cwid <- x$clus.avg.widths
    names(cwid) <- csiz
    print(cwid, ...)
    cat("Individual silhouette widths:\n")
    print(x$si.summary, ...)
    invisible(x)
}


## This was the internal function silhouPlot() in plot.partition() :
plot.silhouette <-
    function(x, nmax.lab = 40, max.strlen = 5,
	     main = NULL, sub = NULL,
	     xlab = expression("Silhouette width " * s[i]),
	     col = "gray", do.col.sort = length(col) > 1,
	     border = 0, cex.names = par("cex.axis"),
	     do.n.k = TRUE, do.clus.stat = TRUE, ...)
{
    if(!is.matrix(x) || ncol(x) != 3)
	stop("No valid silhouette information (#{clusters} =? 1)")
    n <- nrow(x)
    x <- sortSilhouette(x)
    s <- rev(x[, "sil_width"])
    space <- c(0, rev(diff(cli <- x[, "cluster"])))
    space[space != 0] <- 0.5 # gap between clusters
    names <- if(n < nmax.lab)
	substring(rev(rownames(x)), 1, max.strlen)
    if(is.null(main)) {
	main <- "Silhouette plot"
	if(!is.null(cll <- attr(x,"call"))) { # drop initial "silhouette":
	    if(!is.na(charmatch("silhouette", deparse(cll[[1]]))))
		cll[[1]] <- as.name("FF")
	    main <-  paste(main, "of", sub("^FF","", deparse(cll)))
	}
    }
    smry <- summary(x)
    k <- length(nj <- smry$clus.sizes) # k clusters
    if(is.null(sub))
	sub <- paste("Average silhouette width : ",
		     round(smry$avg.width, digits = 2))
    if(do.col.sort && (lc <- length(col)) > 1) {
	if(lc == k)# cluster wise coloring
	    col <- col[cli]
	else if(lc != n)
	    col <- rep(col, length = n)
	col <- rev(col[attr(x, "iOrd")])
    }
    y <- barplot(s, space = space, names = names, xlab = xlab,
		 xlim = c(min(0, min(s)), 1),
		 horiz = TRUE, las = 1, mgp = c(2.5, 1, 0),
		 col = col, border = border, cex.names = cex.names, ...)
    title(main = main, sub = sub, adj = 0)
    if(do.n.k) {
	mtext(paste("n =", n),	adj = 0)
	mtext(substitute(k ~~ "clusters" ~~ C[j], list(k=k)), adj= 1)
    }
    if(do.clus.stat) {
	mtext(expression(paste(j," :  ", n[j]," | ", ave[i %in% Cj] ~~ s[i])),
	      adj = 1.04, line = -1.2)
	y <- rev(y)
	for(j in 1:k) {
	    yj <- mean(y[cli == j])
	    text(1, yj,
                 paste(j,":  ", nj[j]," | ",
                       format(smry$clus.avg.widths[j], digits = 1, nsmall = 2)),
		 xpd = NA, adj = 0.8)
	}
    }
}
