### $Id: plotpart.q,v 1.12 2002/01/23 18:20:33 maechler Exp $
plot.partition <-
function(x, ask = FALSE, which.plots = NULL,
         nmax.lab = 40, max.strlen = 5,
         cor = TRUE, stand = FALSE, lines = 2,
         shade = FALSE, color = FALSE, labels = 0, plotchar = TRUE,
         span = TRUE, xlim = NULL, ylim = NULL, ...)
{
    silhouPlot <- function(x, nmax.lab, max.strlen) {
        if(length(x$silinfo) == 0)
            stop("No silhouette plot available when the number of clusters equals 1." )
        s <- rev(x$silinfo[[1]][, 3])
        space <- c(0, rev(diff(x$silinfo[[1]][, 1])))
        space[space != 0] <- 0.5
        names <- if(length(s) < nmax.lab)
            substring(rev(dimnames(x$silinfo[[1]])[[1]]), 1, max.strlen)
        barplot(s, space = space, names = names,
                xlab = "Silhouette width",
                xlim = c(min(0, min(s)), 1), horiz = TRUE,
                mgp = c(2.5, 1, 0), ...)
        title(main = paste("Silhouette plot of ",
              deparse(x$call)),
              sub = paste("Average silhouette width : ",
              round(x$ silinfo$avg.width, digits = 2)), adj = 0)
    }

    if(is.null(which.plots) && !ask)
        which.plots <- 1:2
    if(ask && is.null(which.plots)) { ## Use `menu' ..
        tmenu <- paste("plot ", ## choices :
                       c("All", "Clusplot", "Silhouette Plot"))
        do.all <- FALSE
        repeat {
            if(!do.all)
                pick <- menu(tmenu, title =
                             "\nMake a plot selection (or 0 to exit):\n") + 1
            switch(pick,
                   return(invisible())# 0 -> exit loop
                   ,
                   do.all <- TRUE# 1 : All
                   ,
                   clusplot(x, cor = cor, stand = stand, lines = lines,
                            shade = shade, color = color, labels = labels,
                            plotchar = plotchar, span = span,
                            xlim = xlim, ylim = ylim, ...)
                   ,
                   silhouPlot(x, nmax.lab, max.strlen)
                   )
            if(do.all) { pick <- pick + 1; do.all <- pick <= length(tmenu) + 1}
        }
    }
    else {
        ask <- prod(par("mfcol")) < length(which.plots) && dev.interactive()
        if(ask) { op <- par(ask = TRUE); on.exit(par(op)) }
        for(i in which.plots)
        switch(i,
               clusplot(x, cor = cor, stand = stand, lines = lines,
                        shade = shade, color = color, labels = labels,
                        plotchar = plotchar, span = span,
                        xlim = xlim, ylim = ylim, ...)
               ,
               silhouPlot(x, nmax.lab, max.strlen)
               )
    }
    invisible()
}

clusplot <- function(x, ...) UseMethod("clusplot")


clusplot.default <-
function(x, clus, diss = FALSE, cor = TRUE, stand = FALSE, lines = 2,
         shade = FALSE, color = FALSE, labels = 0, plotchar = TRUE,
         col.p = "dark green", # was 5 (= shaded col)
         col.txt = col.p,
         span = TRUE, xlim = NULL, ylim = NULL,
         main = paste("CLUSPLOT(", deparse(substitute(x)),")"),
         verbose = getOption("verbose"),
         ...)
{
    if(paste(R.version$major, R.version$minor, sep=".") < 1.5) {
        ## a simplified (add = T) version of R 1.5's cmdscale():
        cmdscale <- function (d, k = 2, add = TRUE, ...) {
            if (any(is.na(d)))
                stop("NA values not allowed in d")
            if (is.null(n <- attr(d, "Size"))) {
                d <- as.matrix(d)
                x <- d^2
                if ((n <- nrow(x)) != ncol(x))
                    stop("Distances must be result of dist or a square matrix")
            }
            else {
                x <- matrix(0, n, n)
                if(add) d0 <- x
                x[row(x) > col(x)] <- d^2
                x <- x + t(x)
                if(add) {
                    d0[row(x) > col(x)] <- d
                    d <- d0 + t(d0)
                }
            }
            storage.mode(x) <- "double"
            x <- .C("dblcen", x=x, as.integer(n), PACKAGE="mva")$x
            if(add) { ## solve the additive constant problem
                i2 <- n + (i <- 1:n)
                Z <- matrix(0, 2*n, 2*n)
                Z[cbind(i2,i)] <- -1
                Z[ i, i2] <- -x
                Z[i2, i2] <- .C("dblcen", x= 2*d, as.integer(n),PACKAGE="mva")$x
                e <- La.eigen(Z,symmetric = FALSE, only.val = TRUE)$values
                add.c <- max(Re(e))
                x <- matrix(double(n*n), n, n)
                non.diag <- row(d) != col(d)
                x[non.diag] <- (d[non.diag] + add.c)^2
            }
            e <- La.eigen(-x/2, symmetric = TRUE)
            ev <- e$values[1:k]
            points <- e$vectors[, 1:k] %*% diag(sqrt(ev), k)
            rn <- if(is.matrix(d)) rownames(d) else names(d)
            dimnames(points) <- list(rn, NULL)
            evalus <- e$values[-n]
            list(points = points, eig = ev, ac = if(add) add.c else 0,
                 GOF = sum(ev)/c(sum(abs(evalus)),
                                 sum(evalus[evalus > 0])))
        }
    }## cmdscale() -- if R version < 1.5

    clas.snijpunt <- function(x, loc, m, n, p)
    {
        if(     loc[n, m] <= x[1, m] && x[1, m] <= loc[p, m]) x[1, ]
        else if(loc[n, m] <= x[2, m] && x[2, m] <= loc[p, m]) x[2, ]
        else NA
    }
    coord.snijp1 <- function(x, gemid)
        x[2, 2] - 2 * x[1, 2] * gemid + x[1, 1] * gemid^2
    coord.snijp2 <- function(x, d2, y)
        ((x[1, 1] * x[2, 2] - x[1, 2]^2) * d2)/y
    coord.snijp3 <- function(xx, y, gemid)
    {
        sy <- sqrt(y)
        sy <- c(sy, -sy)
        cbind(xx[1] + sy,
              xx[2] + gemid*sy)
    }

    ## BEGIN ----

    (main)# eval
    if(is.data.frame(x))
        x <- data.matrix(x)
    if(!is.numeric(x))
        stop("x is not numeric")

    if(diss) {
        if(is.na(min(x)))
            stop(message = "NA-values in x are not allowed.")
        if((data.class(x)) != "dissimilarity") {
            if(is.na(sizeDiss(x))) {
                if((n <- nrow(x)) != ncol(x))
                    stop("Distances must be result of dist or a square matrix.")
                if(all.equal(x, t(x)) != TRUE)
                    stop("the square matrix is not symmetric.")
                labels1 <-
                    if(length(dimnames(x)[[1]]) == 0) 1:nrow(x)
                    else dimnames(x)[[1]]
            }
            else {
                if(!is.vector(x)) {
                    labels1 <- attr(x, "Labels") # possibly NULL
                    x <- as.matrix(x)
                    if((n <- nrow(x)) == ncol(x) &&
                       all.equal(x, t(x)) == TRUE) {
                        labels1 <-
                            if(length(dimnames(x)[[1]]) == 0) 1:nrow(x)
                            else dimnames(x)[[1]]
                    }
                    else {
                        if(is.null(labels1))
                            labels1 <- 1:sizeDiss(x)
                        attr(x, "Size") <- sizeDiss(x)
                    }
                }
                else {
                    attr(x, "Size") <- sizeDiss(x)
                    labels1 <- 1:sizeDiss(x)
                }
            }
        }
        else {
            labels1 <-
                if(length(attr(x, "Labels")) == 0)
                    1:attr(x, "Size")
                else attr(x, "Labels")
        }
        x1 <- cmdscale(x, k = 2, eig = TRUE, add = TRUE)
        if(x1$ac < 0)
            x1 <- cmdscale(x, k = 2, eig = TRUE)
        var.dec <- x1$GOF[2] # always in [0,1]
        x1 <- x1$points
    }
    else { ## Not (diss)
        if(!is.matrix(x)) stop("x is not a data matrix")
        if(is.na(min(x))) { ## any(is.na(x))
            y <- is.na(x)
            y1 <- apply(y, 1, sum)
            y2 <- apply(y, 2, sum)
            if(any(y1 == ncol(x)))
                stop("one or more objects contain only missing values")
            if(any(y2 == nrow(x)))
                stop("one or more variables contain only missing values")
            x <- apply(x, 2, function(x)
                   { x[is.na(x)] <- median(x, na.rm = TRUE); x } )
            cat("Missing values were displaced by the median of the corresponding variable(s)\n")

        }

        labels1 <-
            if(length(dimnames(x)[[1]]) == 0) 1:nrow(x)
            else dimnames(x)[[1]]

        if(ncol(x) == 1) {
            hulp <- rep(0, length(x))
            x1 <- matrix(c(t(x), hulp), ncol = 2)
            var.dec <- 1
        }
        else {
            prim.pr <- princomp(x, scores = TRUE, cor = ncol(x) != 2)
            var.dec <- cumsum(prim.pr$sdev^2/sum(prim.pr$ sdev^2))[2]
            x1 <- prim.pr$scores
            x1 <- cbind(x1[, 1], x1[, 2])
        }
    }

    ## --- The 2D space is setup and points are in x1[,]  (aantal x 2) ---

    clus <- as.vector(clus)
    if(length(clus) != length(x1[, 1]))
        stop("The clustering vector has not the good length")
    clus <- as.factor(clus)
    if(sum(is.na(clus)) != 0)
        stop("NA-values are not allowed in clustering vector")
    if(stand)
        x1 <- scale(x1)

    rangx <- range(x1[, 1])
    rangy <- range(x1[, 2])
    minx <- rangx[1]
    maxx <- rangx[2]
    miny <- rangy[1]
    maxy <- rangy[2]
    levclus <- levels(clus)
    n <- length(levclus) # the number of clusters
    z <- A <- vector("list", n)
    maxima <- loc <- matrix(0, nrow = n, ncol = 2)
    d2 <- verhoud <- numeric(n)
    verhouding <- 0
    ## num1 .. num6 : all used only once -- there are more constants anyway
    num3 <- 90
    num6 <- 70

    for(i in 1:n) { ##-------------  i-th cluster  --------------
	x <- x1[clus == levclus[i],, drop = FALSE ]
        aantal <- nrow(x) # number of observations in cluster [i]
        cov <- var(if(aantal == 1) {
                     if(verbose)
                         cat("cluster",i," has only one observation ..\n")
                     rbind(x, c(0, 0))
                   } else x)
        x.1 <- range(x[, 1])
        y.1 <- range(x[, 2])
        notrank2 <- qr(cov, tol = 0.001)$rank != 2
        if(!span && notrank2) {
            d2[i] <- 1
            if((abs(diff(x.1)) > (diff(rangx)/70)) ||
               (abs(diff(y.1)) > (diff(rangy)/50))) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 +
                          (loc[i, 2] - y.1[1])^2)
                a <- a + 0.05 * a
                num2 <- 40
                if(abs(diff(x.1)) > diff(rangx)/70 ) {
                    ind1 <- which.max(x[,1])
                    ind2 <- which.min(x[,1])
                    q <- atan((x[ind1, 2] - x[ind2, 2])/
                              (x[ind1, 1] - x[ind2, 1]))
                    b <-
                        if(diff(rangy) == 0)
                            1
                        else if(abs(diff(y.1)) > diff(rangy)/50)
                            diff(y.1)/10 ## num1 <- 10
                        else diff(rangy)/num2
                }
                else {
                    b <- if(diff(rangx) == 0) 1 else diff(rangx)/num2
                    q <- pi/2
                }
                D <- diag(c(a^2, b^2))
                R <- rbind(c(cos(q), -sin(q)),
                           c(sin(q),  cos(q)))
                A[[i]] <- (R %*% D) %*% t(R)
            }
            else {
                a <- diff(rangx)/num3
                b <- diff(rangy)/num6
                if(a == 0) a <- 1
                if(b == 0) b <- 1
                A[[i]] <- diag(c(a^2, b^2))
                loc[i, ] <- x[1, ]
            }
            oppervlak <- pi * a * b
        }
        else if(span && notrank2) {
            d2[i] <- 1
            if(sum(x[, 1] != x[1, 1]) != 0 ||
               sum(x[, 2] != x[1, 2]) != 0) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2,
                              y.1[1] + diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 +
                          (loc[i, 2] - y.1[1])^2)
                if(any(x[, 1] != x[1, 1])) {
                    ind1 <- which.max(x[,1])
                    ind2 <- which.min(x[,1])
                    q <- atan((x[ind1, 2] - x[ind2, 2])/
                              (x[ind1, 1] - x[ind2, 1]))
                }
                else {
                    q <- pi/2
                }
                b <- 1e-7
                D <- diag(c(a^2, b^2))
                R <- rbind(c(cos(q), -sin(q)),
                           c(sin(q),  cos(q)))
                A[[i]] <- (R %*% D) %*% t(R)
            }
            else {
                a <- diff(rangx)/num3
                b <- diff(rangy)/num6
                if(a == 0) a <- 1
                if(b == 0) b <- 1
                A[[i]] <- diag(c(a^2, b^2))
                loc[i, ] <- x[1, ]
            }
            oppervlak <- pi * a * b

        }
        else { ## rank2
            if(!span) {
                loc[i, ] <- apply(x, 2, mean)
                d2[i] <- max(mahalanobis(x, loc[i, ], cov))
                ## * (1+ 0.01)^2  --- dropped factor for back-compatibility
            }
            else { ## span and rank2
                if(verbose)
                    cat("span & rank2 : calling \"spannel\" ..\n")
                k <- as.integer(2)
                res <- .Fortran("spannel",
                                aantal,
                                ndep= k,
                                dat = cbind(1., x),
                                sqdist = double(aantal),
                                l1 = double((k+1) ^ 2),
                                double(k),
                                double(k),
                                prob = double(aantal),
                                double(k+1),
                                eps = as.double(0.01),## convergence tol.
                                maxit = as.integer(5000),
                                ierr = as.integer(0),
                                PACKAGE = "cluster")
                if(res$ierr != 0)
                    ## MM : exactmve not available here !
                    cat("Error in Fortran routine for the spanning ellipsoid,",
                        "\n rank problem??\n", sep="")

                cov <- cov.wt(x, res$prob)
                loc[i, ] <- cov$center
                ## NB: cov.wt() in R has extra wt[] scaling; revert here:
                cov <- cov$cov * (1 - sum(cov$wt^2))
                d2[i] <- weighted.mean(res$sqdist, res$prob)

                if(verbose)
                    cat("ellipse( A= (", format(cov[1,]),"*", format(cov[2,2]),
                        "),\n\td2=", format(d2[i]),
                        ", loc[]=", format(loc[i, ]), ")\n")
            }
            A[[i]] <- cov
            ## oppervlak (flam.)  =  area (Engl.)
            oppervlak <- pi * d2[i] * sqrt(cov[1, 1] * cov[2, 2] - cov[1, 2]^2)
        }

        z[[i]] <- ellipsoidPoints(A[[i]], d2[i], loc[i, ], n= 201)
        maxima[i, ] <- z[[i]][201, ]
        rx <- range(z[[i]][, 1])
        ry <- range(z[[i]][, 2])
        minx <- min(minx, rx[1])
        maxx <- max(maxx, rx[2])
        miny <- min(miny, ry[1])
        maxy <- max(maxy, ry[2])
        verhoud[i] <- aantal/oppervlak
        if(verhoud[i] < 1e7)
            verhouding <- verhouding + verhoud[i]
    } ## end for( i-th cluster )

    if(verhouding == 0)
        verhouding <- 1
    ## num4 <- 37 ; num5 <- 3 --- but `41' is another constant
    density <- 3 + (verhoud * 37)/verhouding
    density[density > 41] <- 41
    if (span) {
        if (rangx[1] == rangx[2]) { ## diff(rangx)==0 : x-coords all the same
            minx <- x1[1, 1] - 1
            maxx <- x1[1, 1] + 1
        }
        if (rangy[1] == rangy[2]) { ## diff(rangy)==0 : y-coords all the same
            miny <- x1[1, 2] - 1
            maxy <- x1[1, 2] + 1
        }
    }
    if(!is.null(xlim)) {
        if(xlim[1] < minx) minx <- xlim[1]
        if(xlim[2] > maxx) maxx <- xlim[2]
    }
    if(!is.null(ylim)) {
        if(ylim[1] < miny) miny <- ylim[1]
        if(ylim[2] > maxy) maxy <- ylim[2]
    }

    ## --- Now plotting starts ---

    plot(x1[, 1], x1[, 2], xlim = c(minx, maxx), ylim = c(miny, maxy),
         xlab = "Component 1", ylab = "Component 2",
         main = main,
         type = if(plotchar) "n" else "p", # if(plotchar) add points later
         col = col.p, ...)
    title(sub = paste("These two components explain",
          round(100 * var.dec, digits = 2), "% of the point variability."),
          adj = 0)

    if(color) {
        color1 <- c(2, 4, 6, 3)
        i.verh <- order(verhoud)
        jInd <- if(n > 4) pam(verhoud[i.verh], 4)$clustering else 1:n
        for(i in 1:n) {
            k <- i.verh[i]
            polygon(z[[k]], density = if(shade) density[k] else 0,
                    col = color1[jInd[i]], ...)
        }
    }
    else {
        for(i in 1:n)
            polygon(z[[i]], density = if(shade) density[i] else 0,
                    col = 5, ...)
    }

    ## points after polygon in order to write ON TOP:
    if(plotchar) {
        karakter <- 1:19
        for(i in 1:n) {
            x <- x1[clus == levclus[i],  , drop = FALSE]
            kar <- 1+(i-1) %% 19
            points(x[, 1], x[, 2], pch = karakter[kar], col = col.p, ...)
        }
    }

    if((lines == 1 || lines == 2) && n > 1) {
        ## Draw lines between all pairs of the  n  cluster (centers)
        afstand <- matrix(0, ncol = n, nrow = n)
        for(i in 1:(n - 1)) {
            for(j in (i + 1):n) {
                gemid <- (loc[j, 2] - loc[i, 2])/(loc[j, 1] - loc[i, 1])
                s0 <- coord.snijp1(A[[i]], gemid)
                b0 <- coord.snijp2(A[[i]], d2[i], s0)
                snijp.1 <- coord.snijp3(loc[i,], y=b0, gemid)
                s1 <- coord.snijp1(A[[j]], gemid)
                b1 <- coord.snijp2(A[[j]], d2[j], s1)
                snijp.2 <- coord.snijp3(loc[j,], y=b1, gemid)
                if(loc[i, 1] != loc[j, 1]) {
                    if(loc[i, 1] < loc[j, 1]) {
                        punt.1 <- clas.snijpunt(snijp.1, loc, 1, i, j)
                        punt.2 <- clas.snijpunt(snijp.2, loc, 1, i, j)
                    }
                    else {
                        punt.1 <- clas.snijpunt(snijp.1, loc, 1, j, i)
                        punt.2 <- clas.snijpunt(snijp.2, loc, 1, j, i)
                    }
                }
                else {
                    if(loc[i, 2] < loc[j, 2]) {
                        punt.1 <- clas.snijpunt(snijp.1, loc, 2, i, j)
                        punt.2 <- clas.snijpunt(snijp.2, loc, 2, i, j)
                    }
                    else {
                        punt.1 <- clas.snijpunt(snijp.1, loc, 2, j, i)
                        punt.2 <- clas.snijpunt(snijp.2, loc, 2, j, i)
                    }
                }
                if((punt.1[1] == "NA") || (punt.2[1] == "NA") ||
                   (sqrt((punt.1[1] - loc[i, 1])^2 +
                         (punt.1[2] - loc[i, 2])^2) +
                    sqrt((punt.2[1] - loc[j, 1])^2 +
                         (punt.2[2] - loc[j, 2])^2)) >
                   sqrt((loc[j, 1] - loc[i, 1])^2 +
                        (loc[j, 2] - loc[i, 2])^2))
                {
                    afstand[i, j] <- NA
                }
                else if(lines == 1) {
                    afstand[i, j] <- sqrt((loc[i, 1] - loc[j, 1])^2 +
                                          (loc[i, 2] - loc[j, 2])^2)
                    segments(loc[i, 1], loc[i, 2],
                             loc[j, 1], loc[j, 2], col = 6, ...)
                }
                else { ## lines == 2
                    afstand[i, j] <- sqrt((punt.1[1] - punt.2[1])^2 +
                                          (punt.1[2] - punt.2[2])^2)
                    segments(punt.1[1], punt.1[2],
                             punt.2[1], punt.2[2], col = 6, ...)
                }
            }
        }
        afstand <- t(afstand) + afstand
    }
    else afstand <- NULL

    ## FIXME: The following is *not* elegant..
    if(labels == 1) {
        for(i in 1:n) {
            x1 <- rbind(x1, z[[i]][cumsum(rep(10, 40)), ])
            labels1 <- c(labels1, rep(levclus[i], 40))
        }
        identify(x1[, 1], x1[, 2], labels1, col = col.txt)
    }
    else if(labels == 2) {
        x1 <- rbind(x1, maxima)
        labels1 <- c(labels1, levclus)
        x1[, 1] <- x1[, 1] + (maxx - minx)/130
        x1[, 2] <- x1[, 2] + (maxy - miny)/50
        text(x1, labels = labels1, col = col.txt)
    }
    else if(labels == 3) {
        x1[, 1] <- x1[, 1] + (maxx - minx)/130
        x1[, 2] <- x1[, 2] + (maxy - miny)/50
        text(x1, labels = labels1, col = col.txt)
    }
    else if(labels == 4) {
        maxima[, 1] <- maxima[, 1] + (maxx - minx)/ 130
        maxima[, 2] <- maxima[, 2] + (maxy - miny)/ 50
        text(maxima, labels = levclus, col = col.txt)
    }

    density[density == 41] <- NA
    invisible(list(Distances = afstand, Shading = density))
}

clusplot.partition <- function(x, main = NULL, ...)
{
    if(is.null(main) && !is.null(x$call))
	main <- paste("clusplot(",format(x$call),")", sep="")
    if(length(x$data) != 0 &&
       (!is.na(min(x$data)) || data.class(x) == "clara"))
	clusplot.default(x$data, x$clustering, diss = FALSE, main = main, ...)
    else clusplot.default(x$diss, x$clustering, diss = TRUE, main = main, ...)

}
