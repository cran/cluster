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
              deparse(attr(x, "Call"))),
              sub = paste("Average silhouette width : ",
              round(x$ silinfo$avg.width, digits = 2)), adj = 0)
    }

    if(is.null(which.plots) && !ask)
        which.plots <- 1:2
    if(ask && is.null(which.plots)) { ## Use `menu' ..
        tmenu <- paste("plot:", ## choices :
                       c("All", "Clusplot", "Silhouette Plot"))
        pick <- 0
        while(pick <= length(tmenu) + 1) {
            pick <- menu(tmenu, title =
                         "\nMake a plot selection (or 0 to exit):\n") + 1
            switch(pick,
                   return(invisible())
                   ,
                   clusplot(x, cor = cor, stand = stand, lines = lines,
                            shade = shade, color = color, labels = labels,
                            plotchar = plotchar, span = span,
                            xlim = xlim, ylim = ylim, ...)
                   ,
                   silhouPlot(x, nmax.lab, max.strlen)
                   )
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
         span = TRUE, xlim = NULL, ylim = NULL, ...)
{
    size <- function(d)
    {
        discr <- 1 + 8 * length(d)
        sqrtdiscr <- round(sqrt(discr))
        if(round(sqrtdiscr)^2 != discr) 0 else (1 + sqrtdiscr)/2
    }
    ellipse <- function(A, dist, loc, n = 201)
    {
        ## Return (x,y) points on ellipse boundary
        detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
        yl2 <- A[2, 2] * dist^2
        y <- seq( - sqrt(yl2), sqrt(yl2), leng = n)
        sqrt.discr <- sqrt(detA/A[2, 2]^2 * pmax(0, yl2 - y^2))
        sqrt.discr[c(1, n)] <- 0
        b <- loc[1] + A[1, 2]/A[2, 2] * y
        x1 <- b - sqrt.discr
        x2 <- b + sqrt.discr
        y <- loc[2] + y
        return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
    }

    kleur <- function(n, verhoud, z, dens, col, ...)
    {
        verhoud1 <- order(verhoud)
        if(n <= 4) {
            for(i in 1:n) {
                j <- verhoud1[i]
                polygon(z[[j]], ##not yet density = dens[j],
                        col = col[i], ...)
            }
        }
        else {
## if(exists("pam", mode = "function") == FALSE) {
##   print("Looking for function pam in library(cluster) to compute the color effect for more than 4 clusters." )
##   library(cluster)
## }
            j <- pam(sort(verhoud), 4)$clustering
            for(i in 1:n) {
                q <- verhoud1[i]
                polygon(z[[q]], ##not yet density = dens[q],
                        col = col[j[i]], ...)
            }
        }
    }
    clas.snijpunt <- function(x, loc, m, n, p)
    {
        if(     loc[n, m] <= x[1, m] && x[1, m] <= loc[p, m]) x[1, ]
        else if(loc[n, m] <= x[2, m] && x[2, m] <= loc[p, m]) x[2, ]
        else NA
    }
    plotje <- function(x, ...)
    {
        polygon(x, ##not yet density = 0,
                col = 5, ...)
    }
    notavail <- function(x)
    {
        x[x == "NA"] <- median(x, na.rm = TRUE)
        return(x)
    }
    coord.snijp1 <- function(x, gemid)
    {
        x[2, 2] - 2 * x[1, 2] * gemid + x[1, 1] * gemid^2
    }
    coord.snijp2 <- function(x, dist, y)
    {
        ((x[1, 1] * x[2, 2] - x[1, 2]^2) * dist^2)/y
    }
    coord.snijp3 <- function(x, y, n, gemid)
    {
        matrix(c(x[n, 1] + sqrt(y), x[n, 1] - sqrt(y),
                 x[n, 2] + gemid * sqrt(y),
                 x[n, 2] - gemid * sqrt(y)), ncol = 2)
    }

    ## BEGIN ----

    namx <- deparse(substitute(x))
    if(is.data.frame(x))
        x <- data.matrix(x)
    if(!is.numeric(x))
        stop("x is not numeric")

    labels1 <- NULL
    if(diss) {
        if(is.na(min(x)))
            stop(message = "NA-values in x are not allowed.")
        if((data.class(x)) != "dissimilarity") {
            if((size(x)) == 0) {
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
                    if(length(attr(x, "Labels")) != 0)
                        labels1 <- attr(x, "Labels")
                    x <- as.matrix(x)
                    if((n <- nrow(x)) == ncol(x) &&
                       all.equal(x, t(x)) == TRUE) {
                        labels1 <-
                            if(length(dimnames(x)[[1]]) == 0) 1:nrow(x)
                            else dimnames(x)[[1]]
                    }
                    else {
                        if(is.null(labels1))
                            labels1 <- 1:size(x)
                        attr(x, "Size") <- size(x)
                    }
                }
                else {
                    attr(x, "Size") <- size(x)
                    labels1 <- 1:size(x)
                }
            }
        }
        else {
            labels1 <-
                if(length(attr(x, "Labels")) == 0)
                    1:attr(x, "Size")
                else attr(x, "Labels")
        }
        ##x1 <- cmd(x, k = 2, eig = T, add = T)
        ##if(x1$ac < 0)
        ##	x1 <- cmd(x, k = 2, eig = T)
        x1 <- cmdscale(x, k = 2, eig = TRUE)
        var.dec <- sum(x1$eig)/sum(diag(x1$x))
        if (var.dec < 0) var.dec <- 0
        if (var.dec > 1) var.dec <- 1
        x1 <- x1$points
    }
    else { ## Not (diss)

        if(is.na(min(x))) {
            y <- is.na(x)
            y1 <- apply(y, 1, sum)
            y2 <- apply(y, 2, sum)
            if((sum(y1 == ncol(x)) != 0) && (sum(y2 == nrow(x)) != 0))
                stop("some objects and some variables contain only missing values"
                     )
            if(sum(y1 == nrow(x)) != 0)
                stop("one or more objects contain only missing values")
            if(sum(y2 == nrow(x)) != 0)
                stop("one or more variables contain only missing values")
            print("There were missing values and they were displaced by the median of the corresponding variable(s)"
                  )
            x <- apply(x, 2, notavail)
        }
        if(!is.matrix(x))
            stop("x is not allowed")
        ## ELSE
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
            x1 <- prim.pr$scores

            var.dec <- cumsum(prim.pr$sdev^2/sum(prim.pr$ sdev^2))[2]
            x1 <- cbind(x1[, 1], x1[, 2])
        }
    }
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
    n <- length(levclus)
    z <- A <- as.list(0)
    maxima <- loc <- matrix(0, ncol = 2, nrow = n)
    dist <- verhoud <- as.vector(0)
    verhouding <- 0
    num1 <- 10
    num2 <- 40
    num3 <- 90
    num4 <- 37
    num5 <- 3
    num6 <- 70

    for(i in 1:n) {
        x <- x1[clus == levclus[i], ]
        cov <-
          if(is.vector(x)) {
            x <- matrix(x, ncol = 2, byrow = TRUE)
            var(rbind(x, c(0, 0)))
          }
          else var(x)
        aantal <- nrow(x)
        x.1 <- range(x[, 1])
        y.1 <- range(x[, 2])
        notrank2 <- qr(cov, tol = 0.001)$rank != 2
        if(!span && notrank2) {
            dist[i] <- 1
            if((abs(diff(x.1)) > (diff(rangx)/70)) ||
               (abs(diff(y.1)) > (diff(rangy)/50))) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 + (loc[i, 2] - y.1[1])^2)
                a <- a + 0.05 * a
                if(abs(diff(x.1)) > (diff(rangx)/70)) {
                    ind1 <- (1:aantal)[x[,1]==max(x[,1])][1]
                    ind2 <- (1:aantal)[x[,1]==min(x[,1])][1]
                    q <- atan((x[ind1, 2] - x[ind2, 2])/
                              (x[ind1, 1] - x[ind2, 1]))
                    b <-
                        if(diff(rangy) == 0)
                            1
                        else if(abs(diff(y.1)) > (diff(rangy)/50))
                            diff(y.1)/num1
                        else diff(rangy)/num2
                }
                else {
                    b <- if(diff(rangx) == 0) 1 else diff(rangx)/num2
                    q <- pi/2
                }
                D <- diag(c(a^2, b^2))
                R <- cbind(c(  cos(q), sin(q)),
                           c(- sin(q), cos(q)))
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
            dist[i] <- 1
            if(sum(x[, 1] != x[1, 1]) != 0 ||
               sum(x[, 2] != x[1, 2]) != 0) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2,
                              y.1[1] + diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 +
                          (loc[i, 2] - y.1[1])^2)
                if(sum(x[, 1] != x[1, 1]) != 0) {
                    ind1 <- (1:aantal)[x[,1]==max(x[,1])][1]
                    ind2 <- (1:aantal)[x[,1]==min(x[,1])][1]
                    q <- atan((x[ind1, 2] - x[ind2, 2])/
                              (x[ind1, 1] - x[ind2, 1]))
                }
                else {
                    q <- pi/2
                }
                b <- 1e-7
                D <- diag(c(a^2, b^2))
                R <- cbind(c(  cos(q), sin(q)),
                           c(- sin(q), cos(q)))
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
                dist[i] <- sqrt(max(mahalanobis(x, loc[i, ], cov)))
                dist[i] <- dist[i] + 0.01 * dist[i]
            }
            else { ## span and rank2
                x2 <- cbind(matrix(1, aantal, 1), x)
                l1 <- matrix(0, 3, 3)
                sqdist <- prob <- rep(0, aantal)
                storage.mode(sqdist) <- "double"
                storage.mode(prob) <- "double"
                storage.mode(l1) <- "double"
                storage.mode(x2) <- "double"
                res <- .Fortran("spannel",
                                as.integer(aantal),
                                ndep= as.integer(2),
                                dat = x2,
                                eps = as.double(0.01),
                                sqdist = sqdist,
                                l1,
                                double(2),
                                double(2),
                                prob = prob,
                                double(3),
                                ierr = as.integer(0),
                                PACKAGE = "cluster")
                if(res$ierr != 0)
                    print("Error in Fortran routine computing the MVE-ellipsoid, please use the option exactmve=F"
                          )
                cov <- cov.wt(x, res$prob)$cov
                loc[i, ] <- cov.wt(x, res$prob)$center
                dist[i] <- sqrt(weighted.mean(res$sqdist, res$prob))
            }
            A[[i]] <- cov
            oppervlak <- pi * dist[i]^2 *
                sqrt(cov[1, 1] * cov[2, 2] - cov[1, 2]^2)
        }
        z[[i]] <- ellipse(A[[i]], dist[i], loc[i, ])
        rang <- c(range(z[[i]][, 1]), range(z[[i]][, 2]))
        maxima[i, ] <- z[[i]][201, ]
        minx <- min(minx, rang[1])
        maxx <- max(maxx, rang[2])
        miny <- min(miny, rang[3])
        maxy <- max(maxy, rang[4])
        verhoud[i] <- aantal/oppervlak
        if(verhoud[i] < 1e7)
            verhouding <- verhouding + verhoud[i]
    }
    if(verhouding == 0)
        verhouding <- 1
    density <- (verhoud * num4)/verhouding + num5
    density[density > 41] <- 41
    if (span) {
        if (rangx[1]==rangx[2]) {
            minx <- x1[1, 1] - 1
            maxx <- x1[1, 1] + 1
        }
        if (rangy[1]==rangy[2]) {
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
         main = paste("CLUSPLOT(", namx,")"),
         type = if(plotchar) "n" else "p", # if(plotchar) add points later
         col = col.p, ...)
    title(sub = paste("These two components explain",
          round(100 * var.dec, digits = 2), "% of the point variability."),
          adj = 0)

    color1 <- c(2, 4, 6, 3)

    if(shade && color) {
        kleur(n, verhoud, z, density, color1, ...)
    }
    else if(shade) {
        for(i in 1:n)
            polygon(z[[i]], ##not yet density = density[i],
                    col = 5, ...)
    }
    else if(color) {
        dens <- vector(mode = "numeric", length = n)
        kleur(n, verhoud, z, dens, color1, ...)
    }
    else {
        sapply(z, plotje, ...)
    }

    ## points after polygon in order to write ON TOP:
    if(plotchar) {
        karakter <- c(1:19)
        for(i in 1:n) {
            x <- x1[clus == levclus[i],  , drop = FALSE]
            kar <- 1+(i-1) %% 19
            points(x[, 1], x[, 2], pch = karakter[kar], col = col.p, ...)
        }
    }

    if((lines == 1 || lines == 2) && n > 1) {
        afstand <- matrix(0, ncol = n, nrow = n)
        for(i in 1:(n - 1)) {
            for(j in (i + 1):n) {
                gemid <- (loc[j, 2] - loc[i, 2])/(loc[j, 1] - loc[i, 1])
                s0 <- coord.snijp1(A[[i]], gemid)
                b0 <- coord.snijp2(A[[i]], dist[i], s0)
                snijp.1 <- coord.snijp3(loc, b0, i, gemid)
                s1 <- coord.snijp1(A[[j]], gemid)
                b1 <- coord.snijp2(A[[j]], dist[j], s1)
                snijp.2 <- coord.snijp3(loc, b1, j, gemid)
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
                else {
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

clusplot.partition <- function(x, ...)
{
    if(length(x$data) != 0 &&
       (!is.na(min(x$data)) || data.class(x) == "clara"))
         invisible(clusplot.default(x$data, x$clustering, diss = FALSE, ...))
    else invisible(clusplot.default(x$diss, x$clustering, diss = TRUE, ...))
}




