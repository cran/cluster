"plot.partition" <-
function(x, ask = F, cor = T, stand = F, lines = 2, shade = F, color = F, 
	labels = 0, plotchar = T, span = T, xlim = NULL, ylim = NULL, ...)
{
	choices <- c("All", "Clusplot", "Silhouette Plot")
	choices <- substring(choices, 1, 40)
	tmenu <- paste("plot:", choices)
	pick <- 4
	ask.now <- ask
	z <- NULL
	while(pick <= length(tmenu) + 2) {
		if(ask.now)
			pick <- menu(tmenu, title = 
				"\nMake a plot selection (or 0 to exit):\n") + 
				1
		switch(pick,
			return(invisible(x)),
			ask.now <- F,
			{
				clusplot(x, cor = cor, stand = stand, lines = 
				  lines, shade = shade, color = color, labels
				   = labels, plotchar = plotchar, span = span, 
				  xlim = xlim, ylim = ylim, ...)
			}
			,
			{
				if(length(x$silinfo) == 0)
				  stop(message = 
				    "No silhouette plot available when the number of clusters equals 1."
				    )
				s <- rev(x$silinfo[[1]][, 3])
				space <- c(0, rev(diff(x$silinfo[[1]][, 1])))
				space[space != 0] <- 0.5
				names <- if(length(s) < 40) substring(rev(
				     dimnames(x$silinfo[[1]])[[1]]), 1, 5)
				   else NULL
				barplot(s, space = space, names = names, xlab
				   = "Silhouette width", xlim
				   = c(min(0, min(s)), 1), horiz = T, mgp = c(
				  2.5, 1, 0), ...)
				title(main = paste("Silhouette plot of ", 
				  deparse(attr(x, "Call"))), sub = paste(
				  "Average silhouette width : ", round(x$
				  silinfo$avg.width, digits = 2)), adj = 0)
			}
			)
		if(!ask.now)
			pick <- pick + 1
		if(pick == length(tmenu) + 2)
			ask.now <- ask
	}
	invisible()
}
"clusplot"<-
function(x, ...)
{
	UseMethod("clusplot")
}
"clusplot.default"<-
function(x, clus, diss = F, cor = T, stand = F, lines = 2, shade = F, color = F,
	labels = 0, plotchar = T, span = T, xlim = NULL, ylim = NULL, ...)
{	
	size <- function(d)
	{
		discr <- 1 + 8 * length(d)
		sqrtdiscr <- round(sqrt(discr))
		if(round(sqrtdiscr)^2 != discr)
			return(0)
		(1 + sqrtdiscr)/2
	}
	ellipse <- function(A, dist, loc)
	{
		detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
		ylimit <- sqrt(A[2, 2]) * dist
		y <- seq( - ylimit, ylimit, 0.01 * ylimit)
		sqrt.discr <- sqrt(detA/A[2, 2]^2 * (A[2, 2] * dist^2 - y^2))
		sqrt.discr[c(1, length(sqrt.discr))] <- 0
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
				polygon(z[[j]], density = dens[j], col = col[i],
				  ...)
			}
		}
		else {
			if(exists("pam", mode = "function") == F) {
				print(
				  "Looking for function pam in library(cluster) to compute the color effect for more than 4 clusters."
				  )
				library(cluster)
			}
			j <- pam(sort(verhoud), 4)$clustering
			for(i in 1:n) {
				q <- verhoud1[i]
				polygon(z[[q]], density = dens[q], col = col[j[
				  i]], ...)
			}
		}
	}
	clas.snijpunt <- function(x, loc, m, n, p)
	{
		if((loc[n, m] <= x[1, m]) && (x[1, m] <= loc[p, m])) {
			f <- x[1,  ]
		}
		else {
			if((loc[n, m] <= x[2, m]) && (x[2, m] <= loc[p, m])) {
				f <- x[2,  ]
			}
			else {
				f <- NA
			}
		}
		return(f)
	}
	plotje <- function(x, ...)
	{
		polygon(x, density = 0, col = 5, ...)
	}
	notavail <- function(x)
	{
		x[x == "NA"] <- median(x, na.rm = T)
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
		matrix(c(x[n, 1] + sqrt(y), x[n, 1] - sqrt(y), x[n, 2] + gemid * 
			sqrt(y), x[n, 2] - gemid * sqrt(y)), ncol = 2)
	}
	if(is.data.frame(x)) {
		x <- as.matrix(x)
	}
	if(!is.numeric(x))
		stop(message = "x is not numeric")
	labels1 <- NULL
	if(diss) {
		if(is.na(min(x)))
			stop(message = "NA-values in x are not allowed.")
		if((data.class(x)) != "dissimilarity") {
			if((size(x)) == 0) {
				if((n <- nrow(x)) != ncol(x))
				  stop(message = 
				    "Distances must be result of dist or a square matrix."
				    )
				if(all.equal(x, t(x)) != T)
				  stop(message = 
				    "the square matrix is not symmetric.")
				if(length(dimnames(x)[[1]]) == 0) {
				  labels1 <- c(1:nrow(x))
				}
				else {
				  labels1 <- dimnames(x)[[1]]
				}
			}
			else {
				if(is.vector(x) == F) {
				  if(length(attr(x, "Labels")) != 0)
				    labels1 <- attr(x, "Labels")
				  x <- as.matrix(x)
				  if(((n <- nrow(x)) == ncol(x)) && (all.equal(
				    x, t(x)) == T)) {
				    if(length(dimnames(x)[[1]]) == 0)
				      labels1 <- c(1:nrow(x))
				    else labels1 <- dimnames(x)[[1]]
				  }
				  else {
				    if(is.null(labels1))
				      labels1 <- c(1:size(x))
				    attr(x, "Size") <- size(x)
				  }
				}
				else {
				  attr(x, "Size") <- size(x)
				  labels1 <- c(1:size(x))
				}
			}
		}
		else {
			if(length(attr(x, "Labels")) == 0)
				labels1 <- c(1:attr(x, "Size"))
			else labels1 <- attr(x, "Labels")
		}
		##x1 <- cmd(x, k = 2, eig = T, add = T)
		##if(x1$ac < 0)
		##	x1 <- cmd(x, k = 2, eig = T)
                x1 <- cmdscale(x, k = 2, eig = T)
		var.dec <- sum(x1$eig)/sum(diag(x1$x))
		if (var.dec < 0) var.dec <- 0
		if (var.dec > 1) var.dec <- 1
		x1 <- x1$points
	}
	else {
		if(is.na(min(x))) {
			y <- is.na(x)
			y1 <- apply(y, 1, sum)
			y2 <- apply(y, 2, sum)
			if((sum(y1 == ncol(x)) != 0) && (sum(y2 == nrow(x)) != 
				0))
				stop(message = 
				  "some objects and some variables contain only missing values"
				  )
			if(sum(y1 == nrow(x)) != 0)
				stop(message = 
				  "one or more objects contain only missing values"
				  )
			if(sum(y2 == nrow(x)) != 0)
				stop(message = 
				  "one or more variables contain only missing values"
				  )
			print("There were missing values and they were displaced by the median of the corresponding variable(s)"
				)
			x <- apply(x, 2, notavail)
		}
		if(!is.matrix(x)) {
			stop(message = "x is not allowed")
		}
		else {
			if(length(dimnames(x)[[1]]) == 0) {
				labels1 <- c(1:nrow(x))
			}
			else {
				labels1 <- dimnames(x)[[1]]
			}
			if(ncol(x) == 1) {
				hulp <- rep(0, length(x))
				x1 <- matrix(c(t(x), hulp), ncol = 2)
				var.dec <- 1
			}
			else {
				if(ncol(x) == 2) {
				  prim.pr <- princomp(x, scores = T, cor = F)
				  x1 <- prim.pr$scores
				}
				else {
				  prim.pr <- princomp(x, scores = T, cor = cor)
				  x1 <- prim.pr$scores
				}
				var.dec <- cumsum(prim.pr$sdev^2/sum(prim.pr$
				  sdev^2))[2]
				x1 <- cbind(x1[, 1], x1[, 2])
			}
		}
	}
	clus <- as.vector(clus)
	if(length(clus) != length(x1[, 1]))
		stop(message = "The clustering vector has not the good length")
	clus <- as.factor(clus)
	if(sum(is.na(clus)) != 0)
		stop(message = "NA-values are not allowed in clustering vector"
			)
	if(stand == T) {
		x1 <- scale(x1)
	}
	rangx <- range(x1[, 1])
	rangy <- range(x1[, 2])
	minx <- rangx[1]
	maxx <- rangx[2]
	miny <- rangy[1]
	maxy <- rangy[2]
	levclus <- levels(clus)
	n <- length(levclus)
	z <- A <- as.list(0)
	loc <- matrix(0, ncol = 2, nrow = n)
	dist <- verhoud <- as.vector(0)
	verhouding <- 0
	maxima <- matrix(0, ncol = 2, nrow = n)
	num1 <- 10
	num2 <- 40
	num3 <- 90
	num4 <- 37
	num5 <- 3
	num6 <- 70
	for(i in 1:n) {
		x <- x1[clus == levclus[i],  ]
		if(is.vector(x)) {
			x <- matrix(x, ncol = 2, byrow = T)
			cov <- var(rbind(x, c(0, 0)))
		}
		else {
			cov <- var(x)
		}
		aantal <- nrow(x)
		x.1 <- range(x[, 1])
		y.1 <- range(x[, 2])
		if(span == F && (qr(cov, tol = 0.001)$rank != 2)) {
			dist[i] <- 1
			if((abs(diff(x.1)) > (diff(rangx)/70)) || (abs(diff(y.1
				)) > (diff(rangy)/50))) {
				loc[i,  ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + 
				  diff(y.1)/2)
				a <- sqrt((loc[i, 1] - x.1[1])^2 + (loc[i, 2] - 
				  y.1[1])^2)
				a <- a + 0.050000000000000003 * a
				if(abs(diff(x.1)) > (diff(rangx)/70)) {
				  ind1 <- (1:aantal)[x[,1]==max(x[,1])][1]
				  ind2 <- (1:aantal)[x[,1]==min(x[,1])][1]
				  q <- atan((x[ind1, 2] - x[ind2, 2])/(x[ind1, 
				    1] - x[ind2, 1]))
				  if(diff(rangy) == 0) {
				    b <- 1
				  }
				  else {
				    if(abs(diff(y.1)) > (diff(rangy)/50)) {
				      b <- diff(y.1)/num1
				    }
				    else {
				      b <- diff(rangy)/num2
				    }
				  }
				}
				else {
				  if(diff(rangx) == 0) {
				    b <- 1
				  }
				  else {
				    b <- diff(rangx)/num2
				  }
				  q <- pi/2
				}
				D <- matrix(c(a^2, 0, 0, b^2), ncol = 2)
				R <- matrix(c(cos(q), sin(q),  - sin(q), cos(q)
				  ), ncol = 2)
				A[[i]] <- (R %*% D) %*% t(R)
			}
			else {
				a <- diff(rangx)/num3
				b <- diff(rangy)/num6
				if(a == 0) {
				  a <- 1
				}
				if(b == 0) {
				  b <- 1
				}
				A[[i]] <- matrix(c(a^2, 0, 0, b^2), ncol = 2)
				loc[i,  ] <- x[1,  ]
			}
			oppervlak <- pi * a * b
		}
		else {
			if((span == T && (qr(cov, tol = 0.001)$rank != 2))) {
				dist[i] <- 1
				if((sum(x[, 1] != x[1, 1]) != 0) || (sum(x[, 2] !=
				  x[1, 2]) != 0)) {
				  loc[i,  ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + 
				    diff(y.1)/2)
				  a <- sqrt((loc[i, 1] - x.1[1])^2 + (loc[i, 2] -
				    y.1[1])^2)
				  if(sum(x[, 1] != x[1, 1]) != 0) {
				    ind1 <- (1:aantal)[x[,1]==max(x[,1])][1]
				    ind2 <- (1:aantal)[x[,1]==min(x[,1])][1]
				    q <- atan((x[ind1, 2] - x[ind2, 2])/(x[ind1, 
				      1] - x[ind2, 1]))
				    b <- 9.9999999999999995e-08
				  }
				  else {
				    b <- 9.9999999999999995e-08
				    q <- pi/2
				  }
				  D <- matrix(c(a^2, 0, 0, b^2), ncol = 2)
				  R <- matrix(c(cos(q), sin(q),  - sin(q), cos(
				    q)), ncol = 2)
				  A[[i]] <- (R %*% D) %*% t(R)
				}
				else {
				  a <- diff(rangx)/num3
				  b <- diff(rangy)/num6
				  if(a == 0) {
				    a <- 1
				  }
				  if(b == 0) {
				    b <- 1
				  }
				  A[[i]] <- matrix(c(a^2, 0, 0, b^2), ncol = 2)
				  loc[i,  ] <- x[1,  ]
				}
				oppervlak <- pi * a * b
			}
			else {
				if(span == F) {
				  loc[i,  ] <- apply(x, 2, mean)
				  dist[i] <- sqrt(max(mahalanobis(x, loc[i,  ], 
				    cov)))
				  dist[i] <- dist[i] + 0.01 * dist[i]
				}
				else {
				  x2 <- cbind(matrix(1, aantal, 1), x)
				  l1 <- matrix(0, 3, 3)
				  sqdist <- prob <- rep(0, aantal)
				  storage.mode(sqdist) <- "double"
				  storage.mode(prob) <- "double"
				  storage.mode(l1) <- "double"
				  storage.mode(x2) <- "double"
				  res <- .Fortran("spannel",
				    as.integer(aantal),
				    as.integer(2),
				    x2,
				    as.double(0.01),
				    sqdist = sqdist,
				    l1,
				    double(2),
				    double(2),
				    prob = prob,
				    double(3),
				    stop = as.integer(0))
				  if(res$stop != 0)
				    print(
				      "Error in Fortran routine computing the MVE-ellipsoid, please use the option exactmve=F"
				      )
				  cov <- cov.wt(x, res$prob)$cov
				  loc[i,  ] <- cov.wt(x, res$prob)$center
				  dist[i] <- sqrt(weighted.mean(res$sqdist, res$
				    prob))
				}
				A[[i]] <- cov
				oppervlak <- pi * dist[i]^2 * sqrt(cov[1, 1] * 
				  cov[2, 2] - cov[1, 2]^2)
			}
		}
		z[[i]] <- ellipse(A[[i]], dist[i], loc[i,  ])
		rang <- c(range(z[[i]][, 1]), range(z[[i]][, 2]))
		maxima[i,  ] <- z[[i]][201,  ]
		minx <- min(minx, rang[1])
		maxx <- max(maxx, rang[2])
		miny <- min(miny, rang[3])
		maxy <- max(maxy, rang[4])
		verhoud[i] <- aantal/oppervlak
		if(verhoud[i] < 10000000)
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
	if(is.null(xlim) == F) {
		if(xlim[1] < minx)
			minx <- xlim[1]
		if(xlim[2] > maxx)
			maxx <- xlim[2]
	}
	if(is.null(ylim) == F) {
		if(ylim[1] < miny)
			miny <- ylim[1]
		if(ylim[2] > maxy)
			maxy <- ylim[2]
	}
	if(plotchar == F) {
		plot(x1[, 1], x1[, 2], xlim = c(minx, maxx), ylim = c(miny, 
			maxy), xlab = "Component 1", ylab = "Component 2", col
			 = 5, ...)
	}
	else {
		x <- x1[clus == levclus[1],  , drop = F]
		plot(x[, 1], x[, 2], xlim = c(minx, maxx), ylim = c(miny, maxy),
			xlab = "Component 1", ylab = "Component 2", col = 5, 
			...)
		if(n != 1) {
			karakter <- c(0:18)
			for(i in (2:n)) {
				x <- x1[clus == levclus[i],  , drop = F]
				kar <- i - floor((i - 1)/19) * 19
				points(x[, 1], x[, 2], pch = karakter[kar], col
				   = 5, ...)
			}
		}
	}
	title("CLUSPLOT")
	title(sub = paste("These two components explain", round(100 * var.dec, 
		digits = 2), "% of the point variability."), adj = 0)
	color1 <- c(2, 4, 6, 3)
	if((shade == T) && (color == T)) {
		kleur(n, verhoud, z, density, color1, ...)
	}
	else {
		if(shade == T) {
			for(i in 1:n) {
				polygon(z[[i]], density = density[i], col = 5, 
				  ...)
			}
		}
		else {
			if(color == T) {
				dens <- vector(mode = "numeric", length = n)
				kleur(n, verhoud, z, dens, color1, ...)
			}
			else {
				sapply(z, plotje, ...)
			}
		}
	}
	if(((lines == 1) || (lines == 2)) && (n != 1)) {
		afstand <- matrix(0, ncol = n, nrow = n)
		for(i in 1:(n - 1)) {
			for(j in (i + 1):n) {
				gemid <- (loc[j, 2] - loc[i, 2])/(loc[j, 1] - 
				  loc[i, 1])
				s0 <- coord.snijp1(A[[i]], gemid)
				b0 <- coord.snijp2(A[[i]], dist[i], s0)
				snijp.1 <- coord.snijp3(loc, b0, i, gemid)
				s1 <- coord.snijp1(A[[j]], gemid)
				b1 <- coord.snijp2(A[[j]], dist[j], s1)
				snijp.2 <- coord.snijp3(loc, b1, j, gemid)
				if(loc[i, 1] != loc[j, 1]) {
				  if(loc[i, 1] < loc[j, 1]) {
				    punt.1 <- clas.snijpunt(snijp.1, loc, 1, i, 
				      j)
				    punt.2 <- clas.snijpunt(snijp.2, loc, 1, i, 
				      j)
				  }
				  else {
				    punt.1 <- clas.snijpunt(snijp.1, loc, 1, j, 
				      i)
				    punt.2 <- clas.snijpunt(snijp.2, loc, 1, j, 
				      i)
				  }
				}
				else {
				  if(loc[i, 2] < loc[j, 2]) {
				    punt.1 <- clas.snijpunt(snijp.1, loc, 2, i, 
				      j)
				    punt.2 <- clas.snijpunt(snijp.2, loc, 2, i, 
				      j)
				  }
				  else {
				    punt.1 <- clas.snijpunt(snijp.1, loc, 2, j, 
				      i)
				    punt.2 <- clas.snijpunt(snijp.2, loc, 2, j, 
				      i)
				  }
				}
				if((punt.1[1] == "NA") || (punt.2[1] == "NA")) 
				  {
				  afstand[i, j] <- NA
				}
				else {
				  if((sqrt((punt.1[1] - loc[i, 1])^2 + (punt.1[
				    2] - loc[i, 2])^2) + sqrt((punt.2[1] - loc[
				    j, 1])^2 + (punt.2[2] - loc[j, 2])^2)) > 
				    sqrt((loc[j, 1] - loc[i, 1])^2 + (loc[j, 2] -
				    loc[i, 2])^2)) {
				    afstand[i, j] <- NA
				  }
				  else {
				    if(lines == 1) {
				      afstand[i, j] <- sqrt((loc[i, 1] - loc[j, 
				        1])^2 + (loc[i, 2] - loc[j, 2])^2)
				      segments(loc[i, 1], loc[i, 2], loc[j, 1], 
				        loc[j, 2], col = 5, ...)
				    }
				    else {
				      afstand[i, j] <- sqrt((punt.1[1] - punt.2[
				        1])^2 + (punt.1[2] - punt.2[2])^2)
				      segments(punt.1[1], punt.1[2], punt.2[1], 
				        punt.2[2], col = 5, ...)
				    }
				  }
				}
			}
		}
		afstand <- t(afstand) + afstand
	}
	else afstand <- NULL
	if(labels == 1) {
		for(i in 1:n) {
			x1 <- rbind(x1, z[[i]][cumsum(rep(10, 40)),  ])
			labels1 <- c(labels1, rep(levclus[i], 40))
		}
		identify(x1[, 1], x1[, 2], labels1, col = 5)
	}
	else {
		if(labels == 2) {
			x1 <- rbind(x1, maxima)
			labels1 <- c(labels1, levclus)
			x1[, 1] <- x1[, 1] + (maxx - minx)/130
			x1[, 2] <- x1[, 2] + (maxy - miny)/50
			text(x1, labels = labels1, col = 5)
		}
		else {
			if(labels == 3) {
				x1[, 1] <- x1[, 1] + (maxx - minx)/130
				x1[, 2] <- x1[, 2] + (maxy - miny)/50
				text(x1, labels = labels1, col = 5)
			}
			else {
				if(labels == 4) {
				  maxima[, 1] <- maxima[, 1] + (maxx - minx)/
				    130
				  maxima[, 2] <- maxima[, 2] + (maxy - miny)/50
				  text(maxima, labels = levclus, col = 
				    5)
				}
			}
		}
	}
	density[density == 41] <- NA
	summary. <- list(afstand, density)
	names(summary.) <- c("Distances", "Shading")
	invisible(summary.)
}
"clusplot.partition"<-
function(x, ...)
{
	if(length(x$data) != 0) 
		if(!is.na(min(x$data)))
			invisible(clusplot.default(x$data, x$clustering, 
				diss = F, ...))
		else {
			if(data.class(x) == "clara")
				invisible(clusplot.default(x$data, 
					x$clustering, diss = F, ...))
			else invisible(clusplot.default(x$diss, x$clustering, 
				diss = T, ...))
		}
	else invisible(clusplot.default(x$diss, x$clustering, diss = T, ...))
}


