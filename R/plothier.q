"pltree"<-
function(x, ...)
{
	UseMethod("pltree")
}
"pltree.twins"<-
function(x, ...)
{
	call <- attr(x, "Call")
	labels <- NULL
	if(length(x$order.lab) != 0) {
		names(x$order) <- names(x$order.lab) <- 1:length(x$order)
		labels <- x$order.lab[names(sort(x$order))]
	}
	x <- list(order = x$order, height = sort(x$height), merge = x$merge)
	if(is.null(labels))
		plclust(x, plot = T, ylab = "Height", ...)
	else plclust(x, labels = labels, plot = T, ylab = "Height", 
			...)
	title(main = paste("Clustering tree of ", deparse(call)), adj = 0)
	invisible()
}
"plot.agnes"<-
function(x, ask = F, ...)
{
	choices <- c("All", "Banner", "Clustering Tree")
	choices <- substring(choices, 1, 40)
	tmenu <- paste("plot:", choices)
	pick <- 3
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
				w <- rev(x$height)
				m <- max(x$height)
				w <- rbind(w, m - w)
				barplot(w, xlab = "Height", horiz = T, inside
				   = F, space = 0, axes = F, col = c(0, 2), mgp
				   = c(2.5, 1, 0), ...)
				title(main = paste("Banner of ", deparse(attr(x,
				  "Call"))), sub = paste(
				  "Agglomerative Coefficient = ", round(x$ac, 
				  digits = 2)), adj = 0)
				flrm <- floor(m)
				at.vals <- c(seq(0, flrm, length = 11), m)
				lab.vals <- c(seq(0, flrm, length = 11), round(
				   m, digits = 2))
				axis(1, at = at.vals, labels = lab.vals, ...)
				if(length(x$order) < 35) {
					names <- if (length(x$order.lab) != 0)
					     substring(rev(x$order.lab), 1, 5)
						 else rev(x$order)
					axis(4, at = 0:(length(x$order) - 1), 
					  labels = names, pos = m, mgp = c(3, 
					  1.25, 0), ...)
				} 
			}
			,
			{
				pltree(x, ...)
			}
			)
		if(!ask.now)
			pick <- pick + 1
		if(pick == length(tmenu) + 2)
			ask.now <- ask
	}
	invisible()
}
"plot.diana"<-
function(x, ask = F, ...)
{
	choices <- c("All", "Banner", "Clustering Tree")
	choices <- substring(choices, 1, 40)
	tmenu <- paste("plot:", choices)
	pick <- 3
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
				w <- rev(x$height)
				m <- max(x$height)
				w <- rbind(m - w, w)
				barplot(w, xlab = "Height",  
				  horiz = T, inside = F, space = 0, axes = F, 
				  col = c(2, 0), mgp = c(2.5, 1, 0), ...)
				title(main = paste("Banner of ", deparse(attr(x,
				  "Call"))), sub = paste(
				  "Divisive Coefficient = ", round(x$dc, digits
				   = 2)), adj = 0)
				flrm <- floor(m)
				at.vals <- c(0, seq(0, flrm, length = 11) + 
				  m - flrm)
				lab.vals <- c(round(m, digits = 2), rev(
				  seq(0, flrm, length = 11)))
				axis(1, at = at.vals, labels = lab.vals, ...)
				if(length(x$order) < 35) {
					names <- if (length(x$order.lab) != 0)
					     substring(rev(x$order.lab), 1, 5)
				   		else rev(x$order)
					axis(2, at = 0:(length(x$order) - 1), 
					  labels = names, pos = 0, mgp = c(3, 
					  1.5, 0), ...)
				}
			}
			,
			{
				pltree(x, ...)
			}
			)
		if(!ask.now)
			pick <- pick + 1
		if(pick == length(tmenu) + 2)
			ask.now <- ask
	}
	invisible()
}
"plot.mona" <- 
function(x, ...)
{
	w <- rev(x$step)
	w[w==0] <- max(w)+1
	m <- max(w)
	barplot(rbind(w, m - w), xlab = "Separation step", horiz = T, 
		inside = F, space = 0, axes = F, col = c(2, 0), mgp
		 = c(2.5, 1, 0), ...)
	title(main = paste("Banner of ", deparse(attr(x, "Call"))), adj = 0)
	axis(1, at = 0:m, labels = 0:m, ...)
	if(length(x$order) < 35) {
		names <- if (length(x$order.lab) != 0)
				substring(rev(x$order.lab), 1, 5)
			else rev(x$order)
		axis(2, at = 0:(length(x$order) - 1), labels = names, pos = 0, 
		  mgp = c(3, 1.5, 0), ...)
	}
	names <- rev(x$variable)
	names[rev(x$step) == 0] <- ""
	text(w, 0:(length(x$order) - 2) + 0.5, labels = paste(" ", names), adj
		 = 0, col = 2, ...)
	invisible()
}
