pltree <- function(x, ...) UseMethod("pltree")

pltree.twins <-
  function(x, main = paste("Dendrogram of ", deparse(call)), ...)
{
    call <- attr(x, "Call")
    labels <- NULL
    if(length(x$order.lab) != 0) {
        names(x$order) <- names(x$order.lab) <- 1:length(x$order)
        labels <- x$order.lab[names(sort(x$order))]
    }
    x <- list(order = x$order, height = sort(x$height), merge = x$merge)

 if( sapply(R.version[c("major","minor")], as.numeric) %*% c(10,1) >= 12 ) {
    if(is.null(labels))
         plclust(x,                  main = main, ylab = "Height", ...)
    else plclust(x, labels = labels, main = main, ylab = "Height", ...)
 } else { ## R <= 1.1
    if(is.null(labels))
         plclust(x,                , ylab = "Height", ...)
    else plclust(x, labels = labels, ylab = "Height", ...)
    title(main = main, ...)
 }    
    invisible()
}

## plot.diana() & plot.agnes() are  almost identical;
##  just the bannerplot differs a bit ....

plot.agnes <- function(x, ask = FALSE, which.plots = NULL, 
                       main = NULL,
                       sub = paste("Agglomerative Coefficient = ",
                                   round(x$ac, digits = 2)),
                       adj = 0, nmax.lab = 35, max.strlen = 5, ...)
{
    bannerplot <- function(x, ...)
    {
        w <- rev(x$height)
        m <- max(x$height)
        w <- rbind(w, m - w)
        barplot(w, xlab = "Height", horiz = TRUE, inside = FALSE,
                space = 0, axes = FALSE, col = c(0, 2),
                mgp = c(2.5, 1, 0), ...)
        title(main = main1, sub = sub, adj = adj)
        flrm <- floor(m)
        at.vals <- c(seq(0, flrm, length = 11), m)
        lab.vals<- c(seq(0, flrm, length = 11), round(m, digits = 2))
        axis(1, at = at.vals, labels = lab.vals, ...)
        if(length(x$order) < nmax.lab) {
            names <- if (length(x$order.lab) != 0)
                substring(rev(x$order.lab), 1, max.strlen)
            else rev(x$order)
            axis(4, at = 0:(length(x$order) - 1), 
                 labels = names, pos = m, mgp = c(3, 1.25, 0), ...)
        } 
    }

    if(is.null(main)) {
        ## Different default for banner & pltree:
        cl <- deparse(attr(x, "Call"))
        main1 <- paste("Banner of ", cl)
        main2 <- paste("Dendrogram of ", cl)
    }
    else { # same title for both
        main1 <- .Alias(main)
        main2 <- .Alias(main)
    }
    if(is.null(which.plots)) { ## Use `menu' ..

        choices <- c("All", "Banner", "Clustering Tree")
        choices <- substring(choices, 1, 40)
        tmenu <- paste("plot:", choices)
        pick <- 2
        ask.now <- ask
        while(pick <= length(tmenu) + 2) {
            if(ask.now)
                pick <- menu(tmenu, title = 
                             "\nMake a plot selection (or 0 to exit):\n") + 1
            switch(pick,
                   return(invisible(x)),
                   ask.now <- FALSE,
                   bannerplot(x, ...),
                   pltree(x, main = main2, sub = sub, ...)
                   )
            if(!ask.now)
                pick <- pick + 1
            if(pick == length(tmenu) + 2)
                ask.now <- ask
        }
    }
    else for(i in which.plots)
        switch(i,
               bannerplot(x, ...),
               pltree    (x, main = main2, sub = sub, ...)
               )
    invisible()
}

plot.diana <-
function(x, ask = FALSE, which.plots = NULL,
         main = paste("Banner of ", deparse(attr(x, "Call"))),
         sub  = paste("Divisive Coefficient = ", round(x$dc, digits = 2)),
         adj = 0, nmax.lab = 35, max.strlen = 5, ...)
{
    bannerplot <- function(x, ...)
    {
        w <- rev(x$height)
        m <- max(x$height)
        w <- rbind(m - w, w)
        barplot(w, xlab = "Height", horiz = TRUE, inside = FALSE,
                space = 0, axes = FALSE, col = c(2, 0),
                mgp = c(2.5, 1, 0), ...)
        title(main = main1, sub = sub, adj = adj)
        flrm <- floor(m)
        at.vals <- c(0, seq(0, flrm, length = 11) + m - flrm)
        lab.vals <- c(round(m, digits = 2),
                      rev(seq(0, flrm, length = 11)))
        axis(1, at = at.vals, labels = lab.vals, ...)
        if(length(x$order) < nmax.lab) {
            names <- if (length(x$order.lab) != 0)
                substring(rev(x$order.lab), 1, max.strlen)
            else rev(x$order)
            axis(2, at = 0:(length(x$order) - 1), 
                 labels = names, pos = 0, mgp = c(3, 1.5, 0), ...)
        }
    }

    if(is.null(main)) {
        ## Different default for banner & pltree:
        cl <- deparse(attr(x, "Call"))
        main1 <- paste("Banner of ", cl)
        main2 <- paste("Dendrogram of ", cl)
    }
    else { # same title for both
        main1 <- .Alias(main)
        main2 <- .Alias(main)
    }
    if(is.null(which.plots)) { ## Use `menu' ..
        choices <- c("All", "Banner", "Clustering Tree")
        tmenu <- paste("plot:", choices)
        pick <- 3
        ask.now <- ask
        while(pick <= length(tmenu) + 2) {
            if(ask.now)
                pick <- menu(tmenu, title = 
                             "\nMake a plot selection (or 0 to exit):\n") + 1
            switch(pick,
                   return(invisible(x)),
                   ask.now <- FALSE,
                   bannerplot(x, ...),
                       pltree(x, main = main2, sub = sub, ...)
                   )
            if(!ask.now)
                pick <- pick + 1
            if(pick == length(tmenu) + 2)
                ask.now <- ask
        }
    }
    else for(i in which.plots)
        switch(i,
               bannerplot(x, ...),# i = 1
               pltree    (x, main = main2, sub = sub, ...) # i = 2
               )
    invisible()
}

plot.mona <- function(x, main = paste("Banner of ", deparse(attr(x, "Call"))),
                      col = 2, axes = TRUE, adj = 0,
                      nmax.lab = 35, max.strlen = 5,  ...)
{
    w <- rev(x$step)
    w[w==0] <- max(w)+1
    m <- max(w)
    barplot(rbind(w, m - w), xlab = "Separation step", horiz = TRUE, 
            inside = FALSE, space = 0, axes = FALSE,
            col = c(col, 0), mgp = c(2.5, 1, 0), ...)
    title(main = main, adj = adj, ...)
    if(axes) axis(1, at = 0:m, labels = 0:m, ...)
    if(length(x$order) < nmax.lab) {
        names <- if (length(x$order.lab) != 0)
            substring(rev(x$order.lab), 1, max.strlen)
        else rev(x$order)
        if(axes)
            axis(2, at = 0:(length(x$order) - 1), labels = names, pos = 0, 
                 mgp = c(3, 1.5, 0), las = 1, ...)
    }
    names <- rev(x$variable)
    names[rev(x$step) == 0] <- ""
    text(w, 0:(length(x$order) - 2) + 0.5, labels = paste(" ", names),
         adj = adj, col = col, ...)
    invisible()
}
