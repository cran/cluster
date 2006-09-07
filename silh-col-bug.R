## From: Patrick Connolly <p_connolly@ihug.co.nz>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: Re: [R] Colours in silhouette plots (cluster package)
## Date: Fri, 11 Aug 2006 09:02:26 +1200

## On Thu, 10-Aug-2006 at 12:08PM +0200, Martin Maechler wrote:

## |> Hi Patrick,
## |>
## |> Hmm,
## |>  - why didn't you send this first to the maintainer of 'cluster'
## |>    alone?

## I really thought the problem had more to do with me than with the
## cluster package and someone would set me straignt.  (Sometimes false
## modesty doesn't work.)

## |>  - where is the self-contained reproducible example ?

## Try this:

     data(ruspini)
     pr4 <- pam(ruspini, 4)
     str(si <- silhouette(pr4))
     (ssi <- summary(si))
     plot(si) # silhouette plot
### Works fine.  Now try a vector of colours:

plot(si, nmax= 80, cex.names=0.6, col = c("red", "green", "blue", "purple"))

## The later example does work:
si2 <- silhouette(pr4$clustering, dist(ruspini, "canberra"))
summary(si2) # has small values: "canberra"'s fault
plot(si2, nmax= 80, cex.names=0.6,
     col = c("red", "green", "blue", "purple"))

1
## So, it goes to show my hacking was probably not well-informed.
## Perhaps it does need the attention of the maintainer.

## best

## --
## ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
##    ___    Patrick Connolly
##  {~._.~}          		 Great minds discuss ideas
##  _( Y )_  	  	        Middle minds discuss events
## (:_~*~_:) 	       		 Small minds discuss people
##  (_)-(_)  	                           ..... Anon

## ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

