.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
}

##  for R versions < 1.6:
if(paste(R.version$major, R.version$minor, sep=".") < 1.6)
    stop <- function (..., call. = TRUE)
	.Internal(stop(if (nargs() > 0) paste(..., sep = "")))

## only for ../tests/*.R :  week substitute for R versions < 1.4:
if(!exists("identical", mode = "function"))
    identical <- function(x,y) {## not the *proper* identical(), just a substitute!!
	typeof(x) == typeof(y) &&
	length(x) == length(y) &&
	length(ax <- attributes(x)) == length(ay <- attributes(y)) &&
	all(ax == ay) &&
	all(if(is.list(x) || is.expression(x))
	    sapply(1:length(x), function(i) identical(x[[i]],y[[i]]))
	    else (x == y))
    }

##  for R versions < 1.4:
if(paste(R.version$major, R.version$minor, sep=".") < 1.4)
    warning <- function (..., call. = TRUE)
	.Internal(warning(if (nargs() > 0) paste(..., sep = "")))


## basically for R versions < 1.3:
if(!exists("dev.interactive", mode = "function"))
    dev.interactive <- function ()
    interactive() && .Device %in% c("X11", "GTK", "gnome", "windows", "Macintosh")
