.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
}

## only for ../tests/*.R :  week substitute for R versions < 1.4:
if(!exists("identical", mode = "function"))
    identical <- function(x,y) {## not the *proper* identical(), just a substitute!!
        typeof(x) == typeof(y) &&
        length(x) == length(y) &&
        length(ax <- attributes(x)) == length(ay <- attributes(y)) &&
        all(ax == ay) &&
        all(x == y)
    }

##  for R versions < 1.4:
if(paste(R.version$major, R.version$minor, sep=".") < 1.4)
    warning <- function (..., call. = TRUE) {
        .Internal(warning(if (nargs() > 0) paste(..., sep = "")))
    }

## basically for R versions < 1.3:
if(!exists("dev.interactive", mode = "function"))
    dev.interactive <- function ()
    interactive() && .Device %in% c("X11", "GTK", "gnome", "windows", "Macintosh")
