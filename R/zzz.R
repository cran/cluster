.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
}

##  for R versions < 1.6:
if(paste(R.version$major, R.version$minor, sep=".") < 1.6) {
    stop <- function (..., call. = TRUE)
	.Internal(stop(if (nargs() > 0) paste(..., sep = "")))

    ##	for R versions < 1.5
    if(paste(R.version$major, R.version$minor, sep=".") < 1.5)
	## cheap substitute, used in silhouette.default()
	colSums <- function(x) apply(x, 2, sum)


    ## NOTE: Now (from 1.7.0), require at least 1.4

}# versions < 1.6
