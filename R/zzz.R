.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
}

## no S4 methodology here; speedup :
.noGenerics <- TRUE

## for R versions < 1.8:
if(paste(R.version$major, R.version$minor, sep=".") < 1.8) {

    sQuote <- function(x) {
        if(length(x) == 0) return(character())
        paste("'", x, "'", sep = "")
    }
    dQuote <- function(x) {
        if(length(x) == 0) return(character())
        paste("\"", x, "\"", sep = "")
    }

    ## for R versions < 1.7:
    if(paste(R.version$major, R.version$minor, sep=".") < 1.7) {

        force <- function(x) x

        ## for R versions < 1.6:
        if(paste(R.version$major, R.version$minor, sep=".") < 1.6) {
            stop <- function (..., call. = TRUE)
                .Internal(stop(if (nargs() > 0) paste(..., sep = "")))

            ## for R versions < 1.5
            if(paste(R.version$major, R.version$minor, sep=".") < 1.5)
                ## cheap substitute, used in silhouette.default()
                colSums <- function(x) apply(x, 2, sum)

### NOTE: From cluster 1.7.0, we require at least R 1.4

        } # versions < 1.6

    } # versions < 1.7
} # versions < 1.8
