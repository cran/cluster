## Ensure consistent "diss.." class --- make "namespace-private-global !
dissiCl <- c("dissimilarity", "dist")

if((Rv <- getRversion()) < "3.2.1") {
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
    if(Rv < "3.1.0") {
        anyNA <- function(x) any(is.na(x))
        ## if(Rv < "3.0.0") {
        ##     rep_len <- function(x, length.out) rep(x, length.out=length.out)
        ##     ## if(Rv < "2.15")
        ##     ##     paste0 <- function(...) paste(..., sep = '')
        ## }
    }
}; rm(Rv)

##' Not exported, useful to run CRAN checks faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_CLUSTER_CHECK_EXTRA")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

