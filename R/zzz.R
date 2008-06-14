.onUnload <- function(libpath)
{
    library.dynam.unload("cluster", libpath)
}

## no S4 methodology here; speedup :
.noGenerics <- TRUE

if(paste(R.version$major, R.version$minor, sep=".") < 2.1) {
    ## These are substitutes such that newer code still runs in older R
    gettextf <- function (fmt, ..., domain = NULL)
        sprintf(gettext(fmt, domain = domain), ...)

    gettext <- function (..., domain = NULL) {
        args <- lapply(list(...), as.character)
        ##R 2.1.0: .Internal(gettext(domain, unlist(args)))
        ## Cheap substitute:
        unlist(args)
    }

    ngettext <- function (n, msg1, msg2, domain = NULL) {
        ##R 2.1.0: .Internal(ngettext(n, msg1, msg2, domain))
        ## Cheap substitute:
        if(n == 1) msg1 else msg2
    }
}
