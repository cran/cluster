## For different cluster versions

require(cluster)

if(interactive()) {
    (pkgPath <- .find.package("cluster", verbose = TRUE))
    (verC <- readLines(Dfile <- file.path(pkgPath, "DESCRIPTION"), n = 2)[2])
}

str(d5 <- data.frame(a= c(0, 0, 0,1,0,0, 0,0,1, 0,NA),
                     b= c(NA,0, 1,1,0,1, 0,1,0, 1,0),
                     c= c(0, 1, 1,0,1,NA,1,0,1, 0,NA),
                     d= c(1, 1, 0,1,0,0, 0,0,0, 1,0),
                     e= c(1, NA,0,1,0,0, 0,0,NA,1,1)))
(d0 <- daisy(d5))
(d1 <- daisy(d5, type = list(asymm = 1:5)))
(d2 <- daisy(d5, type = list(symm = 1:2, asymm= 3:5)))
(d2.<- daisy(d5, type = list(     asymm= 3:5)))
stopifnot(identical(c(d2), c(d2.)))

data(flower)
data(agriculture)

## ----------- example(daisy) -----------------------

## Example 1 in ref:
##  Dissimilarities using Euclidean metric and without standardization
(d.agr  <- daisy(agriculture, metric = "euclidean", stand = FALSE))
(d.agr2 <- daisy(agriculture, metric = "manhattan"))


## Example 2 in ref
(dfl0 <- daisy(flower))
stopifnot(identical(c(dfl0),
                    c(daisy(flower, type = list(symm = 1)))) &&
          identical(c(dfl0),
                    c(daisy(flower, type = list(symm = 2)))) &&
          identical(c(dfl0),
                    c(daisy(flower, type = list(symm = 3)))) &&
          identical(c(dfl0),
                    c(daisy(flower, type = list(symm = c(1,3)))))
         )

(dfl1 <- daisy(flower, type = list(asymm = 3)))
(dfl2 <- daisy(flower, type = list(asymm = c(1, 3), ordratio = 7)))
(dfl3 <- daisy(flower, type = list(asymm = 1:3)))

## --- animals
data(animals)
d0 <- daisy(animals)

d1 <- daisy(animals - 1, type=list(asymm=c(2,4)))
(d2 <- daisy(animals - 1, type=list(symm = c(1,3,5,6), asymm=c(2,4))))
stopifnot(c(d1) == c(d2))

d3 <- daisy(2 - animals, type=list(asymm=c(2,4)))
(d4 <- daisy(2 - animals, type=list(symm = c(1,3,5,6), asymm=c(2,4))))
stopifnot(c(d3) == c(d4))

pairs(cbind(d0,d2,d4),
      main = "Animals -- symmetric and asymm. dissimilarities")
