library(cluster)

eh <- ellipsoidhull(cbind(x=1:4, y = 1:4)) #singular
eh

if(paste(R.version$major, R.version$minor, sep=".") >= 1.7)  RNGversion(1.6)
set.seed(157)
for(n in 4:10) { ## n=2 and 3 still differ -- platform dependently!
    cat("n = ",n,"\n")
    x2 <- rnorm(n)
    try(print(ellipsoidhull(cbind(1:n, x2))))
    try(print(ellipsoidhull(cbind(1:n, x2, 4*x2 + rnorm(n)))))
}

x <- rt(100, df = 4)
y <- 100 + 5 * x + rnorm(100)
ellipsoidhull(cbind(x,y))
z <- 10  - 8 * x + y + rnorm(100)
ellipsoidhull(cbind(x,y,z))
