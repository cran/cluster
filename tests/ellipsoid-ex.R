library(cluster)

ellipsoidhull(cbind(x=1:4, y = 1:4)) #singular

set.seed(157)
for(n in 2:10) {
    x2 <- rnorm(n)
    print(ellipsoidhull(cbind(1:n, x2)))
    print(ellipsoidhull(cbind(1:n, x2, 4*x2 + rnorm(n))))
}

x <- rt(100, df = 4)
y <- 100 + 5 * x + rnorm(100)
ellipsoidhull(cbind(x,y))
z <- 10  - 8 * x + y + rnorm(100)
ellipsoidhull(cbind(x,y,z))
