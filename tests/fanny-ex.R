
library(cluster)

###--- An extension of  example(fanny) : -------------------

if(paste(R.version$major, R.version$minor, sep=".") >= 1.7)  RNGversion(1.6)
set.seed(171)
## generate 25 objects, divided into two clusters, and 3 objects lying
## between those clusters.
x <- rbind(cbind(rnorm(10, 0, 0.5), rnorm(10, 0, 0.5)),
           cbind(rnorm(15, 5, 0.5), rnorm(15, 5, 0.5)),
           cbind(rnorm( 3,3.5,0.5), rnorm( 3,3.5,0.5)))

data(ruspini)

.proctime00 <- proc.time()

(fannyx <- fanny(x, 2))
summary(fannyx)
str(fannyx)
## Different platforms differ (even gcc 3.0.1 vs 3.2 on same platform)!
## {70 or 71 iterations}
## ==> No "fanny-ex.Rout.save" is distributed !
## --------------------------------------------
summary(fanny(x,3))# one extra cluster

summary(fanny(ruspini, 3), digits = 9)
summary(fanny(ruspini, 4), digits = 9)# `correct' #{clusters}
summary(fanny(ruspini, 5), digits = 9)

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
