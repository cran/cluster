
library(cluster)

###--- An extension of  example(fanny) : -------------------

set.seed(21)
## generate 10+15 objects in two clusters, plus 3 objects lying
## between those clusters.
x <- rbind(cbind(rnorm(10, 0, 0.5), rnorm(10, 0, 0.5)),
           cbind(rnorm(15, 5, 0.5), rnorm(15, 5, 0.5)),
           cbind(rnorm( 3,3.2,0.5), rnorm( 3,3.2,0.5)))

.proctime00 <- proc.time()

(fannyx <- fanny(x, 2))
summary(fannyx)
str(fannyx)
## Different platforms differ (even gcc 3.0.1 vs 3.2 on same platform)!
## {70 or 71 iterations}
## ==> No "fanny-ex.Rout.save" is distributed !
## --------------------------------------------
summary(fanny(x,3))# one extra cluster

(fanny(x,2, memb.exp = 1.5))
(fanny(x,2, memb.exp = 1.2))
(fanny(x,2, memb.exp = 1.1))
(fanny(x,2, memb.exp = 3))

data(ruspini) # < to run under R 1.9.1
summary(fanny(ruspini, 3), digits = 9)
summary(fanny(ruspini, 4), digits = 9)# `correct' #{clusters}
summary(fanny(ruspini, 5), digits = 9)

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
