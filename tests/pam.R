library(cluster)

(x <- cbind((-3:4)^2, c(0:6,NA), c(1,2,NA,7,NA,8:9,8)))
pam(x,2, metric="manhattan")

## generate 250 objects, divided into 2 clusters.
set.seed(625)
x <- rbind(cbind(rnorm(120, 0,8), rnorm(120, 0,8)),
           cbind(rnorm(130,50,8), rnorm(130,10,8)))

.proctime00 <- proc.time()

summary(px2 <- pam(x, 2))

data(ruspini)
summary(pr4 <- pam(ruspini, 4))
(pr3 <- pam(ruspini, 3))
(pr5 <- pam(ruspini, 5))

data(votes.repub)
summary(pv3 <- pam(votes.repub, 3))
(pv4 <- pam(votes.repub, 4))
(pv6 <- pam(votes.repub, 6))


## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')

