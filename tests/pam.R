library(cluster)

(x <- x0 <- cbind(V1 = (-3:4)^2, V2 = c(0:6,NA), V3 = c(1,2,NA,7,NA,8:9,8)))
(px <- pam(x,2, metric="manhattan"))
stopifnot(identical(x,x0))# DUP=FALSE ..
pd <-  pam(dist(x,"manhattan"), 2)
px2 <- pam(x,2, metric="manhattan", keep.diss=FALSE, keep.data=FALSE)
nms <- c("clustering", "objective", "isolation", "clusinfo", "silinfo")
stopifnot(identical(px[nms], pd[nms]),
          identical(px[nms], px2[nms]),
	  ## and for default dist "euclidean":
	  identical(pam(x,	2)[nms],
		    pam(dist(x),2)[nms])
	  )

if(paste(R.version$major, R.version$minor, sep=".") >= 1.7)  RNGversion(1.6)
set.seed(625)
## generate 250 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(120, 0,8), rnorm(120, 0,8)),
	   cbind(rnorm(130,50,8), rnorm(130,10,8)))

.proctime00 <- proc.time()

summary(px2 <- pam(x, 2))
all.equal(px2[nms], pam(dist(x), 2)[nms], tol = 1e-12) ## TRUE

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

