library(cluster)
## Compare on these:
nms <- c("clustering", "objective", "isolation", "clusinfo", "silinfo")
nm2 <- c("medoids", nms)

(x <- x0 <- cbind(V1 = (-3:4)^2, V2 = c(0:6,NA), V3 = c(1,2,NA,7,NA,8:9,8)))
(px <- pam(x,2, metric="manhattan"))
stopifnot(identical(x,x0))# DUP=FALSE ..
pd <-  pam(dist(x,"manhattan"), 2)
px2 <- pam(x,2, metric="manhattan", keep.diss=FALSE, keep.data=FALSE)

stopifnot(identical(px[nms], pd[nms]),
          identical(px[nms], px2[nms]),
	  ## and for default dist "euclidean":
	  identical(pam(x,	2)[nms],
		    pam(dist(x),2)[nms])
	  )

set.seed(253)
## generate 250 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(120, 0,8), rnorm(120, 0,8)),
	   cbind(rnorm(130,50,8), rnorm(130,10,8)))

.proctime00 <- proc.time()

summary(px2 <- pam(x, 2))
pdx <- pam(dist(x), 2)
all.equal(px2[nms], pdx[nms], tol = 1e-12) ## TRUE
pdxK <- pam(dist(x), 2, keep.diss = TRUE)
stopifnot(identical(pdx[nm2], pdxK[nm2]))

spdx <- silhouette(pdx)
summary(spdx)
spdx
postscript("pam-tst.ps")
if(FALSE)
    plot(spdx)# the silhouette
## is now identical :
plot(pdx)# failed in 1.7.0 -- now only does silhouette

par(mfrow = 2:1)
## new `dist' argument for clusplot():
plot(pdx, dist=dist(x))
## but this should work automagically (via eval()) as well:
plot(pdx)
## or this
clusplot(pdx)

data(ruspini)
summary(pr4 <- pam(ruspini, 4))
(pr3 <- pam(ruspini, 3))
(pr5 <- pam(ruspini, 5))

data(votes.repub)
summary(pv3 <- pam(votes.repub, 3))
(pv4 <- pam(votes.repub, 4))
(pv6 <- pam(votes.repub, 6))

dev.off()

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')

