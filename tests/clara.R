library(cluster)

## generate 1500 objects, divided into 2 clusters.
set.seed(264)
x <- rbind(cbind(rnorm(700, 0,8), rnorm(700, 0,8)),
           cbind(rnorm(800,50,8), rnorm(800,10,8)))
(clara2 <- clara(x, 2))
clara2$clusinfo

clara(x, 2, samples = 50)[c("objective", "medoids", "clusinfo")]

clara(x, 2, samples = 200)[c("objective", "medoids", "clusinfo")]
## Note that this last one is practically identical to the slower  pam() one

x[sample(length(x), 20)] <- NA
clara(x, 2, samples = 50)[c("objective", "medoids", "clusinfo")]
