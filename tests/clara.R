### Note: MM has a version fo this in ../../cluster_tests/clara.R !!

library(cluster)

## generate 1500 objects, divided into 2 clusters.
set.seed(264)
x <- rbind(cbind(rnorm(700, 0,8), rnorm(700, 0,8)),
           cbind(rnorm(800,50,8), rnorm(800,10,8)))

.proctime00 <- proc.time()

summary(clara2 <- clara(x, 2))

clara(x, 2, samples = 50)[c("objective", "medoids", "clusinfo")]

clara(x, 2, samples = 200)[c("objective", "medoids", "clusinfo")]
## Note that this last one is practically identical to the slower  pam() one

x[print(sample(length(x), 20))] <- NA
clara(x, 2, samples = 50)[c("objective", "medoids", "clusinfo")]

###-- Larger example: 2000 objects, divided into 5 clusters.
x5 <- rbind(cbind(rnorm(400, 0,4), rnorm(400, 0,4)),
            cbind(rnorm(400,10,8), rnorm(400,40,6)),
            cbind(rnorm(400,30,4), rnorm(400, 0,4)),
            cbind(rnorm(400,40,4), rnorm(400,20,2)),
            cbind(rnorm(400,50,4), rnorm(400,50,4)))
## plus 1 random dimension
x5 <- cbind(x5, rnorm(nrow(x5)))

clara(x5, 5)
clara(x5, 5, samples = 50)

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')

