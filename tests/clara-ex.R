#### These are *NOT* compared with output in the released version of
###  `cluster'  currently

library(cluster)
if(paste(R.version$major, R.version$minor, sep=".") >= 1.7)  RNGversion(1.6)

data(xclara)
## Try 100 times *different* random samples -- for reliability:
set.seed(421)# (reproducibility)
cl <- matrix(NA,nrow(xclara), 100)
for(i in 1:100)
    cl[,i] <- clara(xclara, 3, rngR = TRUE, keep.data=FALSE, trace=1)$cluster
rcl <- apply(cl,1, range)
## those that are not always in same cluster (5 out of 3000 for this seed):
(iDoubt <- which(rcl[1,] != rcl[2,]))
if(length(iDoubt)) { # (not for all seeds)
  tabD <- apply(cl[iDoubt, , drop = FALSE], 1, table)
  colnames(tabD) <- format(iDoubt)
  t(tabD) # how many times in which clusters
}
