if(FALSE)
     library(cluster, lib.loc="/u/maechler/R/Pkgs/cluster.Rcheck")
data(ruspini)
ru4 <- ruspini[c(1:2,21:22, 45:47),]
dist(ru4, "manhattan")

## This  __still__  shows the "** dysta2(): ... is OUT"  bug:
c4 <- clara(ru4, k=3, met="manhattan", sampsize = 6, trace = 4)

## the result is fine:
clusplot(c4)
