## This came from a bug report on R-help by ge yreyt <tothri2000@yahoo.ca>
## Date: Mon, 9 Jun 2003 16:06:53 -0400 (EDT)
library(cluster)

data(iris)
mdist <- as.dist(1 - cor(t(iris[,1:4])))#dissimlarity
## this is always the same:
hc <- diana(mdist, diss = TRUE, stand = FALSE)

maxk <- 15                # at most 15 clusters
silh.wid <- numeric(maxk)  # myind[k] := the silh.value for k clusters
silh.wid[1] <- NA # 1-cluster: silhouette not defined

for(k in 2:maxk) {
    cat("\n", k,":\n==\n")
    k.gr <- cutree(as.hclust(hc), k = k)
    cat("grouping table: "); print(table(k.gr))
    si <- silhouette(k.gr, mdist)
    cat("silhouette:\n"); print(summary(si))
    silh.wid[k] <- summary(si)$avg.width
    ##      ===
}

# the widths:
silh.wid
#select the number of k clusters with the largest si value :
(myk <- which.min(silh.wid))

postscript(file="silhouette-ex.ps")
## MM:  plot to see how the decision is made
plot(silh.wid, type = 'b', col= "blue", xlab = "k")
axis(1, at=myk, col.axis= "red", font.axis= 2)

dev.off()
