.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
  assign("plclust", plot.hclust,#.Alias
         pos = "package:cluster")
}
