.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
  assign("plclust", .Alias(plot.hclust), pos = "package:cluster")
}
