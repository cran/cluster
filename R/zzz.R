.First.lib <- function(lib, pkg) {
  require(mva)
  library.dynam("cluster", pkg, lib)
}
