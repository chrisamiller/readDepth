
.onLoad <- function(libname, pkgname) {
  if( !require(methods) ) stop("we require methods for package readDepth")
  initRdClass()
  library("doMC")
  registerDoMC()
}
