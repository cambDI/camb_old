.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  packageStartupMessage('Loading camb (chemistry aware model builder)')
  packageStartupMessage('written by Daniel Murrell (dsmurrell@gmail.com) and Isidro Cortes <isidrolauscher@gmail.com>')
  packageStartupMessage('For a complete list of package functions, use ls("package:camb")')
}
