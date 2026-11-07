.onAttach <- function(libname, pkgname) {
  desc <- packageDescription(pkgname, lib.loc = libname)
  if (!is.list(desc))
    desc <- suppressWarnings(packageDescription(pkgname))
  ver <- if (is.list(desc)) desc$Version else getNamespaceVersion(pkgname)
  builddate <- if (is.list(desc) && !is.null(desc$Date)) desc$Date else "unknown"
  packageStartupMessage(paste("\n", pkgname, " (Version ", ver, ", built: ", builddate, ")\n", sep = ""))
  packageStartupMessage("Working directory is ", getwd(), "\n")
}
