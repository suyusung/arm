## =====================================================================
## Regression test for .onAttach  --  issue #25
##
## The bug: .onAttach derived the library location from
## dirname(system.file(package = "arm")), which resolves through
## .libPaths() rather than the library the package was actually loaded
## from, and then called packageDescription(...)$Version without checking
## anything.  When the lookup failed, packageDescription() returned NA --
## a logical, not a list -- and `NA$Version` aborted library(arm) with
##
##   Error: $ operator is invalid for atomic vectors
##
## The crash itself cannot be reproduced here: it needs an R whose
## find.package() does not consult already-loaded namespaces, and R 4.6
## does, which masks the failure.  What can be pinned is the guard -- and
## the requirement that the hook uses the arguments R hands it instead of
## looking the package up again.
## =====================================================================

library(arm)

on_attach <- get(".onAttach", envir = asNamespace("arm"))

## collect startup messages without printing them
catch_msgs <- function(expr) {
  msgs <- character()
  withCallingHandlers(expr, packageStartupMessage = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  msgs
}

## ---- it must take libname/pkgname, not re-derive the location ---------
stopifnot("'.onAttach' no longer takes libname and pkgname" =
            all(c("libname", "pkgname") %in% names(formals(on_attach))))

src <- paste(deparse(body(on_attach)), collapse = "\n")
stopifnot("'.onAttach' resolves the library through system.file() again" =
            !grepl("system.file", src, fixed = TRUE))

## ---- an unreachable library must not abort library(arm) ---------------
empty <- tempfile("nolib"); dir.create(empty)
msgs <- tryCatch(
  catch_msgs(suppressWarnings(on_attach(libname = empty, pkgname = "arm"))),
  error = function(e)
    stop(".onAttach() errored on an unreachable library: ", conditionMessage(e)))
stopifnot("'.onAttach' printed no startup message" =
            any(grepl("Version", msgs, fixed = TRUE)))

## packageDescription() really does return NA here -- this is the input
## that used to be fatal
d <- suppressWarnings(utils::packageDescription("arm", lib.loc = empty))
stopifnot("packageDescription() unexpectedly found arm in an empty library" =
            !is.list(d))

## ---- and with the real library it reports that library's version ------
## (the DESCRIPTION Version, e.g. "1.15-4", not packageVersion()'s "1.15.4")
lib  <- dirname(find.package("arm"))
msgs <- catch_msgs(on_attach(libname = lib, pkgname = "arm"))
stopifnot("'.onAttach' did not report the installed version" =
            any(grepl(utils::packageDescription("arm")$Version, msgs, fixed = TRUE)))
