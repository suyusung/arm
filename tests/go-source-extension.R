## =====================================================================
## Regression tests for go()  --  the source-file extension  --  d9
##
## The bug: go('name') appends ".R", but the help page and the header
## comment in R/go.R said it would run source('name.r').  R scripts are
## ".R", which is what the code does, so the documentation was corrected.
##
## This test pins both sides of that agreement: go() resolves `.R`, and the
## installed help page says `.R` rather than `.r`.
## =====================================================================

library(arm)

## ---- go('name') finds name.R -------------------------------------------
d <- tempfile("go")
dir.create(d)
prog <- file.path(d, "prog")
writeLines("GO_SOURCED_MARKER <- TRUE", paste0(prog, ".R"))
invisible(capture.output(go(prog)))
stopifnot("go('name') did not source name.R (the extension the code appends)" =
            isTRUE(GO_SOURCED_MARKER))

## ---- the help page documents the same extension ------------------------
h <- utils:::.getHelpFile(help("go", package = "arm"))
txt <- paste(capture.output(tools::Rd2txt(h)), collapse = "\n")
stopifnot("man/go.Rd says '.r' while go() appends '.R'" =
            !grepl("myprog.r", txt, fixed = TRUE))
stopifnot("man/go.Rd no longer documents the .R fallback" =
            grepl("myprog.R", txt, fixed = TRUE))
