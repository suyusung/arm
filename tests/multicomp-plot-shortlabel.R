## =====================================================================
## Regression test for the multicomp.plot() short-label help page  --  p5
##
## The bug: the help page promised that label.as.shortlabel = TRUE "takes
## the first 2 characters of the label", but the code calls
## abbreviate(label, minlength = 2), which keeps each label unique (labels
## sharing a two-character prefix would otherwise be indistinguishable).
## The code is the defensible side, so the page was corrected.
##
## The test reads the installed help page and checks the description now
## names abbreviate() and no longer promises "first 2 characters".
## =====================================================================

library(arm)

db <- tools::Rd_db("arm")
txt <- paste(capture.output(tools::Rd2txt(db[["multicomp.plot.Rd"]])), collapse = "\n")
## Rd2txt wraps lines, so collapse the wrapping before matching phrases
txt <- gsub("[[:space:]]+", " ", txt)

stopifnot("multicomp.plot.Rd still says label.as.shortlabel takes the 'first 2 characters'" =
            !grepl("first 2", txt, fixed = TRUE))
stopifnot("multicomp.plot.Rd does not describe the abbreviate() short label" =
            grepl("abbreviate", txt, fixed = TRUE))
