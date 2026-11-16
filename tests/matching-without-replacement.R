## =====================================================================
## Regression test for matching(replace = FALSE) with too few controls
##
## The bug: matching without replacement pairs each treated unit with a
## control not chosen before.  With more treated units than controls the
## pool runs out; matching() then picked index NA, returned NA in
## match.ind for the unmatched treated units, still reported cnts = 1 for
## them, and corrupted `pairs` -- all without a warning.  For
## z = c(1, 1, 1, 0) it returned match.ind = c(4, NA, NA, 1) and
## cnts = c(1, 1, 1, 1), so any weighting by cnts counted three treated
## units as matched to one control.
##
## The fix stops with an error that says why and points to replace = TRUE.
## =====================================================================

library(arm)

res <- tryCatch(matching(z = c(1, 1, 1, 0), score = c(0.1, 0.2, 0.3, 0.4)),
                error = function(e) e)
stopifnot("matching(replace = FALSE) with too few controls did not stop" =
            inherits(res, "error"))
stopifnot("the error does not explain the shortage of controls" =
            grepl("controls", conditionMessage(res)))

## equal numbers still work, and every unit is used exactly once
m <- matching(z = c(1, 1, 0, 0), score = c(0.1, 0.4, 0.35, 0.05))
stopifnot("matching with nt == nc no longer matches everyone" =
            identical(m$match.ind, c(4, 3, 2, 1)),
          identical(m$cnts, c(1, 1, 1, 1)))

## with replacement, a shortage of controls is fine
m <- matching(z = c(1, 1, 1, 0), score = c(0.1, 0.2, 0.3, 0.4), replace = TRUE)
stopifnot("matching(replace = TRUE) with one control changed" =
            identical(m$cnts, c(1, 1, 1, 3)))
