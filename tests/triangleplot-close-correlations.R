## =====================================================================
## Regression test for triangleplot() on close-together correlations  --  p3
##
## The bug: with details = TRUE the legend breaks were rounded to `digits`
## and `digits` was raised until the labels were distinct, but the search
## ran only four times and then stopped with "Too many digits! Try to
## adjust n.col.legend to get better presentation!".  Correlations that
## differ in the third decimal need four digits, so a matrix of them always
## hit the error.
##
## The fix: the search keeps raising `digits` until the rounded breaks are
## distinct (or 15 digits are reached), so close values draw a legend
## instead of an error.
## =====================================================================

library(arm)

## a lower triangle of correlations that agree to 3 decimal places
m <- matrix(NA, 4, 4)
m[lower.tri(m)] <- c(0.5001, 0.5002, 0.5003, 0.5004, 0.5002, 0.5003)
diag(m) <- 1

labels.seen <- new.env()
labels.seen$l <- list()
invisible(suppressMessages(
  trace("text.default", where = asNamespace("graphics"),
        quote(labels.seen$l <- c(labels.seen$l, list(labels))), print = FALSE)))

pdf(NULL)
r <- tryCatch({triangleplot(m); NULL}, error = function(e) conditionMessage(e))
stopifnot("triangleplot() stops with 'Too many digits!' on close-together correlations" =
            is.null(r))

## the color legend actually shows distinct values
lab <- unlist(labels.seen$l)
stopifnot("triangleplot() drew a legend with duplicate labels" =
            length(lab) > 1 && !anyDuplicated(lab))

## details = FALSE still works
r2 <- tryCatch({triangleplot(m, details = FALSE); NULL},
               error = function(e) conditionMessage(e))
stopifnot("triangleplot() fails on close-together correlations with details = FALSE" =
            is.null(r2))

invisible(suppressMessages(untrace("text.default", where = asNamespace("graphics"))))
dev.off()
