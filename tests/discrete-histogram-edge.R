## =====================================================================
## Regression tests for discrete.histogram() edge cases  --  p2
##
## The bugs:
##  (a) a single distinct numeric value left gaps = NA, so bar.width was
##      NA and plot() failed with "need finite 'xlim' values".
##  (b) a factor with an unused level tabulated the empty level, so
##      length(x) != length(prob) and the function stopped with "Length of
##      'x' must be the same as the length of 'prob'".
##  (c) ylab = was passed through ... while the function also passed its
##      own ylab to plot(), giving "formal argument 'ylab' matched by
##      multiple actual arguments" -- the user could not set the label.
##
## The fix: the gaps (and so bar.width) are computed only when there is
## more than one value, a factor is droplevels()-ed before tabulating, and
## a ylab in ... overrides the internal "Probability"/"Count" default.
## =====================================================================

library(arm)

pdf(NULL)   # one null device for all the drawing below

ok <- function(expr)
  tryCatch({expr; TRUE}, error = function(e) FALSE, warning = function(w) FALSE)

## (a) one distinct value
stopifnot("discrete.histogram() errors on a single distinct value" =
            ok(discrete.histogram(rep(5, 10))))
stopifnot("discrete.histogram() errors on a single distinct value with freq = TRUE" =
            ok(discrete.histogram(rep(5, 10), freq = TRUE)))

## (b) a factor with a level that is never observed
f <- factor(c("a", "b", "a", "c"), levels = c("a", "b", "c", "d"))
stopifnot("discrete.histogram() errors on a factor with an unused level" =
            ok(discrete.histogram(f)))
stopifnot("discrete.histogram() errors on a factor with an unused level with freq = TRUE" =
            ok(discrete.histogram(f, freq = TRUE)))

## (c) ylab =
stopifnot("discrete.histogram() errors on ylab=" =
            ok(discrete.histogram(c(1, 2, 3, 1, 2, 1), ylab = "my label")))

## the user's ylab is the one drawn, not the internal default.  The trace is
## on the plot() generic itself: an S3 method cached by an earlier plot()
## call would not be intercepted.
seen <- new.env()
seen$ylab <- NULL
invisible(suppressMessages(
  trace("plot", where = asNamespace("graphics"),
        quote(seen$ylab <- c(seen$ylab, list(...)$ylab)), print = FALSE)))
discrete.histogram(c(1, 2, 3, 1, 2, 1), ylab = "my label")
stopifnot("discrete.histogram() draws its own ylab instead of the user's" =
            identical(seen$ylab, "my label"))

## without ylab= the default label is still Probability (freq = FALSE)
seen$ylab <- NULL
discrete.histogram(c(1, 2, 3, 1, 2, 1))
invisible(suppressMessages(untrace("plot", where = asNamespace("graphics"))))
stopifnot("discrete.histogram() lost its default ylab" =
            identical(seen$ylab, "Probability"))

dev.off()
