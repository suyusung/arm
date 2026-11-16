## =====================================================================
## Regression test: the `Warning` argument of bayesglm()  --  issue #54 (minor)
##
## The bug: `Warning` was an argument of bayesglm() and bayesglm.fit(), it
## was passed along at every call, and the help page promised it "shows the
## error messages for non-convergence and separation" -- but nothing ever
## read it, so Warning = FALSE was silently ignored.
## =====================================================================

library(arm)

## perfectly separated data: the fit runs to the boundary and says so
d <- data.frame(x = 1:10, y = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))

warnings_of <- function(expr) {
  ws <- character()
  withCallingHandlers(expr,
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  ws
}

w.on <- warnings_of(bayesglm(y ~ x, binomial, d, prior.scale = 1e6, Warning = TRUE))
stopifnot("the separated fit no longer warns at all, so this test proves nothing" =
            length(w.on) > 0)

w.off <- warnings_of(bayesglm(y ~ x, binomial, d, prior.scale = 1e6, Warning = FALSE))
stopifnot("Warning = FALSE does not suppress the fit's warnings (issue #54 is back)" =
            length(w.off) == 0L)

## the default is TRUE, so nothing changes for callers who do not pass it
w.def <- warnings_of(bayesglm(y ~ x, binomial, d, prior.scale = 1e6))
stopifnot("the default stopped warning" = length(w.def) > 0)

## a fit that warns is still usable
suppressWarnings(f <- bayesglm(y ~ x, binomial, d, prior.scale = 1e6, Warning = FALSE))
stopifnot("the suppressed fit did not return usable coefficients" =
            length(coef(f)) == 2L)
