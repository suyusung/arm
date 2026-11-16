## =====================================================================
## Regression test: a bad `start` gives a bad error  --  issue #54 (minor)
##
## The bug: bayesglm.fit() compared `start` with 0 before checking its
## length,
##
##   if (length(start) != nvars) {
##       if (start == 0 & length(start) == 1) { ... } else stop(...)
##
## so a `start` of the wrong length that is not a single zero made the
## condition a vector and R reported
##
##   Error in if: the condition has length > 1
##
## instead of naming the length and the coefficients.  glm() accepts a
## single 0 as "start everything at zero" and errors otherwise.
## =====================================================================

library(arm)

set.seed(1)
n <- 100
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))

## the wrong length must say so -- and say what it wanted
e <- tryCatch(bayesglm(y ~ x1 + x2, binomial, d, start = c(0, 0)),
              error = function(e) conditionMessage(e))
stopifnot("a start of the wrong length does not report the length it wanted" =
            is.character(e) && grepl("length of 'start' should equal", e))
stopifnot("the error does not name the coefficients" =
            grepl("x1", e) && grepl("x2", e))

## a non-zero start of the right length still fits
fit <- bayesglm(y ~ x1 + x2, binomial, d, start = c(0, 0, 0))
stopifnot("a start of the correct length should be accepted" =
            length(coef(fit)) == 3L && all(is.finite(coef(fit))))

## the single-zero shorthand still works, as in glm()
fit0 <- bayesglm(y ~ x1 + x2, binomial, d, start = 0)
stopifnot("the single-zero start shorthand stopped working" =
            length(coef(fit0)) == 3L && all(is.finite(coef(fit0))))
