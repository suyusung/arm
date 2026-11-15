## =====================================================================
## Regression tests for standardize()  --  inputs inside a function call
##
## The bug: standardize() rescales each variable of the formula and then
## substitutes the rescaled variable back into the formula.  For
## y ~ log(x) that fits y ~ log(z.x), where z.x = (x - mean(x)) / (2 sd(x))
## is negative for about half the rows: log() returns NaN, lm() drops those
## rows, and the "standardized" model is fitted on a different sample with a
## different meaning -- with only a "NaNs produced" warning.  Below, 64 of
## 100 rows were dropped.
##
## The fix keeps the substitution (y ~ x + I(x^2) still becomes
## z.x + I(z.x^2), a valid model) but compares the number of observations
## of the refit with the original fit's and stops when they differ, naming
## the variable and pointing to `unchanged`.
## =====================================================================

library(arm)

set.seed(5)
d <- data.frame(x = rlnorm(100, 3, 0.5), w = rnorm(100, 10, 2))
d$y <- 2 * log(d$x) + 0.3 * d$w + rnorm(100)
zs <- function(v) (v - mean(v)) / (2 * sd(v))

## ---- log(x): stop instead of silently dropping rows --------------------
m <- lm(y ~ log(x), data = d)
res <- tryCatch(suppressWarnings(standardize(m)), error = function(e) e)
stopifnot("standardize(y ~ log(x)) silently dropped rows instead of stopping" =
            inherits(res, "error"))
stopifnot("the error does not name x and point to 'unchanged'" =
            grepl("rescaling x,", conditionMessage(res), fixed = TRUE),
          grepl("unchanged", conditionMessage(res), fixed = TRUE))

## ---- unchanged = "x" keeps log(x) and every row ------------------------
s <- standardize(m, unchanged = "x")
stopifnot("standardize(y ~ log(x), unchanged = 'x') dropped rows" = nobs(s) == nobs(m),
          "standardize(y ~ log(x), unchanged = 'x') changed the fit" =
            isTRUE(all.equal(unname(coef(s)), unname(coef(m)))))

## ---- I(x^2) is still standardized then squared -------------------------
s2 <- standardize(lm(y ~ x + I(x^2), data = d))
zx <- zs(d$x)
hand <- lm(d$y ~ zx + I(zx^2))
stopifnot("standardize(y ~ x + I(x^2)) no longer gives z.x + I(z.x^2)" =
            identical(names(coef(s2)), c("(Intercept)", "z.x", "I(z.x^2)")),
          isTRUE(all.equal(unname(coef(s2)), unname(coef(hand)))))

## ---- factor() is mapped back, as before --------------------------------
d$g <- rep(1:4, 25)
s3 <- standardize(lm(y ~ w + factor(g), data = d))
zw <- zs(d$w)
hand3 <- lm(d$y ~ zw + factor(d$g))
stopifnot("standardize() with factor() changed" =
            identical(names(coef(s3)), c("(Intercept)", "z.w", paste0("factor(g)", 2:4))),
          isTRUE(all.equal(unname(coef(s3)), unname(coef(hand3)))))
