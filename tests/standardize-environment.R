## =====================================================================
## Regression tests for standardize()  --  where the data are looked up
##
## The bug: standardize() refits the model from its call, and it looked the
## data frame (or the bare variables) up with get() from inside arm, i.e.
## through arm's namespace and then the global environment -- not where the
## model was fitted.  A model fitted inside a function on a local data
## frame `d` was therefore refitted on a *global* `d` if one existed,
## silently, and failed with a confusing error if none did.  In the
## example below the local slope is 2.94; standardize() returned a slope of
## 0.077 for z.x, fitted on the unrelated global data.
##
## The fix evaluates the call, the formula (as the issue #24 fix does) and
## the data in the environment of the model's formula, which is how
## model.frame.lm() finds the data again.
## =====================================================================

library(arm)

## a global data frame that shares the local data's name
set.seed(10)
d <- data.frame(y = rnorm(50), x = rnorm(50, 100, 10))
x <- rnorm(50, 100, 10)

by_hand <- function(y, x) unname(coef(lm(y ~ I((x - mean(x)) / (2 * sd(x))))))

## ---- data = d, fitted inside a function --------------------------------
local_data <- function() {
  set.seed(11)
  d <- data.frame(x = rnorm(50, 5, 1))
  d$y <- 3 * d$x + rnorm(50)
  list(s = coef(standardize(lm(y ~ x, data = d))), hand = by_hand(d$y, d$x))
}
r <- local_data()
stopifnot("standardize() refitted on a global data frame, not the model's data" =
            isTRUE(all.equal(unname(r$s), r$hand)))

## ---- no data argument: bare local vectors ------------------------------
local_vectors <- function() {
  set.seed(12)
  x <- rnorm(50, 5, 1)
  y <- 3 * x + rnorm(50)
  list(s = coef(standardize(lm(y ~ x))), hand = by_hand(y, x))
}
r <- local_vectors()
stopifnot("standardize() used global variables, not the model's variables" =
            isTRUE(all.equal(unname(r$s), r$hand)))

## ---- a formula held in a local object (issue #24 inside a function) ---
local_formula <- function() {
  set.seed(13)
  d <- data.frame(x = rnorm(50, 5, 1))
  d$y <- 3 * d$x + rnorm(50)
  form <- y ~ x
  list(s = coef(standardize(glm(form, data = d))), hand = by_hand(d$y, d$x))
}
r <- local_formula()
stopifnot("standardize() failed on a local formula object fitted on local data" =
            isTRUE(all.equal(unname(r$s), r$hand)))

## ---- standardize() leaves the caller's objects alone -------------------
stopifnot("standardize() modified the global data frame" =
            identical(names(d), c("y", "x")),
          !exists("z.x", envir = globalenv()))
