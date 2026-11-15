## =====================================================================
## Regression test for print.balance(combined = TRUE)
##
## The bug: the combined branch of print.balance() read x$diff.matched.raw,
## a component that balance() never creates (it is x$diff.means.matched).
## cbind() of the raw table with NULL has 5 columns, so selecting columns
## c(4, 9, 5, 10) always failed with "subscript out of bounds":
## print(..., combined = TRUE) could never work.
##
## The fix reads x$diff.means.matched, so the combined table shows the raw
## and matched standardized differences and sd ratios side by side.
## =====================================================================

library(arm)

set.seed(8)
n <- 200
z <- rbinom(n, 1, 0.4)
X <- data.frame(age = rnorm(n, 40 + 5 * z, 10), female = rbinom(n, 1, 0.5),
                inc = rnorm(n, 50, 15))
w <- ifelse(z == 1, 1, runif(n))
invisible(capture.output(B <- balance(X, z, w)))

out <- tryCatch(capture.output(res <- print(B, combined = TRUE)),
                error = function(e) e)
stopifnot("print(balance, combined = TRUE) failed (diff.matched.raw typo)" =
            !inherits(out, "error"))

## print() returns the table it printed: raw and matched std diffs, then
## raw and matched ratios
stopifnot("combined table does not show the matched std diffs" =
            identical(unname(res[, 2]), unname(round(B$diff.means.matched[, 4], 2))))
stopifnot("combined table does not show the raw std diffs" =
            identical(unname(res[, 1]), unname(round(B$diff.means.raw[, 4], 2))))
stopifnot("combined table does not show the matched ratios" =
            identical(unname(res[, 4]), unname(round(B$diff.means.matched[, 5], 2))))
