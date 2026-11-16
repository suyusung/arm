## =====================================================================
## Regression tests for the `ratio` column of balance()  --  d1
##
## The bug: the raw balance table's `ratio` column (column 5 of
## diff.means.raw) was never assigned, so every entry printed as NA ("--")
## in print(balance, combined = TRUE) while the matched table carried a
## number.  The comment above the matched assignment also described the
## opposite ratio from the one the code computes.
##
## The fix computes the raw ratio the same way as the matched one --
## sd(control)/sd(treated), for continuous covariates -- and corrects the
## comment.
## =====================================================================

library(arm)

set.seed(1)
n <- 200
z <- rbinom(n, 1, 0.5)
X <- data.frame(age = rnorm(n, 40, 10), female = rbinom(n, 1, 0.5),
                inc = rnorm(n, 50, 15))
invisible(capture.output(B <- balance(X, z, rep(1, n))))

cont <- names(X)[B$binary == 0]
stopifnot("test setup: expected two continuous covariates" = length(cont) == 2)

## ---- the raw ratio is computed -----------------------------------------
stopifnot("balance() did not compute the raw ratio column" =
            !any(is.na(B$diff.means.raw[cont, "ratio"])))

## ---- and it is sd(control)/sd(treated) ---------------------------------
hand <- sqrt(vapply(cont, function(v) var(X[z == 0, v]) / var(X[z == 1, v]),
                    numeric(1)))
stopifnot("balance() raw ratio is not sd(control)/sd(treated)" =
            isTRUE(all.equal(unname(B$diff.means.raw[cont, "ratio"]), unname(hand))))

## the raw and matched tables must use the same direction
stopifnot("balance() raw and matched ratios use opposite directions" =
            isTRUE(all.equal(unname(B$diff.means.raw[cont, "ratio"]),
                             unname(B$diff.means.matched[cont, "ratio"]))))

## ---- the combined print no longer shows the raw ratio as missing -------
tbl <- capture.output(res <- print(B, combined = TRUE))
stopifnot("print(balance, combined = TRUE) still shows the raw ratio as '--'" =
            !any(is.na(res[cont, 3])))
