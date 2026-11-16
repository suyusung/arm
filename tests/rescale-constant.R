## =====================================================================
## Regression tests for rescale()  --  a constant input  --  d6
##
## The bug: a variable that takes on a single value has sd 0, and rescale()
## divided by 2 sd, so rescale(rep(5, 10)) returned NaN.
##
## The fix centres a constant (subtracting its mean, which is itself), so
## the rescaled vector is 0 everywhere, with NAs preserved.
## =====================================================================

library(arm)

r <- rescale(rep(5, 10))
stopifnot("rescale() of a constant is NaN" = all(is.finite(r)))
stopifnot("rescale() of a constant is not centred at 0" =
            isTRUE(all.equal(unname(r), rep(0, 10))))

stopifnot("rescale() of a constant did not preserve NAs" =
            is.na(rescale(c(5, NA, 5))[2]))
stopifnot("rescale() of a constant factor is NaN" =
            all(is.finite(rescale(factor(rep("a", 4))))))

## a variable with spread is untouched: still (x - mean) / (2 sd)
xx <- c(1, 2, 3, 4, 5)
stopifnot("rescale() of a non-constant changed" =
            isTRUE(all.equal(unname(rescale(xx)),
                             (xx - mean(xx)) / (2 * sd(xx)))))
