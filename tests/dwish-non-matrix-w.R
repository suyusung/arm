## =====================================================================
## Regression tests for dwish()  --  a non-matrix W  --  d8
##
## The bug: dwish() coerces a non-matrix W with `S <- matrix(W)` instead of
## `W <- matrix(W)`, overwriting S and leaving W a plain vector.  The next
## check, nrow(W) != ncol(W), is then NULL != NULL and R stops with
## "argument is of length zero".
##
## dwish() is internal and unused in the package, so there is no
## user-facing call to exercise; this test reaches it directly through
## arm::: and checks the value it returns.
## =====================================================================

library(arm)

## W and S as plain scalars: a 1 x 1 Wishart
d1 <- tryCatch(arm:::dwish(4, 6, 4),
               error = function(e) conditionMessage(e))
stopifnot("dwish() errors when W is not a matrix (S <- matrix(W) typo)" =
            is.numeric(d1) && length(d1) == 1 && is.finite(d1))

## the answer must not depend on whether W was passed as a matrix
stopifnot("dwish() disagrees on scalar and matrix W" =
            isTRUE(all.equal(d1, arm:::dwish(matrix(4), 6, matrix(4)))))

## the 1 x 1 Wishart density, by hand
v <- 6
k <- 1
hand <- 4^(-v / 2) * 4^((v - k - 1) / 2) * exp(-1 / 2 * 1) /
        (gamma((v + 1 - 1) / 2) * 2^(v * k / 2) * pi^(k * (k - 1) / 4))
stopifnot("dwish() returns the wrong value" =
            isTRUE(all.equal(as.numeric(d1), hand)))
