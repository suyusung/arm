## =====================================================================
## Regression tests for balance() on a matrix without column names  --  d3
##
## The bug: balance() reads each covariate with rawdata[, covnames[i]].
## With a matrix that has no colnames, balance() generates generic names
## v1, v2, ... and that character subscript is out of bounds, so the call
## stopped with "subscript out of bounds" even though it advertises
## generic names.
##
## The fix indexes the loop by position, as the rest of the function does.
## =====================================================================

library(arm)

set.seed(5)
n <- 100
M <- cbind(rnorm(n, 40, 10), rnorm(n, 50, 15))
colnames(M) <- NULL
z <- rbinom(n, 1, 0.5)

res <- tryCatch({capture.output(B <- balance(M, z, rep(1, n))); B},
                error = function(e) e)
stopifnot("balance() errors on a matrix without column names" =
            !inherits(res, "error"))
stopifnot("balance() did not generate generic covariate names" =
            identical(res$covnames, c("v1", "v2")))
stopifnot("balance() on an unnamed matrix returned the wrong control mean" =
            isTRUE(all.equal(unname(res$diff.means.raw[, "control"]),
                             unname(apply(M[z == 0, , drop = FALSE], 2, mean)))))
