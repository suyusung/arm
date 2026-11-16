## =====================================================================
## Regression test: bayespolr's `start`  --  issue #54 (minor)
##
## The bug: the help page documents start as c(coefficients, zeta), but the
## value was handed straight to optim(), which works in c(coefficients,
## theta) where the cutpoints are
##
##   zeta <- cumsum(c(theta[1], exp(theta[-1])))
##
## so a start written as documented was read in the wrong parameterisation
## and silently fitted a different model.  Nothing warned, and a zeta that
## was not increasing (impossible, but nothing checked) was accepted.
## =====================================================================

library(arm)

set.seed(3)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- factor(cut(runif(n) + 0.8 * d$x, breaks = c(-Inf, -0.5, 0.5, Inf),
                  labels = c("a", "b", "c")), ordered = TRUE)

fit <- bayespolr(y ~ x, data = d, method = "cloglog")

## ---- a start in the documented form reproduces the fit ------------------
s <- c(coef(fit), fit$zeta)
stopifnot("the documented start is not c(coefficients, zeta)" =
            length(s) == length(coef(fit)) + length(fit$zeta))
fit2 <- bayespolr(y ~ x, data = d, method = "cloglog", start = s)
stopifnot("a start in the documented parameterisation does not reproduce the fit" =
            isTRUE(all.equal(unname(coef(fit2)), unname(coef(fit)), tolerance = 1e-3)) &&
            isTRUE(all.equal(unname(fit2$zeta), unname(fit$zeta), tolerance = 1e-3)))

## ---- a zeta that is not increasing is not a start ----------------------
## (it cannot be one: the cutpoints have to increase)
bad <- c(coef(fit), rev(fit$zeta))
e <- tryCatch(bayespolr(y ~ x, data = d, method = "cloglog", start = bad),
              error = function(e) conditionMessage(e))
stopifnot("a decreasing zeta in `start` is accepted silently" =
            is.character(e) && grepl("increasing", e))

## ---- and the old wrong-length errors still fire -------------------------
e2 <- tryCatch(bayespolr(y ~ x, data = d, method = "cloglog", start = c(1, 2)),
               error = function(e) conditionMessage(e))
stopifnot("a start of the wrong length is not rejected" =
            is.character(e2) && grepl("not of the correct length", e2))
