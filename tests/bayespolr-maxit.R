## =====================================================================
## Regression test: bayespolr's `maxit`  --  issue #54 (minor)
##
## The bug: maxit is documented as the iteration limit of the fit, but it
## reached only the bayesglm.fit() call that produces the starting values.
## The main optimisation runs through optim(), which has its own maxit
## default of 100, so maxit = 1 gave a fully converged fit.
## =====================================================================

library(arm)

set.seed(3)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- factor(cut(runif(n) + 0.8 * d$x, breaks = c(-Inf, -0.5, 0.5, Inf),
                  labels = c("a", "b", "c")), ordered = TRUE)

full <- bayespolr(y ~ x, data = d, method = "cloglog")
one <- suppressWarnings(bayespolr(y ~ x, data = d, method = "cloglog", maxit = 1))

## one iteration cannot be the optimum
stopifnot("maxit = 1 still converges to the optimum, so maxit does not reach the fit" =
            one$deviance > full$deviance * 1.2)

## control passed through `...` still wins over maxit
one.ctl <- suppressWarnings(bayespolr(y ~ x, data = d, method = "cloglog",
                                      maxit = 100, control = list(maxit = 1)))
stopifnot("control = list(maxit = ) in `...` is ignored" =
            one.ctl$deviance > full$deviance * 1.2)

## the default is unchanged: a default fit reaches the optimum
stopifnot("the default fit no longer converges" =
            isTRUE(all.equal(full$deviance, bayespolr(y ~ x, data = d,
                                                      method = "cloglog", maxit = 100)$deviance)))
