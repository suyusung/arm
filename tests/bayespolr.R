## =====================================================================
## Regression tests for bayespolr()
##
## scaled = TRUE (the default) must rescale the prior for each coefficient
## as ?bayespolr describes: prior.scale / range(x) for a two-valued
## predictor and prior.scale / (2 * sd(x)) otherwise.  The block that did
## this was guarded by `xint > 1`, where xint is the column position of the
## intercept (1), so it never ran and every fit used the unscaled prior.
## =====================================================================

library(arm)

set.seed(1)
n <- 400
d <- data.frame(x1 = rbinom(n, 1, 0.5) * 10, x2 = rnorm(n, 0, 0.25))
d$y <- cut(0.3 * d$x1 + 2 * d$x2 + rlogis(n), c(-Inf, 0.5, 2.5, Inf),
           labels = c("L", "M", "H"))

## ---- scaled = TRUE applies the documented per-coefficient scales --------
fit <- bayespolr(y ~ x1 + x2, data = d)
expected <- c(2.5 / diff(range(d$x1)), 2.5 / (2 * sd(d$x2)))
stopifnot("scaled = TRUE did not rescale the prior" =
            isTRUE(all.equal(unname(fit$prior.scale), expected)))

## ---- scaled = FALSE keeps the given scale, expanded per coefficient ------
fit0 <- bayespolr(y ~ x1 + x2, data = d, scaled = FALSE)
stopifnot("scaled = FALSE changed the prior scale" =
            isTRUE(all.equal(unname(fit0$prior.scale), c(2.5, 2.5))))

## ---- and the scaling changes the estimates -------------------------------
stopifnot("scaled had no effect on the estimates" =
            max(abs(coef(fit) - coef(fit0))) > 1e-3)

## ---- an intercept-only model still fits --------------------------------
stopifnot(length(coef(bayespolr(y ~ 1, data = d))) == 0)
