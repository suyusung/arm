## =====================================================================
## Regression tests for bayespolr(): deviance and fitted values
##
## With a flat prior (prior.scale = Inf, prior.df = Inf) and no prior
## counts for the cutpoints, bayespolr() is MASS::polr(), so its
## deviance and fitted probabilities must match polr's.
##
## The deviance used to be 2 * the optimised objective, which also holds
## the log prior and the prior pseudo-counts, and the fitted values
## left out the offset.
## =====================================================================

library(arm)

set.seed(1)
n <- 400
d <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rnorm(n), off = rnorm(n))
d$y <- cut(0.5 * d$x1 + d$x2 + d$off + rlogis(n), c(-Inf, 0, 1.5, Inf),
           labels = c("L", "M", "H"))

ref <- MASS::polr(y ~ x1 + x2 + offset(off), data = d, Hess = TRUE)
fit <- bayespolr(y ~ x1 + x2 + offset(off), data = d, prior.scale = Inf,
                 prior.df = Inf, prior.counts.for.bins = 0)

stopifnot("estimates differ from polr under a flat prior" =
            isTRUE(all.equal(c(coef(fit), fit$zeta), c(coef(ref), ref$zeta),
                             tolerance = 1e-4)))
stopifnot("deviance includes more than the data" =
            isTRUE(all.equal(deviance(fit), deviance(ref), tolerance = 1e-6)))
stopifnot("fitted values ignore the offset" =
            isTRUE(all.equal(unname(fitted(fit)), unname(fitted(ref)),
                             tolerance = 1e-4)))

## ---- without an offset or predictors the fitted values still work -------
fit1 <- bayespolr(y ~ 1, data = d)
stopifnot(isTRUE(all.equal(unname(colMeans(fitted(fit1))),
                           as.numeric(prop.table(table(d$y))),
                           tolerance = 0.01)))
