## =====================================================================
## Regression test: a rank-deficient fit with a flat prior  --  issue #54 (minor)
##
## The bug: with an aliased column (x3 is x1) and a prior flat enough that
## its pseudo-rows carry no weight, the augmented least-squares fit returns
## non-finite coefficients at the first iteration and the loop breaks.  The
## rank-deficiency fix-up after the loop then indexed `coef`, which had
## never been assigned on that path, so the name still resolved to the base
## function coef() and bayesglm() died with
##
##   Error: object of type 'closure' is not subsettable
##
## glm() returns NA for the aliased coefficient instead of erroring.
##
## The fix initialises `coef` (and `dev`) before the loop, so the same
## situation now returns a fit with NA coefficients and the usual warnings.
## Deciding *which* coefficient is aliased, and returning numbers for the
## others as glm does, needs the design to be reduced before fitting; that
## is left for later and is why this test pins the crash, not the estimate.
## =====================================================================

library(arm)

set.seed(1)
n <- 100
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))
d$x3 <- d$x1                     # aliased

fit <- try(bayesglm(y ~ x1 + x2 + x3, binomial, d, prior.scale = 1e8), silent = TRUE)
stopifnot("a rank-deficient fit with a flat prior still errors (issue #54 is back)" =
            !inherits(fit, "try-error"))

## it reports failure rather than inventing numbers
stopifnot("the degenerate fit claims to have converged" = !fit$converged)
stopifnot("the degenerate fit returned finite coefficients" =
            all(is.na(coef(fit))))

## and the same design with the default prior still fits properly: the prior
## is what makes the augmented system non-singular, so nothing changed there
ok <- bayesglm(y ~ x1 + x2 + x3, binomial, d)
stopifnot("the default prior no longer regularises a rank-deficient design" =
            all(is.finite(coef(ok))))
stopifnot("the default-prior fit no longer converges" = isTRUE(ok$converged))
