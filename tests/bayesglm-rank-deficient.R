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
## The fix initialises `coef` (and `dev`) before the loop -- so the rank fix-up
## no longer indexes a function -- and then hands the problem to glm.fit(),
## which detects the aliased columns, returns NA for them and estimates the
## rest.  A flat prior means the prior is not doing anything, so glm's answer
## is the right limit, and it is what this test pins.
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

## it now gives what glm() gives: NA for the aliased coefficient only, and
## numbers for the ones the data can estimate
g <- glm(y ~ x1 + x2 + x3, binomial, d)
stopifnot("the degenerate fit does not fall back to glm's coefficients" =
            isTRUE(all.equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-6)))
stopifnot("the aliased coefficient is not the only NA" =
            identical(is.na(coef(fit)), is.na(coef(g))))
stopifnot("the aliased column is not the missing one" =
            is.na(coef(fit))[["x3"]] && sum(is.na(coef(fit))) == 1L)

## and the same design with the default prior still fits properly: the prior
## is what makes the augmented system non-singular, so nothing changed there
ok <- bayesglm(y ~ x1 + x2 + x3, binomial, d)
stopifnot("the default prior no longer regularises a rank-deficient design" =
            all(is.finite(coef(ok))))
stopifnot("the default-prior fit no longer converges" = isTRUE(ok$converged))
