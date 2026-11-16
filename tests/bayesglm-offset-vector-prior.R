## =====================================================================
## Regression test: bayesglm() with an offset and per-coefficient priors
##
## The bug: when the model has an offset and an intercept, bayesglm() refits
## an intercept-only model to get the null deviance.  That refit passed on
## prior.mean, prior.scale and prior.df, which are documented as vectors of
## length "number of predictors" -- but the intercept-only model has none,
## so a vector prior (length > 1) stopped with "invalid length for
## prior.scale" (or prior.mean / prior.df).  The same call without an offset
## worked.
##
## Identities: with flat priors the null deviance equals glm()'s; and a
## vector prior whose entries are all equal gives the same fit as the scalar.
## =====================================================================

library(arm)

close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-6))
fit <- function(...) tryCatch(bayesglm(...), error = function(e) NULL)

set.seed(5)
n <- 50
x1 <- rnorm(n)
x2 <- rnorm(n)
e <- runif(n, 1, 5)
y <- rpois(n, e * exp(0.2 + 0.4 * x1 - 0.3 * x2))
g <- glm(y ~ x1 + x2 + offset(log(e)), poisson)

b <- fit(y ~ x1 + x2 + offset(log(e)), poisson,
         prior.mean = c(0, 0), prior.scale = c(Inf, Inf), prior.df = c(Inf, Inf),
         prior.scale.for.intercept = Inf, prior.df.for.intercept = Inf)
stopifnot("offset + vector prior: bayesglm stopped in the null-deviance refit" =
            !is.null(b))
stopifnot("offset + vector flat prior: null deviance differs from glm" =
            close(b$null.deviance, g$null.deviance))
stopifnot("offset + vector flat prior: coefficients differ from glm" =
            close(coef(b), coef(g)))

## the offset= argument goes through the same refit
## (called directly: model.frame() cannot see offset= through a wrapper's ...)
b <- tryCatch(bayesglm(y ~ x1 + x2, poisson, offset = log(e),
                       prior.scale = c(Inf, Inf), prior.df = Inf,
                       prior.scale.for.intercept = Inf,
                       prior.df.for.intercept = Inf),
              error = function(e) NULL)
stopifnot("offset= argument + vector prior: null deviance differs from glm" =
            !is.null(b) && close(b$null.deviance, g$null.deviance))

## an all-equal vector prior is the scalar prior
bs <- bayesglm(y ~ x1 + x2 + offset(log(e)), poisson, prior.scale = 2.5)
bv <- fit(y ~ x1 + x2 + offset(log(e)), poisson, prior.scale = c(2.5, 2.5),
          prior.mean = c(0, 0), prior.df = c(1, 1))
stopifnot("offset + vector default prior differs from the scalar prior" =
            !is.null(bv) && close(coef(bv), coef(bs)) &&
              close(bv$null.deviance, bs$null.deviance))
