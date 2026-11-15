## =====================================================================
## Regression test: bayesglm() AIC for a binomial response with trials > 1
##
## The bug: bayesglm.fit() called family$aic(y, n.ok, mu, weights, dev),
## passing the scalar count of non-zero-weight rows where glm.fit() passes
## the per-row number of trials `n` set by family$initialize.  For a
## cbind(successes, failures) response the binomial log-likelihood was then
## computed with the wrong number of trials, so AIC() and logLik() were wrong
## even though the coefficients and deviance were right.
##
## Identity: with flat priors (prior.scale = Inf, prior.df = Inf) bayesglm is
## maximum likelihood, so AIC and logLik must equal glm()'s.
## =====================================================================

library(arm)

flat <- function(...) bayesglm(..., prior.scale = Inf, prior.df = Inf,
                               prior.scale.for.intercept = Inf,
                               prior.df.for.intercept = Inf)

set.seed(1)
n <- 40
x <- rnorm(n)
trials <- sample(3:12, n, TRUE)
s <- rbinom(n, trials, plogis(0.3 + 0.8 * x))
d <- data.frame(s, f = trials - s, x)

g <- glm(cbind(s, f) ~ x, binomial, d)
b <- flat(cbind(s, f) ~ x, binomial, d)
stopifnot("flat-prior binomial coefficients differ from glm" =
            isTRUE(all.equal(coef(b), coef(g), tolerance = 1e-6)))
stopifnot("binomial (trials > 1) AIC differs from glm: wrong n passed to family$aic" =
            isTRUE(all.equal(AIC(b), AIC(g), tolerance = 1e-6)))
stopifnot("binomial (trials > 1) logLik differs from glm" =
            isTRUE(all.equal(as.numeric(logLik(b)), as.numeric(logLik(g)),
                             tolerance = 1e-6)))

## proportion response with weights = trials is the same model
d$p <- d$s / trials
gw <- glm(p ~ x, binomial, d, weights = trials)
bw <- bayesglm(p ~ x, binomial, d, weights = trials,
               prior.scale = Inf, prior.df = Inf,
               prior.scale.for.intercept = Inf, prior.df.for.intercept = Inf)
stopifnot("binomial proportion-with-weights AIC differs from glm" =
            isTRUE(all.equal(AIC(bw), AIC(gw), tolerance = 1e-6)))

## Bernoulli and Poisson were already right; keep them right
y01 <- rbinom(n, 1, plogis(x))
stopifnot("Bernoulli AIC differs from glm" =
            isTRUE(all.equal(AIC(flat(y01 ~ x, binomial)),
                             AIC(glm(y01 ~ x, binomial)), tolerance = 1e-6)))
yp <- rpois(n, exp(0.5 + 0.3 * x))
stopifnot("Poisson AIC differs from glm" =
            isTRUE(all.equal(AIC(flat(yp ~ x, poisson)),
                             AIC(glm(yp ~ x, poisson)), tolerance = 1e-6)))
