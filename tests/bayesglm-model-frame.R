## =====================================================================
## Regression test: bayesglm() builds its design matrix from the model frame
##
## The bug: bayesglm() built the model frame with subset, weights, offset and
## drop.unused.levels, but then built the design matrix from the original
## `data`, and it also deleted na.action from the model.frame() call.  So:
##   * subset = ... gave an "incompatible dimensions" error, because the
##     response had the subset's rows and the design matrix had all rows;
##   * na.action = na.exclude was silently replaced by na.omit (residuals
##     not padded) and na.action = na.fail did not fail;
##   * drop.unused.levels had no effect: a factor with an unused level kept
##     a column for it, so the baseline was the empty level (a flat-prior fit
##     then failed outright, and a default-prior fit gave different estimates
##     from the same data with the level dropped by hand).
##
## Identities: with flat priors bayesglm is maximum likelihood, so it must
## equal glm() on the same subset / na.action / weights / offset; and a factor
## with an unused level must give the same fit as droplevels() of it.
## =====================================================================

library(arm)

close <- function(a, b) !is.null(a) && isTRUE(all.equal(a, b, tolerance = 1e-6))
flat <- list(prior.scale = Inf, prior.df = Inf,
             prior.scale.for.intercept = Inf, prior.df.for.intercept = Inf)
## an error (the unfixed behaviour for subset=) counts as a failed identity
bflat <- function(...) tryCatch(do.call(bayesglm, c(list(...), flat)),
                                 error = function(e) NULL)

set.seed(3)
n <- 40
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)),
                w = runif(n, 0.5, 2), t = rpois(n, 5) + 1)
d$y <- 1 + d$x + rnorm(n)
d$k <- rpois(n, d$t * exp(0.2 + 0.3 * d$x))

## ---- subset ----------------------------------------------------------
b <- bflat(y ~ x, data = d, subset = quote(x > 0))
g <- glm(y ~ x, data = d, subset = x > 0)
stopifnot("subset: flat-prior bayesglm differs from glm on the subset" =
            close(coef(b), coef(g)))
stopifnot("subset: wrong number of fitted values" =
            !is.null(b) && length(fitted(b)) == sum(d$x > 0))

## a subset that empties a factor level drops it (drop.unused.levels = TRUE)
b <- bflat(y ~ g, data = d, subset = quote(g != "c"))
g <- glm(y ~ g, data = d, subset = g != "c")
stopifnot("subset dropping a level: bayesglm differs from glm" =
            close(coef(b), coef(g)))

## ---- na.action -------------------------------------------------------
d2 <- d
d2$x[c(2, 5)] <- NA
b <- bayesglm(y ~ x, data = d2, na.action = na.exclude)
g <- glm(y ~ x, data = d2, na.action = na.exclude)
stopifnot("na.exclude not honoured: residuals not padded to nrow(data)" =
            length(residuals(b)) == n && all(is.na(residuals(b)[c(2, 5)])))
stopifnot("na.exclude not honoured: fitted values not padded" =
            length(fitted(b)) == n)
stopifnot("na.exclude not recorded in the fit" =
            inherits(b$na.action, "exclude"))
stopifnot("na.fail did not fail on missing data" =
            inherits(try(bayesglm(y ~ x, data = d2, na.action = na.fail),
                         silent = TRUE), "try-error"))
## the default (na.omit) is unchanged
b <- bflat(y ~ x, data = d2)
g <- glm(y ~ x, data = d2)
stopifnot("na.omit: flat-prior bayesglm differs from glm" =
            close(coef(b), coef(g)))

## ---- unused factor levels --------------------------------------------
set.seed(7)
m <- 60
d3 <- data.frame(x = rnorm(m),
                 g = factor(sample(c("b", "c", "d"), m, TRUE),
                            levels = c("a", "b", "c", "d")))
d3$y <- 1 + d3$x + (d3$g == "c") - (d3$g == "d") + rnorm(m)
d3$g2 <- droplevels(d3$g)
stopifnot("unused level: flat-prior bayesglm differs from glm" =
            close(unname(coef(bflat(y ~ x + g, data = d3))),
                  unname(coef(glm(y ~ x + g, data = d3)))))
b <- bayesglm(y ~ x + g, data = d3)
b0 <- bayesglm(y ~ x + g2, data = d3)
stopifnot("unused level: drop.unused.levels had no effect (baseline is the empty level)" =
            close(unname(coef(b)), unname(coef(b0))))
nd <- data.frame(x = 0, g = factor(c("b", "c", "d"), levels = levels(d3$g)))
stopifnot("unused level: predict() on newdata fails or disagrees" =
            close(unname(predict(b, nd)),
                  unname(predict(b0, data.frame(x = 0, g2 = nd$g)))))

## ---- weights and offsets still line up with the rows -----------------
b <- bflat(y ~ x, data = d2, weights = quote(w), subset = quote(g != "a"))
g <- glm(y ~ x, data = d2, weights = w, subset = g != "a")
stopifnot("weights + subset + NA: bayesglm differs from glm" =
            close(coef(b), coef(g)))
d2$t[7] <- NA
b <- bflat(k ~ x + offset(log(t)), family = poisson, data = d2,
           subset = quote(g != "b"))
g <- glm(k ~ x + offset(log(t)), family = poisson, data = d2, subset = g != "b")
stopifnot("formula offset + subset + NA: bayesglm differs from glm" =
            !is.null(b) && close(coef(b), coef(g)) && close(deviance(b), deviance(g)))
b <- bflat(k ~ x, family = poisson, data = d2, offset = quote(log(t)),
           subset = quote(g != "b"))
g <- glm(k ~ x, family = poisson, data = d2, offset = log(t), subset = g != "b")
stopifnot("offset argument + subset + NA: bayesglm differs from glm" =
            !is.null(b) && close(coef(b), coef(g)) && close(deviance(b), deviance(g)))
