## =====================================================================
## Regression tests for display(merMod)  --  intercept-only fits  --  d5
##
## The bug: display() copies the fixed effects with coefs[, "coef.est"].
## For a model with a single fixed effect (an intercept-only lmer) that is
## a one-row matrix, and subsetting it drops *both* dimensions, so the
## returned out$coef (and out$se) lose the name "(Intercept)".
##
## The fix restores the coefficient names from rownames(coefs).
## =====================================================================

library(arm)
suppressMessages(library(lme4))

set.seed(2)
g <- factor(rep(1:10, each = 10))
y <- rnorm(100) + as.numeric(g) * 0.1
m0 <- suppressWarnings(lmer(y ~ 1 + (1 | g)))

out <- tryCatch({capture.output(suppressWarnings(res <- display(m0))); res},
                error = function(e) e)
stopifnot("display() failed on an intercept-only lmer" =
            !inherits(out, "error"))
stopifnot("display(intercept-only lmer) dropped the (Intercept) coefficient name" =
            identical(names(out$coef), "(Intercept)"))
stopifnot("display(intercept-only lmer) dropped the (Intercept) se name" =
            identical(names(out$se), "(Intercept)"))

## a model with more than one fixed effect is unchanged
m1 <- lmer(y ~ x1 + (1 | g), data = data.frame(y = y, x1 = rnorm(100), g = g))
out1 <- capture.output(suppressWarnings(res1 <- display(m1)))
stopifnot("display(lmer) lost the fixed-effect names" =
            identical(names(res1$coef), c("(Intercept)", "x1")))
