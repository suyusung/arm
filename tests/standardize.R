## =====================================================================
## Regression tests for standardize()  --  issue #24
##
## The bug: standardize.default() rebuilds the model from object$call, and
## it used call$formula directly.  When the formula had been stored in an
## object -- lm(form, data = d) -- call$formula was the *symbol* `form`
## rather than the formula itself, and standardize() died with
##
##   Error: The object NA must be either numeric or a factor.
##
## It now evaluates the formula (eval(call$formula, parent.frame())), so
## the two ways of writing the same model must agree.
##
## Note that standardize() re-fits from the call plus the data frame it
## finds *by name*, not from the model's stored model frame -- so `d` has
## to stay reachable under its original name.  That is also why the
## comparison below is between two models fitted on the same data frame.
## =====================================================================

library(arm)

set.seed(3)
d <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                g = factor(rep(c("a", "b"), length.out = 60)))
d$y <- d$y + 0.5 * d$x + 0.4 * (d$g == "b")

## the same model written two ways
m_inline <- stats::lm(y ~ x + g, data = d)
form     <- y ~ x + g
m_object <- stats::lm(form, data = d)

s_inline <- standardize(m_inline)
s_object <- standardize(m_object)

stopifnot("standardize() failed on a formula held in an object (#24 is back)" =
            inherits(s_object, "lm"))
stopifnot("inline and stored formulas standardize differently" =
            isTRUE(all.equal(stats::coef(s_inline), stats::coef(s_object))))

## the generated variable names are part of the interface
stopifnot("unexpected coefficient names from standardize()" =
            identical(names(stats::coef(s_inline)),
                      c("(Intercept)", "z.x", "c.g")))

## an unnamed argument list, i.e. plain lm(..., data = d), still resolves
## the data by the name in the call
d2 <- d
m2 <- stats::lm(y ~ x, data = d2)
s2 <- standardize(m2)
stopifnot("standardize() lost the data frame name" =
            identical(names(stats::coef(s2)), c("(Intercept)", "z.x")))

## a standardized response, by request: the response then has mean 0, so the
## intercept is ~0, and the slopes shrink by 2 sd(y)
s3 <- standardize(m_inline, standardize.y = TRUE)
stopifnot("standardize.y = TRUE did not standardize the response" =
            abs(stats::coef(s3)[["(Intercept)"]]) < 1e-8)
stopifnot("standardize.y = TRUE did not rescale the slopes by 2 sd(y)" =
            isTRUE(all.equal(unname(stats::coef(s3)[-1]),
                             unname(stats::coef(s_inline)[-1]) / (2 * stats::sd(d$y)),
                             tolerance = 1e-6)))
