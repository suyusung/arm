## =====================================================================
## Regression test: bayesglm() terms follow keep.order
##
## The bug: bayesglm() builds its design matrix with keep.order = TRUE (the
## default), so for y ~ x:z + x the coefficients are ordered
## (Intercept), x:z, x.  But the terms object stored in the fit came from
## model.frame() with the default keep.order = FALSE, ordered
## (Intercept), x, x:z.  Everything that rebuilds the design matrix from the
## stored terms -- model.matrix(fit), predict(fit, se.fit = TRUE) and
## predict(fit, type = "terms") without newdata -- then multiplied columns
## by the wrong coefficients.
##
## Identities: predict(fit, se.fit = TRUE)$fit equals the stored linear
## predictors; model.matrix(fit) has the coefficients' column order; the
## "terms" contributions plus their constant add up to the linear predictor.
## =====================================================================

library(arm)

close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8))

set.seed(8)
n <- 50
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.5 * d$x + d$x * d$z))

b <- bayesglm(y ~ x:z + x, binomial, d)
stopifnot("keep.order: model.matrix(fit) columns are not in coefficient order" =
            identical(colnames(model.matrix(b)), names(coef(b))))
stopifnot("keep.order: predict(se.fit = TRUE)$fit differs from the linear predictors" =
            close(unname(predict(b, se.fit = TRUE)$fit),
                  unname(b$linear.predictors)))
stopifnot("keep.order: se.fit without newdata differs from se.fit with newdata = data" =
            close(unname(predict(b, se.fit = TRUE)$se.fit),
                  unname(predict(b, newdata = d, se.fit = TRUE)$se.fit)))
tr <- predict(b, type = "terms")
stopifnot("keep.order: type = 'terms' does not add up to the linear predictor" =
            close(unname(rowSums(tr) + attr(tr, "constant")),
                  unname(b$linear.predictors)))

## the ordering really is kept, and keep.order = FALSE is unchanged
stopifnot("keep.order = TRUE did not keep the formula's term order" =
            identical(attr(terms(b), "term.labels"), c("x:z", "x")))
b2 <- bayesglm(y ~ x:z + x, binomial, d, keep.order = FALSE)
stopifnot("keep.order = FALSE: predict(se.fit = TRUE)$fit differs from the linear predictors" =
            close(unname(predict(b2, se.fit = TRUE)$fit),
                  unname(b2$linear.predictors)))
stopifnot("keep.order = FALSE should reorder terms by degree" =
            identical(attr(terms(b2), "term.labels"), c("x", "x:z")))

## a formula given as a character string still works, as it does for glm()
b3 <- bayesglm("y ~ x:z + x", binomial, d)
stopifnot("a character formula should give the same fit as the formula" =
            identical(coef(b3), coef(b)))
