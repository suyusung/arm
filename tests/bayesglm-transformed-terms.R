## =====================================================================
## Regression test: a formula whose term is a transformation  --  found by the
## reverse-dependency check of 1.15-4 (qgcomp's bayes path)
##
## The bug: bayesglm() handed model.matrixBayes() the *model frame*, and that
## function subset it -- keeping only the columns the terms name.  A model frame
## holds the evaluated column of a transformed term, `poly(psi, degree = 3,
## raw = TRUE)`, not the variable `psi` the term refers to, and subsetting drops
## the frame's "terms" attribute, so model.matrix.default tried to rebuild the
## frame from what was left and died with
##
##   Error in eval(predvars, data, env) : object 'psi' not found
##
## qgcomp calls bayesglm() with formulas like this one and a `weights` column
## passed as a bare symbol; that combination is the one that failed.
## =====================================================================

library(arm)

set.seed(42)
d <- data.frame(Ya = rnorm(60), psi = rep(0:5, 10), weights = 1)
d$Ya <- 1 + 0.3 * d$psi + rnorm(60)

## ---- the exact shape qgcomp uses --------------------------------------
f <- as.formula("Ya ~ 1 + poly(psi, degree=3, raw=TRUE)")
b <- bayesglm(f, data = d, weights = weights, x = TRUE)

stopifnot("bayesglm() cannot fit a formula whose term is a transformation" =
            length(coef(b)) == 4L && all(is.finite(coef(b))))
stopifnot("the stored design matrix disagrees with the coefficients" =
            ncol(model.matrix(b)) == length(coef(b)))

## the fit still predicts and still adds up, layer by layer
stopifnot("predict() on the transformed term misbehaves" =
            length(predict(b, newdata = d)) == 60L)
tr <- predict(b, type = "terms")
stopifnot("the term contributions do not add up to the linear predictor" =
            isTRUE(all.equal(unname(rowSums(tr) + attr(tr, "constant")),
                             unname(b$linear.predictors))))

## ---- other transformations, and no weights -----------------------------
b2 <- bayesglm(Ya ~ 1 + log(psi + 1), data = d, weights = weights)
b3 <- bayesglm(Ya ~ 1 + poly(psi, degree = 2, raw = TRUE), data = d)
b4 <- bayesglm(Ya ~ 1 + I(psi^2) + psi, data = d)
stopifnot("log()/I() terms fail too" =
            all(is.finite(coef(b2))) && all(is.finite(coef(b3))) && all(is.finite(coef(b4))))

## ---- a binomial fit with a transformed term ----------------------------
set.seed(7)
d$y <- rbinom(60, 1, plogis(d$psi - 2))
b5 <- bayesglm(y ~ 1 + poly(psi, degree = 2, raw = TRUE), binomial, d)
stopifnot("a binomial fit with a transformed term fails" = all(is.finite(coef(b5))))
stopifnot("the binomial design matrix disagrees with the coefficients" =
            ncol(model.matrix(b5)) == length(coef(b5)))
