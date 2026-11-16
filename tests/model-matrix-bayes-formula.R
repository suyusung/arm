## =====================================================================
## Regression test: what model.matrixBayes(formula, data) returns -- found by
## the reverse-dependency check of 1.15-4 (FindIt's CausalANOVA examples)
##
## model.matrixBayes() copies attr(object, "intercept") onto the terms it builds,
## and a *formula* has no such attribute, so the copy sets it to NULL and the
## matrix comes back WITHOUT an intercept column.  That is an accident of the
## code, but FindIt depends on it: cv.CausalANOVA() calls
## model.matrixBayes(formula, data) and prepends its own intercept,
##
##   cbind(1, X.l[[z]]) %*% coef.cv[[z]]
##
## so a version of arm that adds the intercept to this call makes FindIt's
## examples fail with "non-conformable arguments".
##
## The terms path is the one that must keep the intercept: predict() builds its
## term contributions from object$terms plus the model frame.
## =====================================================================

library(arm)

set.seed(1)
d <- data.frame(y = rnorm(20), x = rnorm(20))

## ---- the formula path: no intercept column -----------------------------
X <- arm:::model.matrixBayes(y ~ x, data = d)
stopifnot("model.matrixBayes(formula, data) grew an intercept column" =
            ncol(X) == 1L && identical(colnames(X), "x"))

## the way FindIt builds the call: a formula assembled inside its own function
build <- function(data) {
  f <- as.formula(paste("y ~", "x"))
  arm:::model.matrixBayes(f, data = data)
}
Xc <- build(d)
stopifnot("the formula-assembled-at-runtime path differs" =
            ncol(Xc) == 1L && identical(colnames(Xc), "x"))

## ---- the terms path: intercept kept ------------------------------------
fit <- bayesglm(y ~ x, data = d)
Xt <- arm:::model.matrixBayes(fit$terms, data = model.frame(fit))
stopifnot("the terms path lost the intercept" =
            ncol(Xt) == 2L && identical(colnames(Xt), c("(Intercept)", "x")))

## ---- no intercept in the formula means no intercept anywhere -----------
Xi <- arm:::model.matrixBayes(y ~ x - 1, data = d)
stopifnot("an intercept appeared for a formula without one" = ncol(Xi) == 1L)

## ---- and the documented path still adds up ------------------------------
tr <- predict(fit, type = "terms", newdata = d)
stopifnot("predict(type = 'terms', newdata =) no longer adds up to the fit" =
            isTRUE(all.equal(unname(rowSums(tr) + attr(tr, "constant")),
                             unname(predict(fit, newdata = d)))))
