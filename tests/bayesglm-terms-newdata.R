## =====================================================================
## Regression test: predict(type = "terms", newdata = )  --  issue #54 (minor)
##
## The bug: predictLM() builds the term-contribution matrix from
## model.matrixBayes(object, ...), and that call dropped the intercept
## column: model.matrixBayes() copies attr(object, "intercept") onto its
## terms, but a *fitted* object has no such attribute (its terms do), so it
## copied NULL and the matrix came back without an intercept.  The term
## labels then had one more entry than the matrix had columns and
##
##   predict(fit, type = "terms", newdata = d)
##
## died with "invalid 'labels'; length 3 should be 1 or 2".  Without
## newdata it worked, because that path reuses the model matrix from the
## fit.
## =====================================================================

library(arm)

set.seed(1)
n <- 100
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))
fit <- bayesglm(y ~ x1 + x2, binomial, d)

tr <- predict(fit, type = "terms", newdata = d)
stopifnot("predict(type = 'terms', newdata =) does not return a matrix" =
            is.matrix(tr))
stopifnot("one column per term is missing" =
            identical(colnames(tr), c("x1", "x2")))

## the contributions plus the constant must add up to the linear predictor,
## which is what "terms" is for
lp <- predict(fit, newdata = d)
stopifnot("the term contributions do not add up to the linear predictor" =
            isTRUE(all.equal(unname(rowSums(tr) + attr(tr, "constant")),
                             unname(lp))))

## and the same on the training data, where it already worked
tr0 <- predict(fit, type = "terms")
stopifnot("the no-newdata path broke" =
            isTRUE(all.equal(unname(rowSums(tr0) + attr(tr0, "constant")),
                             unname(fit$linear.predictors))))
