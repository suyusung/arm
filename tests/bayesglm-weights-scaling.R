## =====================================================================
## Regression test: weights and the scaled prior  --  issue #54 (minor)
##
## The bug: with scaled = TRUE the prior scales were built from spreads
## computed without the weights -- sd(y) for a gaussian response, sd(x) for a
## predictor with more than two values -- and the augmentation used an
## unweighted colMeans().  The dispersion update had the same problem, its
## mean() ran over rows rather than over the weighted sample.  So a weight of
## 2 meant something different from a duplicated row, which is what a weight
## is supposed to mean.
##
## The fix computes all of those spreads with the prior weights.  The
## identity below is the definition: a fit with weights w must equal the fit
## on the data expanded by w.
## =====================================================================

library(arm)

set.seed(11)
n <- 30
x <- rnorm(n)
x[1] <- 6                                   # an outlier, so the weights bite
y <- 1 + 0.5 * x + rnorm(n)
d <- data.frame(x = x, y = y, w = c(200, rep(1, n - 1)))
d2 <- rbind(d[, c("x", "y")], d[rep(1, 199), c("x", "y")])   # row 1 repeated 200x

f1 <- bayesglm(y ~ x, gaussian, d, weights = w)
f2 <- bayesglm(y ~ x, gaussian, d2)

stopifnot("the prior scales still ignore the weights" =
            isTRUE(all.equal(unname(f1$prior.scale), unname(f2$prior.scale),
                             tolerance = 1e-8)))

stopifnot("a weight of 200 is not the same as 200 copies of the row" =
            isTRUE(all.equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-4)))

## unweighted fits are untouched: weights of 1 must change nothing
f3 <- bayesglm(y ~ x, gaussian, d)
f4 <- bayesglm(y ~ x, gaussian, d, weights = rep(1, n))
stopifnot("an unweighted fit changed" =
            isTRUE(all.equal(unname(coef(f3)), unname(coef(f4)), tolerance = 1e-12)))

## and the prior scales still scale with the data, as documented
d$y2 <- d$y * 100
f5 <- bayesglm(y2 ~ x, gaussian, d)
stopifnot("the gaussian prior scale no longer follows the spread of y" =
            isTRUE(all.equal(f5$prior.scale[1], 100 * f3$prior.scale[1], tolerance = 1e-6)))
