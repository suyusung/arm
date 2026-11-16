## =====================================================================
## Regression test: zero prior weights in a Gaussian bayesglm()
##
## The bug: the Gaussian dispersion update in bayesglm.fit() used the
## working response z and weights w restricted to the rows with positive
## weight ("good" rows), but multiplied the full design matrix x by the
## coefficients.  With any zero weight the vectors had different lengths,
## R recycled them (with a "longer object length" warning), and the
## dispersion -- and through it the prior's pull and the coefficients --
## came out wrong.
##
## Identity: an observation with prior weight 0 contributes nothing, so the
## fit must equal the fit with that row dropped.  scaled = FALSE isolates the
## dispersion update from the data-dependent prior scaling (which uses
## sd(x) and sd(y) over all rows, zero-weight rows included).
## =====================================================================

library(arm)

close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-5))

set.seed(4)
n <- 30
x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(n)
w0 <- c(0, rep(1, n - 1))

warned <- FALSE
a <- withCallingHandlers(bayesglm(y ~ x, weights = w0, scaled = FALSE),
                         warning = function(w) {
                           warned <<- TRUE
                           invokeRestart("muffleWarning")
                         })
b <- bayesglm(y[-1] ~ x[-1], scaled = FALSE)

stopifnot("zero weight: Gaussian fit warned (recycling in the dispersion update)" =
            !warned)
stopifnot("zero weight: dispersion differs from dropping the row" =
            close(a$dispersion, b$dispersion))
stopifnot("zero weight: coefficients differ from dropping the row" =
            close(unname(coef(a)), unname(coef(b))))

## several zero weights among unequal positive weights
set.seed(5)
wv <- runif(n, 0.5, 2)
wv[c(3, 10, 17)] <- 0
keep <- wv > 0
a <- bayesglm(y ~ x, weights = wv, scaled = FALSE)
b <- bayesglm(y[keep] ~ x[keep], weights = wv[keep], scaled = FALSE)
stopifnot("several zero weights: fit differs from dropping the rows" =
            close(a$dispersion, b$dispersion) &&
              close(unname(coef(a)), unname(coef(b))))
