## =====================================================================
## Regression test: scaled = TRUE prior scales for a Gaussian bayesglm()
##
## The Rd: with scaled = TRUE the prior scale for a predictor with more than
## two values is prior.scale/(2*sd(x)) (prior.scale/range(x) for a binary
## one), and for a Gaussian response prior.scale is also multiplied by
## 2*sd(y).
##
## The bug: for the Gaussian family bayesglm.fit() saved the prior scales
## after the 2*sd(y) step but before the division by the x scale, and
## restored that copy on every IRLS iteration.  Only the first iteration used
## the documented prior; after that the x scaling was undone, so the prior
## on a coefficient did not depend on the predictor's units.  Rescaling x
## then changed the fit: with a normal prior, coef(x) was 0.105 but
## coef(x/1000)/1000 was 1.4e-7.  Binomial and Poisson fits were not affected.
##
## Identities: coef(y ~ x) equals coef(y ~ I(x/1000))/1000, and coef(10*y) is
## 10*coef(y); with a normal prior the fit equals the closed-form augmented
## least-squares (ridge) solution at the fitted dispersion, with the
## documented prior sd prior.scale*2*sd(y)/(2*sd(x)).
## =====================================================================

library(arm)

close <- function(a, b, tol = 1e-6) isTRUE(all.equal(a, b, tolerance = tol))

set.seed(2)
n <- 30
x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(n)
xk <- x / 1000

## ---- invariance to the units of x ------------------------------------
c1 <- bayesglm(y ~ x, prior.df = Inf, prior.scale = 0.05)
c2 <- bayesglm(y ~ xk, prior.df = Inf, prior.scale = 0.05)
stopifnot("gaussian scaled prior: normal-prior slope not invariant to the units of x" =
            close(coef(c1)[["x"]], coef(c2)[["xk"]] / 1000))
stopifnot("gaussian scaled prior: intercept not invariant to the units of x" =
            close(coef(c1)[[1]], coef(c2)[[1]]))

b1 <- bayesglm(y ~ x)
b2 <- bayesglm(y ~ xk)
stopifnot("gaussian scaled prior: default (Cauchy) slope not invariant to the units of x" =
            close(coef(b1)[["x"]], coef(b2)[["xk"]] / 1000))

## ---- the documented prior scale is the one used ----------------------
stopifnot("gaussian scaled prior: returned prior.scale is not 2.5*2*sd(y)/(2*sd(x))" =
            close(b1$prior.scale[2], 2.5 * 2 * sd(y) / (2 * sd(x))))

## closed form: at the fitted dispersion s2, the normal-prior posterior mode
## is weighted least squares on the data plus one pseudo-observation per
## coefficient (the intercept's is taken at mean(x), as bayesglm does).
c3 <- bayesglm(y ~ x, prior.df = Inf, prior.scale = 0.05,
               prior.df.for.intercept = Inf)
s2 <- c3$dispersion
sd.int <- 10 * 2 * sd(y)
sd.x <- 0.05 * 2 * sd(y) / (2 * sd(x))
Xs <- rbind(cbind(1, x), c(1, mean(x)), c(0, 1))
ws <- c(rep(1, n), sqrt(s2) / sd.int, sqrt(s2) / sd.x)
ridge <- unname(coef(lm.fit(Xs * ws, c(y, 0, 0) * ws)))
stopifnot("gaussian scaled prior: normal-prior fit differs from the augmented least-squares solution" =
            close(unname(coef(c3)), ridge))

## ---- the 2*sd(y) factor keeps the fit equivariant in y ---------------
c4 <- bayesglm(I(10 * y) ~ x, prior.df = Inf, prior.scale = 0.05)
stopifnot("gaussian scaled prior: fit not equivariant to the units of y" =
            close(unname(coef(c4)), 10 * unname(coef(c1))))

## ---- a binary predictor is scaled by its range -----------------------
xb <- rep(c(0, 1), length.out = n)
yb <- 1 + xb + rnorm(n)
e1 <- bayesglm(yb ~ xb, prior.df = Inf, prior.scale = 0.05)
e2 <- bayesglm(yb ~ I(5 * xb), prior.df = Inf, prior.scale = 0.05)
stopifnot("gaussian scaled prior: binary slope not invariant to the range of x" =
            close(coef(e1)[[2]], 5 * coef(e2)[[2]]))
