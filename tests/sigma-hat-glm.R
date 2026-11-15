## =====================================================================
## Regression tests for sigma.hat.glm()
##
## The bug: for any family other than gaussian, sigma.hat.glm() returned
## summary(object)$sigma.  summary.glm() has no "sigma" component, so the
## result was NULL for binomial, poisson, quasi and Gamma fits, and the
## dispersion computed a few lines earlier was thrown away.
##
## sigma.hat() should be sqrt(dispersion) for every family: 1 for
## binomial and poisson, sqrt(estimated dispersion) otherwise.  That is
## also what sim() stores in @sigma for a glm.
## =====================================================================

library(arm)

set.seed(1)
x <- rnorm(50)
yb <- rbinom(50, 1, invlogit(x))
yp <- rpois(50, exp(x))
yg <- rgamma(50, shape = 2, rate = 2 / exp(0.3 * x))

fits <- list(gaussian     = glm(yg ~ x),
             binomial     = glm(yb ~ x, family = binomial),
             poisson      = glm(yp ~ x, family = poisson),
             quasipoisson = glm(yp ~ x, family = quasipoisson),
             Gamma        = glm(yg ~ x, family = Gamma))

for (fit in fits) {
  stopifnot("sigma.hat(glm) is not sqrt(summary(fit)$dispersion)" =
              isTRUE(all.equal(sigma.hat(fit), sqrt(summary(fit)$dispersion))))
  stopifnot("sigma.hat(glm) disagrees with the sigma that sim() draws" =
              isTRUE(all.equal(sigma.hat(fit), sim(fit, 2)@sigma[1])))
}

## binomial and poisson have dispersion fixed at 1
stopifnot("sigma.hat() of a binomial glm is not 1" =
            identical(sigma.hat(fits$binomial), 1))
stopifnot("sigma.hat() of a poisson glm is not 1" =
            identical(sigma.hat(fits$poisson), 1))
