## =====================================================================
## Regression test for sim() on aliased coefficients -- issue #54, s2
##
## Bug (s2): when a model has an aliased (collinear) coefficient, sim(glm)
## kept every coefficient in coefficients(object) and filled the aliased
## one with 0, while sim(lm) let summary()'s coef matrix drop the aliased
## column entirely -- the two methods disagreed on the shape of @coef.
## man/sim.Rd documents @coef as an n.sims x k matrix with k the number of
## coefficients, i.e. the full coefficient vector, so sim(glm)'s shape is
## the documented one and sim(lm) is the side that changed: it now keeps
## the aliased column filled with 0.
## =====================================================================

library(arm)
library(MASS)

set.seed(11)
n <- 50
x <- rnorm(n)
x2 <- 2 * x          # exactly collinear with x: x2 is aliased
y <- rnorm(n)
m <- lm(y ~ x + x2)
g <- glm(y ~ x + x2)

set.seed(11); sm <- sim(m, 5)@coef
set.seed(11); sg <- sim(g, 5)@coef

## ---- s2: lm and glm must agree on the coefficient names ---------------
stopifnot("s2: sim(lm) and sim(glm) disagree on the coefficient names" =
            identical(colnames(sm), colnames(sg)))
stopifnot("s2: sim(lm) dropped the aliased coefficient instead of filling 0" =
            identical(colnames(sm), names(coefficients(m))))
stopifnot("s2: sim(lm) did not fill the aliased coefficient with 0" =
            all(sm[, "x2"] == 0))
stopifnot("s2: sim(lm) zeroed a coefficient that is estimable" =
            all(sm[, c("(Intercept)", "x")] != 0))

## ---- s2: the estimable draws must be unchanged ------------------------
## Rebuild the draws the way sim(lm) does and compare the estimable columns.
summ <- summary(m)
beta.hat <- summ$coef[, 1, drop = FALSE]
sigma.hat <- summ$sigma
N <- summ$df[1] + summ$df[2]
K <- summ$df[1]
set.seed(11)
sig <- numeric(5)
B <- matrix(NA_real_, 5, K)
for (s in 1:5) {
  sig[s] <- sigma.hat * sqrt((N - K) / rchisq(1, N - K))
  B[s, ] <- MASS::mvrnorm(1, beta.hat, summ$cov.unscaled * sig[s]^2)
}
stopifnot("s2: sim(lm)'s estimable draws changed" =
            isTRUE(all.equal(unname(sm[, c("(Intercept)", "x")]), unname(B))))
