## =====================================================================
## Regression tests for sim() index handling at edge sizes
##
## Bug 1: sim() for a glm with n.sims = 1 returned all-zero coefficients.
## MASS::mvrnorm(1, ...) returns a bare named vector, not a 1 x k matrix,
## so dimnames(beta)[[2]] was NULL, the %in% name mask was all FALSE, and
## nothing was copied into the zero-filled result.
##
## Bug 2: sim() for an intercept-only polr (k = 0 regression
## coefficients) used 1:k, which is c(1, 0), so the first cutpoint was
## put into @coef.  @coef must have zero columns and @zeta one column per
## cutpoint.
##
## Both are checked against MASS::mvrnorm() called directly with the same
## seed and the same mean and covariance that sim() uses.
## =====================================================================

library(arm)
library(MASS)

set.seed(1)
x <- rnorm(60)
y <- rbinom(60, 1, invlogit(0.5 + x))
g <- glm(y ~ x, family = binomial)

## the covariance matrix exactly as sim() builds it
summ <- summary(g, correlation = TRUE)
sd.beta <- summ$coef[, 2]
V <- summ$corr * outer(sd.beta, sd.beta)

## ---- glm, n.sims = 1 ---------------------------------------------------
set.seed(42); s1 <- sim(g, 1)@coef
set.seed(42); d1 <- mvrnorm(1, summ$coef[, 1], V)
stopifnot("sim(glm, n.sims = 1) returned all-zero coefficients" =
            any(s1 != 0))
stopifnot("sim(glm, n.sims = 1) is not a 1 x k matrix" =
            identical(dim(s1), c(1L, 2L)))
stopifnot("sim(glm, n.sims = 1) differs from the mvrnorm draw" =
            isTRUE(all.equal(s1[1, ], d1)))

## ---- glm, n.sims > 1 is unchanged -------------------------------------
set.seed(42); s5 <- sim(g, 5)@coef
set.seed(42); d5 <- mvrnorm(5, summ$coef[, 1], V)
stopifnot("sim(glm, n.sims = 5) differs from the mvrnorm draws" =
            isTRUE(all.equal(unname(s5), unname(d5))))

## ---- glm with an aliased coefficient, n.sims = 1 ----------------------
x2 <- 2 * x
ga <- glm(y ~ x + x2, family = binomial)
set.seed(42); sa <- sim(ga, 1)@coef
stopifnot("sim(glm with aliased term, n.sims = 1) lost the estimable draws" =
            all(sa[1, c("(Intercept)", "x")] != 0) && sa[1, "x2"] == 0)

## ---- intercept-only polr ----------------------------------------------
yo <- cut(x, c(-Inf, -0.5, 0.5, Inf))
p0 <- polr(yo ~ 1, Hess = TRUE)
set.seed(7); sp <- sim(p0, 4)
set.seed(7); dp <- mvrnorm(4, p0$zeta, vcov(p0))
stopifnot("sim(polr y ~ 1) put a cutpoint into @coef" =
            ncol(sp@coef) == 0)
stopifnot("sim(polr y ~ 1) @zeta does not have one column per cutpoint" =
            ncol(sp@zeta) == length(p0$zeta))
stopifnot("sim(polr y ~ 1) @zeta differs from the mvrnorm draws" =
            isTRUE(all.equal(unname(sp@zeta), unname(dp))))

## ---- polr with a covariate is unchanged -------------------------------
p1 <- polr(yo ~ y, Hess = TRUE)
set.seed(7); sq <- sim(p1, 4)
set.seed(7); dq <- mvrnorm(4, c(coef(p1), p1$zeta), vcov(p1))
stopifnot("sim(polr) @coef/@zeta split is wrong" =
            isTRUE(all.equal(unname(cbind(sq@coef, sq@zeta)), unname(dq))) &&
            ncol(sq@coef) == 1 && ncol(sq@zeta) == 2)
