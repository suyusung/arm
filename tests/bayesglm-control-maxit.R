## =====================================================================
## Regression test: bayesglm() honours control = list(maxit = )
##
## The bug: bayesglm() always overwrote control$maxit with its own maxit
## argument (default 100), so control = list(maxit = k) was ignored, contrary
## to the Rd ("This can also be controlled by control").
##
## Expected: control$maxit is used when the maxit argument is not supplied;
## an explicitly supplied maxit argument still wins; the default is still 100.
## =====================================================================

library(arm)

set.seed(9)
n <- 100
x <- rnorm(n)
y <- rbinom(n, 1, plogis(0.5 + x))

quiet <- function(expr) suppressWarnings(expr)

b <- quiet(bayesglm(y ~ x, binomial, control = list(maxit = 1)))
stopifnot("control = list(maxit = 1) ignored: more than one IRLS iteration ran" =
            b$iter == 1 && b$control$maxit == 1)

b <- quiet(bayesglm(y ~ x, binomial, maxit = 1))
stopifnot("maxit = 1 argument ignored" = b$iter == 1)

b <- quiet(bayesglm(y ~ x, binomial, control = list(maxit = 1), maxit = 50))
stopifnot("an explicit maxit argument should win over control$maxit" =
            b$iter > 1 && b$control$maxit == 50)

b <- bayesglm(y ~ x, binomial)
stopifnot("default maxit should be 100" =
            b$control$maxit == 100 && b$converged)

## other control settings are passed through alongside maxit
b <- quiet(bayesglm(y ~ x, binomial, control = list(maxit = 2, epsilon = 1e-12)))
stopifnot("control epsilon not passed through with control maxit" =
            b$control$epsilon == 1e-12 && b$iter == 2)
