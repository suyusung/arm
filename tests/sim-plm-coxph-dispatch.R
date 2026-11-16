## =====================================================================
## Regression test for sim() dispatch and covariance scaling -- issue #54, s1
##
## Bug (s1): sim.plm() and sim.coxph() were ordinary, unexported functions,
## not methods of the S4 generic sim(), so sim(<plm>) and sim(<coxph>)
## both errored with "unable to find an inherited method for function
## 'sim'".  Reaching sim.plm() directly (arm:::sim.plm) exposed a second
## defect: vcov(<plm>) is already scaled by the estimated error variance
## (summary()'s standard errors are exactly sqrt(diag(vcov(object)))), yet
## the draw loop multiplied it by sigma[s]^2 as well, counting the
## dispersion twice.  On master the draws' variances are ~2600x vcov(plm)
## (middle diag ratio 2623 for this fixture, with sigma.hat^2 = 2935).
##
## The fixture is the plm package's Grunfeld panel -- small and fixed.  The
## test also verifies the vcov(plm) oracle against summary() before using
## it, so a mistake in the test's own maths fails loudly.
##
## sim.coxph() never had the scaling bug (vcov(coxph) was used as-is); it
## only needed to be dispatched, so its draw covariance is checked too.
## =====================================================================

library(arm)

## plm and survival are neither Depends nor Suggests here, and R CMD check
## runs every tests/*.R, so skip rather than fail where they are missing.
if (!requireNamespace("plm", quietly = TRUE) ||
    !requireNamespace("survival", quietly = TRUE)) {
    cat("plm and/or survival not installed: skipping the sim dispatch test\n")
    quit(save = "no", status = 0)
}
library(plm)
library(survival)

data("Grunfeld", package = "plm")
p <- plm(inv ~ value + capital, data = Grunfeld, model = "within")

## ---- s1: the methods must be reachable through the sim() generic ------
sp <- tryCatch(sim(p, 1000), error = function(e) e)
stopifnot("s1: sim() does not dispatch on a plm object (method unreachable)" =
            !inherits(sp, "error"))

cx <- coxph(Surv(futime, fustat) ~ age + ecog.ps, data = ovarian)
sc <- tryCatch(sim(cx, 1000), error = function(e) e)
stopifnot("s1: sim() does not dispatch on a coxph object (method unreachable)" =
            !inherits(sc, "error"))

## ---- s1: vcov(plm) is the scaled covariance (validate the oracle) -----
stopifnot("s1: the vcov(plm) oracle disagrees with summary(), so the test's maths is wrong" =
            isTRUE(all.equal(unname(sqrt(diag(vcov(p)))),
                             unname(summary(p)$coefficients[, 2]),
                             check.attributes = FALSE, tolerance = 1e-8)))

## ---- s1: sim(plm) draws must have covariance vcov(plm) ----------------
set.seed(21)
sp <- sim(p, 2000)
stopifnot("s1: sim(plm)@coef has the wrong shape" =
            identical(dim(sp@coef), c(2000L, 2L)) &&
            identical(colnames(sp@coef), rownames(summary(p)$coefficients)))

ratio <- diag(var(sp@coef)) / diag(vcov(p))
stopifnot("s1: sim(plm) draws have wildly wrong variances (dispersion double-counted)" =
            all(ratio > 0.7 & ratio < 1.4))
stopifnot("s1: sim(plm) draws are systematically off vcov(plm)" =
            stats::median(ratio) > 0.9 && stats::median(ratio) < 1.1)

## ---- sim(coxph) draws must have covariance vcov(coxph) ----------------
set.seed(22)
sc <- sim(cx, 2000)
ratio.c <- diag(var(sc$coef)) / diag(vcov(cx))
stopifnot("s1: sim(coxph) draws do not match vcov(coxph)" =
            all(ratio.c > 0.7 & ratio.c < 1.4))
