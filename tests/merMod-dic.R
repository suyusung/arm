## =====================================================================
## Regression tests for DIC of merMod fits: extractDIC() and display()
##
## The bug: both computed pD = deviance(refitML(fit)) - (-2 logLik(fit)),
## the gap between two marginal likelihoods, and DIC = deviance + pD.
## That pD is exactly 0 for an ML lmer fit, an arbitrary REML-vs-ML gap
## for a REML fit, and about -110 for the binomial cbpp glmer, where
## deviance() and logLik() are on different scales (DIC was -37).
##
## The fix computes the conditional DIC that BUGS reports, from sim()
## draws: D = -2 log p(y | beta, b, sigma), pD = mean(D) - D(means),
## DIC = D(means) + 2 pD.
##
## The checks below do not reuse the implementation.  D at the fitted
## values is computed by hand from fitted() and sigma() (lmer) or from
## lme4's own conditional deviance fit@resp$aic() (glmer).  For lmer,
## pD is also checked against the trace of the hat matrix, which is the
## effective number of (beta, b) parameters given the variance
## parameters, plus 1 for sigma.  Monte Carlo sd of pD at 2000 draws is
## about 0.15 here, so the tolerances are several sds wide.
## =====================================================================

library(arm)

D.hand <- function(fit) {
  if (isGLMM(fit)) return(fit@resp$aic())
  -2 * sum(dnorm(model.response(model.frame(fit)), fitted(fit), sigma(fit),
                 log = TRUE))
}
dic <- function(fit) { set.seed(1); extractDIC(fit, n.sims = 2000) }

fmREML <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fmML   <- update(fmREML, REML = FALSE)
gb <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
set.seed(3)
g <- rep(1:20, each = 10)
x <- rnorm(200)
yp <- rpois(200, exp(0.5 + rnorm(20, 0, 0.5)[g] + 0.3 * x))
gp <- glmer(yp ~ x + (1 | g), family = poisson)

## ---- pD >= 0 ----------------------------------------------------------
for (fit in list(fmML, fmREML, gb, gp)) {
  d <- dic(fit)
  stopifnot("DIC implies a negative pD against D at the fitted values" =
              (d - D.hand(fit)) / 2 > -0.5)
  stopifnot("extractDIC() pD is missing or negative" =
              isTRUE(attr(d, "pD") >= 0))
}

## ---- (1 | Subject): p_fixed <= pD <= p_fixed + n_groups + 1 ----------
fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
pD1 <- attr(dic(fm1), "pD")
stopifnot("pD for (1 | Subject) is outside [2, 2 + 18 + 1]" =
            isTRUE(pD1 >= 2 && pD1 <= 21))

## ---- glmer: p_fixed <= pD <= p_fixed + n_groups -----------------------
pDb <- attr(dic(gb), "pD")
stopifnot("pD for the cbpp glmer is outside [4, 4 + 15]" =
            isTRUE(pDb >= 4 && pDb <= 19))

## ---- no group effect: pD is about p_fixed + 1 -------------------------
set.seed(4)
xs <- rnorm(200)
ys <- 1 + 2 * xs + rnorm(200)
fs <- suppressMessages(lmer(ys ~ xs + (1 | g)))
stopifnot("fixture is not a singular fit" = getME(fs, "theta") < 1e-4)
pDs <- attr(dic(fs), "pD")
stopifnot("pD for a fit with no group variance is not about 2 + 1" =
            isTRUE(abs(pDs - 3) < 0.5))

## ---- lmer: DIC = D(fitted) + 2 pD and pD = tr(H) + 1 ------------------
for (fit in list(fm1, fmML, fmREML, fs)) {
  d <- dic(fit)
  stopifnot("lmer pD is not trace(hat matrix) + 1" =
              isTRUE(abs(attr(d, "pD") - (sum(hatvalues(fit)) + 1)) < 1))
  stopifnot("lmer DIC is not D(fitted) + 2 (trace(hat matrix) + 1)" =
              abs(d - (D.hand(fit) + 2 * (sum(hatvalues(fit)) + 1))) < 2)
}

## ---- display() reports the same DIC and leaves the RNG alone ----------
invisible(capture.output(o <- display(fm1)))
set.seed(1)
stopifnot("display() DIC differs from extractDIC()" =
            isTRUE(all.equal(o$DIC, unname(c(extractDIC(fm1))))))
set.seed(5); a <- runif(1)
set.seed(5); invisible(capture.output(display(fm1))); b <- runif(1)
stopifnot("display() moved the caller's random number stream" = identical(a, b))
rm(".Random.seed", envir = globalenv())
invisible(capture.output(display(fm1)))
stopifnot("display() created a .Random.seed that was not there" =
            !exists(".Random.seed", envir = globalenv(), inherits = FALSE))
