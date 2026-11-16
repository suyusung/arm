## =====================================================================
## Regression tests for the "Error terms" table of display(merMod)
##
## The bug: for binomial and poisson glmer fits, which have no residual
## scale parameter, display() printed a fabricated "Residual 1.00" SD.
## The helper as.matrix.VarCorr() always appends a Residual row, filling
## a missing scale with 1, and was meant to blank it with
## `if (useScale < 0)`.  But useScale is lme4's useSc flag, which is 0 or
## 1, so that test was never true.
##
## A model without a scale parameter must print no Residual row; lmer
## output must still print the residual sd, equal to sigma(fit).
## =====================================================================

library(arm)

error.terms <- function(fit) {
  txt <- capture.output(display(fit))
  txt[seq(grep("^Error terms:", txt) + 1, grep("^---", txt) - 1)]
}

## ---- binomial glmer ---------------------------------------------------
gb <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
eb <- error.terms(gb)
stopifnot("display(glmer binomial) prints a Residual sd" =
            !any(grepl("Residual", eb)))
stopifnot("display(glmer binomial) lost the herd sd" =
            any(grepl(paste0("herd +\\(Intercept\\) +",
                             fround(attr(VarCorr(gb)$herd, "stddev"), 2)), eb)))

## ---- poisson glmer, correlated random effects --------------------------
set.seed(3)
g <- rep(1:20, each = 10)
x <- rnorm(200)
y <- rpois(200, exp(0.5 + rnorm(20, 0, 0.5)[g] + (0.3 + rnorm(20, 0, 0.3)[g]) * x))
gp <- glmer(y ~ x + (x | g), family = poisson)
ep <- error.terms(gp)
stopifnot("display(glmer poisson) prints a Residual sd" =
            !any(grepl("Residual", ep)))
stopifnot("display(glmer poisson) lost a random-effect row" =
            length(ep) == 3)

## ---- lmer still prints the residual sd ---------------------------------
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
el <- error.terms(fm)
stopifnot("display(lmer) no longer prints the residual sd" =
            any(grepl(paste0("^ *Residual +", fround(sigma(fm), 2)), el)))
stopifnot("display(lmer) error-terms table changed shape" =
            length(el) == 4)
