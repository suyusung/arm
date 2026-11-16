## =====================================================================
## Regression tests for se.coef() / se.ranef() on merMod fits
##
## The bug: when one grouping factor has more than one random-effects
## term, e.g. (1 | Subject) + (0 + Days | Subject), lme4 stores the
## "postVar" attribute of ranef(fit, condVar = TRUE) as a list of arrays
## (one per term) instead of one K x K x J array.  se.coef() and
## se.ranef() took dim() of that list, got NULL, and stopped with
## "'dim' cannot be of length 0".
##
## The standard errors must equal the conditional sds lme4 reports in
## as.data.frame(ranef(fit, condVar = TRUE))$condsd, and fits with one
## term per factor must be unchanged.
## =====================================================================

library(arm)

## condsd from lme4, reshaped to the J x K layout of se.ranef()
condsd <- function(fit, grp) {
  d <- as.data.frame(ranef(fit, condVar = TRUE))
  d <- d[d$grpvar == grp, ]
  re <- ranef(fit)[[grp]]
  out <- sapply(colnames(re), function(term)
    d$condsd[d$term == term][match(rownames(re), d$grp[d$term == term])])
  matrix(out, nrow = nrow(re), dimnames = dimnames(re))
}

## ---- two uncorrelated terms for one factor ----------------------------
fm2 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), sleepstudy)
r2 <- try(se.ranef(fm2), silent = TRUE)
stopifnot("se.ranef() fails when postVar is a list (two terms, one factor)" =
            !inherits(r2, "try-error"))
stopifnot("se.ranef() on a two-term factor differs from lme4's condsd" =
            isTRUE(all.equal(r2$Subject, condsd(fm2, "Subject"))))
c2 <- try(se.coef(fm2), silent = TRUE)
stopifnot("se.coef() fails when postVar is a list (two terms, one factor)" =
            !inherits(c2, "try-error"))
stopifnot("se.coef() on a two-term factor differs from lme4's condsd" =
            isTRUE(all.equal(c2$Subject, condsd(fm2, "Subject"))))

## ---- a 2 x 2 term plus a 1 x 1 term for the same factor ---------------
sleepstudy$Days2 <- sleepstudy$Days^2 / 10
fm3 <- lmer(Reaction ~ Days + (Days | Subject) + (0 + Days2 | Subject),
            sleepstudy, control = lmerControl(check.conv.singular = "ignore"))
stopifnot("se.ranef() on a 2x2 + 1x1 factor differs from lme4's condsd" =
            isTRUE(all.equal(se.ranef(fm3)$Subject, condsd(fm3, "Subject"))))

## ---- single-term fits: unchanged and equal to condsd ------------------
old.se <- function(vars.m) {
  K <- dim(vars.m)[1]; J <- dim(vars.m)[3]
  out <- array(NA, c(J, K))
  for (j in 1:J) out[j, ] <- sqrt(diag(as.matrix(vars.m[, , j])))
  out
}
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
pv1 <- attr(ranef(fm1, condVar = TRUE)$Subject, "postVar")
stopifnot("se.ranef() changed for a correlated single-term fit" =
            isTRUE(all.equal(unname(se.ranef(fm1)$Subject), old.se(pv1))))
stopifnot("se.coef() differs from lme4's condsd for a single-term fit" =
            isTRUE(all.equal(se.coef(fm1)$Subject, condsd(fm1, "Subject"))))

fmP <- lmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin)
sP <- se.ranef(fmP)
for (g in c("plate", "sample")) {
  stopifnot("se.ranef() differs from lme4's condsd for two grouping factors" =
              isTRUE(all.equal(sP[[g]], condsd(fmP, g))))
  stopifnot("se.ranef() changed for a fit with two grouping factors" =
              isTRUE(all.equal(unname(sP[[g]]),
                               old.se(attr(ranef(fmP, condVar = TRUE)[[g]], "postVar")))))
}
