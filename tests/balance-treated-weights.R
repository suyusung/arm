## =====================================================================
## Regression tests for balance()  --  weights on the treated units
##
## The bug: `matched` is documented as a vector of weights "to apply to the
## full dataset", including IPTW weights, and estimand may be ATT, ATC or
## ATE.  But balance() applied the weights only to the controls: the
## treated mean and variance after "matching" were the raw, unweighted
## ones.  For ATE or ATC weights, which reweight the treated units too, the
## reported treated mean, difference and standardized difference were the
## wrong quantities.  For example, with IPTW ATE weights the age column
## below reported a treated mean of 44.27 (the raw mean) instead of the
## weighted 40.19, and an abs. std. diff of 0.417 instead of 0.007.
##
## The fix weights the treated units as well.  When every treated weight is
## 1 -- the ATT case, including the output of matching() -- the result is
## unchanged.
## =====================================================================

library(arm)

set.seed(7)
n <- 400
X <- data.frame(age = rnorm(n, 40, 10), inc = rnorm(n, 50, 15))
z <- rbinom(n, 1, plogis(-0.5 + 0.08 * (X$age - 40)))
ps <- fitted(glm(z ~ age + inc, data = X, family = binomial))

wmean <- function(x, w) sum(w * x) / sum(w)
## weighted variance with the same (unbiased) normalisation as cov.wt()
wvar <- function(x, w) {
  w <- w / sum(w)
  sum(w * (x - sum(w * x))^2) / (1 - sum(w^2))
}

## ---- ATE inverse-probability weights reweight both groups -------------
w <- ifelse(z == 1, 1 / ps, 1 / (1 - ps))
invisible(capture.output(B <- balance(X, z, w, estimand = "ATE")))
for (v in names(X)) {
  x <- X[[v]]
  mt <- wmean(x[z == 1], w[z == 1])
  mc <- wmean(x[z == 0], w[z == 0])
  denom <- sqrt((var(x[z == 1]) + var(x[z == 0])) / 2)
  stopifnot("balance() ignored the weights of the treated units (treat mean)" =
              isTRUE(all.equal(B$diff.means.matched[v, "treat"], mt)))
  stopifnot("balance() control mean is not the weighted control mean" =
              isTRUE(all.equal(B$diff.means.matched[v, "control"], mc)))
  stopifnot("balance() ignored the weights of the treated units (std diff)" =
              isTRUE(all.equal(B$diff.means.matched[v, "abs.std.diff"],
                               abs(mt - mc) / denom)))
  stopifnot("balance() ignored the weights of the treated units (sd ratio)" =
              isTRUE(all.equal(B$diff.means.matched[v, "ratio"],
                               sqrt(wvar(x[z == 0], w[z == 0]) /
                                      wvar(x[z == 1], w[z == 1])))))
}

## ---- ATT matching (treated weights all 1) is unchanged ----------------
m <- matching(z = z, score = ps, replace = TRUE)
invisible(capture.output(A <- balance(X, z, m$cnts, estimand = "ATT")))
for (v in names(X)) {
  x <- X[[v]]
  wc <- m$cnts[z == 0]
  mt <- mean(x[z == 1])
  mc <- wmean(x[z == 0], wc)
  stopifnot("ATT balance with unit treated weights changed" =
              isTRUE(all.equal(unname(A$diff.means.matched[v, ]),
                               c(mt, mc, mt - mc,
                                 abs(mt - mc) / sd(x[z == 1]),
                                 sqrt(wvar(x[z == 0], wc) / var(x[z == 1]))),
                               tolerance = 1e-12)))
}
