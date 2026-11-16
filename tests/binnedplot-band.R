## =====================================================================
## Regression test for the binned-residual band  --  issue #53
##
## The bug: binned.resids() computed each bin's +-2 SE band as
## 2*sd(y[bin])/sqrt(n), the empirical standard deviation of the residuals in
## that bin.  In a low-probability bin where every outcome is 0 the residuals
## are all -p, so they barely vary, the sd is near 0 and the band collapses
## onto the axis -- a correct model then shows 40%+ of its bins outside the
## band instead of about 5%.  Gelman & Hill (2007), section 5.6 use the
## model-based 2*sqrt(p*(1-p)/n) with p the bin's mean expected value.
##
## The fix: that model-based band is now the default, and the empirical one is
## kept behind se = "empirical".  The model-based band assumes x is on the
## probability scale and y is the response residual y - p, so binned.resids()
## warns when x falls outside [0, 1].
## =====================================================================

library(arm)

## ---- the default band is the model-based one ---------------------------
set.seed(42)
n <- 900
p <- runif(n, 0.01, 0.05)          # low-probability bins: where it collapsed
r <- rbinom(n, 1, p) - p
nb <- floor(sqrt(n))
b <- binned.resids(p, r, nb)$binned

stopifnot("the default band is not 2*sqrt(p*(1-p)/n) (issue #53 is back)" =
            isTRUE(all.equal(unname(b[, "2se"]),
                             2 * sqrt(b[, "xbar"] * (1 - b[, "xbar"]) / b[, "n"]))))

stopifnot("the band collapses in low-probability bins (issue #53 is back)" =
            all(b[, "2se"] > 1e-3))

## a correct model should keep roughly 5% of its bins outside the band; the
## collapsing one put 40% outside
stopifnot("too many bins outside the default band for a true model" =
            mean(abs(b[, "ybar"]) > b[, "2se"]) < 0.15)

## ---- se = "empirical" still gives the old band -------------------------
be <- binned.resids(p, r, nb, se = "empirical")$binned
emp <- 2 * sapply(split(seq_along(p), cut(p, binned.resids(p, r, nb)$xbreaks)),
                  function(i) sd(r[i]) / sqrt(length(i)))

stopifnot("se = 'empirical' is not 2*sd(y)/sqrt(n)" =
            isTRUE(all.equal(sort(unname(be[, "2se"])), sort(unname(emp[!is.na(emp)])))))

stopifnot("the empirical band no longer collapses where it used to" =
            min(be[, "2se"]) < 1e-3)
stopifnot("the empirical band is no longer the one that misbehaves" =
            mean(abs(be[, "ybar"]) > be[, "2se"]) > 0.25)

## ---- the two bands agree when p is mid-range ---------------------------
set.seed(7)
p2 <- runif(900, 0.3, 0.7)
r2 <- rbinom(900, 1, p2) - p2
b2 <- binned.resids(p2, r2, 30)$binned
e2 <- binned.resids(p2, r2, 30, se = "empirical")$binned
stopifnot("the bands disagree where they used to agree" =
            max(abs(b2[, "2se"] - e2[, "2se"])) < 0.05)

## ---- x outside [0, 1] is not a probability scale -----------------------
w <- tryCatch({ binned.resids(seq(0, 3, length.out = 50), rnorm(50), 5); NULL },
              warning = function(w) conditionMessage(w))
stopifnot("no warning when x is outside [0, 1] with the model band" =
            !is.null(w) && grepl("response", w))
w2 <- tryCatch({ binned.resids(runif(50), rnorm(50), 5); NULL },
               warning = function(w) conditionMessage(w))
stopifnot("a warning was raised for a valid probability-scale x" = is.null(w2))

## ---- a bad `se` is an error, not a silent wrong answer -----------------
stopifnot("an unknown se should be an error" =
            inherits(tryCatch(binned.resids(p, r, 10, se = "model-based"),
                              error = function(e) e), "error"))

## ---- both bands plot ---------------------------------------------------
grDevices::pdf(tempfile(fileext = ".pdf"))
binnedplot(p, r)
binnedplot(p, r, se = "empirical")
grDevices::dev.off()
