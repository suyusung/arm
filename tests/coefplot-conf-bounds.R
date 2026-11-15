## =====================================================================
## Regression tests for coefplot() with lower.conf.bounds/upper.conf.bounds
##
## The bugs: coefplot.default() took its axis range from
## coefs.l/coefs.h = coefs -/+ CI*sds.
##  (a) With conf bounds but no sds, coefs.l was never defined, so the
##      documented call coefplot(coefs, lower.conf.bounds=, upper.conf.bounds=)
##      failed with "object 'coefs.l' not found".
##  (b) With sds plus conf bounds, the segments drawn are the bounds but the
##      range was still coefs -/+ 2*sds, so bounds outside it were clipped
##      silently.
## The fix takes the range from the interval endpoints that are drawn.
## =====================================================================

library(arm)

b <- c(0.2, 1.4, 2.3)
s <- c(0.1, 0.2, 0.2)
lo <- c(-3, 1.0, 2.0)
hi <- c(0.5, 6, 2.6)

## run a plot on a null device and return par("usr") afterwards
usr_after <- function(...) {
  pdf(NULL)
  on.exit(dev.off())
  res <- tryCatch({coefplot(...); par("usr")},
                  error = function(e) conditionMessage(e))
  res
}
## the visible coefficient axis: x for vertical plots, y for horizontal ones
covers <- function(usr, vertical) {
  if (is.character(usr)) return(FALSE)
  r <- if (vertical) usr[1:2] else usr[3:4]
  r[1] <= min(lo) && r[2] >= max(hi)
}

for (vertical in c(TRUE, FALSE)) {
  u <- usr_after(b, lower.conf.bounds = lo, upper.conf.bounds = hi,
                 vertical = vertical)
  stopifnot("conf bounds without sds give an error ('coefs.l' not found)" =
              !is.character(u))
  stopifnot("axis range without sds does not cover the conf bounds" =
              covers(u, vertical))

  u <- usr_after(b, s, lower.conf.bounds = lo, upper.conf.bounds = hi,
                 vertical = vertical)
  stopifnot("conf bounds outside coefs +/- 2 sds are clipped" =
              covers(u, vertical))
}

## without conf bounds the range is still coefs +/- CI*sds
for (CI in 1:2) {
  u <- usr_after(b, s, CI = CI)
  r <- range(b - CI * s, b + CI * s)
  stopifnot("range without conf bounds is no longer coefs +/- CI*sds" =
              isTRUE(all.equal(u[1:2], r + c(-1, 1) * 0.04 * diff(r))))
}
