## =====================================================================
## Regression test for coefplot() with a CI other than 1 or 2  --  p1
##
## The bug: coefplot.default() drew a fixed +/-1 sd bar for every CI other
## than 2 (vertical = TRUE) and nothing at all for a CI other than 1 or 2
## (vertical = FALSE), while the axis range was coefs +/- CI*sds.  A call
## coefplot(b, s, CI = 3) therefore showed a 3-sd axis under a 1-sd bar, or
## an empty plot when vertical = FALSE.
##
## The fix: the outer interval drawn is coefs +/- CI*sds, so any CI plots
## that interval over the +/-1 sd inner bar.  CI = 1 and CI = 2 keep their
## documented behaviour.
## =====================================================================

library(arm)

segs <- new.env()
segs$x <- NULL
invisible(suppressMessages(
  trace("segments", where = asNamespace("graphics"),
        quote(segs$x <- rbind(segs$x, cbind(x0 = x0, y0 = y0, x1 = x1, y1 = y1))),
        print = FALSE)))

b <- c(0.2, 1.4, 2.3)
s <- c(0.1, 0.2, 0.2)

## the coefficient-axis endpoints of the segments drawn by one coefplot() call
draw <- function(CI, vertical, add = FALSE) {
  segs$x <- NULL
  pdf(NULL)
  on.exit(dev.off())
  coefplot(b, s, CI = CI, vertical = vertical)
  if (add) {
    segs$x <- NULL
    coefplot(b, s, CI = CI, vertical = vertical, add = TRUE)
  }
  if (is.null(segs$x)) return(numeric(0))
  if (vertical) c(segs$x[, "x0"], segs$x[, "x1"]) else c(segs$x[, "y0"], segs$x[, "y1"])
}

## CI = 3: the +/- 3 sd bar over the +/- 1 sd bar, in both orientations
for (vertical in c(TRUE, FALSE)) {
  x <- draw(3, vertical)
  stopifnot("coefplot(CI = 3) draws no interval when vertical = FALSE" =
              length(x) > 0)
  stopifnot("coefplot(CI = 3) draws only the +/- 1 sd bar" =
              length(x) == 4 * length(b))
  stopifnot("coefplot(CI = 3) does not draw the coefs +/- 3 sds interval" =
              isTRUE(all.equal(range(x), range(b - 3 * s, b + 3 * s))))
}

## the same interval is drawn on top when add = TRUE
for (vertical in c(TRUE, FALSE)) {
  x <- draw(3, vertical, add = TRUE)
  stopifnot("coefplot(CI = 3, add = TRUE) does not draw the +/- 3 sd interval" =
              length(x) == 4 * length(b) &&
                isTRUE(all.equal(range(x), range(b - 3 * s, b + 3 * s))))
}

## the documented behaviour is unchanged
x1 <- draw(1, TRUE)
stopifnot("CI = 1 no longer draws the +/- 1 sd interval only" =
            length(x1) == 2 * length(b) &&
              isTRUE(all.equal(range(x1), range(b - s, b + s))))
x2 <- draw(2, TRUE)
stopifnot("CI = 2 no longer draws +/- 2 sds over +/- 1 sd" =
            length(x2) == 4 * length(b) &&
              isTRUE(all.equal(range(x2), range(b - 2 * s, b + 2 * s))))

invisible(suppressMessages(untrace("segments", where = asNamespace("graphics"))))
