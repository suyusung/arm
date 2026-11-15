## =====================================================================
## Regression test for coefplot() on a "bugs" object  --  pch.pts
##
## The bug: the bugs method has a documented pch.pts argument, but every
## points() call inside it hard-coded pch=20, so pch.pts was ignored.
## The test traces graphics::points.default to record the pch actually
## drawn, for CI = 1 and 2, vertical and horizontal, with and without add.
## =====================================================================

library(arm)

## a minimal bugs object: coefplot only reads $summary quantiles
sm <- cbind("2.5%" = c(-1, 0, 1), "25%" = c(-.5, .5, 1.5), "50%" = c(0, 1, 2),
            "75%" = c(.5, 1.5, 2.5), "97.5%" = c(1, 2, 3))
rownames(sm) <- c("alpha", "beta", "gamma")
fake <- structure(list(summary = sm), class = "bugs")

seen <- new.env()
invisible(suppressMessages(
  trace("points.default", where = asNamespace("graphics"),
        quote(seen$pch <- c(seen$pch, list(...)$pch)), print = FALSE)))

pdf(NULL)
for (CI in 1:2) for (vertical in c(TRUE, FALSE)) {
  seen$pch <- NULL
  coefplot(fake, CI = CI, vertical = vertical, pch.pts = 17)
  coefplot(fake, CI = CI, vertical = vertical, pch.pts = 17, add = TRUE)
  stopifnot("coefplot(<bugs>) ignores pch.pts (pch=20 hard-coded)" =
              length(seen$pch) == 2 && all(seen$pch == 17))
}
## the default is still a solid dot
seen$pch <- NULL
coefplot(fake)
stopifnot("default pch for coefplot(<bugs>) is no longer 20" =
            identical(seen$pch, 20))
invisible(dev.off())
suppressMessages(untrace("points.default", where = asNamespace("graphics")))
