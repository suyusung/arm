## =====================================================================
## Regression test for binned.resids()$xbreaks
##
## The bug: the help page says binned.resids() returns xbreaks, "the
## breakpoints of the bins", but the function set xbreaks <- NULL and never
## filled it, so it was always NULL.  It now returns the breaks passed to
## cut(): nclass + 1 values, bin i being the interval (xbreaks[i], xbreaks[i+1]].
## =====================================================================

library(arm)

check <- function(x, nclass) {
  y <- stats::rnorm(length(x))
  res <- binned.resids(x, y, nclass)
  br <- res$xbreaks
  b <- res$binned
  stopifnot("binned.resids()$xbreaks is NULL" = !is.null(br))
  stopifnot("xbreaks does not have one more value than there are bins" =
              length(br) == nrow(b) + 1)
  stopifnot("xbreaks is not increasing" = all(diff(br) > 0))
  ## every bin's observed range lies inside its interval
  stopifnot("a bin's x values fall outside (xbreaks[i], xbreaks[i+1]]" =
              all(br[-length(br)] < b[, "x.lo"] & b[, "x.hi"] <= br[-1]))
  ## and the breaks reproduce the bin counts, so no x is outside its bin
  stopifnot("xbreaks does not reproduce the bin counts" =
              all(tabulate(findInterval(x, br, left.open = TRUE), nrow(b)) ==
                    b[, "n"]))
}

set.seed(1)
check(stats::runif(400), 20)
check(stats::rnorm(57), 10)
## heavy ties: some requested bins merge, xbreaks must follow
check(round(stats::runif(300) * 5), 17)
check(stats::rnorm(9), 1)

## Equal scores at multiple boundaries stay together. The first group of
## repeated scores has no lower neighbor.
x <- c(0.5, 1, 0, 0.5, 1, 0, 0.5, 0, 1, 0.5, 0, 1)
y <- seq_along(x) - x
z <- binned.resids(x, y, 6, se = "empirical")
by.score <- split(y, x)
stopifnot(identical(z$xbreaks, c(-Inf, 0.25, 0.75, Inf)),
          identical(unname(z$binned[, "n"]), c(4, 4, 4)),
          isTRUE(all.equal(unname(z$binned[, "ybar"]), unname(vapply(by.score, mean, 0)))),
          isTRUE(all.equal(unname(z$binned[, "2se"]), unname(vapply(by.score, sd, 0)))))

## A constant prediction with a single requested class leaves an empty second bin.
empty <- suppressWarnings(binned.resids(rep(0.5, 4), 1:4, 1))
stopifnot(identical(empty$binned[, "n"], c(4, 0)),
          is.nan(empty$binned[2, "ybar"]))
