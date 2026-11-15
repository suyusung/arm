## =====================================================================
## Regression tests for plot.balance()  --  labels with which.covs
##
## The bug: with which.covs = "binary" or "cont", plot.balance() keeps only
## the matching covariates and subsets `covnames` accordingly, but a
## user-supplied `longcovnames` was not subset.  The y-axis then used the
## first K long names, so the labels sat next to the wrong points (with
## binary covariates 2 and 4, the plot labelled them with the names of
## covariates 1 and 2).
##
## The fix subsets `longcovnames` the same way as `covnames`.  The labels
## passed to axis(2) are captured with trace() on a null device.
## =====================================================================

library(arm)

set.seed(9)
n <- 200
z <- rbinom(n, 1, 0.4)
X <- data.frame(age = rnorm(n), female = rbinom(n, 1, 0.5),
                income = rnorm(n), married = rbinom(n, 1, 0.5))
invisible(capture.output(B <- balance(X, z, rep(1, n))))
long <- c("Age (years)", "Female", "Income ($1000)", "Married")
stopifnot("test setup: expected covariates 2 and 4 to be binary" =
            identical(B$binary, c(0, 1, 0, 1)))

captured <- new.env()
trace("axis", where = asNamespace("arm"), print = FALSE,
      tracer = quote(if (side == 2) assign("labels", labels, envir = captured)))
grDevices::pdf(NULL)
y_labels <- function(...) {
  invisible(capture.output(plot(B, ...)))
  get("labels", envir = captured)
}

## the plot runs from top to bottom in the original order, so axis(2)
## receives the labels in reverse
stopifnot("plot.balance(which.covs = 'binary') put long names on the wrong points" =
            identical(y_labels(which.covs = "binary", longcovnames = long),
                      rev(long[B$binary == 1])))
stopifnot("plot.balance(which.covs = 'cont') put long names on the wrong points" =
            identical(y_labels(which.covs = "cont", longcovnames = long),
                      rev(long[B$binary == 0])))
stopifnot("plot.balance(which.covs = 'mixed') long names changed" =
            identical(y_labels(which.covs = "mixed", longcovnames = long), rev(long)))
stopifnot("plot.balance(which.covs = 'binary') default labels changed" =
            identical(y_labels(which.covs = "binary"),
                      rev(names(X)[B$binary == 1])))

grDevices::dev.off()
untrace("axis", where = asNamespace("arm"))
