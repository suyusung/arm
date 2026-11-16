## =====================================================================
## Regression tests for plot.balance()  --  the ... argument  --  d2
##
## The bug: plot.balance() accepts ... ("other plot options may be passed
## to this function", man/balance.Rd) and forwards them to plot() only on
## the plot = FALSE branch.  With the default plot = TRUE the extra
## arguments were silently dropped.
##
## `type` is already supplied by plot.balance(), so forwarding a `type`
## through ... must collide -- "formal argument \"type\" matched by
## multiple actual arguments".  If ... is dropped there is no error at all,
## which is what the test below detects.
##
## The fix forwards ... on the plotting branch too.
## =====================================================================

library(arm)

set.seed(9)
n <- 100
z <- rbinom(n, 1, 0.4)
X <- data.frame(age = rnorm(n), female = rbinom(n, 1, 0.5), income = rnorm(n))
invisible(capture.output(B <- balance(X, z, rep(1, n))))

bad <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

grDevices::pdf(NULL)
forwarded <- bad(capture.output(plot(B, type = "l")))
grDevices::dev.off()

stopifnot("plot.balance(plot = TRUE) dropped its ... argument" = forwarded)
