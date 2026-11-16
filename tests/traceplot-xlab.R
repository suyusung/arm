## =====================================================================
## Regression test for the traceplot() x axis label  --  p4
##
## The bug: the x axis was labelled "iteration", but the values plotted are
## the kept draws (1:n.keep), not the iteration numbers of the chain -- with
## burn-in and thinning the two differ, so the label read as if it were the
## iteration number.
##
## The fix: the axis is labelled "Kept draws", which is exactly what the
## positions are.
## =====================================================================

library(arm)

sm <- array(rnorm(50 * 2 * 2), dim = c(50, 2, 2),
            dimnames = list(NULL, c("chain1", "chain2"), c("alpha", "beta")))
fake <- structure(list(n.chains = 2L, n.keep = 50L, sims.array = sm),
                  class = "bugs")

seen <- new.env()
seen$xlab <- NULL
invisible(suppressMessages(
  trace("plot", where = asNamespace("graphics"),
        quote(seen$xlab <- c(seen$xlab, list(...)$xlab)), print = FALSE)))

pdf(NULL)
traceplot(fake, ask = FALSE)
dev.off()

invisible(suppressMessages(untrace("plot", where = asNamespace("graphics"))))

stopifnot("traceplot() labels the kept draws as 'iteration'" =
            length(seen$xlab) == 2 && all(seen$xlab == "Kept draws"))
