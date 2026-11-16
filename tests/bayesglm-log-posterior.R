## =====================================================================
## Regression test: the printed unnormalized log posterior  --  issue #54 (minor)
##
## The bug: with print.unnormalized.log.posterior = TRUE the code read
##
##   logprior <- sum(dt(coefs.hat, prior.df, prior.mean, log = TRUE))
##   xb <- invlogit(x %*% coefs.hat)
##
## dt()'s third argument is the non-centrality, not the location, and it has
## no scale argument at all, so the prior mean was not the prior mean and the
## prior scale was ignored; and invlogit() was used whatever the link, so a
## probit fit printed a log likelihood that did not match logLik().  The log
## likelihood also ignored the binomial totals (the weights).
## =====================================================================

library(arm)

set.seed(1)
n <- 200
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))

printed <- function(fit) {
  txt <- paste(capture.output(fit <- eval(fit)), collapse = " ")
  c(logprior = as.numeric(sub(".*log prior: ([^,]+),.*", "\\1", txt)),
    loglik  = as.numeric(sub(".*log likelihood: ([^,]+),.*", "\\1", txt)),
    post    = as.numeric(sub(".*unnormalized log posterior: ([^,]+).*", "\\1", txt)))
}

## a probit fit: the printed log likelihood is the fit's own log likelihood
fp <- bayesglm(y ~ x1 + x2, binomial(link = "probit"), d,
               print.unnormalized.log.posterior = TRUE)
p <- printed(quote(bayesglm(y ~ x1 + x2, binomial(link = "probit"), d,
                            print.unnormalized.log.posterior = TRUE)))

stopifnot("the printed log likelihood is not the fit's logLik (wrong link?)" =
            isTRUE(all.equal(unname(p["loglik"]), as.numeric(logLik(fp)),
                             tolerance = 1e-6)))

## the printed log prior is the scaled-t prior of the fit, up to nothing
lp <- sum(dt((coef(fp) - fp$prior.mean) / fp$prior.scale, fp$prior.df, log = TRUE) -
            log(fp$prior.scale))
stopifnot("the printed log prior is not the fit's prior density" =
            isTRUE(all.equal(unname(p["logprior"]), lp, tolerance = 1e-6)))

stopifnot("the printed posterior is not likelihood + prior" =
            isTRUE(all.equal(unname(p["post"]), unname(p["loglik"] + p["logprior"]),
                             tolerance = 1e-6)))

## and the same for the default logit link
fl <- bayesglm(y ~ x1 + x2, binomial, d, print.unnormalized.log.posterior = TRUE)
pl <- printed(quote(bayesglm(y ~ x1 + x2, binomial, d,
                             print.unnormalized.log.posterior = TRUE)))
stopifnot("the printed log likelihood is wrong for the logit link too" =
            isTRUE(all.equal(unname(pl["loglik"]), as.numeric(logLik(fl)),
                             tolerance = 1e-6)))
