## =====================================================================
## Regression test for sim() on merMod  --  issue #26 / pull request #27
##
## R CMD check runs every *.R file in tests/ and fails if any of them
## errors, so these need no test framework and add no dependency.
##
## The bug: arm::sim() built its draw matrix from the sparse Cholesky
## factor getME(fit, "L"), which carries a fill-reducing permutation, and
## never undid that permutation; it also transposed Lambda.  Fixed effects
## came out right, so the bug survived for years, but the random-effect
## draws had variances far from the model's.
##
## Choosing the fixture is the whole difficulty.  The bug is invisible on
## most models: with only (1 | g) terms Lambda is diagonal and t(Lambda)
## == Lambda, and with a trivial permutation there is nothing to undo.  A
## small crossed fixture (20 x 10) has both properties and still shows
## only ~10% error -- it passes against the broken code, so it is useless
## here.  The fixture below is the issue's own crossed, unbalanced design
## (cohorts of 30/100/100 pupils, (1 | pupil) + (1 + school_type | test)),
## where the broken code puts 53% of the random-effect variances outside
## +/-40% while the fixed code puts none outside +/-5%.
##
## The two guards below are asserted, not trusted: a fixture that quietly
## stops exercising the bug is worse than no test at all.
## =====================================================================

library(arm)
library(lme4)
library(Matrix)

## ---- fixture: crossed, unbalanced, with correlated slopes ------------
set.seed(2026)
n_pupil <- c(A = 30, B = 100, C = 100)
n_test  <- 36
cohort  <- rep(c("A", "B", "C"), each = n_test / 3)

pupils <- data.frame(cohort = rep(names(n_pupil), n_pupil))
pupils$pupil <- factor(seq_len(nrow(pupils)))
pupils$school_type <- factor(
  sample(c("public", "private", "charter"), nrow(pupils), TRUE, prob = c(.7, .2, .1)),
  levels = c("public", "private", "charter"))
pupils$sex <- factor(sample(c("F", "M"), nrow(pupils), TRUE))
pupils$ses <- stats::rnorm(nrow(pupils))
pupils$pupil_int <- stats::rnorm(nrow(pupils), 0, 0.1)

tests <- data.frame(test = factor(seq_len(n_test)), test_id = seq_len(n_test),
                    cohort = cohort)

d <- merge(pupils, tests, by = "cohort")
d <- d[sample(nrow(d), nrow(d) * 0.5), ]          # each pupil sits ~half the tests
d$b     <- stats::rnorm(nrow(d), 0, 0.2)
d$score <- 0.5 + d$pupil_int + d$b + 0.03 * d$ses + stats::rnorm(nrow(d), 0, 0.15)

fit <- suppressMessages(lme4::lmer(
  score ~ school_type + sex + ses + (1 | pupil) + (1 + school_type | test), d))

## ---- the fixture must actually exercise the bug ----------------------
stopifnot("fixture has an identity permutation, so it cannot catch #26" =
            any(lme4::getME(fit, "L")@perm !=
                  seq_along(lme4::getME(fit, "L")@perm) - 1L))

Lm <- as.matrix(lme4::getME(fit, "Lambda"))
stopifnot("fixture's Lambda is diagonal, so the transpose bug cannot show" =
            any(Lm[row(Lm) != col(Lm)] != 0))

## ---- the exact covariance the draws should have ----------------------
## Computed independently of the code under test, from Zt, Lambdat, X and
## sigma.  Vb is the *marginal* posterior covariance of the random effects
## b, in the natural (unpermuted) order.
##
## Careful: Var(u) is s2 times the [1:q, 1:q] block of the INVERSE of M,
## which is not s2 * solve(M[1:q, 1:q]) -- a block of an inverse is not the
## inverse of a block, and here the difference is a factor of ~1.6.  Get
## this wrong and the test fails against correct code.
q  <- lme4::getME(fit, "q")
p  <- lme4::getME(fit, "p")
Zt <- lme4::getME(fit, "Zt")
X  <- lme4::getME(fit, "X")
Lt <- lme4::getME(fit, "Lambdat")
s2 <- stats::sigma(fit)^2
Lambda <- lme4::getME(fit, "Lambda")

ZL <- Matrix::t(Zt) %*% Matrix::t(Lt)
M  <- rbind(
  cbind(Matrix::crossprod(ZL) + Matrix::Diagonal(ncol(ZL)), Matrix::crossprod(ZL, X)),
  cbind(Matrix::crossprod(X, ZL),                            Matrix::crossprod(X)))
Mi <- s2 * Matrix::solve(M)
Vb <- as.matrix(Lambda %*% Mi[seq_len(q), seq_len(q)] %*% Matrix::t(Lambda))

## sanity: the fixed-effect block of the same inverse must equal vcov(fit).
## If this fails the test's own maths is wrong, not arm's.
stopifnot("oracle disagrees with vcov(fit), so the test's maths is wrong" =
            max(abs(as.matrix(Mi[q + seq_len(p), q + seq_len(p)]) -
                      as.matrix(stats::vcov(fit)))) < 1e-8)

## ---- draw, reassemble (n_sims x q) in natural order, and compare -----
n_sims <- 1000
set.seed(7)
sims <- arm::sim(fit, n_sims)

parts <- lapply(names(lme4::getME(fit, "cnms")), function(g) {
  a <- sims@ranef[[g]]                       # (n_sims, n_groups, n_coefficients)
  n_g <- dim(a)[2]; n_c <- dim(a)[3]
  ## the b vector is level-major, coefficient-minor
  out <- matrix(NA_real_, n_sims, n_g * n_c)
  for (k in seq_len(n_c)) out[, seq(k, by = n_c, length.out = n_g)] <- a[, , k]
  out
})
D <- do.call(cbind, parts)
stopifnot("reassembled draw matrix has the wrong width" = ncol(D) == q)

ratio <- diag(stats::cov(D)) / diag(Vb)

stopifnot("sim() draws have wildly wrong variances (issue #26 is back)" =
            all(ratio > 0.6 & ratio < 1.6))
stopifnot("sim() draws are systematically off" =
            stats::median(ratio) > 0.9 && stats::median(ratio) < 1.1)

## the draws must also be centred where ranef() is -- a variance ratio of 1
## cannot detect a mean assigned to the wrong group
mu <- colMeans(D)
se <- apply(D, 2, stats::sd) / sqrt(n_sims)
stopifnot("sim() draws are not centred on getME(fit, 'b')" =
            max(abs(mu - lme4::getME(fit, "b")@x) / (6 * se + 1e-4)) < 1)

## and the fixed effects must still match vcov(fit), which they always did
stopifnot("sim() fixed-effect draws no longer match vcov(fit)" =
            max(abs(stats::var(sims@fixef) - as.matrix(stats::vcov(fit)))) <
              0.25 * max(abs(stats::vcov(fit))))
