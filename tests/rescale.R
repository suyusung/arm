## =====================================================================
## Regression tests for rescale()  --  issue #19
##
## The bug: rescale() turned a non-numeric input into numbers with
## as.numeric(factor(x)), which in R gives 1..k, not 0..k-1.  For a
## two-level factor with binary.inputs = "-0.5,0.5" that produced 0.5/1.5
## instead of -0.5/0.5.  The fix maps a binary factor onto 0/1 first.
##
## The same bug survived for *numeric* inputs: "-0.5,0.5" simply subtracted
## 0.5 from x, which is only correct when the variable is already coded 0/1,
## so rescale(c(1, 2), "-0.5,0.5") gave 0.5, 1.5.  It now maps the two
## observed values onto -0.5 and 0.5, so both codings agree.
##
## Also pinned: an unrecognised binary.inputs must be an error.  It used to
## fall through every branch and return NULL invisibly for a binary variable.
## =====================================================================

library(arm)

## ---- a two-level factor is a 0/1 variable, not a 1/2 variable ---------
f <- gl(2, 1)
stopifnot("binary factor not mapped to -0.5/0.5 (issue #19 is back)" =
            isTRUE(all.equal(rescale(f, binary.inputs = "-0.5,0.5"), c(-0.5, 0.5))))
stopifnot("'0/1' did not map a binary factor onto 0/1" =
            isTRUE(all.equal(rescale(f, binary.inputs = "0/1"), c(0, 1))))
stopifnot("'center' did not centre a binary factor" =
            isTRUE(all.equal(rescale(f, binary.inputs = "center"), c(-0.5, 0.5))))

## ---- a numeric 0/1 variable behaves the same --------------------------
stopifnot("numeric 0/1 did not map to -0.5/0.5" =
            isTRUE(all.equal(rescale(c(0, 1), binary.inputs = "-0.5,0.5"),
                             c(-0.5, 0.5))))
stopifnot("'0/1' did not normalise a numeric binary variable" =
            isTRUE(all.equal(rescale(c(2, 5), binary.inputs = "0/1"), c(0, 1))))

## ---- the two levels are mapped, whatever they happen to be coded ------
stopifnot("numeric 1/2 did not map to -0.5/0.5 (#19 for numeric input)" =
            isTRUE(all.equal(rescale(c(1, 2), binary.inputs = "-0.5,0.5"),
                             c(-0.5, 0.5))))
stopifnot("'0/1' is not invariant to the coding" =
            isTRUE(all.equal(rescale(c(7, 9), binary.inputs = "0/1"), c(0, 1))))
stopifnot("a factor and the equivalent numeric disagree" =
            isTRUE(all.equal(rescale(factor(c("a", "b")), binary.inputs = "-0.5,0.5"),
                             rescale(c(1, 2), binary.inputs = "-0.5,0.5"))))
## centring and full scaling were already coding-invariant
stopifnot("'center' is not invariant to the coding" =
            isTRUE(all.equal(rescale(c(1, 2), binary.inputs = "center"),
                             rescale(c(0, 1), binary.inputs = "center"))))
stopifnot("'full' is not invariant to the coding" =
            isTRUE(all.equal(rescale(c(1, 2), binary.inputs = "full"),
                             rescale(c(0, 1), binary.inputs = "full"))))

## ---- NAs pass through; the scale is taken from the observed values ----
f2 <- factor(rep(c("a", "b"), 5))
f2[1] <- NA
r <- rescale(f2, binary.inputs = "-0.5,0.5")
stopifnot("NA did not stay NA" = is.na(r[1]))
stopifnot("NA disturbed the other values" =
            isTRUE(all.equal(sort(unique(r[!is.na(r)])), c(-0.5, 0.5))))

## ---- deterministic: the same input must always give the same answer ---
## (issue #21 reported rescale() returning different results on different runs)
set.seed(1)
x <- stats::rnorm(50)
stopifnot("rescale() is not deterministic" =
            identical(rescale(x), rescale(x)))
stopifnot("rescale() depends on the RNG state for a fixed input" =
            identical(rescale(x), {set.seed(99); rescale(x)}))

## ---- a continuous variable is scaled by two standard deviations --------
xx <- c(1, 2, 3, 4, 5)
stopifnot("continuous variable not scaled by two standard deviations" =
            isTRUE(all.equal(rescale(xx), (xx - mean(xx)) / (2 * stats::sd(xx)))))

## ---- a factor with more than two levels ---------------------------------
f3 <- factor(c("a", "b", "b", "c", "c", "c"))
r3 <- rescale(f3)
stopifnot("multi-level factor should be scaled, not centred" =
            length(unique(r3)) == 3)

## ---- an unrecognised binary.inputs is an error, not a silent NULL ------
## For a binary variable, a value outside the four documented options used to
## fall through every branch and return NULL invisibly; for a continuous one
## it was ignored.  Both are now an error.
bad <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot("a bad binary.inputs gave no error for a binary numeric" =
            bad(rescale(c(0, 1), binary.inputs = "bogus")))
stopifnot("a bad binary.inputs gave no error for a factor" =
            bad(rescale(gl(2, 1), binary.inputs = "bogus")))
stopifnot("a bad binary.inputs gave no error for a continuous variable" =
            bad(rescale(c(1, 2, 3, 4), binary.inputs = "bogus")))
stopifnot("a non-character binary.inputs gave no error" =
            bad(rescale(gl(2, 1), binary.inputs = 1)))
stopifnot("a length-2 binary.inputs gave no error" =
            bad(rescale(gl(2, 1), binary.inputs = c("center", "full"))))
stopifnot("a NULL binary.inputs gave no error" =
            bad(rescale(gl(2, 1), binary.inputs = NULL)))

## the message has to name the offending value, or it is not actionable
msg <- tryCatch(rescale(gl(2, 1), binary.inputs = "bogus"),
                error = function(e) conditionMessage(e))
stopifnot("the error does not name the offending value" =
            grepl("bogus", msg, fixed = TRUE))
stopifnot("the error does not list the allowed values" =
            all(vapply(c("center", "0/1", "-0.5,0.5", "full"),
                       function(v) grepl(v, msg, fixed = TRUE), logical(1))))

## and the four documented values all still work
stopifnot("a documented binary.inputs was rejected" =
            all(vapply(c("center", "0/1", "-0.5,0.5", "full"),
                       function(v) !bad(rescale(gl(2, 1), binary.inputs = v)),
                       logical(1))))
