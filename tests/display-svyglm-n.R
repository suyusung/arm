## =====================================================================
## Regression tests for display() on svyglm fits
##
## The bug: display() reported n as summary(fit)$df[1] + summary(fit)$df[2].
## For a svyglm the second element is the residual design degrees of
## freedom, which survey counts in primary sampling units, not
## observations.  For the two-stage apiclus2 design (40 districts, 126
## schools) display() printed "n = 40".
##
## n must be the number of observations used in the fit.  The checks use
## the row count of the data, not nobs(), so they share no code with the
## fix.
## =====================================================================

library(arm)
library(survey)

data(api)
dclus2 <- svydesign(id = ~dnum + snum, weights = ~pw, data = apiclus2)
g <- svyglm(api00 ~ ell + meals, design = dclus2)
invisible(capture.output(o <- display(g)))
stopifnot("display(svyglm) n is not the number of observations" =
            o$n == nrow(apiclus2))
stopifnot("display(svyglm) k is not the number of coefficients" =
            o$k == 3)

## rows dropped for missing values are not counted
apin <- apiclus2
apin$ell[1:5] <- NA
gn <- svyglm(api00 ~ ell + meals,
             design = svydesign(id = ~dnum + snum, weights = ~pw, data = apin))
invisible(capture.output(on <- display(gn)))
stopifnot("display(svyglm) n counts rows dropped for NA" =
            on$n == sum(complete.cases(apin[, c("api00", "ell", "meals")])))

## a subdomain counts only its own observations
ge <- svyglm(api00 ~ ell + meals, design = subset(dclus2, stype == "E"))
invisible(capture.output(oe <- display(ge)))
stopifnot("display(svyglm) n on a subset design is not the subset size" =
            oe$n == sum(apiclus2$stype == "E"))

## the printed line agrees with the returned value
txt <- capture.output(display(g))
stopifnot("display(svyglm) prints the wrong n" =
            any(grepl(paste0("n = ", nrow(apiclus2), ", k = 3"), txt, fixed = TRUE)))
