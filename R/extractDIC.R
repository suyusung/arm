

extractDIC <- function(fit,...){
  UseMethod("extractDIC")
}


extractDIC.merMod <- function(fit, n.sims = 1000, ...){
        r <- dic.merMod(fit, n.sims)
        DIC <- r[["DIC"]]
        names(DIC) <- "DIC"
        attr(DIC, "pD") <- r[["pD"]]
        return(DIC)
}


# Conditional DIC as BUGS computes it, from sim() draws of the fixed
# effects, random effects and sigma, with the variance parameters held at
# their estimates.  D = -2 log p(y | beta, b, sigma), pD = mean of D over
# the draws minus D at the posterior means, and DIC = D(means) + 2 pD.
# A non-NULL seed is used for the draws and the caller's RNG state is
# restored afterwards, so display() does not move the user's stream.
dic.merMod <- function(fit, n.sims, seed = NULL){
    if (!is.null(seed)) {
        genv <- globalenv()
        if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
            old.seed <- get(".Random.seed", envir = genv, inherits = FALSE)
            on.exit(assign(".Random.seed", old.seed, envir = genv))
        } else {
            on.exit(rm(".Random.seed", envir = genv))
        }
        set.seed(seed)
    }
    is.glmm <- getME(fit, "devcomp")$dims[["GLMM"]] != 0L
    fam <- family(fit)
    if (is.glmm && !(fam$family %in% c("binomial", "poisson"))) {
        warning("DIC is computed only for lmer fits and binomial or poisson glmer fits")
        return(c(DIC = NA_real_, pD = NA_real_))
    }
    y <- fit@resp$y
    w <- fit@resp$weights
    D <- function(mu, sigma){
        if (is.glmm) {
            # -2 log p(y | mu) weighted, evaluated columnwise on the
            # (n obs x n.sims) matrix of draws -- reproduces exactly what
            # family$aic() computes per draw (binomial folds the weights
            # into the number of trials when any n > 1), vectorized.
            if (fam$family == "binomial") {
                n <- fit@resp$n
                m <- if (any(n > 1)) n else w
                coef <- ifelse(m > 0, w / m, 0)
                d <- dbinom(round(m * y), round(m), mu, log = TRUE)
            } else {
                coef <- w
                d <- dpois(y, mu, log = TRUE)
            }
            # density functions drop the matrix dims of mu -- restore the
            # (n obs x n.sims) shape, as the gaussian branch below does
            d <- matrix(d, nrow = length(y))
            d[coef == 0, ] <- 0    # zero-weight rows: 0 * -Inf would be NaN
            -2 * colSums(d * coef)
        } else {
            sd <- outer(1 / sqrt(w), sigma)
            -2 * colSums(matrix(dnorm(y, mu, sd, log = TRUE), nrow = length(y)))
        }
    }
    s <- sim(fit, n.sims)
    s.mean <- new("sim.merMod",
                  fixef = matrix(colMeans(s@fixef), 1,
                                 dimnames = list(NULL, colnames(s@fixef))),
                  ranef = lapply(s@ranef, function(a)
                      array(apply(a, c(2, 3), mean), c(1, dim(a)[2:3]),
                            c(list(NULL), dimnames(a)[2:3]))),
                  sigma = mean(s@sigma))
    Dbar <- mean(D(fitted(s, fit), s@sigma))
    Dhat <- D(fitted(s.mean, fit), s.mean@sigma)
    pD <- Dbar - Dhat
    return(c(DIC = Dhat + 2 * pD, pD = pD))
}




#
#extractAIC.mer <- function(fit,...){
##     REML <- fit@dims["REML"]
##    llik <- logLik(fit, REML)
##    AIC <- AIC(llik)
##    names(AIC) <- "AIC"
##    return(AIC)
#    L <- logLik(refitML(fit))
#    edf <- attr(L,"df")
#    out <- c(edf,-2*L + k*edf)
#    return(out)
#} 
