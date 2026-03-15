setMethod("sim", signature(object = "lm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1,drop = FALSE]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, rownames(beta.hat))
    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }

    ans <- new("sim",
                coef = beta,
                sigma = sigma)
    return (ans)
    }
)



setMethod("sim", signature(object = "glm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object, correlation=TRUE, dispersion = object$dispersion)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    beta.hat <- coef[,1,drop=FALSE]
    sd.beta <- coef[,2,drop=FALSE]
    corr.beta <- summ$corr
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
    #beta <- array (NA, c(n.sims,k))
#    dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])
#    for (s in 1:n.sims){
#      beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta)
#    }
    beta <- MASS::mvrnorm (n.sims, beta.hat, V.beta)
    # Added by Masanao
    beta2 <- array (0, c(n.sims,length(coefficients(object))))
    dimnames(beta2) <- list (NULL, names(coefficients(object)))
    beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta
    # Added by Masanao
    sigma <- rep (sqrt(summ$dispersion), n.sims)

    ans <- new("sim",
                coef = beta2,
                sigma = sigma)
    return(ans)
    }
)


setMethod("sim", signature(object = "polr"),
  function(object, n.sims = 100) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for this function. Please install it.")
    }

    # Extract coefficients and thresholds
    coefs <- coef(object)
    zeta <- object$zeta
    
    # Number of regression coefficients
    k <- length(coefs)
    
    # Variance-covariance matrix of all parameters (coefficients + thresholds)
    Sigma <- vcov(object)
    
    # Draw parameters from multivariate normal distribution
    parameters <- MASS::mvrnorm(n = n.sims, mu = c(coefs, zeta), Sigma = Sigma)
    
    # If only one simulation, ensure 'parameters' has dimension (1, nparams)
    if (n.sims == 1) {
      parameters <- matrix(parameters, nrow = 1)
    }
    
    # Create new "sim.polr" object with coefficients and thresholds separated
    ans <- new("sim.polr",
               coef = parameters[, 1:k, drop = FALSE],
               zeta = parameters[, (k + 1):ncol(parameters), drop = FALSE]
    )
    
    return(ans)
  }
)

sim.coxph <- function(object, n.sims = 100) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for this function to work. Please install it.")
  }
  
  # Check object class
  if (!inherits(object, "coxph")) {
    stop("Input object is not a coxph model.")
  }
  
  # Extract coef estimates
  beta.hat <- coef(object)
  if (is.null(beta.hat)) stop("Could not extract coefficients from the coxph object.")

  # Extract variance-covariance matrix of coefficient estimates
  V.beta <- tryCatch(vcov(object),
                     error = function(e) stop("Could not extract variance-covariance matrix from coxph object."))

  k <- length(beta.hat)
  beta.sim <- matrix(NA, nrow = n.sims, ncol = k)
  colnames(beta.sim) <- names(beta.hat)
  
  for (i in seq_len(n.sims)) {
    beta.sim[i, ] <- MASS::mvrnorm(1, mu = beta.hat, Sigma = V.beta)
  }
  
  # For consistency with sim.plm, return a list of simulated coefs
  ans <- list(coef = beta.sim)
  class(ans) <- "sim.coxph"
  
  return(ans)
}

sim.plm <- function(object, n.sims = 100) {
  # Load required package
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for this function to work. Please install it.")
  }
  
  # Extract model frame (data used in the fitted model)
  mf <- tryCatch(model.frame(object),
                 error = function(e) stop("Cannot extract model frame from the object."))
  
  # Extract coefficients and their standard errors
  summ <- tryCatch(summary(object),
                   error = function(e) stop("Cannot compute summary for the model object."))
  
  # Try to get coefficients table reliably
  coefmat <- tryCatch(
    {
      # Most models have coef matrix under summ$coefficients or summ$coef
      if (!is.null(summ$coefficients)) {
        summ$coefficients[, 1:2, drop = FALSE]
      } else if (!is.null(summ$coef)) {
        summ$coef[, 1:2, drop = FALSE]
      } else {
        stop("Coefficient matrix not found in summary(object).")
      }
    },
    error = function(e) stop("Error extracting coefficient matrix from summary: ", e$message)
  )
  dimnames(coefmat)[[2]] <- c("coef.est", "coef.sd")
  
  # Number of observations (rows) and parameters
  n <- nrow(mf)
  k <- nrow(coefmat)
  
  # Estimate residual standard deviation
  # Use deviance or residual variance if available, otherwise fallback
  sigma.hat <- tryCatch({
    dev <- deviance(object)
    if(is.null(dev)) stop("deviance() returned NULL")
    sqrt(dev / (n - k))
  }, error = function(e) {
    # fallback: try sigma method
    if("sigma" %in% methods(class = class(object))) {
      sigma(object)
    } else {
      # fallback: use residuals to estimate sigma
      res <- residuals(object)
      if (is.null(res)) stop("Cannot estimate residual standard deviation")
      sqrt(sum(res^2) / (n - k))
    }
  })
  
  # Unscale covariance matrix of coefficients if necessary
  # Sometimes vcov(object) already includes sigma^2, so unscale it:
  Vbeta_raw <- tryCatch(vcov(object),
                        error = function(e) stop("Failed to get vcov of the object"))
  
  # Check if Vbeta_raw is scaled by sigma.hat^2 (heuristic)
  # If max diagonal is more than 100 times sigma.hat^2, probably not scaled
  diag_vcov <- diag(Vbeta_raw)
  if(all(diag_vcov > 0) && max(diag_vcov) > 100 * sigma.hat^2) {
    # Assume vcov returns unscaled covariance of coefficients: multiply by sigma^2
    V.beta <- Vbeta_raw * sigma.hat^2
  } else {
    # Assume vcov gives scaled covariance (multiplied by sigma^2), so unscale it
    V.beta <- Vbeta_raw
  }
  
  # Preallocate outputs
  sigma <- numeric(n.sims)
  beta <- matrix(NA, nrow = n.sims, ncol = k)
  colnames(beta) <- rownames(coefmat)
  
  for (s in seq_len(n.sims)) {
    # Draw sigma from scaled inverse-chi-squared distribution
    sigma[s] <- sigma.hat * sqrt( (n - k) / rchisq(1, df = n - k) )
    
    # Draw beta conditional on sigma
    beta[s, ] <- MASS::mvrnorm(1, mu = coefmat[, "coef.est"], Sigma = V.beta * sigma[s]^2)
  }
  
  ans <- list(coef = beta, sigma = sigma)
  class(ans) <- "sim"
  return(ans)
}



#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100)
#    {
#    #object <- summary(object)
##    if (lapply(object@bVar,sum)<=0|sum(unlist(lapply(object@bVar, is.na)))>0){
##        object@call$control <- list(usePQL=TRUE)
##        object <- lmer(object@call$formula)
#    #}
#    #sc <- attr (VarCorr (object), "sc")
#    # simulate unmodeled coefficients
#
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "unmodeled"
#    }
#    # simulate coefficients within groups
#    #coef <- ranef (object)
#    #estimate <- ranef(object, postVar=TRUE)
#    #vars <- object@bVar
#    #beta.bygroup <- vars
#
#    sc <- attr (VarCorr (object), "sc")
#    coef <- ranef(object, postVar=TRUE)
#    beta.bygroup <- c(coef)
#    n.groupings <- length (coef)
#    for (m in 1:n.groupings){
#      #vars.m <- vars[[m]]
#      vars.m <- attr (coef[[m]], "postVar")
#      K <- dim(vars.m)[1]
#      J <- dim(vars.m)[3]
#      beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#      bhat <- coef[[m]]
#      for (j in 1:J){
#        V.beta <- untriangle(vars.m[,,j])#*sc^2
#        beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#      }
#      dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#    }
#    betas <- c (beta.unmodeled, beta.bygroup)
#    return (betas)
#    }
#)

#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100, ranef=TRUE)
#    {
#    # simulate unmodeled coefficients
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "fixef"#"unmodeled"
#      coef <- beta.unmodeled
#    }
#    if(ranef){
#      # simulate coefficients within groups
#      sc <- attr (VarCorr (object), "sc")  # scale
#      #coef <- ranef (object)
#      #estimate <- ranef(object, postVar=TRUE)
#      coef <- ranef(object, postVar=TRUE)
#      beta.bygroup <- coef
#      n.groupings <- length (coef)
#      for (m in 1:n.groupings){
#        bhat <- as.matrix(coef[[m]]) # to suit the use of mvrnorm
#        vars.m <- attr (coef[[m]], "postVar")
#        K <- dim(vars.m)[1]
#        J <- dim(vars.m)[3]
#        beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#        for (j in 1:J){
#          V.beta <- .untriangle(vars.m[,,j])#*sc^2
#          beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#        }
#        dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#      }
#      coef <- c (beta.unmodeled, beta.bygroup)
#      }
#    return (coef)
#    }
#)
