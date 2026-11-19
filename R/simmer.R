# simulations of sigma, fixef, and ranef drawn from a posterior
# under a flat prior and conditioned on estimate of ranef covar
setMethod("sim", signature(object = "merMod"),
          function(object, n.sims=100)
{
  # information is conditional on hyperparameters
  # information is of [ranef, fixef]
  getInverseInformationLeftFactor <- function(regression) {
    Lz  <- getME(regression, "L");
    Rzx <- getME(regression, "RZX");
    Rx  <- getME(regression, "RX");

    # upper left, lower right, and lower left blocks of left-factor of inverse
    #solveFunc <- getMethod("solve", signature(a = "CHMfactor", b = "diagonalMatrix"));
    #Rz.inv  <- t(solveFunc(Lz, Diagonal(Lz@Dim[1]), "L"));
    Rz.inv  <- t(solve(Lz, Diagonal(Lz@Dim[1]), system = "L"))
    Rx.inv  <- solve(Rx)
    Rzx.inv <- -Rz.inv %*% Rzx %*% Rx.inv

    # Undo the fill-reducing permutation (rows), via a sparse permutation matrix
    perm <- Lz@perm + 1L
    Pt <- t(sparseMatrix(i = seq_along(perm), j = perm, x = 1))
    Rz.inv  <- Pt %*% Rz.inv
    Rzx.inv <- Pt %*% Rzx.inv

    # this is me figuring some stuff out. new lmer doesn't permute Zt apparently
    #
    #Lz.tmp <- as(Lz, "sparseMatrix");
    #P.chol <- as(Lz@perm + 1, "pMatrix");
    #Zt <- getME(regression, "Zt");
    #W <- Diagonal(numObs, regression@resp$sqrtXwt);
    ## P.ranef <- getRanefPerm(regression);
    #Lambdat <- getME(regression, "Lambdat") # t(P.ranef) %*% getME(regression, "Lambdat") %*% P.ranef;
    #A <- Lambdat %*% Zt;
    #C <- A %*% W;
    #L.hyp <- Cholesky(tcrossprod(P.chol %*% C), Imult = 1, LDL = FALSE, perm = FALSE);
    #L.hyp@perm <- Lz@perm;
    #L.hyp@type[1] <- 2L;
    #browser();

    #P.ranef <- getRanefPerm(model);
    #Lambda <- P.ranef %*% getRanefChol(model) %*% t(P.ranef);
    Lambda <- getME(regression, "Lambda");

    return(list(ul = Lambda %*% Rz.inv, ur = Lambda %*% Rzx.inv, lr = Rx.inv));
  }

  regression <- object;
  devcomp <- getME(regression, "devcomp");
  dims <- devcomp$dims;

  if (dims[["NLMM"]] != 0L) stop("sim not yet implemented for nlmms");

  numObs    <- dims[["n"]];
  numRanef  <- dims[["q"]];
  numFixef  <- dims[["p"]];
  numLevels <- dims[["reTrms"]];

  isLinearMixedModel <- dims[["GLMM"]] == 0L && dims[["NLMM"]] == 0L;
  numEffects <- numRanef + numFixef;
  numDoF <- numObs - numFixef;

  # pertain to simulations that we do all as a single vector
  ranefRange <- 1:numRanef;
  fixefRange <- numRanef + 1:numFixef;

  # stuff used to rearrange ranef into usable form
  groupsPerUniqueFactor <- lapply(regression@flist, levels);
  factorPerLevel <- attr(regression@flist, "assign");

  coefficientNamesPerLevel <- regression@cnms;
  numCoefficientsPerLevel <- as.numeric(sapply(coefficientNamesPerLevel, length));
  numGroupsPerLevel <- as.numeric(sapply(groupsPerUniqueFactor[factorPerLevel], length));
  numRanefsPerLevel <- numCoefficientsPerLevel * numGroupsPerLevel;

  # storage for sims
  simulatedSD    <- if (isLinearMixedModel) { rep(NA, n.sims); } else { NA };
  simulatedRanef <- vector("list", numLevels);
  names(simulatedRanef) <- names(regression@cnms);
  for (i in 1:numLevels) {
    simulatedRanef[[i]] <- array(NA, c(n.sims, numGroupsPerLevel[i], numCoefficientsPerLevel[i]),
                                 list(NULL, groupsPerUniqueFactor[[factorPerLevel[i]]], coefficientNamesPerLevel[[i]]));
  }

  simulatedFixef <- matrix(NA, n.sims, numFixef,
                           dimnames = list(NULL, names(fixef(regression))));


  # "b" are the rotated random effects, i.e. what ranef() returns in
  # a rearranged format.
  effectsMean <- c(getME(regression, "b")@x, getME(regression, "beta"));
  effectsCovLeftFactor <- getInverseInformationLeftFactor(regression);

  # The heavy work below is batched, but the DRAWS stay one-at-a-time:
  # rgamma() is rejection-based and consumes a variable number of uniforms,
  # so drawing all gammas up front would shift the shared Mersenne-Twister
  # stream relative to the old interleaved loop (rgamma(1) then numEffects
  # rnorms, per draw).  Keeping the loop preserves bit-identical draws --
  # what display() documents for its fixed-seed DIC -- while the products
  # run as one BLAS matmul and each ranef array is built once with array()
  # instead of arr[i,,] <- copying the whole array on every draw (quadratic
  # in n.sims, and the heavy GC that pushed the display() example over R
  # CMD check's 5s NOTE threshold).  Same level-major, coefficient-minor
  # layout as before.
  sphericalEffects <- matrix(0, n.sims, numEffects)
  if (isLinearMixedModel) {
    for (i in seq_len(n.sims)) {
      simulatedSD[i] <- sqrt(1 / rgamma(1, 0.5 * numDoF,
                                        0.5 * devcomp$cmp[["pwrss"]]))
      sphericalEffects[i, ] <- rnorm(numEffects, 0, simulatedSD[i])
    }
  } else {
    for (i in seq_len(n.sims)) {
      sphericalEffects[i, ] <- rnorm(numEffects)
    }
  }

  leftFactor <- matrix(0, numEffects, numEffects);
  leftFactor[ranefRange, ranefRange] <- as.matrix(effectsCovLeftFactor$ul);
  leftFactor[ranefRange, fixefRange] <- as.matrix(effectsCovLeftFactor$ur);
  leftFactor[fixefRange, fixefRange] <- as.matrix(effectsCovLeftFactor$lr);

  simulatedEffects <- sphericalEffects %*% t(leftFactor) +
                        matrix(effectsMean, n.sims, numEffects, byrow = TRUE);

  simulatedFixef[] <- simulatedEffects[, fixefRange, drop = FALSE];

  ranefStart <- cumsum(c(1L, numRanefsPerLevel[-numLevels]));
  for (k in 1:numLevels) {
    cols <- (ranefStart[k] - 1) + seq_len(numRanefsPerLevel[k]);
    slice <- simulatedEffects[, cols, drop = FALSE];
    # slice rows are draws, columns are (group, coefficient) level-major,
    # coefficient-minor -- the same layout matrix(v, ncol = n_c, byrow = TRUE)
    # gave per draw; an array built straight from the slice would come out
    # coefficient-major, so build (draws, coefficient, group) and transpose.
    simulatedRanef[[k]] <- aperm(
      array(slice, c(n.sims, numCoefficientsPerLevel[k], numGroupsPerLevel[k])),
      c(1, 3, 2));
    dimnames(simulatedRanef[[k]]) <-
      list(NULL, groupsPerUniqueFactor[[factorPerLevel[k]]],
           coefficientNamesPerLevel[[k]]);
  }

  ans <- new("sim.merMod",
             "fixef" = simulatedFixef,
             "ranef" = simulatedRanef,
             "sigma" = simulatedSD);
  return(ans);
});
