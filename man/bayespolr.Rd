\name{bayespolr}
\docType{class}
% Classes
\alias{bayespolr-class}
% functions
\alias{bayespolr}
% display methods
\alias{print,bayespolr-method}
\alias{show,bayespolr-method}

\title{Bayesian Ordered Logistic or Probit Regression}

\description{Bayesian functions for ordered logistic or probit modeling
        with independent normal, t, or Cauchy prior distribution
        for the coefficients.
}
\usage{
bayespolr(formula, data, weights, start,  ...,
    subset, na.action, contrasts = NULL,
    Hess = TRUE, model = TRUE,
    method = c("logistic", "probit", "cloglog", "cauchit"),
    drop.unused.levels=TRUE,
    prior.mean = 0,
    prior.scale = 2.5,
    prior.df = 1,
    prior.counts.for.bins = NULL,
    min.prior.scale=1e-12,
    scaled = TRUE,
    maxit = 100,
    print.unnormalized.log.posterior = FALSE)
}

\arguments{
  \item{formula}{a formula expression as for regression models, of the form
          \code{response ~ predictors}. The response should be a factor
          (preferably an ordered factor), which will be interpreted as
          an ordinal response, with levels ordered as in the factor.  A
          proportional odds model will be fitted.  The model must have
          an intercept: attempts to remove one will lead to a warning
          and be ignored. An offset may be used.  See the documentation
          of \code{formula} for other details.}

  \item{data}{an optional data frame in which to interpret the variables
          occurring in \code{formula}.}

  \item{weights}{optional case weights in fitting.  Default to 1.}

  \item{start}{initial values for the parameters.  This is in the format
          \code{c(coefficients, zeta)}}

  \item{\dots}{additional arguments to be passed to \code{optim}, most often a
          \code{control} argument.}

  \item{subset}{expression saying which subset of the rows of the data should
           be used in the fit.  All observations are included by
          default.}

  \item{na.action}{a function to filter missing data.}

  \item{contrasts}{a list of contrasts to be used for some or all of the
          factors appearing as variables in the model formula.}

  \item{Hess}{logical for whether the Hessian (the observed information
          matrix) should be returned.}

  \item{model}{logical for whether the model matrix should be returned.}

  \item{method}{logistic or probit or complementary log-log or cauchit
          (corresponding to a Cauchy latent variable and only available
          in R >= 2.1.0).}

  \item{drop.unused.levels}{default \code{TRUE}, if \code{FALSE}, it interpolates the
        intermediate values if the data have integer levels.}

  \item{prior.mean}{prior mean for the coefficients: default is 0.
    Can be a vector of length equal to the number of predictors
    (not counting the intercepts). If it is a scalar, it is
    expanded to the length of this vector.}

  \item{prior.scale}{prior scale for the coefficients: default is 2.5.
    Can be a vector of length equal to the number of predictors
    (not counting the intercepts). If it is a scalar, it is expanded
    to the length of this vector.}

  \item{prior.df}{for t distribution: default is 1 (Cauchy).
    Set to \code{Inf} to get normal prior distributions. Can
    be a vector of length equal to the number of predictors (not counting
    the intercepts). If it is a scalar, it is expanded to the length of this
    vector.}

  \item{prior.counts.for.bins}{default is \code{NULL}, which will augment the data by
    giving each cut point a \code{1/levels(y)}. To use a noninformative prior, assign
    prior.counts.for.bins = 0. If it is a scalar, it is expanded to the number
    of levels of y.}

  \item{min.prior.scale}{Minimum prior scale for the coefficients: default is 1e-12.}


  \item{scaled}{if \code{scaled = TRUE}, then the prior distribution is rescaled.
    Can be a vector of length equal to the number of cutpoints
    (intercepts). If it is a scalar, it is expanded to the length of
    this vector.}

  \item{maxit}{integer giving the maximal number of IWLS iterations, default is 100. This can also be controlled by \code{control}.}

  \item{print.unnormalized.log.posterior}{display the unnormalized log posterior
    likelihood for bayesglm fit, default=\code{FALSE}}
}
\details{
The program is a simple alteration of \code{\link[MASS]{polr}} in
\code{VR} version 7.2-31 that augments the
loglikelihood with the log of the t prior distributions for the
coefficients.

We use Student-t prior distributions for the coefficients. The prior
distributions for the intercepts (the cutpoints) are set so they apply
to the value when all predictors are set to their mean values.

If scaled=TRUE, the scales for the prior distributions of the
coefficients are determined as follows: For a predictor with only one
value, we just use \code{prior.scale}. For a predictor with two
values, we use prior.scale/range(x).
For a predictor with more than two values, we use prior.scale/(2*sd(x)).
}

\value{
See \code{polr} for details.

\item{prior.mean}{prior means for the cofficients.}
\item{prior.scale}{prior scales for the cofficients.}
\item{prior.df}{prior dfs for the cofficients.}
\item{prior.counts.for.bins}{prior counts for the cutpoints.}
}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn};
        Maria Grazia Pittau \email{grazia@stat.columbia.edu}
}
\seealso{\code{\link{bayesglm}},
        \code{\link[MASS]{polr}}
}

\examples{
    M1 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    display (M1)

    M2 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
        prior.scale=Inf, prior.df=Inf) # Same as M1
    display (M2)

    M3 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    display (M3)

    M4 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
        prior.scale=2.5, prior.df=1)  # Same as M3
    display (M4)

    M5 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
        prior.scale=2.5, prior.df=7)
    display (M5)

    M6 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
        prior.scale=2.5, prior.df=Inf)
    display (M6)

    # Assign priors
    M7 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
        prior.mean=rep(0,6), prior.scale=rep(2.5,6), prior.df=c(1,1,1,7,7,7))
    display (M7)


    #### Another example
    y <- factor (rep (1:10,1:10))
    x <- rnorm (length(y))
    x <- x - mean(x)

    M8 <- polr (y ~ x)
    display (M8)

    M9 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=0)
    display (M9) # same as M1

    M10 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=10000)
    display (M10)


    #### Another example

    y <- factor (rep (1:3,1:3))
    x <- rnorm (length(y))
    x <- x - mean(x)

    M11 <- polr (y ~ x)
    display (M11)

    M12 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=0)
    display (M12) # same as M1

    M13 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=1)
    display (M13)

    M14 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=10)
    display (M14)

}

\keyword{models}
\keyword{methods}
\keyword{regression}
