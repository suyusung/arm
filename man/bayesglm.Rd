\name{bayesglm}
\docType{class}
% Classes
\alias{bayesglm-class}
%\alias{bayesglm.h-class}
% functions
\alias{bayesglm}
\alias{bayesglm.fit}
% display methods
\alias{print,bayesglm-method}
%\alias{print,bayesglm.h-method}
\alias{show,bayesglm-method}
%\alias{show,bayesglm.h-method}
\alias{predict.bayesglm}

\title{Bayesian generalized linear models.}

\description{Bayesian functions for generalized linear modeling
        with independent normal, t, or Cauchy prior distribution
        for the coefficients.}

\usage{
bayesglm (formula, family = gaussian, data,
    weights, subset, na.action,
    start = NULL, etastart, mustart,
    offset, control = list(...),
    model = TRUE, method = "glm.fit",
    x = FALSE, y = TRUE, contrasts = NULL,
    drop.unused.levels = TRUE,
    prior.mean = 0,
    prior.scale = NULL,
    prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale=1e-12,
    scaled = TRUE, keep.order=TRUE,
    drop.baseline=TRUE,
    maxit=100, 
    print.unnormalized.log.posterior=FALSE,
    Warning=TRUE,...)

bayesglm.fit (x, y, weights = rep(1, nobs),
    start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), family = gaussian(),
    control = list(), intercept = TRUE,
    prior.mean = 0,
    prior.scale = NULL,
    prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale=1e-12, scaled = TRUE,
    print.unnormalized.log.posterior=FALSE, Warning=TRUE)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{family}{a description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function.  (See \code{\link{family}} for details of
    family functions.)}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glm} is called.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{starting values for the parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  This should be \code{NULL} or a numeric vector of
    length either one or equal to the number of cases.
    One or more \code{\link{offset}} terms can be included in the
    formula instead or as well, and if both are specified their sum is
    used.  See \code{\link{model.offset}}.}

  \item{control}{a list of parameters for controlling the fitting
    process.  See the documentation for \code{\link{glm.control}}
    for details.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{method}{the method to be used in fitting the model.
    The default method \code{"glm.fit"} uses iteratively reweighted
    least squares (IWLS).  The only current alternative is
    \code{"model.frame"} which returns the model frame and does no fitting.}

  \item{x, y}{For \code{glm}:
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.

    For \code{glm.fit}: \code{x} is a design matrix of dimension \code{n
      * p}, and \code{y} is a vector of observations of length \code{n}.
  }

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{drop.unused.levels}{default TRUE, if FALSE, it interpolates the
        intermediate values if the data have integer levels.}

  \item{intercept}{logical. Should an intercept be included in the
    \emph{null} model?}

  \item{prior.mean}{prior mean for the coefficients: default is 0. Can be a vector
    of length equal to the number of predictors
    (not counting the intercept, if any). If it is a scalar, it is
    expanded to the length of this vector.}

  \item{prior.scale}{prior scale for the coefficients: default is NULL; if is NULL, for
    a logit model, prior.scale is 2.5; for a probit model, prior scale is 2.5*1.6.
    Can be a vector of length equal to the number of predictors
    (not counting the intercept, if any). If it is a scalar, it is
    expanded to the length of this vector.}

  \item{prior.df}{prior degrees of freedom for the coefficients.
    For t distribution: default is 1 (Cauchy). Set to Inf to
    get normal prior distributions. Can be a vector of length equal to
    the number of predictors (not counting the intercept, if any).
    If it is a scalar, it is expanded to the length of this vector.}

  \item{prior.mean.for.intercept}{prior mean for the intercept: default
    is 0. See \sQuote{Details}.}

  \item{prior.scale.for.intercept}{prior scale for the intercept: default is NULL; for
    a logit model, prior scale for intercept is 10;
    for probit model, prior scale for intercept is rescaled as 10*1.6.}

  \item{prior.df.for.intercept}{prior degrees of freedom for the intercept: default is 1.}

  \item{min.prior.scale}{Minimum prior scale for the coefficients: default is 1e-12.}

  \item{scaled}{scaled=TRUE, the scales for the prior distributions of the coefficients
  are determined as follows: For a predictor with only one value,
  we just use prior.scale. For a predictor with two values,
  we use prior.scale/range(x). For a predictor with more than two values,
  we use prior.scale/(2*sd(x)). If the response is Gaussian, prior.scale is also
  multiplied by 2 * sd(y). Default is TRUE}

  \item{keep.order}{a logical value indicating whether the terms should
    keep their positions. If \code{FALSE} the terms are reordered so
    that main effects come first, followed by the interactions,
    all second-order, all third-order and so on.  Effects of a given
    order are kept in the order specified. Default is TRUE.}

  \item{drop.baseline}{Drop the base level of categorical x's, default is TRUE.}

  \item{maxit}{integer giving the maximal number of IWLS iterations, default is 100. This can also be controlled by \code{control}.}
  \item{print.unnormalized.log.posterior}{display the unnormalized log posterior likelihood for bayesglm, default=FALSE}
  \item{Warning}{default is TRUE, which will show the error messages of not convergence and separation.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The program is a simple alteration of \code{glm()} that uses an approximate EM
algorithm to update the betas at each step using an augmented regression
to represent the prior information.

We use Student-t prior distributions for the coefficients. The prior
distribution for the constant term is set so it applies to the value
when all predictors are set to their mean values.

If scaled=TRUE, the scales for the prior distributions of the
coefficients are determined as follows: For a predictor with only one
value, we just use prior.scale. For a predictor with two values, we use
prior.scale/range(x). For a predictor with more than two values, we use
prior.scale/(2*sd(x)).

We include all the \code{glm()} arguments but we haven't tested that all the
options (e.g., \code{offsets}, \code{contrasts},
\code{deviance} for the null model) all work.

The new arguments here are: \code{prior.mean}, \code{prior.scale},
\code{prior.scale.for.intercept}, \code{prior.df}, \code{prior.df.for.intercept}and
\code{scaled}.
}

\value{See \code{\link[stats]{glm}} for details.

\item{prior.mean}{prior means for the coefficients and the intercept.}
\item{prior.scale}{prior scales for the coefficients}
\item{prior.df}{prior dfs for the coefficients.}
\item{prior.scale.for.intercept}{prior scale for the intercept}
\item{prior.df.for.intercept}{prior df for the intercept}

}

\references{Andrew Gelman, Aleks Jakulin, Maria Grazia Pittau and Yu-Sung Su. (2009).
            \dQuote{A Weakly Informative Default Prior Distribution For
            Logistic And Other Regression Models.}
            \emph{The Annals of Applied Statistics} 2 (4): 1360--1383.
            \url{http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf}
}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn};
        Daniel Lee \email{bearlee@alum.mit.edu};
        Aleks Jakulin \email{Jakulin@stat.columbia.edu}}

\seealso{
    \code{\link[stats]{glm}},
    \code{\link{bayespolr}}
}

\examples{
  n <- 100
  x1 <- rnorm (n)
  x2 <- rbinom (n, 1, .5)
  b0 <- 1
  b1 <- 1.5
  b2 <- 2
  y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))

  M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M1)

  M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=Inf, prior.df=Inf)
  display (M2)  # just a test:  this should be identical to classical logit

  M3 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))
    # default Cauchy prior with scale 2.5
  display (M3)

  M4 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.df=1)
    # Same as M3, explicitly specifying Cauchy prior with scale 2.5
  display (M4)

  M5 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.df=7)   # t_7 prior with scale 2.5
  display (M5)

  M6 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.df=Inf)  # normal prior with scale 2.5
  display (M6)

# Create separation:  set y=1 whenever x2=1
# Now it should blow up without the prior!

  y <- ifelse (x2==1, 1, y)

  M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M1)

  M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=Inf, prior.scale.for.intercept=Inf) # Same as M1
  display (M2)

  M3 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M3)

  M4 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.scale.for.intercept=10)  # Same as M3
  display (M4)

  M5 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.df=7)
  display (M5)

  M6 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"),
    prior.scale=2.5, prior.df=Inf)
  display (M6)

  # bayesglm with gaussian family (bayes lm)
  sigma <- 5
  y2 <- rnorm (n, b0+b1*x1+b2*x2, sigma)
  M7 <- bayesglm (y2 ~ x1 + x2, prior.scale=Inf, prior.df=Inf)
  display (M7)


  # bayesglm with categorical variables
  z1 <- trunc(runif(n, 4, 9))
  levels(factor(z1))
  z2 <- trunc(runif(n, 15, 19))
  levels(factor(z2))

  ## drop the base level (R default)
  M8 <- bayesglm (y ~ x1 + factor(z1) + factor(z2),
    family=binomial(link="logit"), prior.scale=2.5, prior.df=Inf)
  display (M8)

  ## keep all levels with the intercept, keep the variable order
  M9 <- bayesglm (y ~ x1 + x1:x2 + factor(z1) + x2 + factor(z2),
    family=binomial(link="logit"),
    prior.mean=rep(0,12),
    prior.scale=rep(2.5,12),
    prior.df=rep(Inf,12),
    prior.mean.for.intercept=0,
    prior.scale.for.intercept=10,
    prior.df.for.intercept=1,
    drop.baseline=FALSE, keep.order=TRUE)
  display (M9)

  ## keep all levels without the intercept
  M10 <- bayesglm (y ~ x1 + factor(z1) + x1:x2 + factor(z2)-1,
    family=binomial(link="logit"),
    prior.mean=rep(0,11),
    prior.scale=rep(2.5,11),
    prior.df=rep(Inf,11),
    drop.baseline=FALSE)
  display (M10)

}

\keyword{models}
\keyword{methods}
\keyword{regression}
