\name{se.coef}
%\docType{genericFunction}
\alias{se.coef}
\alias{se.coef,lm-method}
\alias{se.coef,glm-method}
\alias{se.coef,merMod-method}
\alias{se.fixef}
\alias{se.ranef}

\title{Extract Standard Errors of Model Coefficients}

\description{
   These functions extract standard errors of model coefficients
   from objects returned by modeling functions.
}
\usage{
se.coef (object, \dots)
se.fixef (object)
se.ranef (object)

\S4method{se.coef}{lm}(object)
\S4method{se.coef}{glm}(object)
\S4method{se.coef}{merMod}(object)
}

\arguments{
  \item{object}{object of \code{lm}, \code{glm} and \code{merMod} fit}
  \item{\dots}{other arguments}
}
\value{
  \code{se.coef} gives lists of standard errors for \code{coef},
  \code{se.fixef} gives a vector of standard errors for \code{fixef} and
  \code{se.ranef} gives a list of standard errors for \code{ranef}.
}
\details{
  \code{se.coef} extracts standard errors from objects
    returned by modeling functions.
  \code{se.fixef} extracts standard errors of the fixed effects
    from objects returned by lmer and glmer functions.
  \code{se.ranef} extracts standard errors of the random effects
    from objects returned by lmer and glmer functions.
}

\seealso{
  \code{\link{display}},
  \code{\link{coef}},
  \code{\link{sigma.hat}},
}

\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
        }

\examples{
# Here's a simple example of a model of the form, y = a + bx + error,
# with 10 observations in each of 10 groups, and with both the
# intercept and the slope varying by group.  First we set up the model and data.

   group <- rep(1:10, rep(10,10))
   mu.a <- 0
   sigma.a <- 2
   mu.b <- 3
   sigma.b <- 4
   rho <- 0
   Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b,
                    rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
   sigma.y <- 1
   ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
   a <- ab[,1]
   b <- ab[,2]
#
   x <- rnorm (100)
   y1 <- rnorm (100, a[group] + b[group]*x, sigma.y)
   y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))

#  lm fit
   M1 <- lm (y1 ~ x)
   se.coef (M1)

#  glm fit
   M2 <- glm (y2 ~ x)
   se.coef (M2)

#  lmer fit
   M3 <- lmer (y1 ~ x + (1 + x |group))
   se.coef (M3)
   se.fixef (M3)
   se.ranef (M3)

#  glmer fit
   M4 <- glmer (y2 ~ 1 + (0 + x |group), family=binomial(link="logit"))
   se.coef (M4)
   se.fixef (M4)
   se.ranef (M4)
}
\keyword{manip}
\keyword{methods}
\keyword{models}
