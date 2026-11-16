\name{extractDIC}
%\docType{genericFunction}
\alias{extractDIC}
\alias{extractDIC.merMod}
\alias{extractAIC.merMod}

\title{
Extract AIC and DIC from a \sQuote{mer} model
}

\description{
 Computes the (generalized) Akaike *A*n *I*nformation *C*riterion 
 and *D*eviance *I*nformation *C*riterion for a mer model.
}
\usage{
extractDIC(fit,\dots)
\method{extractDIC}{merMod}(fit, n.sims = 1000, \dots)

%\method{extractAIC}{merMod}(fit,\dots)
}

\arguments{
  \item{fit}{a fitted \code{merMod} model, usually the result of \code{lmer}.}
  \item{n.sims}{number of \code{\link{sim}} draws used to compute DIC.}
  \item{\dots}{further arguments (currently unused).}
}

\details{
  DIC is the conditional DIC that BUGS reports.  \code{sim(fit, n.sims)}
  draws the fixed effects, the random effects and, for \code{lmer}, the
  residual sd, with the variance parameters held at their estimates.  With
  \eqn{D(\theta) = -2 \log p(y \mid \beta, b, \sigma)}{D(theta) = -2 log p(y | beta, b, sigma)},
  \eqn{\bar{D}}{Dbar} is the mean of \eqn{D} over the draws,
  \eqn{\hat{D}}{Dhat} is \eqn{D} at the means of the draws,
  \eqn{p_D = \bar{D} - \hat{D}}{pD = Dbar - Dhat}, and
  \eqn{DIC = \hat{D} + 2 p_D}{DIC = Dhat + 2 pD}.  The value is random;
  call \code{set.seed} first to reproduce it.  Only \code{lmer} fits and
  binomial or poisson \code{glmer} fits are supported; for other families
  the result is \code{NA} with a warning.
}

\value{
  The DIC, named \code{"DIC"}, with the effective number of parameters
  \eqn{p_D}{pD} as attribute \code{"pD"}.
}


\author{
  Andrew Gelman \email{gelman@stat.columbia.edu};
  Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}


\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
extractAIC(fm1)
extractDIC(fm1)
}

\keyword{manip}
\keyword{methods}
