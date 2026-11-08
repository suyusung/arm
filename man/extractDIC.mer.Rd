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
\method{extractDIC}{merMod}(fit,\dots)

%\method{extractAIC}{merMod}(fit,\dots)
}

\arguments{
  \item{fit}{a fitted \code{merMod} model, usually the result of \code{lmer}.}
  \item{\dots}{further arguments (currently unused).}
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
