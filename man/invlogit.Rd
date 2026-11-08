\name{invlogit}
\alias{invlogit}
\alias{logit}

\title{Logit and Inverse-Logit Functions}

\description{
  The inverse-logit (logistic) transformation, which maps continuous values
  onto the interval (0, 1), and the logit transformation, which maps
  probabilities back onto the real line.
}
\usage{
logit(x)
invlogit(x)
}

\arguments{
  \item{x}{A vector of continuous values}
}

\details{
        The inverse-logit function is defined as
        \eqn{\mathrm{logit}^{-1}(x) = e^x/(1+e^x)}{logit^-1(x) = e^x/(1+e^x)},
        which transforms continuous values to the range (0, 1).  This is
        necessary because probabilities must lie between 0 and 1, and it maps
        from the linear predictor to the probability scale.
}
\value{
  For \code{invlogit}, a vector of probabilities.  For \code{logit}, a vector
  of log-odds.
}
\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.}

\author{Andrew Gelman \email{gelman@stat.columbia.edu}, 
        M.Grazia Pittau \email{grazia@stat.columbia.edu}
        }


\examples{
 n <- 100
 x1 <- rnorm (n)
 x2 <- rbinom (n, 1, .5)
 b0 <- 1
 b1 <- 1.5
 b2 <- 2
 Inv.logit <- invlogit(b0+b1*x1+b2*x2)
 plot(b0+b1*x1+b2*x2, Inv.logit)
}

\keyword{models}
 
