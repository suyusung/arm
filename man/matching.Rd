\name{matching}
\alias{matching}

\title{Single Nearest Neighbor Matching}

\description{
Function for matching treatment units to control units using a propensity score.
}

\usage{
matching(z, score, replace=FALSE)
}

\arguments{
  \item{z}{vector of indicators for treatment or control.}
  
  \item{score}{vector of the propensity scores in the same order as z.}
  
  \item{replace}{whether the control units could be reused for matching, 
      default is \code{FALSE}.} 
}

\details{The function matches each treatment unit in turn to the control unit
    (not previously chosen) with the closest propensity score.}

\value{
A list with components
\describe{
  \item{match.ind}{a vector of indices giving the unit that each unit is
    matched to; 0 means matched to nothing.}
  \item{cnts}{the number of times each unit will be used in any subsequent
    analysis (1 for each treated unit, and, for each control unit, the number
    of treated units it is matched to).}
  \item{matches}{for \code{replace = FALSE}, a matrix recording which treated
    observations were matched to which controls.}
  \item{pairs}{for \code{replace = TRUE}, an indicator for each pair.}
}
}
 
\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.}

\author{Jennifer Hill \email{jh1030@columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}

\seealso{\code{\link{balance}}
        }

\examples{
# matching first
data(lalonde)
attach(lalonde)
fit <- glm(treat ~ re74 + re75 + age + factor(educ) + 
            black + hisp + married + nodegr + u74 + u75, 
            family=binomial(link="logit"))
pscores <- predict(fit, type="response")
matches <- matching(z=lalonde$treat, score=pscores)
matched <- matches$cnts

# balance check!
b.stats <- balance(lalonde, treat, matched)
print(b.stats)
plot(b.stats)
}

\keyword{models}
\keyword{methods}
