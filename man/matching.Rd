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

\details{The function matches each treatment unit to the control unit with the
    closest propensity score.  With \code{replace = FALSE} the treatment units are
    matched in turn and a control unit already chosen is not available again.
    With \code{replace = TRUE} a control unit can be matched to several treatment
    units, and ties are broken at random, so set a seed to reproduce the result.}

\value{
A list with components
\describe{
  \item{match.ind}{for \code{replace = FALSE}, a vector of length
    \code{length(z)} giving the index of the unit that each unit is matched to;
    0 means matched to nothing.  For \code{replace = TRUE}, a vector of length
    \code{2 * sum(z)}: the indices of the treated units followed by the indices
    of the control units they are matched to, in the same order.}
  \item{cnts}{the number of times each unit will be used in any subsequent
    analysis (1 for each treated unit, and, for each control unit, the number
    of treated units it is matched to).}
  \item{matches}{for \code{replace = TRUE} only, a 0/1 matrix with one row per
    treated unit and one column per control unit (each in the order they appear
    in \code{z}), recording which control each treated unit was matched to.}
  \item{pairs}{for \code{replace = FALSE} only, a vector of length
    \code{length(z)} giving the pair number of each matched unit; \code{NA}
    for unmatched control units.}
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
