\name{matching}
\alias{matching}

\title{Single Nearest Neighborhood Matching}

\description{
Function for processing matching with propensity score 
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

\details{Function for matching each treatment unit in turn 
    the control unit (not previously chosen) with the 
    closest propensity score }

\value{
The function returns a vector of indices that the corresponding 
unit is matched to. 0 means matched to nothing.}
 
\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.}

\author{Jeniffer Hill \email{jh1030@columbia.edu};
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
