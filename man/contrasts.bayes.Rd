\name{contrast.bayes}
\alias{contr.bayes.ordered}
\alias{contr.bayes.unordered}

\title{Contrast Matrices}

\description{
  Return a matrix of contrasts used in \code{\link{bayesglm}}.
}

\usage{
contr.bayes.unordered(n, base = 1, contrasts = TRUE)
contr.bayes.ordered (n, scores = 1:n, contrasts = TRUE)
}

\arguments{
  \item{n}{a vector of levels for a factor, or the number of levels.}
  \item{base}{an integer specifying which group is considered the baseline
          group. Ignored if \code{contrasts} is \code{FALSE}.}
  \item{contrasts}{a logical indicating whether contrasts should be computed.}
  \item{scores}{the set of values over which orthogonal polynomials are to be
          computed.}
}

\details{
   These functions are adapted from \code{contr.treatment} and \code{contr.poly}
   in \code{\link{stats}} package.  The purpose for these functions are to keep
   the baseline levels of categorical variables and thus to suit the use of 
   \code{\link{bayesglm}}.
   
   \code{contr.bayes.unordered} is equivalent to \code{contr.treatment} whereas
   \code{contr.bayes.ordered} is equivalent to \code{contr.poly}.
   
}


\author{Yu-Sung Su \email{suyusung@tsinghua.edu.cn}}

\seealso{
  \code{\link{C}},
  \code{\link{contr.helmert}},
  \code{\link{contr.poly}},
  \code{\link{contr.sum}},
  \code{\link{contr.treatment}};
  \code{\link{glm}},
  \code{\link{aov}},
  \code{\link{lm}},
  \code{\link{bayesglm}}.
}

\examples{
 cat.var <- rep(1:3, 5)
 dim(contr.bayes.unordered(cat.var))
 # 15*15 baseline level kept!
 dim(contr.treatment(cat.var))
 # 15*14
}
\keyword{design}
\keyword{regression}
\keyword{array}
\keyword{manip}
