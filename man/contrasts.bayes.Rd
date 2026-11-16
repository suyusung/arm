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
  \item{base}{an integer between 1 and the number of levels.  It is only
          checked to be in range (when \code{contrasts} is \code{TRUE}); it does
          not change the result, since no level is dropped.}
  \item{contrasts}{a logical.  For \code{contr.bayes.unordered} it does not
          change the result.  For \code{contr.bayes.ordered}, \code{FALSE}
          replaces the first (constant) column by a column of ones.}
  \item{scores}{the set of values over which orthogonal polynomials are to be
          computed.}
}

\details{
   These functions are adapted from \code{contr.treatment} and \code{contr.poly}
   in the \code{stats} package.  Their purpose is to keep the baseline levels
   of categorical variables, which suits the use of 
   \code{\link{bayesglm}}.
   
   \code{contr.bayes.unordered} returns the n by n identity matrix: one
   indicator column for every level, where \code{contr.treatment} drops the
   column of the baseline level.  \code{contr.bayes.ordered} returns n
   orthonormal polynomial columns: a constant column named \code{"^0"} followed
   by the n - 1 columns of \code{contr.poly}.
   
}


\author{Yu-Sung Su \email{suyusung@tsinghua.edu.cn}}

\seealso{
  \code{\link{C}},
  \code{\link{contr.helmert}},
  \code{\link{contr.poly}},
  \code{\link{contr.sum}},
  \code{\link{contr.treatment}},
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
