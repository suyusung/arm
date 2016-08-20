\name{rescale}
\alias{rescale}

\title{Function for Standardizing by Centering and Dividing by 2 sd's}
\description{
This function standardizes a variable by centering and dividing by 2 sd's with 
exceptions for binary variables.
}

\usage{
rescale(x, binary.inputs="center")
}

\arguments{
  \item{x}{a vector}
  \item{binary.inputs}{options for standardizing binary variables, default is \code{center}; \code{0/1} keeps original scale;
    \code{-0.5,0.5} rescales 0 as -0.5 and 1 as 0.5; \code{center} substracts the mean; and 
    \code{full} substracts the mean and divids by 2 sd.}
}

\value{
 the standardized vector
}

\references{Andrew Gelman. (2008). 
    \dQuote{Scaling regression inputs by dividing by two standard deviations}. 
    \emph{Statistics in Medicine} 27: 2865--2873.
    \url{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}
}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}

\seealso{\code{\link{standardize}}
    }

\examples{
  # Set up the fake data
  n <- 100
  x <- rnorm (n, 2, 1)
  x1 <- rnorm (n)
  x1 <- (x1-mean(x1))/(2*sd(x1))   # standardization
  x2 <- rbinom (n, 1, .5)
  b0 <- 1
  b1 <- 1.5
  b2 <- 2
  y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))
  rescale(x, "full")
  rescale(y, "center")
}
\keyword{manip}
