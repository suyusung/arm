\name{fround}
\alias{fround}
\alias{pfround}

\title{Formatting the Rounding of Numbers}

\description{
    \code{fround} rounds the values in its first argument to the specified
    number of decimal places and returns them as character strings.

    \code{pfround} does the same, but prints the result without the
    surrounding quotes that R normally puts around character output.

}

\usage{
fround(x, digits)
pfround(x, digits)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{digits}{integer indicating the precision to be used.}
}


\author{Andrew Gelman \email{gelman@stat.columbia.edu}; 
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
        }

\seealso{
    \code{\link{round}}
        }

\examples{
    x <- rnorm(1)
    fround(x, digits=2)
    pfround(x, digits=2)
}

\keyword{manip}
\keyword{print}
