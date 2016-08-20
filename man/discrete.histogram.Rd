\name{discrete.histogram}
\alias{discrete.histogram}
\alias{discrete.hist}

\title{Histogram for Discrete Distributions}

\description{Creates a prettier histogram for discrete distributions}

\usage{
discrete.histogram (x, prob, prob2=NULL, prob3=NULL,
    xlab="x", xaxs.label=NULL, yaxs.label=NULL, bar.width=NULL,
    freq=FALSE, prob.col="blue", prob2.col="red", prob3.col="gray", ...)
}

\arguments{
  \item{x}{The vector of x's}
  \item{prob}{The probabilities for the x's}
  \item{prob2}{A second vector of probabilities of the x's}
  \item{prob3}{A third vector of probabilities of the x's}
  \item{xlab}{Label for the x axis}
  \item{xaxs.label}{Label for the x's}
  \item{yaxs.label}{Label for the y axis}
  \item{bar.width}{Width of the bars}
  \item{freq}{If TRUE, shows a frequency histogram as opposed to probability.}
  \item{prob.col}{The color of the first set of histogram bars.}
  \item{prob2.col}{The color of the second set of histogram bars.}
  \item{prob3.col}{The color of the third set of histogram bars.}
  \item{...}{Additional arguments passed to function \code{plot}}
  }

\details{This function displays a histogram for discrete
  probability distributions.
}


\examples{
a <- c(3,4,0,0,5,1,1,1,1,0)
discrete.histogram (a)

x <- c(0,1,3,4,5)
p <- c(.3,.4,.1,.1,.1)
discrete.histogram (x,p)

x <- c(0,1,3,4,5)
y <- c(3,4,1,1,1)
discrete.histogram (x,y)
}


\keyword{dplot}
