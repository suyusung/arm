\name{residual.plot}
\alias{residual.plot}

\title{residual plot for the observed values}
\description{
  Plots the residual of observed variable.
}
\usage{
residual.plot(Expected, Residuals, sigma, main = deparse(substitute(Expected)), 
  col.pts = "blue", col.ctr = "red", col.sgm = "black", cex = 0.5, 
  gray.scale = FALSE, xlab = "Predicted", ylab = "Residuals", ...)
}

\arguments{
  \item{Expected}{ Expected value. }
  \item{Residuals}{ Residual value. }
  \item{sigma}{ Standard error. }
  \item{main}{ main for the plot.  See \code{plot} for detail.}
  \item{col.pts}{ Color of the points. }
  \item{col.ctr}{ Color of the line at zero. }
  \item{col.sgm}{ Color of standard error line. }
  \item{cex}{ A numerical value giving the amount by which plotting text
          and symbols should be magnified relative to the default.  See par for detail. }
  \item{gray.scale}{ If \code{TRUE}, makes the plot into black and white.  This option overwrites the color specification.  Default is FALSE. }
  \item{xlab}{ Label for x axis. }
  \item{ylab}{ Label for y axis. }
  \item{\dots}{ Additional parameters passed to \code{plot} function. }
}

\value{
  Plot to visualize pattern of residulal value for the expected value.
}

\author{
  Masanao Yajima \email{yajima@stat.columbia.edu},
  M.Grazia Pittau \email{grazia@stat.columbia.edu}
}

\examples{
old.par <- par(no.readonly = TRUE)

x <- rnorm(100)
y <- rnorm(100)
fit <- lm(y~x)
y.hat <- fitted(fit)
u <- resid(fit)
sigma <- sigma.hat(fit)
residual.plot(y.hat, u, sigma)

par(old.par)
}


\keyword{hplot}
