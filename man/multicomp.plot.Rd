\name{multicomp.plot}
\alias{multicomp.plot}
\alias{mcplot}

\title{Multiple Comparison Plot}

\description{
  Plots significant difference of simulated array.
}

\usage{
multicomp.plot(object, alpha = 0.05, main = "Multiple Comparison Plot",
  label = NULL, shortlabel = NULL, show.pvalue = FALSE,
  label.as.shortlabel = FALSE, label.on.which.axis = 3,
  col.low = "lightsteelblue", col.same = "white", col.high = "lightslateblue",
  vertical.line = TRUE, horizontal.line = FALSE,
  vertical.line.lty = 1, horizontal.line.lty = 1, mar=c(3.5,3.5,3.5,3.5))
}

\arguments{
  \item{object}{Simulated array of coefficients, columns being different variables
      and rows being simulated result.}
  \item{alpha}{Level of significance to compare.}
  \item{main}{Main label.}
  \item{label}{Labels for simulated parameters.}
  \item{shortlabel}{Short labels to put into the plot.}
  \item{show.pvalue}{Default is FALSE, if set to TRUE replaces short label
      with Bayesian p value. }
  \item{label.as.shortlabel}{Default is FALSE, if set to TRUE takes first 2 character
      of label and use it as short label.}
  \item{label.on.which.axis}{default is the 3rd (top) axis.}
  \item{col.low}{Color of significantly low coefficients.}
  \item{col.same}{Color of not significant difference.}
  \item{col.high}{Color of significantly high coefficients.}
  \item{vertical.line}{Default is TRUE, if set to FALSE does not draw vertical line.}
  \item{horizontal.line}{Default is FALSE, if set to TRUE draws horizontal line.}
  \item{vertical.line.lty}{Line type of vertical line.}
  \item{horizontal.line.lty}{Line type of horizontal line.}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)}
          which gives the number of lines of margin to be specified on
          the four sides of the plot. The default is \code{c(3.5,3.5,3.5,3.5)}.}

}

\value{
  \item{pvalue}{Array of Bayesian p value.}
  \item{significant}{Array of significance.}
}

\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.
}


\author{
  Masanao Yajima \email{yajima@stat.columbia.edu},
  Andrew Gelman \email{gelman@stat.columbia.edu}
}

\seealso{
  \code{\link{coefplot}}
}

\examples{
old.par <- par(no.readonly = TRUE)

# example 1
simulation.array <- data.frame(coef1=rnorm(100,10,2), coef2=rnorm(100,5,2),
                      coef3=rnorm(100,0,1), coef4=rnorm(100,-5,3),
                      coef5=rnorm(100,-2,1))
short.lab <- c("c01", "c02", "c03", "c04", "c05")
multicomp.plot(simulation.array[,1:4], label.as.shortlabel=TRUE)

# wraper for multicomp.plot
mcplot(simulation.array, shortlabel = short.lab)

# example 2
data(lalonde)
M1 <- lm(re78 ~ treat + re74 + re75 + age + educ + u74 + u75, data=lalonde)
M1.sim <- sim(M1)
lm.sim <- coef(M1.sim)[,-1]
multicomp.plot(lm.sim, label.as.shortlabel=TRUE, label.on.which.axis=2)

par(old.par)
}
\keyword{hplot}
