\name{binnedplot}
\alias{binnedplot}
\alias{binned.resids}

\title{Binned Residual Plot}
\description{
  A function that plots averages of y versus averages of x and can be 
  useful to plot residuals for logistic regression.
}
\usage{
 binnedplot(x ,y, nclass=NULL, 
    xlab="Expected Values", ylab="Average residual", 
    main="Binned residual plot", 
    cex.pts=0.8, col.pts=1, col.int="gray", 
    se = c("model", "empirical"), ...)
}

\arguments{
  \item{x}{The expected values from the logistic regression.}
  \item{y}{The residual values from logistic regression (observed values
            minus expected values).}
  \item{nclass}{Number of categories (bins), based on their fitted values, into
                  which the data are divided.  Default \code{NULL}, in which case
                  \code{nclass = floor(sqrt(length(x)))} if n >= 100,
                  \code{nclass = 10} if 10 < n < 100, and
                  \code{nclass = floor(n/2)} if n <= 10.}
  \item{xlab}{a label for the x axis, default is "Expected Values".}
  \item{ylab}{a label for the y axis, default is "Average residual".}       
  \item{main}{a main title for the plot, default is "Binned residual plot". 
    See also \code{title}.}
  \item{cex.pts}{The size of points, default=0.8.}
  \item{col.pts}{color of points, default is black.}
  \item{col.int}{color of intervals, default is gray.}
  \item{se}{which standard error the plotted bounds use.  \code{"model"},
    the default, uses \code{2*sqrt(p*(1-p)/n)} with \code{p} the bin's mean
    expected value, as in Gelman and Hill (2007, section 5.6);
    \code{"empirical"} uses the empirical standard deviation of the residuals
    in the bin, as earlier versions of the package used.  See Details.}
  \item{...}{Graphical parameters to be passed to methods}
}

\details{
In logistic regression, as with linear regression, the residuals can be defined as 
observed minus expected values. The data are discrete and so are the residuals. 
As a result, plots of raw residuals from logistic regression are generally not useful.
The binned residuals plot instead, after dividing the data into categories (bins) based
on their fitted values, plots the average residual versus the average fitted
value for each bin.

The bounds drawn by \code{binnedplot} are, by default, the model-based ones of
Gelman and Hill (2007, section 5.6), \code{2*sqrt(p*(1-p)/n)}, where \code{p} is
the bin's mean expected value and \code{n} the number of observations in the bin.
That formula assumes \code{x} is on the probability scale and that \code{y} is the
response residual \code{y - p}, so call it with \code{predict(fit, type =
"response")} and \code{resid(fit, type = "response")}; on the link scale, or with
deviance residuals, the bounds mean nothing (\code{binned.resids} warns when
\code{x} falls outside \code{[0, 1]}).

With \code{se = "empirical"} the bounds use the empirical standard deviation of
the residuals in the bin, \code{2*sd(y)/sqrt(n)}, as earlier versions always did.
It is unreliable when bins hold few events: in a bin where every outcome is 0 the
residuals are all \code{-p}, so their sd is near zero and the band collapses onto
the axis, and a bin with a single event gives a band that is far too wide.
}
\value{
\code{binnedplot} is called for its side effect: a plot in which the gray lines
indicate plus and minus 2 standard-error bounds, within which one would expect
about 95\% of the binned residuals to fall if the model were actually true.
By default the standard error is the model-based \code{sqrt(p*(1-p)/n)}; see
Details, and the \code{se} argument for the empirical alternative.

\code{binned.resids} returns a list with components \code{binned}, a matrix of
bin summaries, and \code{xbreaks}, the breakpoints of the bins.
}
\references{ Andrew Gelman and Jennifer Hill, 
             Data Analysis Using Regression and Multilevel/Hierarchical Models, 
             Cambridge University Press, 2006.}

\author{M. Grazia Pittau \email{grazia@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}
        
\note{There is typically some arbitrariness in choosing the number of bins: 
     each bin should contain enough points so that the averaged residuals are not too noisy,
     but it helps to have also many bins so as to see more local patterns in the residuals (see 
     Gelman and Hill, Data Analysis Using Regression and Multilevel/Hierarchical Models, p. 97).}

\seealso{\code{\link{par}}, \code{\link{plot}}}

\examples{
old.par <- par(no.readonly = TRUE)
 data(lalonde)
 attach(lalonde)
 fit <- glm(treat ~ re74 + re75 + educ + black + hisp + married 
               + nodegr + u74 + u75, family=binomial(link="logit"))
 x <- predict(fit, type = "response")
 y <- resid(fit, type = "response")
 binnedplot(x, y)
 binnedplot(x, y, se = "empirical")   # the empirical bounds
par(old.par)
}
\keyword{dplot}
