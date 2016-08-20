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
    cex.pts=0.8, col.pts=1, col.int="gray", ...)
}

\arguments{
  \item{x}{The expected values from the logistic regression.}
  \item{y}{The residuals values from logistic regression (observed values
            minus expected values).}
  \item{nclass}{Number of categories (bins) based on their fitted values in which
                  the data are divided. Default=NULL and will take the value of nclass
                  according to the $n$ such that if $n >=100$, 
                  nclass=floor(sqrt(length(x))); if $10<n<100$, nclass=10;
                  if $n<10$, nclass=floor(n/2).}
  \item{xlab}{a label for the x axis, default is "Expected Values".}
  \item{ylab}{a label for the y axis, default is "Average residual".}       
  \item{main}{a main title for the plot, default is "Binned residual plot". 
    See also \code{title}.}
  \item{cex.pts}{The size of points, default=0.8.}
  \item{col.pts}{color of points, default is black}
  \item{col.int}{color of intervals, default is gray}
  \item{...}{Graphical parameters to be passed to methods}
}

\details{
In logistic regression, as with linear regression, the residuals can be defined as 
observed minus expected values. The data are discrete and so are the residuals. 
As a result, plots of raw residuals from logistic regression are generally not useful.
The binned residuals plot instead, after dividing the data into categories (bins) based
on their fitted values, plots the average residual versus the average fitted
value for each bin.
}
\value{
A plot in which the gray lines indicate plus and minus 2 standard-error
bounds, within which one would expect about 95\% of the binned residuals to fall,
if the model were actually true. 
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
     Gelman and Hill, Data Analysis Using Regression and Multilevel/Hierarchical Models, pag 97).}

\seealso{\code{\link{par}}, \code{\link{plot}}}

\examples{
old.par <- par(no.readonly = TRUE)
 data(lalonde)
 attach(lalonde)
 fit <- glm(treat ~ re74 + re75 + educ + black + hisp + married 
               + nodegr + u74 + u75, family=binomial(link="logit"))
 x <- predict(fit)
 y <- resid(fit)
 binnedplot(x,y)
par(old.par)
}
\keyword{dplot}
