\name{coefplot}
%\docType{genericFunction}
\alias{coefplot}
\alias{coefplot.default}
\alias{coefplot,numeric-method}
\alias{coefplot,lm-method}
\alias{coefplot,glm-method}
\alias{coefplot,bugs-method}
\alias{coefplot,polr-method}


\title{Generic Function for Making Coefficient Plot}
\description{
 Functions that plot the coefficients plus and minus 1 and 2 sd 
 from a lm, glm, bugs, and polr fits.
}
\usage{
coefplot(object,\dots)

\method{coefplot}{default}(coefs, sds, CI=2, 
                lower.conf.bounds, upper.conf.bounds,
                varnames=NULL, vertical=TRUE, 
                v.axis=TRUE, h.axis=TRUE,
                cex.var=0.8, cex.pts=0.9, 
                col.pts=1, pch.pts=20, var.las=2,
                main=NULL, xlab=NULL, ylab=NULL, mar=c(1,3,5.1,2),
                plot=TRUE, add=FALSE, offset=.1, \dots)

\S4method{coefplot}{bugs}(object, var.idx=NULL, varnames=NULL, 
            CI=1, vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE, 
            cex.var=0.8, cex.pts=0.9, 
            col.pts=1, pch.pts=20, var.las=2, 
            main=NULL, xlab=NULL, ylab=NULL, 
            plot=TRUE, add=FALSE, offset=.1,
            mar=c(1,3,5.1,2), \dots)

\S4method{coefplot}{numeric}(object, \dots)
\S4method{coefplot}{lm}(object, varnames=NULL, intercept=FALSE, \dots)
\S4method{coefplot}{glm}(object, varnames=NULL, intercept=FALSE, \dots)
\S4method{coefplot}{polr}(object, varnames=NULL, \dots)
}

\arguments{
  \item{object}{fitted objects-lm, glm, bugs and polr, or a vector of coefficients.}
  \item{...}{further arguments passed to or from other methods.}
  \item{coefs}{a vector of coefficients.}
  \item{sds}{a vector of sds of coefficients.}
  \item{CI}{confidence interval, default is 2, which will plot plus and minus 2 sds or
    95\% CI. If CI=1, plot plus and minus 1 sds or 50\% CI instead.}
  \item{lower.conf.bounds}{lower bounds of confidence intervals.}
  \item{upper.conf.bounds}{upper bounds of confidence intervals.}
  \item{varnames}{a vector of variable names, default is NULL, which will use
    the names of variables; if specified, the length of varnames must be equal to
    the length of predictors, including the intercept.}
  \item{vertical}{orientation of the plot, default is TRUE which will plot
    variable names in the 2nd axis. If FALSE, plot variable names in the first
    axis instead.}
  \item{v.axis}{default is TRUE, which shows the bottom axis--axis(1).}
  \item{h.axis}{default is TRUE, which shows the left axis--axis(2).}
  \item{cex.var}{The fontsize of the varible names, default=0.8.}
  \item{cex.pts}{The size of data points, default=0.9.}
  \item{col.pts}{color of points and segments, default is black.}
  \item{pch.pts}{symbol of points, default is solid dot.}
  \item{var.las}{the orientation of variable names against the axis, default is 2. 
    see the usage of \code{las} in \code{\link{par}}.}
  \item{main}{The main title (on top) using font and size (character
          expansion) \code{par("font.main")} and color \code{par("col.main")}.}
  \item{xlab}{X axis label using font and character expansion
          \code{par("font.lab")} and color \code{par("col.lab")}.}
  \item{ylab}{Y axis label, same font attributes as \code{xlab}.}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)}
          which gives the number of lines of margin to be specified on
          the four sides of the plot. The default is \code{c(1,3,5.1,2)}.}
  \item{plot}{default is TRUE, plot the estimates.}
  \item{add}{if add=TRUE, plot over the existing plot. default is FALSE.}
  \item{offset}{add extra spaces to separate from the existing dots. default is 0.1.}
%  \item{lower.bound}{default is -Inf.}
  \item{var.idx}{the index of the variables of a bugs object, default
        is NULL which will plot all the variables.}
  \item{intercept}{If TRUE will plot intercept, default=FALSE to get better presentation.}
}

\details{
  This function plots coefficients from bugs, lm, glm and polr
  with 1 sd and 2 sd interval bars.
}

\value{
   Plot of the coefficients from a bugs, lm or glm fit. You can add the intercept, 
   the variable names and the display the result of the fitted model. 
}

\references{Andrew Gelman and Jennifer Hill, 
            Data Analysis Using Regression and Multilevel/Hierarchical Models, 
            Cambridge University Press, 2006.}
            
\author{Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}
        
\seealso{\code{\link{display}},
         \code{\link[graphics]{par}},
         \code{\link[stats]{lm}},
         \code{\link[stats]{glm}},
         \code{\link{bayesglm}},
         \code{\link[graphics]{plot}}
        }
\examples{
old.par <- par(no.readonly = TRUE)

 y1 <- rnorm(1000,50,23)
 y2 <- rbinom(1000,1,prob=0.72)
 x1 <- rnorm(1000,50,2) 
 x2 <- rbinom(1000,1,prob=0.63) 
 x3 <- rpois(1000, 2) 
 x4 <- runif(1000,40,100) 
 x5 <- rbeta(1000,2,2) 
 
 longnames <- c("a long name01","a long name02","a long name03",
                "a long name04","a long name05")

 fit1 <- lm(y1 ~ x1 + x2 + x3 + x4 + x5)
 fit2 <- glm(y2 ~ x1 + x2 + x3 + x4 + x5, 
            family=binomial(link="logit"))
 op <- par()
 # plot 1
 par (mfrow=c(2,2))
 coefplot(fit1)
 coefplot(fit2, col.pts="blue")
 
 # plot 2
 longnames <- c("(Intercept)", longnames) 
 coefplot(fit1, longnames, intercept=TRUE, CI=1)
 
 # plot 3
 coefplot(fit2, vertical=FALSE, var.las=1, frame.plot=TRUE)
 
 # plot 4: comparison to show bayesglm works better than glm
 n <- 100
 x1 <- rnorm (n)
 x2 <- rbinom (n, 1, .5)
 b0 <- 1
 b1 <- 1.5
 b2 <- 2
 y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))
 y <- ifelse (x2==1, 1, y)
 x1 <- rescale(x1)
 x2 <- rescale(x2, "center")
 
 M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
       display (M1)
 M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))
       display (M2)

#=================== 
#    stacked plot
#===================
  coefplot(M2, xlim=c(-1,5), intercept=TRUE)
  coefplot(M1, add=TRUE, col.pts="red")
  
#==================== 
# arrayed plot       
#====================
  par(mfrow=c(1,2))
  x.scale <- c(0, 7.5) # fix x.scale for comparison
  coefplot(M1, xlim=x.scale, main="glm", intercept=TRUE)
  coefplot(M2, xlim=x.scale, main="bayesglm", intercept=TRUE)

# plot 5: the ordered logit model from polr
 M3 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
 coefplot(M3, main="polr")
   
 M4 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
 coefplot(M4, main="bayespolr", add=TRUE, col.pts="red")

## plot 6: plot bugs & lmer
# par <- op
# M5 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
# M5.sim <- mcsamp(M5)
# coefplot(M5.sim, var.idx=5:22, CI=1, ylim=c(18,1), main="lmer model")


# plot 7: plot coefficients & sds vectors
 coef.vect <- c(0.2, 1.4, 2.3, 0.5)
 sd.vect <- c(0.12, 0.24, 0.23, 0.15)
 longnames <- c("var1", "var2", "var3", "var4")
 coefplot (coef.vect, sd.vect, varnames=longnames, main="Regression Estimates")
 coefplot (coef.vect, sd.vect, varnames=longnames, vertical=FALSE, 
    var.las=1, main="Regression Estimates")
    
par(old.par)
}

\keyword{hplot}
\keyword{dplot}
\keyword{methods}
\keyword{manip}
