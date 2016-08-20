\name{balance}
\docType{class}
% Classes
\alias{balance-class}
% Function
\alias{balance}
% display methods
\alias{print.balance}
\alias{plot.balance}

\title{Functions to compute the balance statistics}
\description{
  This function computes the balance statistics before and after matching.
}
\usage{
balance(rawdata, matched, pscore.fit, factor=TRUE)

\method{print}{balance}(x, \dots, digits = 2)

\method{plot}{balance}(x, longcovnames = NULL, 
    main = "Standardized Difference in Means",
    v.axis=TRUE, cex.main = 1, 
    cex.vars = 0.8, cex.pts = 0.8, 
    mar=c(0,3,5.1,2), plot=TRUE, \dots)
}

\arguments{
  \item{rawdata}{data before using \code{matching} function, see the example below.}
  \item{matched}{matched data using \code{matching} function, see the example below.}
  \item{pscore.fit}{glm.fit object to get propensity scores.}
  \item{factor}{default is \code{TRUE} which will display the factorized categorical variables. In
        a situation where no equal levels of factorized categorical variables is observed,
        use factor=FALSE to proceed.}
  \item{x}{an object return by the balance function.}
  \item{digits}{minimal number of \emph{significant} digits, default is 2.}
  \item{longcovnames}{long covariate names.  If not provided, plot will 
            use covariate variable name by default}
  \item{main}{The main title (on top) using font and size (character
          expansion) \code{par("font.main")} and color \code{par("col.main")};
          default title is \code{Standardized Difference in Means}.}
  \item{v.axis}{default is \code{TRUE}, which shows the top axis--axis(3).}
  \item{cex.main}{font size of main title}
  \item{cex.vars}{font size of variabel names}
  \item{cex.pts}{point size of the estimates}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)}
          which gives the number of lines of margin to be specified on
          the four sides of the plot. The default is \code{c(0,3,5.1,2)}.}
  \item{plot}{default is \code{TRUE}, which will plot the plot.}
  \item{\dots}{other plot options may be passed to this function}
}
\details{
  This function plots the balance statistics before and after matching.  
  The open circle dots represent the unmatched balance statistics.  The 
  solid dots represent the matched balance statistics. The closer the value
  of the estimates to the zero, the better the treated and control groups are 
  balanced after matching.
}

\note{ 
 The function does not work with predictors that contain factor(x), log(x) or all
 other data transformation. Create new objects for these variables.  Attach
 them into the original dataset before doing the matching procedure.
}

\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press. (Chapter 10)}

\author{Jennifer Hill \email{jh1030@columbia.edu}; 
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
        }

\seealso{\code{\link{matching}},
        \code{\link{par}}
        }
        
\examples{
# matching first
old.par <- par(no.readonly = TRUE)
data(lalonde)
attach(lalonde)
fit <- glm(treat ~ re74 + re75 + age + factor(educ) + 
            black + hisp + married + nodegr + u74 + u75, 
            family=binomial(link="logit"))
pscores <- predict(fit, type="link")
matches <- matching(z=lalonde$treat, score=pscores)
matched <- lalonde[matches$matched,]

# balance check
b.stats <- balance(lalonde, matched, fit)
print(b.stats)
plot(b.stats)
par(old.par)
}


\keyword{methods}
\keyword{manip}
\keyword{hplot}
\keyword{dplot}
