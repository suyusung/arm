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
balance(rawdata, treat, matched, estimand="ATT")

\method{print}{balance}(x, \dots, combined = FALSE, digits = 2)

\method{plot}{balance}(x, longcovnames=NULL, which.covs="mixed",
    v.axis=TRUE, cex.main=1, cex.vars=1, cex.pts=1,
    mar=c(4, 3, 5.1, 2), plot=TRUE, x.max = NULL, \ldots)
}

\arguments{
  \item{rawdata}{The full covariate dataset}
  \item{treat}{the vector of treatment assignments for the full dataset}
  \item{matched}{vector of weights to apply to the full dataset to create the
   restructured data: for matching without replacement these will all be 0's and 1's; 
   for one-to-one matching with replacement these will all be non-negative
   integers; for IPTW or more complicated matching methods these could be any non-negative numbers}
  \item{estimand}{can either be \code{ATT}, \code{ATC}, or \code{ATE}, default is \code{ATT}}
  \item{x}{an object return by the balance function.}
  \item{combined}{default is \code{FALSE}}
  \item{digits}{minimal number of \emph{significant} digits, default is 2.}
  \item{longcovnames}{long covariate names.  If not provided, plot will 
            use covariate variable name by default}
  \item{which.covs}{\code{mixed} then it plots all as std diffs;
    \code{binary} it only plots binary and as abs unstd diffs;
    \code{cont} it only plots non-binary and as abs std diffs}
  \item{v.axis}{default is \code{TRUE}, which shows the top axis--axis(3).}
  \item{cex.main}{font size of main title}
  \item{cex.vars}{font size of variabel names}
  \item{cex.pts}{point size of the estimates}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)}
          which gives the number of lines of margin to be specified on
          the four sides of the plot. The default is \code{c(0,3,5.1,2)}.}
  \item{plot}{default is \code{TRUE}, which will plot the plot.}
  \item{x.max}{set the max of the \code{xlim}, default is \code{NULL}}
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

\author{Jennifer Hill \email{jennifer.hill@nyu.edu}; 
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
matched <- matches$cnts

# balance check
b.stats <- balance(lalonde, treat, matched, estimand = "ATT")
print(b.stats)
plot(b.stats)
par(old.par)
}


\keyword{methods}
\keyword{manip}
\keyword{hplot}
\keyword{dplot}
