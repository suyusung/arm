\name{lalonde}
\alias{lalonde}
\docType{data}
\title{Lalonde Dataset}
\description{
  Dataset used by Dehejia and Wahba (1999) to evaluate propensity score matching.
}
\usage{data(lalonde)}
\format{
  A data frame with 445 observations on the following 12 variables.
  \describe{
    \item{age}{age in years.}
    \item{educ}{years of schooling.}
    \item{black}{indicator variable for blacks.}
    \item{hisp}{indicator variable for Hispanics.}
    \item{married}{indicator variable for martial status.}
    \item{nodegr}{indicator variable for high school diploma.}
    \item{re74}{real earnings in 1974.}
    \item{re75}{real earnings in 1975.}
    \item{re78}{real earnings in 1978.}
    \item{u74}{indicator variable for earnings in 1974 being zero.}
    \item{u75}{indicator variable for earnings in 1975 being zero.}
    \item{treat}{an indicator variable for treatment status.}
  }
}
\details{
  Two demos are provided which use this dataset.  The first,
  \code{DehejiaWahba}, replicates one of the models from Dehejia and
  Wahba (1999).  The second demo, \code{AbadieImbens}, replicates the
  models produced by Abadie and Imbens
  \url{http://scholar.harvard.edu/imbens/scholar_software/matching-estimators}.
  Many of these models are found to produce good balance for the Lalonde
  data.
  }
\references{
  Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in
  Non-Experimental Studies: Re-Evaluating the
  Evaluation of Training Programs.'' \emph{Journal of the American Statistical
  Association} 94 (448): 1053-1062.

  LaLonde, Robert. 1986.  ``Evaluating the Econometric Evaluations of
  Training Programs.''  \emph{American Economic Review} 76:604-620.
}

\note{This documentation is adapted from \code{Matching} package.}

\seealso{\code{\link{matching}},
        \code{\link[Matching]{GenMatch}}
            \code{\link{balance}}
            }

\examples{
data(lalonde)
}
\keyword{datasets}
