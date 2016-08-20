\name{triangleplot}
\alias{triangleplot}

\title{Triangle Plot}

\description{
Function for making a triangle plot from a square matrix
}

\usage{
triangleplot (x, y=NULL, cutpts=NULL, details=TRUE, 
            n.col.legend=5, cex.col=0.7, 
            cex.var=0.9, digits=1, color=FALSE)
}

\arguments{
  \item{x}{a square matrix.}
  \item{y}{a vector of names that corresponds to each element of the square matrix x.}
  \item{cutpts}{a vector of cutting points for color legend, default is \code{NULL}. 
        The function will decide the cutting points if cutpts is not assigned.}
  \item{details}{show more than one digits correlaton values. Default 
              is \code{TRUE}. \code{FALSE} is suggested to get readable output.}
  \item{n.col.legend}{number of legend for the color thermometer}
  \item{cex.col}{font size of the color thermometer.}
  \item{cex.var}{font size of the variable names.}
  \item{digits}{number of digits shown in the text of the color theromoeter.}
  \item{color}{color of the plot, default is FALSE, which uses gray scale.}
}

\details{
The function makes a triangle plot from a square matrix, e.g., the correlation plot, see 
\code{\link{corrplot}}.  If a square matrix contains missing values, the cells of missing values
will be marked \code{x}.
}


\author{
Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}
        
\seealso{\code{\link{corrplot}},
        \code{\link[graphics]{par}}
}

\examples{
old.par <- par(no.readonly = TRUE)

 # create a square matrix
 x <- matrix(runif(1600, 0, 1), 40, 40)
 
 # fig 1
 triangleplot(x)
 
 # fig 2 assign cutting points
 triangleplot(x, cutpts=c(0,0.25,0.5,0.75,1), digits=2)
 
 # fig 3 if x contains missing value
 x[12,13] <- x[13,12] <- NA
 x[25,27] <- x[27,25] <- NA
 triangleplot(x)
 
par(old.par)

#
#library(RColorBrewer)
#cormat <-  cor(iris[,-5])
#triangleplot2(cormat,color = brewer.pal( 5, "RdBu" ),
# n.col.legend=5, cex.col=0.7, cex.var=0.5)

}

\keyword{dplot}
