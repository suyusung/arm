\name{GO}
\docType{class}
% Classes
\alias{GO-class}
% functions
\alias{go}
\alias{G}

\title{Function to Recall Last Source File}

\description{
    A function that is like \code{source()}, but that recalls the last source file names by default.
}
\usage{
go(\dots, add=FALSE, timer=FALSE)
}
\arguments{
  \item{\dots}{list of filenames as character strings.}
  \item{add}{add these names to the current list; if \code{FALSE}, replace the list instead.}
  \item{timer}{time the execution of \code{go()}.}
}



\author{
Jouni Kerman \email{jouni@kerman.com} \email{kerman@stat.columbia.edu}
}


\examples{
    go('myprog')            # will run source('myprog.r')
    go()                    # will run source('myprog.r') again
    go('somelib',add=TRUE)  # will run source('myprog.r') and source('somelib.r')
    go('myprog','somelib')  # same as above
    go('mytest')            # will run source('mytest') only
    go()                    # runs source('mytest') again
    G                       # short cut to call go()
 }

\keyword{methods}
\keyword{manip}
