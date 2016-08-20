\name{traceplot}
%\docType{genericFunction}
\alias{traceplot}
\alias{traceplot.default}
\alias{traceplot,mcmc.list-method}
\alias{traceplot,bugs-method}


\title{Trace plot of \sQuote{bugs} object}

\usage{
\S4method{traceplot}{bugs}( x, mfrow = c( 1, 1 ), varname = NULL,
  match.head = TRUE, ask = TRUE,
  col = rainbow( x$n.chains ),
  lty = 1, lwd = 1, \dots)
}

\arguments{
  \item{x}{A bugs object}
  \item{mfrow}{graphical parameter (see \code{par})}
  \item{varname}{vector of variable names to plot}
  \item{match.head}{ matches the variable names by the beginning of the variable names in bugs object}
  \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see
          \code{par(ask=.)}.}
  \item{col}{graphical parameter (see \code{par})}
  \item{lty}{graphical parameter (see \code{par})}
  \item{lwd}{graphical parameter (see \code{par})}
  \item{\dots}{further graphical parameters}
}


\description{
Displays a plot of iterations \emph{vs.} sampled values for each variable
in the chain, with a separate plot per variable.
}

\author{
  Masanao Yajima \email{yajima@stat.columbia.edu}.
  Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}

\seealso{
   \code{\link[coda]{densplot}}, \code{\link[coda]{plot.mcmc}},
   \code{\link[coda]{traceplot}}
}

\keyword{hplot}
