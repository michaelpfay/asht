\name{bfControl}
\alias{bfControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algorithm control arguments for Behrens-Fisher test
}
\description{
Usually these arguments do not need to be changed.
}
\usage{
bfControl(calcmethod = c("int", "mc"), epsilon = 10^(-8), nmc = 10^5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calcmethod}{
either 'int' for numeric integration (default), or 'mc' for Monte Carlo estimation}
  \item{epsilon}{
small value input into \code{\link{pbf}} or \code{\link{qbf}}
}
  \item{nmc}{
number of Monte Carlo replications used when calcmethod='mc'
}
}
\details{
When calcmethod='int' does numeric integration. This calls \code{\link{qbf}} which uses the \code{epsilon} argument. 

When calcmethod='mc' does Monte Carlo estimation of p-value and confidence interval. Uses \code{nmc} as the number of replicates.
}
\value{
A list with elements:
\item{calcmethod}{either 'int' or 'mc'}
\item{epsilon}{small positive number}
\item{nmc}{number of Monte Carlo replications used when calcmethod='mc'}
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
