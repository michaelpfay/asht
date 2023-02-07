\name{methodRuleWMW}
\alias{methodRuleWMW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to pick the method for \code{wmwTest} given the data and \code{exact} argument.
}
\description{
Inputs x,y,exact, and chooseLimit, and outputs the method, one of 'asymptotic', 'exact.ce', or 'exact.mc'
}
\usage{
methodRuleWMW(x, y, exact, chooseLimit = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of numeric responses from group 1
}
  \item{y}{
vector of numeric responses from group 2
}
  \item{exact}{
logical, should exact methods be used?
}
  \item{chooseLimit}{
boundary for choosing between methods
}
}
\details{
Let J=choose(m+n,n) where length(x)=m and length(y)=n. If exact=NULL then return 'exact.ce' if J is less than or equal to chooseLimit, otherwise return 'asymptotic'. If exact=FALSE return 'asymptotic'. If exact=TRUE then 
return 'exact.ce' if J is less than or equal to chooseLimit, otherwise return 'exact.mc'.

}
\value{ a character vector with one element, either:
\item{asymptotic}{for using normal approximations}
\item{exact.ce}{for using exact methods with complete enumeration}
\item{exact.mc}{for a Monte Carlo implementation of the exact version of the test}

}
%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
methodRuleWMW(rnorm(7),rnorm(7),exact=NULL)
methodRuleWMW(rnorm(7),rnorm(8),exact=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
