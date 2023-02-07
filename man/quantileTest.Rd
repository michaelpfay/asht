\name{quantileTest}
\alias{quantileTest}
\alias{quantileTest.default}
\alias{quantileTest.ordered}
\alias{medianTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests and Confidence Intervals about a Quantile.
}
\description{
The ath quantile of a distribution is the value, q, such that 
F(q-) <= a <= F(q), where F(x)=Pr[X <= x]. 
These are exact tests and confidence intervals on independent observations that do not any assumptions on the distribution, F. For example, the tests are exact when data are discrete or continuous, and when the distribution is non-symmetric.
}
\usage{
\method{quantileTest}{ordered}(x,\dots)

\method{quantileTest}{default}(x, q = 0, prob = 0.5, 
   alternative = c("two.sided", "less", "greater"), 
   conf.level = 0.95, \dots)

medianTest(x, m=0, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of numeric, integer or ordered factor values
}
  \item{q}{
null quantile for test
}
\item{m}{ null median for test }
  \item{prob}{ quantile
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
  \item{conf.level}{
confidence level of the interval
}
\item{\dots}{further arguments to be passed to or from methods.}
}
\details{
A test on the quantile. The \code{medianTest} is just a wrapper function to call \code{quantileTest} with \code{prob}=0.5. 

Ordinal factors may be used. The calculations just use as.numeric(x) for the factors, then return the character associated with that value. Estimates that are between two ordered factors, say "C" and "D", return the character "C/D". 
}
\value{
A list of class 'htest'. 
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{}
\author{
Michael P. Fay
}
%\note{}


\seealso{
\code{\link{signTest}}
}
\examples{
## For Poisson(mean=2.5) the median is 2
x<-rpois(20,2.5)
medianTest(x)
x<-ordered(c(rep("A",10),rep("B",60),rep("C",30)),levels=c("A","B","C"))
xnum<-as.numeric(x)
quantileTest(xnum,q=2,prob=0.705)
quantileTest(x,q=2,prob=0.705)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}