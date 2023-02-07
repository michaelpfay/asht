\name{abcnonHtest}
\alias{abcnonHtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric ABC (Approximate Bootstrap Confidence) intervals.
}
\description{
A hypothesis testing function using the nonparametric ABC intervals.
}
\usage{
abcnonHtest(x, tt, nullValue = NULL, conf.level = 0.95, 
   alternative = c("two.sided", "less", "greater"), epsilon = 0.001, minp = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{the data. Must be either a vector, or a matrix whose rows are
the observations} 
\item{tt}{function defining the parameter in the resampling form
\code{tt(p,x)}, where \code{p} is the vector of proportions and \code{x}
is the data} 
  \item{nullValue}{
null value of the parameter for the two-sided hypothesis test, or boundary of null parameter space for one-sided ones
}
  \item{conf.level}{
confidence level for interval
}
  \item{alternative}{
	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
\item{epsilon}{optional argument specifying step size for finite
difference calculations} 
  \item{minp}{
minimum p-value (used in uniroot search to give a bound, toe two.sided alternatives actual minimum is 2*minp)
}
}
\details{
Calculates the nonparametric ABC confidence interval of DiCiccio and Efron (1992). See also Efron and Tibshirani (1993). 

The p-values are calculated by solving for confidence limit that just touches the \code{nullValue}. If it is outside of the range (minp, 1-minp) for one-sided p-values, then it is set to minp.
If it is outside the range (2*minp, 1- 2*minp) for two-sided p-values, then it is set to 2*minp.
}
\value{
A value of class "htest" containing the following components:
\item{p.value}{p-value for test defined by alternative and nullValue}
\item{estimate}{estimate of the parameter, calculated using \code{x} and the \code{tt} function }
\item{conf.int}{confidence interval for the parameter associated with \code{tt} }
\item{null.value}{ the \code{nullValue} (or null boundary) for the hypothesis test}
\item{alternative}{ a character string describing the alternative hypothesis }
\item{method}{a character string describing the kind of test}
\item{data.name}{a character string giving the name of the data
and the function}
}
\references{
DiCiccio, T and Efron, B (1992). More accurate confidence intervals in exponential families. Biometrika 79: 231-245.

Efron, B and Tibshirani, RJ (1993). An introduction to the bootstrap. Chapman and Hall: New York.
}
\author{
the function is modification of \code{abcnon} in the \code{bootstrap} R package, originally written by Rob Tibshirani, modifications by M.P. Fay
}
%\note{}
\seealso{
See also \code{\link[bootstrap]{abcnon}}.
}
\examples{
# compute abc intervals for the mean
x <- c(2,4,12,4,6,3,5,7,6)
theta <- function(p,x) {sum(p*x)/sum(p)}
## smallest p-value is 2*minp for two-sided alternatives
abcnonHtest(x, theta, nullValue=0)  
## test null at 95% confidence limit is like just barely
## rejecting at the two-sided 5% level, so p-value is 0.05
abcnonHtest(x, theta, nullValue=4.072772)  
# compute abc intervals for the correlation
set.seed(1)
x <- matrix(rnorm(20),ncol=2)
theta <- function(p, x)
{
    x1m <- sum(p * x[, 1])/sum(p)
    x2m <- sum(p * x[, 2])/sum(p)
    num <- sum(p * (x[, 1] - x1m) * (x[, 2] - x2m))
    den <- sqrt(sum(p * (x[, 1] - x1m)^2) *
              sum(p * (x[, 2] - x2m)^2))
    return(num/den)
}
abcnonHtest(x, theta) 
## compare with 
\dontrun{
library(bootstrap)
abcnon(x, theta, alpha=c(.025,.975))$limits[,"abc"]}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
