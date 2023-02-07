\name{signTest}
\alias{signTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact Sign Test with Confidence Intervals
}
\description{
Uses \code{link{binom.exact}} or \code{\link{mcnemarExactDP}} to create sign test with confidence intervals on different parameters.
Mid-p versions are available for some parameterizations (see details).
}
\usage{
signTest(x, stat=c("cd","cpp","ud"), nullparm=NULL, 
   alternative=c("two.sided","less","greater"), conf.level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector
}
  \item{stat}{
statistic for estimates and confidence intervals,
"cd"= conditional difference: proportion positive - proportion negative, 
"cpp"= conditional proportion positive, and 
"ud"= unconditional difference: proportion positive-proportion negative
(conditional proportions are out of non-zero values, unconditional are out of all values)
}
\item{nullparm}{
null parameter value associated with \code{stat}, NULL value defaults to the exact  sign test (i.e., \code{stat="cd"} and 
code{stat="ud"} gives 0, and \code{stat="cpp"} gives 0.5).
}
\item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
\item{conf.level}{
confidence level of the interval
}
  \item{\dots}{
arguments passed to \code{\link{binom.exact}}
or \code{\link{mcnemarExactDP}}
}
}
\details{

The sign test is a conditional test, conditioning on the total number of non-zero observations testing that the proportion positive is different (or less, or greater) than the proportion negative.  When the responses are differences in paired binary observations this 
is the same as a McNemar test. 

This function gives estimates and confidence intervals compatible with the exact sign test for three different parameterizations.
Let n.pos, n.neg,
%=\code{sum(sign(x)==1)} 
and 
n.nonzero
%=\code{sum(sign(x)!=0)}
be the number of positive, negative, and non-zero observations respectively out of 
n=\code{length(x)}.
The conditional proportion positive are n.pos/n.nonzero, and the unconditional proportion positive are n.pos/n. Similarly, 
the conditional proportion negative are n.neg/n.nonzero and the unconditional proportion negative are n.neg/n. 
When \code{stat='cd'} the parameterization is the conditional difference in proportions (pos-neg),
and when \code{stat='ud'} the parameterization is the uncondtional difference in proportions (pos-neg).
The third parameterization is \code{stat='cpp'} the conditional proportion positive. 
The argument \code{nullparm} gives the null value of the test when \code{alternative='two.sided'}. When \code{nullparm=NULL}, this gives 
the traditional sign test, where \code{nullparm=0} for \code{stat='cd'} and \code{stat='ud'} and 
\code{nullparm=0.5} for \code{stat='cpp'}.  

Conditioning on m=n.nonzero, Y is binomial with parameters m and beta. So when \code{stat='cpp'} the parameter we are estimating is beta,
and when  \code{stat='cd'} the parameter we are estimating is beta - (1-beta) = 2*beta-1. 
We use \code{\link{binom.exact}} to do the p-value and confidence interval calculations. Thus, \code{midp} versions and different two-sided methods (given by \code{tsmethod}) can be calculated. 


Unconditionally, we treat M (the number non-zero) as a random variable, and assume M is binomial with parameters n and theta.  
When  \code{stat='ud'} the parameter we are estimating is delta = theta*(2*beta-1), which is the unconditional difference: (proportion positive out of the total) - (proprtion negative out of the total). We use \code{\link{mcnemarExactDP}} to do the the p-value and confidence interval calculations. The methods associated with that function are described in Fay and Lumbard (2020). As of now, when \code{stat='ud'} a \code{midp} version is not available, and the only two-sided method available is a 'central' one, meaning the error for the 95\% confidence interval is bounded by 2.5\% on each side. 

}
\value{
A list of class 'htest' (use \code{str} to see elements)
\item{statistic}{vector of number of positive, negative, zero, and non-zero}
\item{estimate}{vector of estimates related to \code{stat} argument}
\item{p.value}{p.value associated with alternative}
\item{conf.int}{confidence interval}
\item{null.value}{null parameter value}
\item{alternative}{alternative}
\item{method}{description of method}
\item{data.name}{name of x argument}
}
\references{
Fay MP, and Lumbard, K (2020). Confidence Intervals for Difference in Proportions for Matched Pairs Compatible with Exact McNemar's or Sign Tests. (unpublished manuscript).
}
\author{
Michael P. Fay
}
\note{
The sign test can be interpreted as a test that the median is zero assuming continuous data. If you want to test on the median without making continuity assumptions use \code{\link{medianTest}}.

Previous versions of \code{signTest} had \code{stat='pos-neg'} and \code{stat='prop pos'}, which are 
now referred to as \code{stat='cd'} and \code{stat='cpp'}, respectively. The old names give a warning, but may be removed in future versions.  
}
%\seealso{%\code{\link{help}}, ~~~}
\examples{
x<-c(rep(-1,10),rep(0,60),rep(1,30))
signTest(x, stat='cd')
signTest(x, stat='cpp')
signTest(x, stat='ud')
# sample median is zero, 
# and not surprisingly the median test 
# properly gives a large p-value
medianTest(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}