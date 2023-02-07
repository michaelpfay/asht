\name{var1Test}
\alias{var1Test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One Sample Test of Normal Variance
}
\description{
Give tests and confidence intervals on the variance of a sample from a normal distribution. 
}
\usage{
var1Test(x, nullVar = 1, 
   alternative = c("two.sided", "less", "greater"), 
   conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector
}
  \item{nullVar}{
null variance, or variance on the boundary between the null and alternative hypotheses
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
  \item{conf.level}{
confidence level of the interval
}
}
\details{
Tests derived from normality assumption.
}
\value{
A list of class 'htest' (use str to see elements)
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Dudewicz, EJ and Mishra, SN (1988) Modern Mathematical Statistics. Wiley. (Section 9.6).
}
\author{
Michael P. Fay
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
var1Test(rnorm(25))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
