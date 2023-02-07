\name{cvTest}
\alias{cvTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coefficient of Variation Test
}
\description{
One-sample coefficient of variation tests and confidence intervals based on either normal or lognormal assumptions.
}
\usage{
cvTest(x, nullCV = 1, 
alternative = c("two.sided", "less", "greater"), 
conf.level = 0.95, distn = c("normal", "lognormal"), 
CVmax = 10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector
}
  \item{nullCV}{
null coefficient of variation, or CV on boundary between null and alternative hypotheses
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
  \item{conf.level}{
confidence level of the interval
}
  \item{distn}{
assumed distribution
}
  \item{CVmax}{
maximum coefficient of variation used in uniroot CI searches when distn='normal'
}
}
%\details{}
\value{
A list of class 'htest'
\item{statistic}{mean}
\item{parameter}{stadard deviation}
\item{estimate}{estimate of coefficient of variation: sd(x)/mean(x) for distn='normal', and sqrt(exp(var(log(x)))-1) for distn='lognormal'}
\item{p.value}{p.value associated with alternative}
\item{conf.int}{confidence interval}
\item{null.value}{null CV}
\item{alternative}{alternative}
\item{method}{description of method}
}
\references{
Koopmans, Owen, Rosenblatt (1964) "Confidence intervals for the coefficient of variation for the normal and log normal distributions" Biometrika 25-32.
}
\author{
Michael P. Fay
}
%\note{}

%% ~other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~See Also as \code{\link{help}}, ~~~}
\examples{
cvTest(rnorm(25,mean=3,sd=.2),distn="normal")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
