\name{wsrTest}
\alias{wsrTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact Wilcoxon Signed Rank Test
}
\description{
Calculates the exact Wilcoxon signed rank test (using Pratt's method if there are zero values). Gives exact matching confidence intervals based on repeated calls to \code{\link{wilcoxsign_test}}, and gives associated Hodges-Lehmann estimator of center of the symmetric distribution of the difference.
}
\usage{
wsrTest(x, y = NULL, conf.int = TRUE, conf.level = 0.95, 
   mu = 0, alternative = c("two.sided", "less", "greater"),
   digits = NULL, tieDigits=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector, either the difference (if y=NULL) or the first of the paired responses (so difference is x-y).
}
  \item{y}{
second of paired differences. If NULL assumes x is the vector of paired differences.
}
  \item{conf.int}{
logica, calculate confidence interval on median of differences
}
  \item{conf.level}{
confidence level
}
  \item{mu}{
null median difference
}
  \item{alternative}{
alternative hypothesis
}
  \item{digits}{
number of digits for accuracy of confidence intervals, results are accurate to round(ci,cidigits). If digits=NULL picks about 4 digits if the range of the differences is 0 to 1, with similar accuracy as the range changes (see details). 
}
\item{tieDigits}{
number of digits to round x and y, values closer than that number of digits are treated as tied. This is to avoid rankings based on computer error.
}
}
\details{
The Wilcoxon signed rank test tests the null hypothesis of 
whether a set of values (x values, if y=NULL)
or differences (x-y, if y!=NULL) are symmetric about \code{mu}. 

This function calculates the exact Wilcoxon signed rank test using the Pratt method if there are zeros. In other words, rank the differences equal to zero together with the absolute value of the differences, but then permute the signs of only the non-zero ranks. The p-values are calculated using \code{\link{wilcoxsign_test}}, this function is just a wrapper to get confidence intervals.


When \code{conf.int=TRUE}, we get an estimator of the center of the symmetric distribution of the differences
based on the shift value where the one-sided p-values are equal (or the middle of the range if there are many values where they are equal). This type of estimator is called a Hodges-Lehmann estimator (see for example, Hodges and Lehmann, 1983).
The upper confidence limit when alternative='less' is  the
smallest shift value that gives a one-sided (alternative='less') p-value that is  
less than alpha=1-conf.level. Analogously, the lower confidence limit when alternative='greater' is the largest shift value that gives a one-sided (alternative='greater') p-value that is less than alpha. When alternative='two.sided' 
the confidence interval is the union of the two one-sided intervals 
each with level 1-alpha/2 (where alpha=1-conf.level).
Under the symmetry assumption, the center of a symmetric distribution is its median, pseudo-median, and mean.

}
\value{
An object of class 'htest', list with elements:
\item{estimate}{estimator of median difference}
\item{p.value}{p.value associated with alternative}
\item{conf.int}{confidence interval}
\item{null.value}{null median difference}
\item{alternative}{alternative}
\item{method}{description of method}
}
\references{
Pratt, JW (1959). Remarks on zeros and ties in the Wilcoxon signed rank procedures. JASA 54(287) 655-667.

Hodges, JL, and Lehmann, EL (1983). Hodges-Lehmann Estimators.
In Encyclopedian of Statistics, Volume 3. Editors S. Kotz and NL Johnson. Wiley: New York.
}
\author{
Michael P. Fay
}
\note{
The estimator and confidence interval here are different than the ones used in \code{\link[stats]{wilcox.test}}
(with paired=TRUE and exact=TRUE).

}

%% ~Make other sections like Warning  \section{Warning }{....} ~

\seealso{
\code{\link{wilcoxsign_test}}
}
\examples{
wsrTest((-3:8))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
