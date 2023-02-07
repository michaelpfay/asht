\name{asht-package}
\alias{asht-package}
\alias{asht}
\docType{package}
\title{
Applied Statistical Hypothesis Tests
}
\description{
Test and confidence intervals for some applied statistical hypothesis tests.
}
\details{
\tabular{ll}{
Package: \tab asht\cr
Type: \tab Package\cr
Version: \tab 0.9.8\cr
Date: \tab 2022-09-15\cr
License: \tab GPL-3\cr
}
A collection of statistical hypothesis tests, with a focus on non-asymptotic tests. Some tests are \code{\link{medianTest}} for exact tests and
confidence intervals about a median, \code{\link{quantileTest}} which generalizes \code{\link{medianTest}} for other quantiles besides the median,  \code{\link{signTest}} to run the exact sign test, \code{\link{bfTest}} to run the Behrens-Fisher test,
\code{\link{abcnonHtest}} to calculate ABC intervals and tests,
\code{\link{wmwTest}} to run the Wilcoxon-Mann-Whitney test (i.e., Wilcoxon rank sum test, or Mann-Whitney U test) and calculate confidence intervals on the Mann-Whitney parameter. In rare cases, the function only gives a confidence interval and and estimate and does not test a specific hypothesis (see \code{\link{prevSeSp}} which estimates prevalence accounting for sensitivity and specificity).
}
\author{
Michael P. Fay

Maintainer: Michael P. Fay <mfay@niaid.nih.gov>
}
%\references{}
\keyword{ package }
\keyword{htest}
%\seealso{~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~}
%\examples{~~ simple examples of the most important functions ~~}
