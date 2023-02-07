% edited from t.test.Rd
\name{bfTest}
\alias{bfTest}
\alias{bfTest.default}
\alias{bfTest.formula}
\title{Behrens-Fisher Test}
\description{
Tests for a difference in means from two normally distributed variates with possibly different variances. 
}
\usage{
bfTest(x, \dots)

\method{bfTest}{default}(x, y,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, conf.level = 0.95, control=bfControl(), \dots)

\method{bfTest}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{mu}{a number indicating the true value of the 
    difference in means }
  \item{conf.level}{confidence level of the interval.}
  \item{control}{a list of arguments used for determining the calculation algorithm, see \code{\link{bfControl}} }
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric variable giving the data values and \code{rhs} a factor
    with two levels giving the corresponding groups.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  Fisher (1935) developed a fiducial 
  solution to the two-sample difference in means problem
  with normally distributed data with different variances. That has become known as the Behrens-Fisher solution.  Robinson (1976) showed through extensive simulations, that the Behrens-Fisher solution is valid (i.e., the test gives  
  type I error rate less than the significance level, and its confidence intervals on the difference in means have coverage at least as large as the nominal confidence level). 
  
 
  
  
  The following are the same as with the usual t-test in \code{\link{t.test}}.
  \code{alternative = "greater"} is the alternative that \code{x} has a
  larger mean than \code{y}. Missing values are silently removed. If the input data are effectively constant an error is generated.
  
  
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the t-statistic.}
  \item{parameter}{R = \code{atan(SEMx/SEMy)} used in Behrens-Fisher distribution,
  where SEMx=std error of the mean of x, see \code{\link{pbf}}, but not used in calculation for this function
}
\item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the difference in means 
    (mean.x-mean.y) 
    appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated means}
  \item{null.value}{the specified hypothesized value of the mean
    difference}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string describing the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Fisher, RA (1935). The fiducial argument in statistical inference.
Annals of Eugenics. 6, 391-398.

Robinson, G (1976). Properties of Students t and of the Behrens-Fisher solution to the two means problem. The Annals of
Statistics 4, 963-971 (Corr: 1982, p. 321).

}

\seealso{
The more common solution for this problem is Welch's t-test (the default in \code{\link{t.test}}). Welch's t-test does not guarantee that the type I error rate is less than the significance level, but it appears to work well in most cases.  
}
\examples{
## Classical example: Student's sleep data
## Traditional interface
with(sleep, bfTest(extra[group == 1], extra[group == 2]))
## Formula interface
bfTest(extra ~ group, data = sleep)
## Results are simular to Welch's t-test, 
## but a little more conservative
t.test(extra~group,data=sleep)
}
\keyword{htest}
