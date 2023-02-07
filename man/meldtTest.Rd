% edited from t.test.Rd
\name{meldtTest}
\alias{meldtTest}
\title{
Meld t Test
}
\description{
Tests for a difference in parameters, when the parameter estimates are independent and both have t distributions.
}
\usage{
meldtTest(x, y, alternative = c("two.sided", "less", "greater"), delta = 0, 
    conf.level = 0.95, control = bfControl(), \dots)

}
\arguments{
  \item{x}{a list from the first group with objects: estimate (estimate of parameter),
  stderr (standard error of the estimate), and df 
  (degrees of freedom associated with t distribution)}
  \item{y}{a list from the second group with objects: estimate, stderr, and df}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{delta}{a number indicating the null hypothesis value of the 
    difference in parameters when  \code{alternative="two.sided"}. See details for one-sided hypotheses}
  \item{conf.level}{confidence level of the interval.}
  \item{control}{a list of arguments used for determining the calculation algorithm, see \code{\link{bfControl}} }
    \item{\dots}{further arguments to be passed to or from methods (currently not used)}
}
\details{
Suppose x$estimate and y$estimate estimate the parameters xParm and yParm. Let Delta=yParm-xParm. This function tests hypotheses of the form, 
\itemize{
\item alternative="two.sided" tests H0: Delta=delta versus H1: Delta != delta
\item alternative="less" tests H0: Delta >= delta versus H1: Delta< delta
\item alternative="greater" tests H0: Delta <= delta versus H1: Delta> delta
}   

The test uses the theory of melding (Fay, Proschan and Brittain, 2015). The idea is to use confidence distribution random variables (CD-RVs). It is easiest to understand the melding confidence intervals by looking at the Monte Carlo implementation. Let nmc be the number of Monte Carlo replicates, then the simulated CD-RV associated with x are
Bx = x$estimate + x$stderr * rt(nmc,df=x$df). Similarly define By. Then the 95 percent melded confidence interval for Delta=yParm-xParm is estimated by
quantile(By-Bx, probs=c(0.025,0.975)). 


When the estimates are means from normal distributions, then the meldtTest reduces to the Behrens-Fisher solution (see \code{\link{bfTest}}).  

Only one of \code{x$stderr} or \code{y$stderr} may be zero.
    
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the t-statistic.}
  \item{parameter}{R = \code{atan(x$stderr/y$stderr)} used in Behrens-Fisher distribution, see \code{\link{pbf}}
}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the difference in means appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{means and difference in means estimates}
  \item{null.value}{the specified hypothesized value of the
    difference in parameters}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string describing the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\section{Warning}{
  If the two estimates are not independent, this function may give invalid p-values and confidence intervals!
}

\references{
Fay, MP, Proschan, MA, Brittain, E (2015). Combining One-sample confidence procedures for inference in the two-sample case. Biometrics. 71: 146-156.
}
\author{
Michael P. Fay
}
\seealso{\code{\link{bfTest}} and \code{\link{pbf}}}
\examples{
## Classical example: Student's sleep data
## Compare to bfTest
xValues<- sleep$extra[sleep$group==1]
yValues<- sleep$extra[sleep$group==2]


x<-list(estimate=mean(xValues),
    stderr=sd(xValues)/sqrt(length(xValues)),
    df=length(xValues)-1)
y<-list(estimate=mean(yValues),
    stderr=sd(yValues)/sqrt(length(yValues)),
    df=length(yValues)-1)
bfTest(xValues,yValues)
# by convention the meldtTest does mean(y)-mean(x)
meldtTest(x,y)
meldtTest(y,x)
}
\keyword{htest}
