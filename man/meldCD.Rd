\name{meldCD}
\alias{meldCD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meld Two Confidence Distributions
}
\description{
Melding is a very general way of combining two independent confidence interval proceedures to create a confidence interval on a function of the two associated parameters (e.g., difference or ratio). 
}
\usage{
meldCD(H1, H2, nullparm = NULL, parmtype = c("difference", "ratio", "oddsratio"), 
    conf.level = 0.95, alternative = c("two.sided", "less", "greater"), 
    estimate = c("median", "mean"), lim = c(-Inf, Inf), parmGrid = NULL,
    nmc = 1e5, ngrid = 1e4, calcmethod = "int", epsilon=1e-8, utol=1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H1}{
a function representing the confidence distribution for parameter 1 (see details)
}
  \item{H2}{
a function representing the confidence distribution for parameter 2
}
  \item{nullparm}{
null parameter value for the parameter defined by parmtype
}
  \item{parmtype}{
parameter type, 'difference' gives parm2-parm1, 'ratio' gives parm2/parm1 (for 'oddsratio' see details).
}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{conf.level}{confidence level of the interval.}
  \item{estimate}{
type of estimate derived from each confidence distribution, either 'median' or 'mean' 
}
  \item{lim}{
a vector with limits on the parameters (both parameters should have the same limits)
}
  \item{parmGrid}{
a vector of a grid of possible values of the parameter, if NULL one is produced based on the lim argument
}
  \item{nmc}{
number of Monte Carlo replications, used if calcmethod='mc'
}
  \item{ngrid}{
minimum number of elements in the parameter grid, used if parmGrid=NULL 
}
  \item{calcmethod}{
calculation method, either 'int' (numeric integration) or 'mc' (Monte Carlo)
}
\item{epsilon}{small value for warning check, we want the minimum of the CD over the parameter grid to be greater than epsilon, and the maximum to be less than 1-epsilon}
\item{utol}{small value for passing to tol option in uniroot for confidence interval calculations}
}
\details{
For continuous responses, a confidence distribution (CD) is like a frequentist posterior distribution. 
We represent the CDs as cumulative distribution functions in the parameter space. 
The CD gets its name because it is created from the confidence interval process.
If (L,U) is the 1-alpha confidence interval for group 1, then H1(L) = alpha/2 
and H1(U)=1-alpha/2. Typically, the the CDs can be formulated as one-sided 
(alternative='greater') p-value functions, 
or 1-p for alternative='less', where the main function argument is  the boundary value on the parameter space between the null and alternative. See binomial example below.


The median of the CD can be used as an estimate of the parameter. 

We want inferences on a function of the parameters, say g(parm1, parm2), where when 
\itemize{
\item parmtype="difference" then g(parm1,parm2)=parm2-parm1
\item parmtype="ratio" then g(parm1,parm2)=parm2/parm1
\item parmtype="oddsratio" then g(parm1,parm2)=(parm2*(1-parm1))/(parm1*(1-parm2)).
}   
The function g(parm1, parm2) must be increasing in parm2 and decreasing in parm1, 
so for example normal CDs (or any with a range -Inf to Inf) are not allowed for parmtype='ratio'.
The \code{lim} argument checks to see if the parmtype is allowed.

Let T1 and T2 be simulated independent random variables associated with the CDs H1 and H2. 
Then to get a two-sided 1-alpha confidence interval on 
g(parm1,parm2) we can use \code{quantile(g(T1,T2),probs=c(alpha/2,1-alpha/2))}.
This is basically how it works when \code{calcmethod='mc'}. When \code{calcmethod='int'}
then numeric integration is used. 

For discrete responses, to ensure validity of the resulting confidence intervals, 
each group uses either a lower or upper CD, depending on the one-sided alternative. 
Thus, confidence intervals for two-sided alternatives cannot be calculated in one call 
to the \code{meldCD} for discrete data. See Fay, Proschan, and Brittain (2015) and the example. 

}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the mean appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{vector of parameter estimates for each group and using the parmtype,
     uses the median of the CDs for estimates}
  \item{null.value}{the specified hypothesized value of the
    difference in parameters}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string describing the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Fay, MP, Proschan, MA, Brittain, E (2015). Combining One-sample confidence procedures for inference in the two-sample case. Biometrics. 71: 146-156.
}
\author{
Michael P. Fay
}
%\note{}
\section{Warning }{The function has not been tested for discrete confidence distributions.
Note most confidence distributions for discrete data are not discrete CDs because the parameters are 
continuous.}
\seealso{ \code{\link{meldtTest}} and \code{\link{binomMeld.test}} for special cases.}
\examples{
x1<-4
n1<-11
x2<- 13
n2<-24


#  we use the upper and lower CDs
# this is needed for discrete data to ensure valid intervals
H1L<-function(theta){ pbeta(theta,x1,n1-x1+1)}
# Note, this is just a p-value function that inputs the null boundary value: 
binom.test(x1,n1,p=.4,alternative="greater")$p.value
H1L(.4)
H1U<-function(theta){ pbeta(theta,x1+1,n1-x1)}
# Note, but this is just a function for 1-p that inputs the null boundary value: 
1-binom.test(x1,n1,p=.4,alternative="less")$p.value
H1U(.4)
H2L<-function(theta){ pbeta(theta,x2,n2-x2+1)}
H2U<-function(theta){ pbeta(theta,x2+1,n2-x2)}

meldCD(H1U,H2L, lim=c(0,1),conf.level=0.975,alternative="greater")
meldCD(H1L,H2U, lim=c(0,1),conf.level=0.975,alternative="less")
# notice that the estimates are different than the usual 
# difference in sample proportions
require(exact2x2)
binomMeld.test(x1,n1,x2,n2, conf.level=0.975, alternative="greater")
# compare to two-.sided from 
binomMeld.test(x1,n1,x2,n2, conf.level=0.95, alternative="two.sided")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
