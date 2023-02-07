\name{prevSeSp}
\alias{prevSeSp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate prevalence with confidence interval accounting for sensitivity and specificity}
\description{
Using the method of Lang and Reiczigel (2014), estimate prevalence and get a confidence interval 
adjusting for the sensitivity and specificity (including accounting for the variability of the sensitivity and specificity estimates).
}
\usage{
prevSeSp(AP, nP, Se, nSe, Sp, nSp, conf.level = 0.95, neg.to.zero=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AP}{
apparent prevalence (proportion positive by test)
}
  \item{nP}{
number tested for AP
}
  \item{Se}{
estimated sensitivity (true positive rate)
}
  \item{nSe}{
number of positive controls used to estimate sensitivity
}
  \item{Sp}{
estimated specificity (1- false positive rate)
}
  \item{nSp}{
number of negative controls used to estimate specificity
}
  \item{conf.level}{
confidence level
}
\item{neg.to.zero}{
logical, should negative prevalence estimates and lower confidence limits be set to zero?
}
}
\details{
When measuring the prevalence of some disease in a population, it is useful to adjust for the fact that the test for the disease may not be perfect. 
We adjust the apparent prevalence (the proportion of people tested positive) for the sensitivity (true positive rate: proportion of the population that
has the disease that tests positive) and the specificity (1-false positive rate: proportion of the  population that do not have the disease that tests negative).
So if the true prevalence is \eqn{\theta} and the true sensitivity and specificity are \eqn{Se} and \eqn{Sp}, then the expected value of the apparent prevalence is 
the sum of the expected proportion of true positive results and the expected proportion of false positive results:  
\deqn{AP = \theta Se + (1-Sp) (1-\theta).} Plugging in the estimates (and using the same notation for the estimates as the true values) 
and solving for \eqn{\theta} we get the estimate of prevalence of 
\deqn{\theta = \frac{AP - (1-Sp)}{Se -(1-Sp)}.}

Lang and Reiczigel (2014) developed an approximate confidence interval for the prevalence that not only adjusts for the sensitivity and specificity, but also 
adjusts for the fact that the sensitivity is estimated from a sample of true positive individuals (\code{nSe}) and
the specificity is estimate from a sample of true negative individuals (\code{nSp}).


If the estimated false positive rate (1-specificity) is larger than the apparent prevalence, the prevalence estimate will be negative. This occurs because we observe
a smaller proportion of positive results than we would expect from a population known not to have the disease. The lower confidence limit can also be negative because of the variability in the specificity estimate. The default with \code{neg.to.zero=TRUE} sets those negative estimates and lower confidence limits to zero. 


The Lang-Reiczigel method uses an idea discussed in Agresti and Coull (1998) to get approximate confidence intervals. 
For 95\% confidence intervals, the idea is similar to adding 2 positive and 2 negative individuals to the apparent prevalence results, 
and adding 1 positive and 1 negative individual to the sensitivity and specificity test results, then using asymptotic normality. Simulations 
in Lang and Reiczigel (2014) show the method works well for true sensitivity and specificity each in ranges from 70\% to over 90\%.

}
\value{
   A list with class \code{"htest"} containing the following components:
  %\item{p.value}{the p-value of the test.}
  \item{estimate}{the adjusted prevalence estimate, adjusted for sensitivity and specificity}
  \item{statistic}{the estimated sensitivity given by \code{Se}}
  \item{parameter}{the estimated specificity given by \code{Sp}}
  \item{conf.int}{a confidence interval for the prevalence.}
  %\item{null.value}{}
  %\item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{the character string describing the output.}
  \item{data.name}{a character string giving the unadjusted prevalence value and the sample size used to estimate it (nP).}
}
\references{
Agresti, A., Coull, B.A., 1998. Approximate is better than 'exact'for interval estimation of binomial proportions. Am. Stat. 52,119-126.

Lang, Z. and Reiczigel, J., 2014. Confidence limits for prevalence of disease adjusted for estimated sensitivity and specificity. Preventive veterinary medicine, 113(1), pp.13-22.

}
\author{
Michael P. Fay
}
\note{There is a typo in equation 4 of Lang and Reiczigel (2014), the \eqn{(1+\hat{P})^2} should be \eqn{(1-\hat{P})^2}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{truePrev} in package \CRANpkg{prevalence} for Bayesian methods for this problem (but this requires JAGS (Just Another Gibbs Sampler), a separate software that can be called from R if it is installed on the user's system.)
}
\examples{
# Example 1 of Lang and Reiczigel, 2014
# 95% CI should be 0.349, 0.372
prevSeSp(AP=4060/11284,nP=11284,Se=178/179,nSe=179,Sp=358/359, nSp=359)

# Example 2 of Lang and Reiczigel, 2014
# 95% CI should be  0, 0.053
prevSeSp(AP=51/2971,nP=2971,Se=32/33,nSe=33,Sp=20/20, nSp=20)

# Example 3 of Lang and Reiczigel, 2014
# 95% CI should be 0 and 0.147
prevSeSp(AP=0.06,nP=11862,Se=0.80,nSe=10,Sp=1, nSp=12)

# Example 4 of Lang and Reiczigel, 2014
# 95% CI should be 0.58 to 0.87
prevSeSp(AP=259/509,nP=509,Se=84/127,nSe=127,Sp=96/109, nSp=109)
# 95% CI should be 0.037 to 0.195
prevSeSp(AP=51/509,nP=509,Se=23/41,nSe=41,Sp=187/195, nSp=195)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
