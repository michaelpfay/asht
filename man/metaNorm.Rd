\name{metaNorm}
\alias{metaNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta analysis of normally distributed parameters with assumed known variance
}
\description{
Performs either a random effects meta analysis (Paule-Mandel method or Dersimonian-Laird method) or a fixed effects meta analysis. 
}
\usage{
metaNorm(y, s2, method = c("PM", "DL", "fixed"), df = NULL, nullparm = 0,
    conf.level = 0.95, alternative = c("two.sided", "less", "greater"),  
    niter = 100, epsilon = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
vector of parameter estimates
}
  \item{s2}{
vector of variances of parameter estimates
}
  \item{method}{
either "PM" (Paule-Mandel random effects method), "DL" (Dersimonian-Laird random effects method) or "fixed" (fixed effects method)
}
  \item{df}{
degrees of freedom, NULL gives either df=k-1 (method="PM"), df=Inf (method="DL" or "fixed")
}
  \item{conf.level}{
confidence level 
}
  \item{alternative}{
type of alternative hypothesis
}
  \item{nullparm}{
null value of the parameter for calculating the p-value
}
  \item{niter}{
maximum number of iterations for method="PM" 
}
  \item{epsilon}{
small number for determining convergence of Paule-Mandel method.
}
}
\details{
Assume you have a vector of treatment effect estimates from K studies (y), together with variance estimates (s2). Assume that 
y[i] is distributed normal with mean theta[i] and variance s2[i], and assume the theta[i] (the latent treatment effect for the ith study) is normally distributed with mean theta and variance tau2 (tau^2). Assume independence between studies. 

We are interested in estimating the weighted average of the theta[i]. If tau2 is known, then an efficient estimator weighs each study proportional to the inverse of its variance, w[i] = 1/(tau2 + s2[i]). We can either assume tau2=0, and we have a fixed effects model (in other words, the treatment effect is constant across all the studies), or estimate tau2. The method for estimating tau2 either uses a simple method of moments estimator of Dersimonian and Laird (1986), or an iterative method of moments estimator of Paule and Mandel (1982). Dersimonian and Kacker (2007) give the details. 

For the Paule-Mandel estimator, to account for the fact that we are estimating tau2, we default to using a t-distribution with K-1 degrees of freedom (for motivation see Brittain, Fay and Follmann, 2012, Supplement, Section 3).

}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{a vector of [1] the estimator of tau2 and [2] the t-statistic (or Z-statistic)}
  \item{parameter}{degrees of freedom of the t-distribution (df=Inf gives a normal distribution)}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval}
  \item{estimate}{a vector of [1] the estimated weighted means and [2] the estimated standard error of the weighted means}
  \item{null.value}{the specified hypothesized value of the weighted means}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string describing the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Brittain, Fay, and Follmann (2012) A valid formulation of the analysis of noninferiority trials under random effects meta-analysis. Biostatistics 13(4): 637-649.


Dersimonian, R and Kacker, R (2007) Random-effects model for meta-analysis of clinical trials: an update. Contemporary Clinical Trials 28:105-144. 

Dersimonian, R and Laird, N. (1986). Meta-analysis in clinical trials. Controled Clinical Trials. 7:177-187.

Paule, RC and Mandel, J (1982). Consensus values and weighting factors. J Res Natl Bur Stand 87: 377-385.


}
\author{
Michael P. Fay
}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \CRANpkg{meta}   package on CRAN
}
\examples{
# Data from Table III of Teo et al, BMJ 303:1499-1503
# Effects of intravenous magnesium in suspected acute myocardial
# infarction: overview of randomised trials
# xt/nt = deaths/total in treatment group (magnesium)
# xc/nc = deaths/total in control group
xt<-c(1,9,2,1,10,1,1)
nt<-c(40,135,200,48,150,59,25)
xc<-c(2,23,7,1,8,9,3)
nc<-c(36,135,200,46,148,56,23)

rt<- xt/nt
rc<- xc/nc
logOR<- log(rt*(1-rc)/(rc*(1-rt)))
varLogOR<- 1/(nt*rt*(1-rt)) + 1/(nc*rc*(1-rc))

# Compare weighted mean and std err to Table 4 of Dersimonian and Kacker, 2007
metaNorm(logOR,varLogOR,method="PM")
metaNorm(logOR,varLogOR,method="DL")
metaNorm(logOR,varLogOR,method="fixed")
# Compare tau values to Table 3 of Dersimonian and Kacker, 2007
sqrt( metaNorm(logOR,varLogOR,method="PM")$statistic["tau squared"] )
sqrt( metaNorm(logOR,varLogOR,method="DL")$statistic["tau squared"] )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
