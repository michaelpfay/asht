\name{pbf}
\alias{pbf}
\alias{qbf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Behrens-Fisher distribution
}
\description{
Used in \code{\link{bfTest}}. Distribution of the t-statistic used in 
Welch's t-test. The distribution depends on the sample sizes of the two groups, and the ratio of the two standard errors of the means.
}
\usage{
pbf(q, n1, n2, R = NULL, s1 = NULL, s2 = NULL, epsilon = 10^(-8))
qbf(p,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
vector of quantiles
}
\item{p}{
vector of probabilities
}
  \item{n1}{
sample size in group 1}
  \item{n2}{
sample size in group 2}
  \item{R}{
a function of the ratio of the two standard errors of the mean. Specifically, 
\code{atan((s1/sqrt(n1))/(s2/sqrt(n2)))}
  }
  \item{s1}{
sample standard deviation in group 1
}
  \item{s2}{
sample standard deviation in group 2
}
  \item{epsilon}{
a small positive number used to avoid computer errors}
}
\details{
The user must supply either \code{s1} and \code{s2} or supply \code{R}.
Suppose m1 and m2 are the means of the two groups and D is the true difference in means. Then the Behrens-Fisher test statistic is 
\code{T=(m1-m2-D)/sqrt(s1/sqrt(n1)+ s2/sqrt(n2))}. The value T can be equivalently written as \code{T=T1*sin(R)+T2*cos(R)}, where T1 and T2 are t random variables with n1-1 and n2-1 degrees of freedom. The cumulative distribution of T is found by numeric integration.

We rewrite \eqn{Pr[T<=q]} as \eqn{ \int_{-Inf}^{Inf} Pr[T2<= (q-u*sin(R))/cos(R) | T1=u] Pr[T1=u] du}.


}
\value{
pbf gives the distribution function, qbf givds the quantile function.
}
\references{
Kim, S-H, and Cohen, AS (1996). Table of percentage points of the Behrens-Fisher distribution. Journal of Statistical Computation and Simulation. 55(3) 181-187.
}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{
# See Table 1 from Kim and Cohen, 1996
# at v1=n1-1=8 and v2=n2-1=12 with 45 degrees = 45*pi/180 radians
# for 0.95th percentile
# Table gives: 1.83496
qbf(0.95,9,13,45*pi/180)
# check Inf degrees of freedom, should give qnorm value
qbf(.95,Inf,Inf,45*pi/180)
qnorm(.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
