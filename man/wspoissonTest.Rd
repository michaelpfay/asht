\name{wspoissonTest}
\alias{wspoissonTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test and Confidence Intervals on Weighted Sum of Poissons
}
\description{
The test is not as important as the confidence intervals, which are often used for directly standardized rates. The default uses the gamma method by fay and Feuer (1997), which by all simulations appears to retain nominal coverage for any set of parameters or weights. There is a mid-p-like version that is less conservative.    
}
\usage{
wspoissonTest(x, w, nullValue = NULL, 
   alternative = c("two.sided", "less", "greater"), 
   conf.level = 0.95, midp = FALSE, nmc = 0,  
   wmtype = c("max", "mean", "minmaxavg", "tcz"), 
   mult = 1, unirootTolFactor=10^(-6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of counts (each assumed Poisson with a different parameter)
}
  \item{w}{
a vector of weights.
}
  \item{nullValue}{
a null hypothesis value of the weighted sum of the Poisson means, if NULL no test is done.
}
  \item{alternative}{
type of alternative hypothesis 
}
  \item{conf.level}{
confidence level
}
  \item{midp}{
logical, should the mid-p confidence distribution method be used}
  \item{nmc}{
Caculation method when midp=TRUE. If nmc=0 (default) does calculations that are very accurate using uniroot. If nmc>0 does Monte Carlo simulations. The Monte Carlo simulations are not needed for general use.
}
  \item{wmtype}{
type of modification for the gamma confidence interval, 'max' is the original gamma method that adds \code{max(w)} to \code{sum(x*w)} for the upper interval, 'mean' adds \code{mean(w)}, 'minmaxavg' adds \code{mean(c(min(w),max(w))}, 'tcz' does a modification of Tiwari, Clegg, and Zou (2006). 
}
  \item{mult}{
a factor to multiply the estimate and confidence intervals by, to give rates per \code{mult}
}
\item{unirootTolFactor}{
tol factor used in uniroot for calculating when midp=TRUE and nmc=0. Value multiplies by a value close to the quantile of interest in confidence interval, so that if the standardized rates are very small (e.g., 0.00001 before using mult) then the uniroot tol will be unirootTolFactor times that.
}
}
\details{
Fay and Feuer (1997) developed the gamma method (wmtype='max') for calculating confidence intervals on directly standardized rates. The assumptions is that the k by 1 vector of counts, x, are Poisson with an unknown k by 1 vector of means, theta. There are standardizing weights, w. We are interested in \code{sum(theta*w)}.

For age-standardization, \code{x} is the vector of counts of the event for each of the k age groups. The weights are \code{n.standard/(n.x *sum(n.standard)}, where \code{n.x[i]} is the  person-years associated \code{x[i]} and \code{n.standard[i]} is  person-years fro the standard population associated with the ith  age group.

Since the gamma method is conservative, Tiwari, Clegg, and Zou (2006) proposed a modification (wmtype='tcz') and also explored (wmtype='mean'). 

Ng, Filardo, and Zheng (2008) studied these and other methods (for example, wmtype='minmaxavg') through extensive simulations. They showed that the gamma method (wmtype='max') was the only method that maintained at least nominal  coverage in all the simulations. But that method is conservative. 

Fay and Kim (2017) proposed the mid-p gamma method. It appears less conservative, while appearing to retain the nominal coverage in almost all simulations. It is calculated by numeric calculations using uniroot. 
}
\value{
a list of class \code{htest}, containing:
\item{statistic}{k=length(x)}
\item{parameter}{a vector with sample variance of the calibrated weights (so sum(w)=k), and mult (only if mult !=1)}
\item{p.value}{p-value, set to NA if null.value=NULL}
\item{conf.int}{confidence interval on true directly standardized rate, sum(theta*w) }
\item{estimate}{ directly standardized rate, sum(x*w) }
\item{null.value}{null hypothesis value for true DSR}
\item{alternative}{alternative hypothesis}
\item{method}{description of method}
\item{data.name}{desciption of data}
}
\references{
Fay and Feuer (1997). "Confidence intervals for directly standardized rates: a method based on the gamma distribution."
 Statistics in Medicine. 16: 791-801.

Fay and Kim (2017). "Confidence intervals for directly standardized rates using mid-p gamma intervals." 
Biometrical Journal. 59(2): 377-387.

Ng, Filardo, and Zheng (2008). "Confidenc interval estimating procedures for standardized incidence rates." Computational Statistics and Data Analysis 52: 3501-3516.

Tiwari, Clegg, and Zou (2006). "Efficient interval estimation for age-adjusted cancer rates." Statistical Methods in Medical Research. 15: 547-569.

}
\author{Michael P. Fay}
%\note{}

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
## birth data on Down's syndrome from Michigan, 1950-1964
## see Table II  of Fay and Feuer (1997)
##xfive= counts for mothers who have had 5 or more children
## nfive and ntotal are number of live births 
xfive<-c(0,8,63,112,262,295)
nfive<-c(327,30666,123419,149919,104088,34392)
ntotal<-c(319933,931318,786511,488235,237863,61313)
## use mult =10^5 to give rates per 100,000
## gamma method of Fay and Feuer (1997) is default
wspoissonTest(xfive,ntotal/(nfive*sum(ntotal)),mult=10^5)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

