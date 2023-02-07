\name{tukeyWelsch}
\alias{tukeyWelsch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tukey-Welsch Pairwise Tests
}
\description{
Calculate pairwise comparisons between groups levels using step down correction for multiple testing.
}
\usage{
tukeyWelsch(y, g, method = c("aov", "kw", "sr", "user"), 
    pvalfunc = NULL, padjfunc = padjTW, maxnTest=10^4,nTestMessage=FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response vector
}
  \item{g}{
grouping vector or factor
}
  \item{method}{
type of method for tests, one of 'aov' (ANOVA which is a t-test for the pairwise comparisons) 'kw' (Kruskal-Wallis test, which is a Wilcoxon-Mann-Whitney test for the pairwise comparisons), 
'sr' (studentized range test), or 'user' (user supplied function, see details).
}
  \item{pvalfunc}{
function to test for effects and return a p-value. Used if method='user'. See details. 
}
  \item{padjfunc}{
function that takes the unadjusted p-value vector from a stage, and returns the adjusted p-value vector for that stage (see details)
}
\item{maxnTest}{maximum number of tests, if the number of tests is larger than \code{maxnTest} then
gives an error}
\item{nTestMessage}{
logical, print a message at the start of calculations telling how many tests will be calculated
}
  \item{\dots}{
additional arguments to pass to XXX (if method='aov'), YYY (if method='kw')
or pvalfunc (if method='user')
}
}
\details{

This function does a k-sample test (either one-way ANOVA [method='aov'] or Kruskal-Wallis test [method='kw']) on the responses when the g vector has k levels. Then it does all the pairwise comparisons (either t.tests [method='aov'] or Wilcoxon-Mann-Whitney tests [method='kw']) giving multiple 
comparison adjusted p-values. The adjustment uses a step-down method, that is different from (and potentially more powerful than) the single step procedures in \code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}}. The method is described in Einot and Gabriel (1975) [for the anova case] and Campbell and Skillings (1985) for the Kruskal-Wallis case. See also Hochberg and Tamhane (1987, p. 111 for 'aov' case, and p. 247-248 for the 'kw' case). 


Here are the details. First, the k-sample test is done, where the type of test is determined  by the \code{method}. The function repeats that type of test k-1 times, leaving out a different level of the group each time. These are k-1 tests, each having k-1 levels. This process repeats itself (i.e., do choose(k,k-2) tests each having k-2 levels, then do choose(k,k-3) tests each having k-3 levels, etc) until we get to the choose(k,2) pairwise tests. Reject at level aj = 1- (1-alpha)^(j/k), for all tests where there are j groups, for j=2,..,k-2 and at level aj=alpha for j=k-1 and k. These adjusted significance levels are known as the Tukey-Welch (see Hochberg and Tamhane, p. 111) or Ryan (see Einot and Gabriel, 1975) levels. Then we only reject each pairwise comparison, if we reject at all null hypotheses that contain that pair.

We convert this procedure into adjusted p-values, by finding the lowest alpha level such that each pairwise comparison would be rejected, that is its adjusted p-value. The \code{padjfunc} is a function that takes the unadjusted p-values 
and gives the adjustment for each level by itself. For example, the default uses the Tukey-Welch adjusted significance levels, and the function solves for alpha as a function of aj (i.e.,  inputs unadjP and returns either  1-(1-unadjP^(k/j) for j=2,3,...k-2 or unadjP for j=k-1 or k)). Then taking the individual level adjusted p-values, we define the step-down adjusted p-value for each pairwise comparison as the maximum of all the individual level adjusted p-values for each hypothesis that contains the pair as part of its groups tested. 

When k=3, this method gives an adjusted p-value for each pairwise comparison that is the maximum of the k-sample test p-value and the unadjusted p-value for the two-sample test using that pair of levels.


When method='user' the function uses the \code{pvalfunc} function to test for p-values. The function must input y and g and output the p-value for the j-sample test, where j is the number of levels present in g. So the function must be defined when j=2,3,...,k. 

}
\value{
An object of class 'tukeyWelsch', a list with elements:
\item{fullResults}{a list of all the intermediate p-values (unadjusted and adjusted). Not printed by default}
\item{method}{description of method}
\item{data.name}{description of input data}
\item{ksample.pvalue}{p-value for k-sample test}
\item{pairwise.pvalues}{vector of adjusted p-values for pairwise comparisons}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Campbell and Skillings (1985) JASA 998-1003.

Einot and Gabriel (1975) JASA 574-583.

Hochberg, Y and Tamhane, AC (1987) Multiple Comparison Procedures. Wiley: New York.
}
\author{
Michael P. Fay
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{pairwise.wilcox.test}} and   \code{\link{pairwise.t.test}}
}
\examples{
##
createData<-function(n,props,shifts,ry=rnorm){
  k<-length(props)
  if (round(sum(props),8)!=1) stop("sum of props must be 1")
  props<- props/sum(props)
  if (length(shifts)!=k) stop("length of shifts must equal length of props")
  g<-rep(1:k,as.vector(rmultinom(1,n,props)))
  y<-ry(n)
  for (i in 1:k){
    y[g==i]<-y[g==i]+shifts[i]
    }
  list(y=y,g=g)
}
set.seed(1)
d<-createData(100,c(.2,.3,.2,.3),c(0,0,0,1))
tukeyWelsch(d$y,factor(d$g),method="kw")
tukeyWelsch(d$y,factor(d$g),method="aov")
tukeyWelsch(d$y,factor(d$g),method="sr")
TukeyHSD(aov(d$y~factor(d$g)))[[1]][,"p adj"]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
