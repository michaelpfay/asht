\name{wmwTest}
\alias{wmwTest}
\alias{wmwTest.default}
\alias{wmwTest.formula}
\alias{wmwTest.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilcoxon-Mann-Whitney test with Confidence Interval on Mann-Whitney Parameter
}
\description{
 The \code{wmwTest} function calculates the 
Wilcoxon-Mann-Whitney test (normal approximation, exact complete enumeration, and exact Mante Carlo implementation) together with confidence intervals on the Mann-Whitney parameter, Pr[ X<Y] + 0.5 Pr[X=Y].
}
\usage{
wmwTest(x, \dots)

\method{wmwTest}{default}(x, y, alternative = c("two.sided", "less", "greater"), 
   phiNull = 0.5, exact = NULL, correct = TRUE, conf.int = TRUE, conf.level = 0.95, 
   latentContinuous = FALSE, method = NULL, methodRule = methodRuleWMW, 
   tsmethod = c("central", "abs"), control = wmwControl(),\dots)

\method{wmwTest}{formula}(formula, data, subset, na.action, \dots)

\method{wmwTest}{matrix}(x,\dots)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a (non-empty) numeric vector of data values from group 1, or a contingency table matrix with 2 rows with the top row representing group 1 and the bottom row group 2}
  \item{y}{an optional (non-empty) numeric vector of data values from group 2}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{phiNull}{
null hypothesis value for the Mann-Whitney parameter, Pr[X<Y]+0.5*Pr[X=Y].
Defaults to 0.5.
}
\item{exact}{ logical, should exact test be calculated? (see method) }
  \item{correct}{a logical indicating whether to apply continuity correction in the normal approximation for the p-value and confidence interval (when method='asymptotic')
}
\item{conf.int}{logical, should confidence intervals be calculated?}
  \item{conf.level}{confidence level of the interval.}
\item{latentContinuous}{logical, should estimates and confidence intervals be presented as 
latent continuous parameters? (see details)}
\item{method}{ character defining method, one of 'asymptotic', 'exact.ce' (exact by complete enumeration), 'exact.mc' (exact by Monte Carlo approximation). NULL value defaults to result of methodRule}
\item{methodRule}{ function that inputs x,y, and exact and outputs a method, see \code{\link{methodRuleWMW}}}
\item{tsmethod}{two-sided method, either  'central' (double the one-sided p-values) or 'abs' (test statistic uses absolute value of difference in phi estimate and phiNull)} 
\item{control}{a list of arguments for control of algorithms and output, see \code{\link{wmwControl}}}
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
The function \code{wmwTest} evaluates the Wilcoxon-Mann-Whitney test (also called the Mann-Whitney U test or the Wilcoxon rank sum test). The WMW test is a permutation two-sample rank test, and the test may be evaluated under many different sets of assumptions (Fay and Proschan, 2010). The least restrictive set of assumptions tests the null hypothesis that the two distributions of the two samples are equal versus the alternative that they are different. Unfortunately, with only those assumptions, we cannot get confidence intervals on the Mann-Whitney parameter, phi= Pr[ X<Y] + 0.5 Pr[X=Y]. In order to get confidence intervals on phi, we need additional assumptions, and for this function we use the proportional odds assumption. This assumption can be interpreted as saying that there exists some unknown monotonic transformation of the responses that leads to a location shift in a logistic distribution. This can work for discrete data (i.e., with ties allowed) if we interpret discrete responses as a grouping of some underlying latent continuous response. The proportional odds assumption is less restrictive that the assumption used in \code{wilcox.test}, which assumes a location shift on the unknown continuous distribution of the untransformed data. 


In summary, the two-sided p-value can be interpreted as testing the null that the two 
distributions are equal, and the confidence intervals on the Mann-Whitney parameter 
are intrepreted under the proportional odds assumption.
In general the confidence intervals are compatible with the associated p-values,
for details see Fay and Malinovsky (2018).

There is a choice of three methods. When \code{method='asymptotic'}, the test is implemented using a normal approximation, with (\code{correct=TRUE}) or without 
(\code{correct=FALSE}) a continuity correction. The resulting p-values should 
match \code{\link{wilcox.test}} (when 
\code{paired=FALSE}  and \code{exact=FALSE}). When \code{method='exact.ce'}, the test is implemented using complete enumeration of all permutations, and hence is only tractible for very small sample sizes (less than 10 in each group). When \code{method='exact.mc'}, the test is implemented using Monte Carlo with B=10^4 replications (change B with \code{control=controlWMW(nMC=B)}). As B gets larger the p-value approaches the exact one. (See 'note' section, sometimes the \code{method='exact.mc'} will not work.)    

The  \code{tsmethod='central'} gives two-sided p-value that is equal to 
\code{min(1,min(2*pless,2*pgreater))}. 
Alternatively,  \code{tsmethod='abs'} gives the two-sided method, 
which is based on the test statistic  |phi - phiNull|.
Under the proportional odds assumption,  \code{tsmethod='central'} allows us to interpret p.value/2 as one-sided p-values (this is not allowed using \code{tsmethod='abs'}). 
With continuous data, the p-values will be the same, but with ties they can be different.

 
From the two groups x (or top row of contringency table, or first factor in rhs of formula) and y (or bottom row of contingency table, or second factor in rhs of formula) the Mann-Whitney parameter represents Pr[X<Y]+0.5Pr[X=Y]. It is also the area under the curve of an ROC curve (see Hanley and McNeil, 1982). The confidence interval when \code{method='asymptotic'} generalizes the Method 5 of Newcombe (2006), which was a 
score-type modification of the interval of Hanley and McNeil (1982). The generalization is that the confidence interval 
adjusts for ties and allows a continuity correction (see examples below).


The \code{methodRule} function allows automatic choice of the method of calculation based on the data and the \code{exact} argument. 

When the data are discrete, we can treat the data as if they are a grouping of some underlying continuous responses. Using the proportional odds assumption, we can then translate the Mann-Whitney parameter on the observed discrete data into the Mann-Whitney parameter on the latent continuous data (when \code{latentContinuous=TRUE} and using the default 
\code{control=controlWMW(latentOutput='mw')}). You can also translate the results into the proportional odds parameter on the latent continuous responses (when \code{latentContinuous=TRUE} and using  
\code{control=controlWMW(latentOutput='po')}). Translation is done with \code{\link{latentTransform}}.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{U statistic estimate of the Mann-Whitney parameter.}
  \item{parameter}{tie factor}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the Mann-Whitney parameter 
  appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated difference in means}
  \item{null.value}{the specified hypothesized value of the mean
    difference}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string describing the test.}
  \item{data.name}{a character string giving the name(s) of the data.}

}
\references{
Fay, MP and Malinovsky, Y (2018). Confidence Intervals of the Mann-Whitney Parameter that are Compatible with the Wilcoxon-Mann-Whitney Test. Statistics in Medicine:
DOI: 10.1002/sim.7890.

Fay, MP and Proschan MA (2010). Wilcoxon-Mann-Whitney of t-test? On assumptions for hypothesis tests and multiple interpretations of decision rules. Statistics Surveys 4:1-39.

Hanley, JA, and McNeil, BJ (1982). The Meaning and Use of the Area under a Receiver Operating Characteristic (ROC) Curve. Radiology 143: 29-36.

Newcombe, Robert G. (2006). Confidence intervals for an effect size measure based on the Mann-Whitney statistic. Part 2: asymptotic methods and evaluation. Statistics in medicine 25(4): 559-573.
}
%\author{%%  ~~who you are~~}
\note{
The \code{method='exact.mc'} can sometimes fail. The issue is that for some Monte Carlo simulations the one-sided p-value function is not monotonic, even in for data sets where the one-sided p-value would be monotonic if we could do complete enumeration. In this case, the confidence limit will be set to \code{NA} and a warning will suggest trying \code{method='asymptotic'}
or \code{method='exact.ce'} if feasible. Here is an example where that occurs: 
\code{set.seed(1);} \code{g<- c(rep(0,6),1,rep(0,4),1,rep(0,3),1,1,0,1,1,0,rep(1,5));}
\code{y<-1:26;} \code{wmwTest(y~g,exact=TRUE)}.
}

\section{Warning }{
The algorithm for calculating the confidence interval when \code{tsmethod='abs'} is not guaranteed to give the correct value. It is possible to skip over a value. For more accurate results 
increase \code{control=wmwControl(rcheckgrid)}
and \code{control=wmwControl(ncheckgrid)}
} 
\seealso{
See \code{\link{wilcox.test}} for either exact p-value or the same asymptotic p-value and confidence interval on location shift under the shift assumption. 

See \code{\link{wilcox_test}} for exact p-value and exact confidence interval on location shift.
}
\examples{
# data from Table 1 of Hanley and McNeil (also given in Table 1 of Newcombe, 2006)
HMdata<-matrix(c(33,3,6,2,6,2,11,11,2,33),nrow=2,dimnames=
   list(c("Normal","Abnormal"),
   c("Definitely Normal",
   "Probably Normal",
   "Questionable",
   "Probably Abnormal",
   "Definitely Abnormal")))
HMdata
# to match Newcombe (2006, Table 1, Method 5) exactly 
# use correct=FALSE and RemoveTeAdjustment=TRUE
wmwTest(HMdata, correct=FALSE, RemoveTieAdjustment=TRUE)
# generally smaller intervals with closer to nominal coverage with 
# tie adjustment and continuity correction
wmwTest(HMdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
