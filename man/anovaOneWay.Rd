\name{anovaOneWay}
\alias{anovaOneWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Way ANOVA
}
\description{
Do one-way ANOVA with estimates and confidence intervals on 
parameters. The parameter is called tau.sq and is the 
weighted sum of the square of the difference between 
the true means and the weighted average of the true means. 
Allows var.equal=FALSE using the Brown-Forsythe method 
that generalizes Welch's t-test to the k-sample problem.
}
\usage{
anovaOneWay(y, g, var.equal = TRUE, nullValue = 0, 
    parm =c("ICC","varb"), conf.level = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
numeric vector of responses
}
  \item{g}{
group membership vector (may be numeric, character, or factor) 
}
  \item{var.equal}{
logical, are the variances for all groups be equal? TRUE gives usual anova, FALSE gives Brown-Forsythe method.
}
  \item{nullValue}{
null value of tau.square (between group variance) or tau.sq/sigma.sq (must be 0 now)
}
  \item{parm}{
type of parameter, either 'ICC' (the parameter that R square estimates for this problem) or 
'varb' (the between group variance).
}
  \item{conf.level}{
confidence level for the confidence interval. Default is 0.90 so that when the p-value<0.05, the two-sided confidence interval will exclude 0.
}
}
\details{
The typical way to get the p-value for a one-way anova is 
\code{anova(lm(y~g))}. This function was written to add two new features. 

First, using the method of Brown and Forsythe (1974a), the function allows for non-equal variances between the groups. This is one generalization of Welch's t-test to the one-way ANOVA case. Brown and Forsythe (1974b) give simulations showing 
 that the type I error rate is close to the nominal (under the nomrality assumption with different variances).

 Second, the function gives confidence intervals on either 'ICC' or 'varb'.
The 'varb' (the between-group variance) is \code{sum((na/n)*(ua-u)^2)} where na is a vector of length k giving the sample size in each group, n is the total sample size, and ua is a vector of the k means in the groups, and u is the overall mean. Let varw be the within-group variance, then ICC=varb/(varb+varw).
ICC is the intraclass correlation coefficient, and in this situation it is the 
parameter that the R square is estimating.

}
\value{
A object of class 'htest'.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Brown and Forsythe (1974a). Biometrics 30:719-724.

Brown and Forsythe (1974b). Technometrics 16: 129-132.
}
\author{
Michael P. Fay
}
\note{

Note also that it is possible to get a 90 pct confidence interval for varb that is (0,0). This occurs when the group means are much closer to each other than they would be expected to be by chance, given the observed variability between observations within the groups. 


%If the variances are very different between the groups, it is not clear how %well the confidence interval will work. 

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
require(datasets)
library(asht)
ChickWeightTime20<-ChickWeight[ChickWeight$Time==20,]

anovaOneWay(1:10,c(rep(1,4),rep(2,6)))
anova(lm(weight~Diet,data=ChickWeightTime20))
t.test(ChickWeightTime20$weight[ChickWeightTime20$Diet==1],
       ChickWeightTime20$weight[ChickWeightTime20$Diet==2],
       var.equal=FALSE)
anovaOneWay(ChickWeightTime20$weight, ChickWeightTime20$Diet,       
     var.equal=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
