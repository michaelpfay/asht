\name{WprevSeSp}
\alias{WprevSeSp}
\title{Weighted prevalence inferences adjusted for sensitivity and specificity}
\description{
Prevalence inferences from weighted survey data adjusted for
sensitivity and specificity estimates as well as their variability.
}
\usage{
WprevSeSp(
  x,
  n,
  w,
  cn,
  mn,
  cp,
  mp,
  method = c("binomial", "poisson"),
  conf.level = 0.95,
  nmc = 1e+05,
  seed = 49201
)
}
\arguments{
\item{x}{integer vector of apparent positive tests for each group}

\item{n}{integer vector of number of tests for each group}

\item{w}{numeric vector of weights for each group (must sum to 1)}

\item{cn}{number of positive tests for negative controls}

\item{mn}{number of negative controls tested}

\item{cp}{number of positive tests for positive controls}

\item{mp}{number of positive controls tested}

\item{method}{either "binomial" or "poisson"}

\item{conf.level}{confidence level of the interval}

\item{nmc}{number of Monte Carlo replications}

\item{seed}{seed for random number generation}
}
\value{
A list with class "htest" containing the following components:
\item{estimate}{the adjusted prevalence estimate, adjusted for sensitivity and specificity}
\item{statistic}{the estimated sensitivity given by cp / mp}
\item{parameter}{the estimated specificity given by 1 - cn/ mn}
\item{conf.int}{a confidence interval for the prevalence}
\item{data.name}{a character string giving the unadjusted prevalence value}
\item{method}{the character string describing the output}
}
\details{
This function is for estimating prevalence from a population using a complex 
survey with an imperfect assay. The survey data are partitioned into K strata 
(or groups), and the ith stratum is weighted by w[i] and we  observe a 
proportion x[i]/n[i] that 
test positive on the assay. The weights are required to sum to 1, so that
the apparent prevalence (i.e., the prevalence 
assuming the assay is perfect) is a weighted average of the proportions of the 
strata (i.e., \code{sum(w*(x/n))}). The assay's sensitivity and specificity is  
measured by two populations, a positive control population gives a sensitivity 
estimate (\code{cp/mp}), and a negative control population gives a 
specificity estimate (\code{1 - cn/mn}). The adjusted prevalence estimate (adjusted 
for sensitivity and specificity) is the standard one described by Rogan and 
Gladen (1978). The \code{WprevSeSP} gives confidence intervals on the prevalence
that account for the sampling variability in the survey 
(by a multinomial approximation), as well as for the varibility in the 
sensitivity and specificity estimates, but it does not account for the variability 
of the weights. 

There are two methods, the 'binomial' method is based on the method of Korn and 
Graubard (1998) which is developed using a modification of exact binomial 
intervals, and the 'poisson' method which is based on the gamma method of Fay and 
Feuer (1997) which is developed assuming a weighted sum of Poisson variates. 
Simulations show that the 'poisson' method had greater than nominal coverage in all
simulated scenarios at the cost of conservativeness, while the 'binomial' method
was less conservative but had less than nominal coverage in some scenarios.

For details see Bayer, Fay and Graubard (2023). 
}
\references{
Bayer, D.M., Fay, M.P., and Graubard, B.I. (2023). ``Confidence intervals for
prevalence estimates from complex surveys with imperfect assays'' 	
arXiv:2205.13494.

Fay MP, Feuer EJ. Confidence intervals for directly standardized rates: a method 
based on the gamma distribution. Statistics in Medicine 1997; 16(7): 791–801.

Korn EL, Graubard BI. Confidence intervals for proportions with small expected number of positive counts estimated from
survey data. Survey Methodology 1998; 24: 193–201.

Rogan WJ, Gladen B. Estimating prevalence from the results of a screening test. Am J Epidemiol 1978; 107(1): 71–76.
}
\examples{
example_data_WprevSeSp <- list(
  x = c(
    53, 47, 63, 50, 54, 54, 57, 51, 66, 51, 52, 48, 37, 44, 59,
    55, 50, 58, 52, 54, 41, 45, 49, 54, 37, 53, 57, 58, 55, 55, 56,
    42, 58, 47, 49, 63, 54, 54, 54, 41, 43, 56, 44, 49, 47, 45, 62,
    53, 54, 47
  ),
  n = c(
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200
  ),
  w = c(
    0.0205962892844504, 0.0204062236737538, 0.0203843096373626,
    0.0202785701233134, 0.0202617051778543, 0.0202138087214499, 0.0201972974884707,
    0.0201818190015587, 0.0201631543739836, 0.0201560795402158, 0.0201555234250465,
    0.0201461978246263, 0.0201342022821394, 0.0201264004067009, 0.0201167314250592,
    0.0201015081093692, 0.0201003484427457, 0.0201002680000886, 0.0200817537259523,
    0.0200573433887284, 0.0200443907258367, 0.0200358187073312, 0.0200349749335002,
    0.0200264994605187, 0.0200112846914561, 0.020006219121804, 0.0199975642569458,
    0.0199649774153205, 0.0199614929059539, 0.0199426355876479, 0.0199334287088002,
    0.0199298633246975, 0.0199150015155486, 0.0199063452368827, 0.0198920051366782,
    0.0198877425787182, 0.0198679831412633, 0.0198500844815989, 0.0198381388412286,
    0.0198348595904904, 0.0198348180141822, 0.0198174510243331, 0.0197922036364436,
    0.0197821574067888, 0.0197204417557631, 0.0197004976818864, 0.019682896458092,
    0.019649677766428, 0.0196158425485035, 0.019563169292488
  ),
  cn = 77,
  cp = 58,
  mn = 300,
  mp = 60
)


WprevSeSp(
  method = "binomial",
  x = example_data_WprevSeSp$x,
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

WprevSeSp(
  method = "poisson",
  x = example_data_WprevSeSp$x,
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

}
\keyword{htest}