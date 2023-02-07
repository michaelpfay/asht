\name{wmwControl}
\alias{wmwControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Arguments passed to wmwTest.
}
\description{
Creates a list of arguments that are used for algorithm 
control and output control.
}
\usage{
wmwControl(nMC = 10^4, epsilon = 10^(-8), 
   digits = 10, latentOutput = c("mw", "po"),
   removeTieAdjustment = FALSE, ncheckgrid=100,
   rcheckgrid=0.1, Vmethod="LAPH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nMC}{
number of Monte Carlo replications
}
  \item{epsilon}{
small number, used in uniroot for limits. Rarely needs to be changed.
}
  \item{digits}{
number of digits for rounding. This is needed so that true ties are treated as ties by the computer, otherwise there can be non-trival errors in the p.values 
}
  \item{latentOutput}{either 'mw' (Mann-Whitney parameter) or 'po' (proportional odds parameter).
Describe output of parameter, estimate, and confidence intervals when \code{latentContinuous=TRUE}.
}
\item{removeTieAdjustment}{logical, should tie adjustment be removed? If there are no ties TRUE and FALSE give the same answer. \code{removeTieAdjustment=TRUE} is not recommended, only use TRUE when used with \code{correct=FALSE} to reproduce method 5 of Newcombe (2006), see example.
}
\item{ncheckgrid}{number of elements in the grid for checking 
after the uniroot call when \code{tsmethod='abs'}
}
\item{rcheckgrid}{range of how far back or forward to look in the checking after the uniroot call when \code{tsmethod='abs'}
}
\item{Vmethod}{character to determine the variance function for the asymptotic method. 'LAPH' gives the combination Lehmann alternative and proportional hazards method, and 'PO' gives the proportional odds method \code{\link{Vpo}}. The methods are nearly identical, but the 'PO' method takes longer, so the default is 'LAPH'. 
}
}
%\details{%% }
\value{
A list with values of the same name as the arguments of the function.
}
%\references{}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%%  \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
