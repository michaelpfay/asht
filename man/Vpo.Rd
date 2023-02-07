\name{Vpo}
\alias{Vpo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance for estimated Mann-Whitney parameter under proportional odds.
}
\description{
A function to calculate the variance of the estimated Mann-Whitney parameter under the proportional odds model.
}
\usage{
Vpo(PHI, tf, ny, nx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PHI}{
Mann-Whitney parameter
}
  \item{tf}{
tie factor
}
  \item{ny}{
sample size for the y variable
}
  \item{nx}{
sample size for the x variable  
}
}
\details{
Uses uniroot to find the variance numerically. The approximation using the LAPH method (Lehmann alternative model combined with the proportional hazards model) is much faster and is the default.  
}
\value{
the variance of the estimated Mann-Whitney parameter
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
Vpo(.7,.8,5,7)
## Compare to LAPH method
VLAPH<- function(PHI,tf=tiefactor,ny=n.y,nx=n.x){
        tf*(PHI*(1-PHI)/(ny*nx))*
          (1+((ny+nx-2)/2)*
             ((1-PHI)/(2-PHI)+PHI/(1+PHI)))
      }   
VLAPH(.7,.8,5,7)      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
