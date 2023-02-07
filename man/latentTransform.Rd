\name{latentTransform}
\alias{latentTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform  Mann-Whiteny parameter to latent Mann-Whitney parameter
}
\description{
Take data and Mann-Whitney parameter values and transform them to latent continuous Mann-Whitney parameters.
}
\usage{
latentTransform(x, y, phiValues, output = c("mw", "po"), epsilon = 10^(-6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector of responses from group 1
}
  \item{y}{
numeric vector of responses from group 2
}
  \item{phiValues}{
vector of Mann-Whitney parameters to be transformed
}
  \item{output}{
either 'mw'  (to outpout latnet continuous Mann-Whitney parameter) or 'po' (to output proportional odds parameter)
}
  \item{epsilon}{
small value for limits of proportional odds parameter, smallest is epsilon and largest is 1/epsilon
}
}
\details{
Calculates the empirical distribution from both groups and uses that and the sample sizes to transform the \code{phiValues} to their associated latent continuous Mann-Whitney parameter. Extreme values may not be transformed by this method and are replaced by 0 or 1, where those extreme values are data dependent (and depend on \code{epsilon}). See the code and its comments for details.
}
\value{
a vector of latent continuous Mann-Whitney parameters (when output='mw') or proportional odds parameters (when output='po')
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
\keyword{internal}% use one of  RShowDoc("KEYWORDS")

