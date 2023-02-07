\name{ama1c1cpg}
\alias{ama1c1cpg}
\docType{data}
\title{
Three arm phase 1 malaria vaccine trial
}
\description{
Growth inhibition responses from a three arm vaccine trial (Mullen, et al, 2008).
}
\usage{data("ama1c1cpg")}
\format{
  A data frame with 58 observations on the following 2 variables.
  \describe{
    \item{\code{vacc}}{a factor representing the three arms of the trial. The levels are: \code{20ug+CPG} \code{80ug} \code{80ug+CPG}}
    \item{\code{resp}}{a numeric vector giving the response: day 70 sera percent in vitro growth inhibition of the 3D7 malaria parasite.}
  }
}
%\details{}
%\source{}
\references{
Mullen, GE, Ellis, RD, Miura, K, Malkin, E, Nolan, C, Han, M, Fay, MP, Saul, A, Zhu, D, Rausch, K, Moretz, S, Shou, H, Long, CA, Miller, LH, Treanor, J. 2008. Phase 1 trail of ama1-c1/alhydrogel plus cpg 7909: an asexual blood-stage vaccine for plasmodium falciparum malaria. PLoS ONE. 
3(8):32940.
}
\examples{
data(ama1c1cpg)
## maybe str(ama1c1cpg) ; plot(ama1c1cpg) ...
}
\keyword{datasets}
