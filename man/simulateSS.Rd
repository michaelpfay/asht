\name{simulateSS}
\alias{simulateSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate sample sizes
}
\description{
A function that simulates sample sizes in an efficient manner. Inputs two functions: (1) a decision function that returns 1=reject, or 0=fail to reject, and (2) a data generating function.
}
\usage{
simulateSS(decFunc, dataGenFunc, nstart = 100, numBatches = 100, repsPerBatch = 100, 
   power = 0.3, alpha = 0.025, nrepeatSwitch = 3, printSteps = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{decFunc}{
decision function, inputs data from dataGenFunc and outputs 1 (reject) or 0 (fail to reject).
}
  \item{dataGenFunc}{
data generating function, inputs a sample size and outputs simulated data object. Class of the data must match input for decFunc. 
}
  \item{nstart}{
starting sample size value
}
  \item{numBatches}{
number of batches (must be at least 5), default=100
}
  \item{repsPerBatch}{
number of replications per batch (must be at least 10), default=100
}
  \item{power}{
power desired
}
  \item{alpha}{
one-sided alpha level, used for estimating power from batches by normal approximation
}
  \item{nrepeatSwitch}{
one of 2,3,4 or 5. default=3. If nrepeatSwitch batch estimates of sample size are the same in a row, then switch to an up-and-down method (adding or subtracting 1 to sample size).
}
  \item{printSteps}{
logical, print intermediate steps of algorithm?
}
}
\details{
This is an algorithm proposed in Fay and Brittain (2022, Chapter 20). Here are the details of the algorithm. For step 1, we pick a starting sample size, say $N_1$, and the number of replications within a batch, $m$,
and the total number of batches, $b_{tot}$.
We simulate $m$ data sets with sample size $N_1$, and get the proportion of rejections, say $P_1$.
Then we use a normal approximation to estimate the target sample size, say $N_{norm}$. In step 2, we replicate $m$ data sets with sample size $N_2 = N_{norm}$
to get the associated proportion of rejections, say $P_2$. We repeat 2 more batches with $N_3=N_{norm}/2$ and $N_4=2 N_{norm}$,
to get proportions $P_3$, and $P_4$. Then in step 3, we use isotonic regression (which forces monotonicity, power to be non-decreasing with sample size) on the 4 observed pairs ($(N_1,P_1),\ldots,(N_4,P_4)$), and linear interpolation to get our best estimate of
the sample size at the target power, $N_{target}$. We use that estimate of $N_{target}$ for our sample size for the next
batch of simulations. This idea is of using the best estimate of the target for the next iteration is studied in
Wu (1985, see Section 3).   Step 4 is iterative, for the $i$th batch we repeat the isotonic regression, except now with $N_i$ estimated from the  first $(i-1)$ observation pairs.  We repeat step 4 until either the number of batches is $b_{tot}$,
or the current sample size estimate is the same as the last nrepeatSwitch-1 estimates, in which case we switch to an up-and-down-like method.
For each iteration of the up-and-down-like method, if the current proportion of rejections from the last batch of $m$ replicates is greater than the target power, then subtract 1 from the
current sample size estimate, otherwise add 1. Continue with that up-and-down-like method until we reach the number of batches equal to $b_{tot}$. The up-and-down-like method was added because sometimes the algorithm would get stuck in too large of a sample size estimate.
}
\value{
A list with elements:
\item{N}{vector of estimated sample sizes at the end of each batch}
\item{P}{vector of power estimates at the end of each batch}
\item{Nstar}{estimate of sample size, not necessarily an integer}
\item{Nestimate}{integer estimate of sample size equal to ceiling(Nstar)}
}
\references{
Fay, M.P. and Brittain, E.H. (2022). Statistical Hypothesis Testing in Context. Cambridge University Press. New York. 

Wu, CJ (1985). Efficient sequential designs with binary data. Journal of the American Statistical Association. 19: 1085-1098.
}
\author{
Michael P. Fay
}
%\note{%%  ~~further notes~~}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
# simple example to show method
# simulate 2-sample t-test power
# for this simple case, better to just use power.t.test
power.t.test(delta=.5,sig.level=0.025,power=.8,
   type="two.sample",alternative="one.sided")
decFunc<-function(d){
  ifelse(t.test(d$y1,d$y2,alternative="less")$p.value<=0.025,1,0)
}
dataGenFunc<-function(n){
  list(y1=rnorm(n,0),y2=rnorm(n,.5))
}
# for example use on 20 batches with 20 per batch
set.seed(1)
simulateSS(decFunc,dataGenFunc,nstart=100,numBatches=100,repsPerBatch=100,
   power=0.80, alpha=0.025,printSteps=FALSE)
}
\keyword{htest}