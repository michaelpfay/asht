meldtTest<-function (x, y, alternative = c("two.sided", "less", "greater"), 
                     delta = 0, conf.level = 0.95, control = bfControl(), ...) 
{
  # mostly copied from bfTest, but modified to work for estimates and their stderrs and dfs
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) 
    stop("'delta' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                               conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  
  est1<-x$estimate
  se1<-x$stderr
  df1<-x$df
  est2<-y$estimate
  se2<-y$stderr
  df2<-y$df
  calcmethod <- control$calcmethod
  epsilon <- control$epsilon
  nsim <- control$nmc
  Rval<- atan(se2/se1)  
  if (df1 < 1 | df2 < 1) 
    stop("df1 and df2 must be at least 1")
  stderr <- sqrt(se1^2 + se2^2)
  alternative <- match.arg(alternative)
  alpha <- 1 - conf.level
  if (calcmethod == "mc") {
    R2 <- (est2 + se2 * rt(nsim, df2))
    R1 <- (est1 + se1 * rt(nsim, df1))
    D <- R2 - R1
    if (alternative == "two.sided") {
      CI <- quantile(D, probs = c(alpha/2, 1 - alpha/2))
    }
    else if (alternative == "greater") {
      CI <- c(quantile(D, probs = alpha), Inf)
    }
    else if (alternative == "less") {
      CI <- c(-Inf, quantile(D, probs = 1 - alpha))
    }
    p.L <- length(D[D <= delta])/nsim
  }
  else if (calcmethod == "int") {

    ifunc <- function(x, b = delta) {
      pt((b + est1 - est2 + se1 * x)/se2, 
         df2) * dt(x, df1)
    }
    p.L <- integrate(ifunc, -Inf, Inf)$value
    if (alternative == "two.sided") {
      cbf <- qbf(1 - alpha/2, n1=df1+1, n2=df2+1, R=Rval, 
                 epsilon = epsilon)
      CI <- (est2 - est1) + c(-1, 1) * cbf * stderr
    }
    else if (alternative == "greater") {
      cbf <- qbf(1 - alpha, n1=df1+1, n2=df2+1, R=Rval, epsilon = epsilon)
      lower <- (est2 - est1) - cbf * stderr
      CI <- c(lower, Inf)
    }
    else if (alternative == "less") {
      cbf <- qbf(1 - alpha, n1=df1+1, n2=df2+1, R=Rval, epsilon = epsilon)
      upper <- (est2 - est1) + cbf * stderr
      CI <- c(-Inf, upper)
    }
  }
  pval <- switch(alternative, two.sided = min(2 * p.L, 2 * 
                                                (1 - p.L)), less = 1 - p.L, greater = p.L)
  tstat <- (est2 - est1 - delta)/stderr
  #df <- paste("[",df1,df2,"]",sep=",")
  difference<- est2-est1
  estimate <- c(est1, est2, difference)
  names(estimate) <- c("x$estimate", "y$estimate","diff in estimates")
  names(tstat) <- "t"
  names(Rval) <- "R"
  names(delta) <- "difference in parameters"
  attr(CI, "conf.level") <- conf.level
  method <- "Melding on t Distributions (Generalized Behrens-Fisher)"
  if (calcmethod == "mc") {
    method <- paste(method, "(Monte Carlo implementation, nmc=", 
                    nsim, ")")
  }
  rval <- list(statistic = tstat, parameter =Rval, p.value = pval, 
               conf.int = CI, estimate = estimate, null.value = delta, 
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
