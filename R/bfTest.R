## functions needed for Behrens-Fisher test

## cumulative distribution function for B-F statistic
pbf<-function(q,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8)){
  if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
    R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  } else if (is.null(s1) & is.null(s2)){
    # do nothing R<-R
  } else stop("supply one of R or (s1 and s2)")
  
  sinR<-sin(R)
  cosR<-cos(R)
  
  # if sinR=0 then cosR=1 and D=T2
  # and if cosR=0 then sinR=1 and D=T1
  if (abs(sinR)<epsilon){
    out<-pt(q,n2-1)
  } else if (abs(cosR)<epsilon){
    out<-pt(q,n1-1)
  } else if (cosR>0 & cosR<1 & sinR>0 & sinR<1){
    ## usual case       
    ifunc<-function(u,Y){
      pt((Y-u*sinR)/cosR,n2-1)*dt(u,n1-1)
    }
    nq<-length(q)
    out<-rep(NA,nq)
    for (i in 1:nq){
      out[i]<-integrate(ifunc,-Inf,Inf,Y=q[i])$value
    }
  }
  
  out
}
#pbf(c(.9,.95,.975),12,6,R=pi/2-0.00001)
#pt(c(.9,.95,.975),11)
#z<-pbf(c(.345,.45,.6,7),12,6,s1=10,s2=1)
#z


qbf<-function(p,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8)){
  if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
    R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  } else if (is.null(s1) & is.null(s2)){
    # do nothing R<-R
  } else stop("supply one of R or (s1 and s2)")
  
  ## if R=0 then cos(R)=1 and D ~ T2
  if (abs(R)<epsilon){
    out<-qt(p,n2-1)
  } else if (abs(R-pi/2)<epsilon){
    out<-qt(p,n1-1)
  } else {
    rootfunc<-function(q,P){
      P-pbf(q,n1,n2,R=R, epsilon=epsilon)
    }
    np<-length(p)
    out<-rep(NA,np)
    INTERVAL<-c(qt(epsilon,1),-qt(epsilon,1))
    for (i in 1:np){
      out[i]<-uniroot(rootfunc,interval=INTERVAL,P=p[i])$root
    }
    out
  }
  out
}

#pbf(2.91999,2,2,R=0)
#qbf(c(.95),5,8,R=pi/4)
#qbf(c(.05),4,15,R=pi/8)
bfTest<-function (x, ...){
  UseMethod("bfTest")
} 
  
bfTest.formula<-function (formula, data, subset, na.action, ...) 
{
  ### copied from t.test.formula, only needed to change do.call line
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("bfTest", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L) 
    names(y$estimate) <- paste("mean in group", levels(g))
  y
}

bfControl<-function(calcmethod=c("int","mc"),epsilon=10^(-8), nmc=10^5){
   calcmethod<-match.arg(calcmethod)
   if (nmc<2) stop("nmc must get greater than 2, for quite accurate answers use 10^5 or more")
   if (epsilon<=0 | epsilon>1) stop("epsilon should be a very small positive number")   
   list(calcmethod=calcmethod,epsilon=epsilon, nmc=nmc)
}

bfTest.default<-function(x, y, 
                        alternative = c("two.sided", "less", "greater"), 
  mu = 0, conf.level = 0.95, control=bfControl(), 
  ...){
  ## copy first few lines from t.test.default
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                               conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  dname <- paste0("x=",deparse(substitute(x)), "and y=", deparse(substitute(y)))
  ## t.test.default only tests on non-missing data, do the same
  yok <- !is.na(y)
  y <- y[yok]
  xok <- !is.na(x)
  x <- x[xok]
  
  calcmethod<-control$calcmethod
  epsilon<-control$epsilon
  nsim<-control$nmc
  
  
  u2<-mean(x)
  s2<-sqrt(var(x))
  n2<-length(x)
  u1<-mean(y)
  s1<-sqrt(var(y))
  n1<-length(y)
  Rval<- atan((s2/sqrt(n2))/(s1/sqrt(n1)))  

  if (n1<2 | n2<2) stop("n1 and n2 must be at least 2")
  stderr<- sqrt(s1^2/n1 + s2^2/n2 )
  
  alternative<-match.arg(alternative)
  alpha<-1-conf.level   
  
  
  if (calcmethod=="mc"){
    R2<- (u2 + (s2/sqrt(n2))*rt(nsim,n2-1) )  
    R1<- (u1 + (s1/sqrt(n1))*rt(nsim,n1-1) )
    D<- R2-R1
    if (alternative=="two.sided"){
      CI<-quantile(D,probs=c(alpha/2,1-alpha/2))
    } else if (alternative=="greater"){
      CI<- c(quantile(D,probs=alpha),Inf)
    } else if (alternative=="less"){
      CI<- c(-Inf, quantile(D,probs=1-alpha))
    }
    p.L<- length(D[D<=mu])/nsim
  } else if (calcmethod=="int"){
    ifunc<- function(x,b=mu){
      pt( (b+u1-u2+(s1/sqrt(n1))*x)/(s2/sqrt(n2)),n2-1) * dt(x,n1-1)
    }
    p.L<-integrate(ifunc,-Inf,Inf)$value
    if (alternative=="two.sided"){
      cbf<-qbf(1-alpha/2,n1,n2,s1=s1,s2=s2,epsilon=epsilon)
      CI<- (u2-u1) + c(-1,1)*cbf*stderr
    } else if (alternative=="greater"){
      cbf<-qbf(1-alpha,n1,n2,s1=s1,s2=s2,epsilon=epsilon)
      lower<-(u2-u1)-cbf*stderr
      CI<-c(lower,Inf)
    } else if (alternative=="less"){
      cbf<-qbf(1-alpha,n1,n2,s1=s1,s2=s2,epsilon=epsilon)
      upper<-(u2-u1)+cbf*stderr
      CI<-c(-Inf,upper)
    }
  }
  pval<- switch(alternative, two.sided=min(2*p.L,2*(1-p.L)),
                   less=1-p.L, greater=p.L)
  
  tstat <- (u2 - u1 - mu)/stderr
  df<-NA
  names(Rval)<-"R"
  estimate <- c(u2, u1)
  names(estimate) <- c("mean of x", "mean of y")
  
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- "difference in means"
  attr(CI, "conf.level") <- conf.level
  
  method<-"Behrens-Fisher test"
  if (calcmethod=="mc"){
    method<-paste(method,"(Monte Carlo implementation, nmc=",nsim,")")
  }
  
  rval <- list(statistic = tstat, parameter = Rval, p.value = pval, 
               conf.int = CI, estimate = estimate, null.value = mu, 
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
