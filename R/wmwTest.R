## wmwTest functions
## For details, see Fay and Malinovsky (Statistics in Medicine, 2018;37:3991-4006)




wmwTestExact<-function(x,y, alternative=c("two.sided","less","greater"), 
                          tsmethod=c("central","abs"), phiNull=0.5, conf.int=TRUE, conf.level=0.95, 
                          nMC=0, digits=10, ncheckgrid=100,rcheckgrid=0.1){
  alternative<-match.arg(alternative)
  tsmethod<-match.arg(tsmethod)

  
  
  m<-length(x)
  n<-length(y)
  N<-m+n
  # create zobs=vector of group membership
  # ordered by the responses
  # in Fay and Malinovsky (2018, Section 7) this vector
  # is w, in this function it is zobs
  u<-c(x,y)
  zobs<-c(rep(0,m),rep(1,n))
  o<-order(u)
  u<-u[o]
  zobs<-zobs[o]
  
  # Create Zmat, each row is a possible allocation of the groups
  # in Fay and Mainovsky, the jth row is w_j one of the permutations 
  if (nMC==0){
    # complete set of all permutations of groups
    Zmat<-chooseMatrix(N,n)
  } else {
    # create Monte Carlo set of permutations of the groups
    Zmat<-matrix(NA,nMC+1,N)
    Zmat[1,]<-zobs
    for (i in 1:nMC){
      Zmat[i+1,]<-sample(zobs,replace=FALSE)
    }
  }
  J<- nrow(Zmat)
  Nstar<- t(apply(Zmat,1,cumsum))
  Mstar<- t(apply(1-Zmat,1,cumsum))
  # define reverse cumsum, cumsum from end towards the beginning
  rcumsum<-function(x){
    rev(cumsum(rev(x)))
  }
  #rcumsum(1:4)
  Nmat<- t(apply(Zmat,1,rcumsum))
  Mmat<- t(apply(1-Zmat,1,rcumsum))
  
  Rmid<-matrix(rep(rank(u),J),J,N,byrow=TRUE)
  Wmid<- apply(Zmat*Rmid,1,sum)
  # Calcuate Phi for each row of Zmat
  Phimid<- (1/(m*n))*(Wmid - (n*(n+1))/2)
  Phimidobs<- (1/(m*n))*( sum(zobs*rank(u)) - (n*(n+1))/2) 
  Phimid<-round(Phimid,digits)
  Phimidobs<-round(Phimidobs,digits)
  # changes with Version 0.9.4, need PhiMin and PhiMax
  # for bounds when there is ties
  PhiMin<- min(Phimid)
  PhiMax<- max(Phimid)
  
  sumlog<-function(x){ sum( log(x) ) }
  
  # Calculate pi^{PH}(phi), see equation 16 Fay and Malinovsky (2018)
  calcPiPH<-function(phi,m,n,Mmat,Nmat){
    lnumerator<- lfactorial(m)+lfactorial(n)+n*log(1-phi)+m*log(phi)
    denomMat<- phi*Mmat + (1-phi)*Nmat
    ldenom<- apply(denomMat,1,sumlog)
    pijPH<- exp(lnumerator - ldenom)
    pijPH
  }
  # Calculate pi^{LA}(phi), see equation 17 Fay and Malinovsky (2018)
  calcPiLA<-function(phi,m,n,Mstar,Nstar){
    lnumerator<- lfactorial(m)+lfactorial(n)+n*log(phi)+m*log(1-phi)
    denomMat<- (1-phi)*Mstar + phi*Nstar
    ldenom<- apply(denomMat,1,sumlog)
    pijPH<- exp(lnumerator - ldenom)
    pijPH
  }
  
  # equation 18 Fay and Malinovsky
  calcPiLAPH<- function(phi,m,n,Mmat,Nmat,Mstar,Nstar){
    0.5*calcPiLA(phi,m,n,Mstar,Nstar)+ 0.5*calcPiPH(phi,m,n,Mmat,Nmat)
  }
  

  
  calcPval<-function(phi,m,n,Mmat,Nmat,Mstar,Nstar,
                     Phimid,Phimidobs, alternative,tsmethod, doMC=FALSE){
    Pi<-calcPiLAPH(phi,m,n,Mmat,Nmat,Mstar,Nstar)
    if (doMC){ Pi<- Pi/sum(Pi) }
    if (alternative=="less"){
      pval<-sum(Pi[Phimid<=Phimidobs])
    } else if (alternative=="greater"){
      pval<-sum(Pi[Phimid>=Phimidobs])
    } else if (alternative=="two.sided"){
      if(tsmethod=="central"){
        pval.less<- sum(Pi[Phimid<=Phimidobs])
        pval.greater<- sum(Pi[Phimid>=Phimidobs])
        pval<- min(1,2*min(pval.less,pval.greater))
      } else if (tsmethod=="abs"){
        pval<- sum(Pi[ round(abs(Phimid-phi),digits)>=round(abs(Phimidobs-phi),digits)])
      }
    }
    pval
  }
  
  calcCI<-function(m,n,Mmat,Nmat,Mstar,Nstar,
                   Phimid,Phimidobs,alternative,tsmethod, conf.level, eps=10^(-digits), doMC=FALSE){
    rootfunc<-function(phi,Alpha, Alternative){
      calcPval(phi,m,n,Mmat,Nmat,Mstar,Nstar,
               Phimid, Phimidobs, Alternative,tsmethod, doMC) - Alpha
    }
    alpha<- 1-conf.level
    if (alternative=="two.sided" & tsmethod=="central"){ alpha<- alpha/2 }
    ci<-c(0,1)
    # changes with Version 0.9.4 change Phimidobs<1 to Phimidobs< PhiMax
    if ((alternative=="less" | (alternative=="two.sided" & tsmethod=="central")) & Phimidobs<PhiMax){
       #ci[2]<- uniroot(rootfunc,c(eps,1-eps),tol=eps,Alpha=alpha, Alternative="less")$root
       # since for some Monte Carlo simulations the uniroot fails
       # set up try function to give warning and give a NA for the estimate
       ciRoot2<- try(uniroot(rootfunc,c(eps,1-eps),tol=eps,Alpha=alpha, 
                             Alternative="less")$root, silent = TRUE)
       # Sept 16, 2022: use 'inherits(x,"classname")' instead of 'class(x)=="classname")
       # as recommended by check --as-cran
       if (inherits(ciRoot2,"try-error")){
         ci[2]<-NA
         warning("uniroot failed, upper confidence limit set to NA, try method='asymptotic' or if the sample size is small enough method='exact.ce'  ")
       } else {
         ci[2]<-ciRoot2
       }
    } 
    # changes with Version 0.9.4 change Phimidobs>0 to Phimidobs> PhiMin
    if ((alternative=="greater" | (alternative=="two.sided" & tsmethod=="central")) & Phimidobs>PhiMin){
      # see above comments about Monte Carlo simulation uniroot failure
      ciRoot1<- try(uniroot(rootfunc,c(eps,1-eps),tol=eps,Alpha=alpha, 
                        Alternative="greater")$root, silent=TRUE)
      # Sept 16, 2022: use 'inherits(x,"classname")' instead of 'class(x)=="classname")
      # as recommended by check --as-cran
      if (inherits(ciRoot1,"try-error")){
         ci[1]<-NA
         warning("uniroot failed, lower confidence limit set to NA, try method='asymptotic' or if the sample size is small enough method='exact.ce'  ")
       } else {
         ci[1]<-ciRoot1
       }
      
    } 
    if (alternative=="two.sided" & tsmethod=="abs"){
      ## for tsmethod='abs' the pvalue function GENERALLY (but not for sure)
      ##  increases from 0 to Phimidobs then decreases from Phimidobs to 1
      ## but the p-value function can have jumps and non-monotonic parts
      ## it can be like the Fisher-Irwin p-value 
      ## (see Fig 1 of the supplement to Fay, 2010, Biostatistics)
      ## 
      getupper<-getlower<-TRUE
      pvallo<- calcPval(eps,m,n,Mmat,Nmat,Mstar,Nstar,
                        Phimid, Phimidobs, alternative,tsmethod, doMC)
      pvalhi<-calcPval(1-eps,m,n,Mmat,Nmat,Mstar,Nstar,
                       Phimid, Phimidobs, alternative,tsmethod, doMC)
      
      
      if (Phimidobs>1-eps){
        # p-value for phi=1 is around 1, so upper CL should remain 1
        getupper<-FALSE
        pvalPhimidobs<- pvalhi
      } else if (Phimidobs<eps){
        # p-value for phi=0 is around 1, so lower CL should remain 0
        getlower<-FALSE
        pvalPhimidobs<- pvallo
      } else {
        pvalPhimidobs<- calcPval(Phimidobs,m,n,Mmat,Nmat,Mstar,Nstar,
                                 Phimid, Phimidobs, alternative,tsmethod, doMC)
        if (pvallo>alpha){
          # p-value for phi=eps is greater than alpha, so lower cL=0
          getlower<-FALSE
        }
        if (pvalhi>alpha){
          # p-value for phi=1-eps is greater than alpha, so upper cL=1
          getupper<-FALSE
        }
      }
      
      if (getlower){
        # first use uniroot to find a root, then check below that root 
        # to see if there are any more
        cilo<-uniroot(rootfunc,c(eps,Phimidobs),f.lower=pvallo-alpha,f.upper=pvalPhimidobs-alpha,
                       tol=eps, Alpha=alpha, Alternative=alternative)$root
        # check every point on checkgrid to see if it has a pvalue greater than alpha
        checkgrid<- seq(from=max(eps,cilo-rcheckgrid),to=cilo,length.out=ncheckgrid)
        pvals<-rep(alpha,ncheckgrid)
        # set pvals=alpha, then replace all but the last one (at cilo)
        for (i in 1:(ncheckgrid-1)){
          pvals[i]<-calcPval(checkgrid[i],m,n,Mmat,Nmat,Mstar,Nstar,
                             Phimid, Phimidobs, alternative,tsmethod, doMC)
        }
        ## remember cilo=checkgrid[ncheckgrid] so if none of the p-vals are >=alpha
        ## then ci[1]=cilo
        ci[1]<- min(checkgrid[pvals>=alpha])
      }
      if (getupper){
        # first use uniroot to find a root, then check below that root 
        # to see if there are any more
        cihi<-uniroot(rootfunc,c(Phimidobs,1-eps),f.lower=pvalPhimidobs-alpha,f.upper=pvalhi-alpha,
                      tol=eps, Alpha=alpha, Alternative=alternative)$root
        # check every point on checkgrid to see if it has a pvalue greater than alpha
        checkgrid<- seq(from=cihi,to=min(1-eps,cihi+rcheckgrid),length.out=ncheckgrid)
        pvals<-rep(alpha,ncheckgrid)
        # set pvals=alpha, then replace all but the last one (at cilo)
        for (i in 2:ncheckgrid){
          pvals[i]<-calcPval(checkgrid[i],m,n,Mmat,Nmat,Mstar,Nstar,
                             Phimid, Phimidobs, alternative,tsmethod, doMC)
        }
        ## remember cilo=checkgrid[ncheckgrid] so if none of the p-vals are >=alpha
        ## then ci[1]=cilo
        ci[2]<- max(checkgrid[pvals>=alpha])
      }
    }
    ci
  }
  
  
  p.value<-  calcPval(phiNull,m,n,Mmat,Nmat,Mstar,Nstar,
                                Phimid,Phimidobs, alternative,tsmethod, doMC=(nMC>0))
  if (conf.int){
    ci<-  calcCI(m,n,Mmat,Nmat,Mstar,Nstar,
                       Phimid,Phimidobs,alternative,tsmethod, conf.level, eps=10^(-digits), doMC=(nMC>0))
  } else {
    ci<-c(NA,NA)
  }
  METHOD<-"exact Wilcoxon-Man-Whitney test"
  if (nMC>0) METHOD<-paste0(METHOD," (Monte Carlo with nMC=",nMC,")")

  list(p.value=p.value,estimate=Phimidobs,conf.int=ci, METHOD=METHOD)
}





#set.seed(13)
#x<-rnorm(7)
#y<-rnorm(7)+.5
#t0<-proc.time()
#wmwTestExact(x,y,nMC=0)
#t1<-proc.time()
#t1-t0
#wmwTestExact(x,y,nMC=10^4)
#t2<-proc.time()
#t2-t1
#wmwTestAsymptotic(x,y)
#t3<-proc.time()
#t3-t2
#x<-c(8,7,2,5,9)
#y<-c(9,9,9,8)
#wmwTestExact(x,y,nMC=0)
#wmwTestExact(x,y,nMC=10^5)


Vpo<- function(PHI,tf,ny,nx){
  findu<-function(p){
    rootfunc<-function(u,P=p){
      P - ((exp(u))/(-1+exp(u))^2)*(-1+exp(u) +log(1+exp(-u))-log(1+exp(u)))
    }
    if (p==0.5){
      U<- 0
      return(U)
    } else if (p>0.5){
      lo<- 10^(-5)
      hi<- 100
    } else {
      hi<- -10^(-5)
      lo<- -100
    }
    U<-uniroot(rootfunc,lower=lo,upper=hi,tol=10^(-3))$root
    U
  }
  
  u<-findu(PHI)
  if (PHI==0.5){
    Q<- 1/3 
  } else {
    Q<-2*exp(2*u)*(log(1+exp(-u))-log(1+exp(u))+sinh(u))/(-1+exp(u))^3
  }  
  
  V<- tf*(PHI*(1-PHI) + (nx+ny-2)*(Q-PHI^2))/(nx*ny)
  V
}





wmwTestAsymptotic<-function(x,y, phiNull=0.5, conf.level=0.95, alternative="two.sided",
                            correct=TRUE, epsilon=10^(-8),RemoveTieAdjustment=FALSE,Vmethod="LAPH"){
  r <- rank(c(y, x))
  n.y <- as.double(length(y))
  n.x<- as.double(length(x))
  n<-n.y+n.x
  sum.y<- sum(r[1:n.y])
  phi<- (sum.y - n.y*(n.y+1)/2 )/(n.y*n.x)
  NTIES <- table(r)
  tiefactor<- 1 -   sum(NTIES^3 - NTIES)/(n*(n+1)*(n-1))
  
  if (RemoveTieAdjustment){
    if (tiefactor<1) warning("tie factor removed, but there are ties in the data.")
    tiefactor<-1
  } 
  METHOD<- "Wilcoxon-Mann-Whitney test"
  
  if (tiefactor==0){
    # all values in both groups are equal
    ci<-c(0,1)
    p.value<-1
  } else {
    if (Vmethod=="LAPH"){
      V<- function(PHI,tf=tiefactor,ny=n.y,nx=n.x){
        tf*(PHI*(1-PHI)/(ny*nx))*
          (1+((ny+nx-2)/2)*
             ((1-PHI)/(2-PHI)+PHI/(1+PHI)))
      }      
    } else {
      V<- function(PHI,tf=tiefactor,ny=n.y,nx=n.x){
        Vpo(PHI,tf,ny,nx)
      }
    }

    
    CORRECTIONLESS<-CORRECTIONGR<-0
    if (correct) {
      # Changed to have separate correction for less and greater
      # then make two.sided p equal to min of 2*one-sided p.values 
      # OLD METHOD in wilcox.test:
      # CORRECTION <- switch(alternative, two.sided = sign(phi-phiNull) * 
      #                      0.5, greater = 0.5, less = -0.5)
      # CORRECTION<- CORRECTION/(n.y*n.x)
      ## We are using phi instead of the sum of the ranks in one group
      ## so we need to divide the correction by n.y*n.x
      CORRECTIONLESS<- -.5/(n.y*n.x)
      CORRECTIONGR<- 0.5/(n.y*n.x)
      METHOD <- paste(METHOD, "with continuity correction")
    }
    z0Less <- (phi-phiNull-CORRECTIONLESS)/sqrt(V(phiNull))
    z0Gr<- (phi-phiNull-CORRECTIONGR)/sqrt(V(phiNull))
    PVAL <- switch(alternative, less = pnorm(z0Less), 
                   greater = pnorm(z0Gr, lower.tail = FALSE), 
                   two.sided = 2 * min(pnorm(z0Less), pnorm(z0Gr, lower.tail = FALSE)))
    ## Now do confidence interval on phi
    alpha <- 1 - conf.level
    wfunc <- function(phiNull, zq, correction=0) {
      (phi-phiNull-correction)/sqrt(V(phiNull)) - zq
    }
    phimin<- epsilon
    phimax<- 1- epsilon
    root <- function(zq, corr) {
      f.lower <- wfunc(phimin, zq, corr)
      if (f.lower <= 0) 
        return(phimin)
      f.upper <- wfunc(phimax, zq, corr)
      if (f.upper >= 0) 
        return(phimax)
      uniroot(wfunc, c(phimin, phimax), f.lower = f.lower, 
              f.upper = f.upper, tol = epsilon, zq = zq, correction=corr)$root
    }
    ci <- switch(alternative, two.sided = {
      l <- ifelse(phi==0,0,root(zq = qnorm(alpha/2, lower.tail = FALSE), 
                                corr=CORRECTIONGR))
      u <- ifelse(phi==1,1,root(zq = qnorm(alpha/2), 
                                corr=CORRECTIONLESS))
      c(l, u)
    }, greater = {
      l <- ifelse(phi==0,0,root(zq = qnorm(alpha, lower.tail = FALSE), corr=CORRECTIONGR))
      c(l, 1)
    }, less = {
      u <- ifelse(phi==1,1,root(zq = qnorm(alpha),corr=CORRECTIONLESS))
      c(0, u)
    })
  }

  list(p.value=PVAL,estimate=phi,conf.int=ci, tiefactor=tiefactor,METHOD=METHOD)
}

#set.seed(12)
#x<-rpois(16,13)
#y<-rpois(27,14)
#wmwTestAsymptotic(x,y)
#wmwTestAsymptotic(x,y,Vmethod="PO")

#wmwTestAsymptotic(x,y,alternative="two.sided",conf.level=1-w$p.value)
#wmwTestAsymptotic(x,y,alternative="greater",phiNull=w$conf.int[1])
#wmwTestAsymptotic(x,y,alternative="less",phiNull=w$conf.int[2])
#wilcox.test(x,y)
#1-(1/(16*27))*(183)


wmwTest <-
  function (x, ...) {
    UseMethod("wmwTest")  
  }

wmwTest.matrix <-
  function(x,...){
    if (nrow(x)!=2L){
      stop("matrix must have 2 rows representing 2 groups")
    }
    k<-ncol(x)
    X<-rep(1:k,times=x[1,])
    Y<-rep(1:k,times=x[2,])
    rowNames<- dimnames(x)[[1]]
    if (!is.null(rowNames)){
      DNAME<- paste(rowNames[1],"and",rowNames[2])
    } else {
      DNAME<- "x (top row) and y (bottom row)"
    }
    output <- wmwTest.default(X,Y,...)
    output$data.name<-DNAME
    output
  }

wmwTest.formula <-
  function (formula, data, subset, na.action, ...) 
  {
    ## copied from wilcox.test.formula 
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
    y <- do.call("wmwTest", c(DATA, list(...)))
    y$data.name <- DNAME
    y
  }



methodRuleWMW<-function(x,y,exact,chooseLimit=5000){
  m<-length(x)
  n<-length(y)
  J<-choose(m+n,n)
  if (is.null(exact)){
    if (J<=chooseLimit){
      rule<-"exact.ce"
    } else {
      rule<-"asymptotic"
    }
  } else if (exact){
    if (J>chooseLimit){
      rule<-"exact.mc"
    } else {
      rule<-"exact.ce"
    }
  } else {
    rule<-"asymptotic"
  }
  rule
}


latentTransform<-function(x,y,phiValues, output=c("mw","po"), epsilon=10^(-6)){
  output<-match.arg(output)
  getH<-function(x,y){
    H<-ecdf(c(x,y))
    u<- sort(unique(c(x,y)))
    H<-H(u)
    H
  }
  # H=cdf of all responses
  H<-getH(x,y)
  m<-length(x)
  n<-length(y)
  transPO_to_MW<-function(theta){
    # Translate proportional odds parameter, theta,
    # to the Mann-Whitney parameter
    # given H
    # H=empirical cdf for all individuals
    # set H = (m*F + n*G)/(m+n)
    # and use the proportional odds relationship 
    # between F and G: F/(1-F) = theta*G/(1-G) to get 
    # G in terms of F:     G = F/(theta + (1-theta)*F)
    # solve for  F=distn of X1,...Xm
    getF<-function(m,n,theta,H){
      # use quadratic formula
      A<- m*(1-theta)
      B<- m*theta - (m+n)*H*(1-theta)+n
      C<- -(m+n)*H*theta
      # quadratic formula, positve sign give proper cdf for F
      F<- (-B + sqrt(B^2 - 4*A*C))/(2*A)
      F
    }
    F<-getF(m,n,theta,H)
    # use algebra: from F/(1-F) = theta*G/(1-G) to get....
    G<- F/(theta + (1-theta)*F)
    getPhiFG<-function(F,G){
      k<-length(F)
      f<- F - c(0,F[-k])
      g<- G - c(0,G[-k])
      sum(F*g) - 0.5 * sum(g*f)
    }
    phi<- getPhiFG(F,G)
    phi
  }
  
  rootfunc<-function(theta,phiTarget){
    phiTarget - transPO_to_MW(theta)
  }
  minPhi<- transPO_to_MW(epsilon)
  maxPhi<- transPO_to_MW(1/epsilon)
  II<-    (phiValues<minPhi | phiValues>maxPhi)  
  if (any(II)) warning(paste("latent continuous transform for phiGrouped<",prettyNum(minPhi), 
                             " set to 0, and phiGrouped>", prettyNum(maxPhi)," set to 1.",
                             "To change those boundaries, change epsilon."))
  
  modPhiValues<-phiValues  
  modPhiValues[phiValues<minPhi]<-0
  modPhiValues[phiValues>maxPhi]<-1
  
  I<-  !(modPhiValues %in% c(0,.5,1) )
  phiValues<- phiValues[I]
  if (length(phiValues)>0){
    thetaValues<-rep(NA,length(phiValues))
    for (i in 1:length(phiValues)){
      thetaValues[i]<-uniroot(rootfunc,c(epsilon,1/epsilon),phiTarget=phiValues[i])$root
    }
    if (output=="po"){
      out<-thetaValues
    } else {
      getPhiTheta<-function(theta){
        phi<- theta*(-1+theta-log(theta))/(-1+theta)^2
        phi[theta==0]<- 0
        phi[theta==1]<-0.5
        phi[theta==Inf]<-1
        phi
      }
      out<-rep(NA,length(phiValues))
      for (i in 1:length(phiValues)){
        out[i]<-getPhiTheta(thetaValues[i])
      }
    }
    outPhi<- modPhiValues
    outPhi[I]<-out
  } else { outPhi<- modPhiValues }
  outPhi
}

#x<-c(rep(1,3),rep(2,5),rep(3,4))
#y<-c(rep(1,1),rep(2,6),rep(3,12))
#w<-wmwTestAsymptotic(x,y)
#latentTransform(x,y,c(w$estimate,w$conf.int,.5), output=c("mw"), epsilon=10^(-6))



wmwControl<-function(nMC=10^4, epsilon=10^(-8),digits=10, latentOutput=c("mw","po"), 
                     removeTieAdjustment=FALSE, ncheckgrid=100,rcheckgrid=0.1, Vmethod="LAPH"){
  latentOutput<-match.arg(latentOutput)
  list(nMC=nMC,epsilon=epsilon,digits=digits, latentOutput=latentOutput, 
       removeTieAdjustment=removeTieAdjustment, ncheckgrid=ncheckgrid,
       rcheckgrid=rcheckgrid, Vmethod=Vmethod)
}

wmwTest.default <-
  function(x,y,
           alternative=c("two.sided","less","greater"),
           phiNull=0.5, exact=NULL, correct=TRUE, 
           conf.int=TRUE, conf.level=0.95, 
           latentContinuous=FALSE, 
           method=NULL,
           methodRule=methodRuleWMW,
           tsmethod=c("central","abs"),
           control=wmwControl(),...){
    
    
    alternative<-match.arg(alternative)
    tsmethod<-match.arg(tsmethod)
   
    
    if (!((length(phiNull) == 1L) && is.finite(phiNull) && 
          (phiNull > 0) && (phiNull < 1))) 
      stop("'phiNull' must be a single number between 0 and 1")
    if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
          (conf.level > 0) && (conf.level < 1))) 
      stop("'conf.level' must be a single number between 0 and 1")
    if (!is.numeric(x)) 
      stop("'x' must be numeric")
    if (!is.null(y)) {
      if (!is.numeric(y)) 
        stop("'y' must be numeric")
      DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    
    
    
    
    if (any(is.na(y)) | any(is.na(x))){
      warning("NAs removed before calculations")
      y<- y[!is.na(y)]
      x<- x[!is.na(x)]
    }
    # check length
    if (length(y)==0 | length(x)==0){ stop("length of y or x equals 0")}
    
    if (is.null(method)){
      method<-methodRule(x,y,exact)
    }
    # use match.arg so that if method="asy" it will work
    method<-match.arg(method,c("asymptotic","exact.ce","exact.mc"))
    if (method=="asymptotic"){
      wout<-wmwTestAsymptotic(x,y, phiNull, conf.level, alternative,
                         correct, epsilon=control$epsilon,
                         RemoveTieAdjustment=control$removeTieAdjustment, Vmethod=control$Vmethod)
    } else if (method=="exact.ce"){
      wout<-wmwTestExact(x,y, alternative, tsmethod, phiNull, conf.int, conf.level, 
                         nMC=0, digits=control$digits, ncheckgrid=control$ncheckgrid,
                         rcheckgrid=control$rcheckgrid)
    } else if (method=="exact.mc"){
      wout<-wmwTestExact(x,y, alternative, tsmethod, phiNull, conf.int, conf.level, 
                         nMC=control$nMC, digits=control$digits, ncheckgrid=control$ncheckgrid,
                         rcheckgrid=control$rcheckgrid)
    } else stop("method must be either 'asymptotic', 'exact.ce' or 'exact.mc' ")
    
    p.value<-wout$p.value
    cint<- wout$conf.int
    phi<-wout$estimate
    # tiefactor is NULL for exact methods
    tiefactor<-wout$tiefactor
    METHOD<-wout$METHOD
    
    attr(cint, "conf.level") <- conf.level    
    
    names(phiNull)<- "Mann-Whitney parameter"
    names(phi)<-"Mann-Whitney estimate"
    if (!is.null(tiefactor)){ names(tiefactor)<-"tie factor" }
    est<- stat<- phi
    if (latentContinuous){
      # change estimate and CI to latent MW parameter from prop odds model
      if (!is.na(cint[1])){
        ltout<- latentTransform(x,y,c(est,phiNull,cint), control$latentOutput, control$epsilon)
        est<- ltout[1]
        phiNull<-ltout[2]
        cint<-ltout[3:4]
      } else {
        ltout<- latentTransform(x,y,c(est,phiNull), control$latentOutput, control$epsilon)
        est<- ltout[1]
        phiNull<-ltout[2]
        cint<-ltout[3:4]
      }
      if (control$latentOutput=="mw"){
        names(phiNull)<- "latent continuous Mann-Whitney parameter"
        names(est)<-"latent continuous Mann-Whitney estimate"
      } else if (control$latentOutput=="po"){
        names(phiNull)<- "latent proportional odds parameter"
        names(est)<-"latent proportional odds estimate"
      }
    }
    if (phiNull==0.5 & alternative=="two.sided"){
      phiNull<-NULL
      alternative<-"two distributions are not equal"
      METHOD<-paste0(METHOD,"\n (confidence interval requires proportional odds assumption, but test does not)")
    } else {
      METHOD<-paste0(METHOD," (test and confidence interval require proportional odds assumption)")
    }
    output <- list(statistic = stat, estimate=est, parameter = tiefactor, p.value = as.numeric(p.value), 
                   null.value = phiNull, alternative = alternative, method = METHOD, 
                   data.name = DNAME, conf.int=cint)
    class(output)<-"htest"
    output
  }


#wmwTest(c(12,3,2,4,5,6),c(12,14,15,2,15,7))
