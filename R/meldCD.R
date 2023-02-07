


meldCD <-
  function(H1,
           H2,
           nullparm = NULL,
           parmtype = c("difference", "ratio", "oddsratio"),
           conf.level = 0.95,
           alternative = c("two.sided", "less", "greater"),
           estimate = c("median", "mean"),
           lim = c(-Inf, Inf),
           parmGrid = NULL,
           nmc = 1e5,
           ngrid = 1e4,
           calcmethod = "int", epsilon=1e-8, utol=1e-8) {
    
    ## CHANGES TO MAKE: Look into default parmGrid when using difference in Poisson rates
    ##      where the rates are very small, say 1e-4.
    ##      Look at the way I defined parmGrid in 
    #       /PUBLISHED/2015 Fay Proschan Brittain/R/poisson2sampleDiff/poisson2sampleDiff.R 
    
    alternative <- match.arg(alternative)
    parmtype <- match.arg(parmtype)
    estimate <- match.arg(estimate)
    
    
    

    
    
    
    if (!missing(conf.level) &&
        (length(conf.level) != 1 || !is.finite(conf.level) ||
         conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    dname <-
      paste(deparse(substitute(H1)), "and", deparse(substitute(H2)))
    
    
    
    if (min(lim) < 0 &
        (parmtype == "ratio" |
         parmtype == "oddsratio"))
      stop("cannot have min(lim)<0 when parmtype='ratio' or 'oddsratio' ")
    if (max(lim) > 1 &
        (parmtype == "oddsratio"))
      stop("cannot have max(lim)>1 when parmtype='oddsratio' ")
    
    if (parmtype == "difference") {
      g <- function(x1, x2) {
        x2 - x1
      }
      if (is.null(nullparm))
        nullparm <- 0
      names(nullparm) <- "parm2 - parm1"
    } else if (parmtype == "ratio") {
      g <- function(x1, x2) {
        x2 / x1
      }
      if (is.null(nullparm))
        nullparm <- 1
      names(nullparm) <- "parm2/parm1"
    } else if (parmtype == "oddsratio") {
      g <- function(x1, x2) {
        (x2 * (1 - x1)) / (x1 * (1 - x2))
      }
      if (is.null(nullparm))
        nullparm <- 1
      names(nullparm) <- "(parm2(1-parm1))/(parm1(1-parm2))"
    }
    
    glim <- c(g(lim[2], lim[1]), g(lim[1], lim[2]))
    
    ## define parmGrid
    if (is.null(parmGrid)) {
      if (lim[1] >= 0 & lim[2] <= 1) {
        parmGrid <- seq(from = lim[1],
                     to = lim[2],
                     length.out = ngrid)
      } else if (lim[1] == 0 & lim[2] == Inf) {
        Power2 <- ceiling(log2(ngrid))
        parmGrid <- exact2x2::power2gridRatio(power2 = Power2)
        parmGrid <- parmGrid[parmGrid < Inf]
      } else if (lim[1] == -Inf & lim[2] == Inf) {
        Power2 <- ceiling(log2(ngrid / 2))
        posparmGrid <- exact2x2::power2gridRatio(power2 = Power2)
        parmGrid <- c(-rev(posparmGrid), posparmGrid)
        parmGrid <- parmGrid[parmGrid < Inf & parmGrid > -Inf]
      } else {
        stop(
          "no default method to create parmGrid for this 'lim', enter vector of parmGrid within range of 'lim' "
        )
      }
    }
    
    alpha <- 1 - conf.level
    F1 <- H1(parmGrid)
    if (min(F1)>epsilon){ warning(paste("H1(min(parmGrid))=",min(F1),"(should be approximately 0) results could be questionable")) }
    if (max(F1)<1-epsilon){ warning(paste("H1(max(parmGrid))=",max(F1),"(should be approximately 1) results could be questionable")) }
    
    F2 <- H2(parmGrid)
    if (min(F2)>epsilon){ warning(paste("H2(min(parmGrid))=",min(F2),"(should be approximately 0) results could be questionable")) }
    if (max(F2)<1-epsilon){ warning(paste("H2(max(parmGrid))=",max(F2),"(should be approximately 1) results could be questionable")) }
    
    if (calcmethod == "mc") {
      U1 <- runif(nmc)
      U2 <- runif(nmc)
      x1 <- approx(x = F1, y = parmGrid, xout = U1)$y
      x2 <- approx(x = F2, y = parmGrid, xout = U2)$y
      if (estimate == "median") {
        est1 <- approx(x = F1, y = parmGrid, xout = 0.5)$y
        est2 <- approx(x = F2, y = parmGrid, xout = 0.5)$y
      } else if (estimate == "mean") {
        f1 <- diff(c(0, F1))
        est1 <- sum(parmGrid * f1)
        f2 <- diff(c(0, F2))
        est2 <- sum(parmGrid * f2)
      }
      estVector <- g(x1, x2)
      nev <- length(estVector)
      if (alternative == "two.sided") {
        pless <- (length(estVector[estVector >= nullparm]) + 1) / (nev + 1)
        pgr <-
          (length(estVector[estVector <= nullparm]) + 1) / (nev + 1)
        ci <-
          quantile(estVector, probs = c(alpha / 2, 1 - alpha / 2))
        p.value <- min(1, 2 * min(pless, pgr))
      } else if (alternative == "less") {
        pless <- (length(estVector[estVector >= nullparm]) + 1) / (nev + 1)
        ci <- c(glim[1], quantile(estVector, probs = c(1 - alpha)))
        p.value <- pless
      }   else if (alternative == "greater") {
        pgr <- (length(estVector[estVector <= nullparm]) + 1) / (nev + 1)
        ci <- c(quantile(estVector, probs = c(alpha)), glim[2])
        p.value <- pgr
      }
    }  else if (calcmethod == "int") {
      npg <- length(parmGrid)
      xmid <- 0.5 * parmGrid[-1] + 0.5 * parmGrid[-npg]
      h1 <- H1(parmGrid[-1]) - H1(parmGrid[-length(parmGrid)])
      h1 <- h1 / sum(h1)
      h2 <- H2(parmGrid[-1]) - H2(parmGrid[-length(parmGrid)])
      h2 <- h2 / sum(h2)
      root1 <- function(x, q = 0.5) {
        H1(x) - q
      }
      root2 <- function(x, q = 0.5) {
        H2(x) - q
      }
      if (parmtype == "difference") {
        # if H2(x)=0 for x<lim[1] and H2(x)=1 for x>lim[2]
        # then we can define it simply as
        #      Hg <- function(x) {
        #          sum(H2(xmid + x) * h1)
        #      }
        # But in case it is not defined that way, write it like this
        Hg <- function(x) {
          xeval<- xmid+x
          if (all(xeval<lim[1])){
            out<- 0
          } else if (all(xeval>lim[2])){
            out<- 1
          } else {
            H2eval<-rep(NA,length(h1))
            H2eval[xeval<lim[1]]<- 0
            H2eval[xeval>lim[2]]<-1
            I<- xeval>=lim[1] & xeval<=lim[2]
            H2eval[I]<- H2(xeval[I]) 
            out<-sum(H2eval * h1)            
          }
          out
        }
      } else if (parmtype == "ratio") {
        # easy method for when H2(x)=0 if x<lim[1] and H2(x)=1 if x>lim[2]
        #Hg <- function(x) {
        #  sum(H2(xmid * x) * h1)
        #}
        Hg <- function(x) {
          xeval<- xmid*x
          if (all(xeval<lim[1])){
            out<- 0
          } else if (all(xeval>lim[2])){
            out<- 1
          } else {
            H2eval<-rep(NA,length(h1))
            H2eval[xeval<lim[1]]<- 0
            H2eval[xeval>lim[2]]<-1
            I<- xeval>=lim[1] & xeval<=lim[2]
            H2eval[I]<- H2(xeval[I]) 
            out<-sum(H2eval * h1)            
          }
          out
        }
        
      } else if (parmtype == "oddsratio") {
        #Hg <- function(x) {
        #  sum(H2(xmid * x / (1 - xmid + xmid * x)) * h1)
        #}
        Hg <- function(x) {
          xeval<- xmid * x / (1 - xmid + xmid * x)
          if (all(xeval<lim[1])){
            out<- 0
          } else if (all(xeval>lim[2])){
            out<- 1
          } else {
            H2eval<-rep(NA,length(h1))
            H2eval[xeval<lim[1]]<- 0
            H2eval[xeval>lim[2]]<-1
            I<- xeval>=lim[1] & xeval<=lim[2]
            H2eval[I]<- H2(xeval[I]) 
            out<-sum(H2eval * h1)            
          }
          out
        }
      }
      
      rootg <- function(x, q = 0.5) {
        Hg(x) - q
      }
      
      
      if (estimate == "median") {
        if (-sign(root1(min(parmGrid)))+sign(root1(max(parmGrid)))==2){
          est1 <-  uniroot(root1, interval = range(parmGrid), q = 0.5, tol=utol)$root
       } else {
          warning("cannot calculate median of H1")
          est1<-NA
        }
        if (-sign(root2(min(parmGrid)))+sign(root2(max(parmGrid)))==2){
          est2 <-  uniroot(root2, interval = range(parmGrid), q = 0.5, tol=utol)$root
        } else {
          if (H2(min(parmGrid))==0.5){
            warning("median of CD for group 2 may not be unique, set to minimum parmGrid value")
            est2<-min(parmGrid)
          } else           if (H2(max(parmGrid))==0.5){
            warning("median of CD for group 2 may not be unique, set to maximum parmGrid value")
            est2<-max(parmGrid)
          } else {
            warning("median of CD for group 2 not found, consider using calcmethod='mc' or explicitly specifying parmGrid or increasing its range")
            est2<-NA
          }
        }
      } else if (estimate == "mean") {
        est1 <- sum(h1 * xmid)
        est2 <- sum(h2 * xmid)
      }
      pgr <- Hg(nullparm)
      pless <- 1 - pgr
      
      # get limits to search for roots of rootg
      GLIM<-c(g(max(xmid),min(xmid)),g(min(xmid),max(xmid)))
      
      if (alternative == "two.sided") {
        cilo <-  uniroot(rootg, interval = GLIM, q = alpha / 2, tol=utol)$root
        cihi <-
          uniroot(rootg,
                  interval = GLIM,
                  q = 1 - alpha / 2, tol=utol)$root
        ci <- c(cilo, cihi)
        p.value <- min(1, 2 * min(pless, pgr))
      } else if (alternative == "less") {
        cihi <- uniroot(rootg,
                        interval = GLIM,
                        q = 1 - alpha, tol=utol)$root
        ci <- c(glim[1], cihi)
        p.value <- pless
      }   else if (alternative == "greater") {
        cilo <- uniroot(rootg, interval = GLIM, q = alpha, tol=utol)$root
        ci <- c(cilo, glim[2])
        p.value <- pgr
      }
    }
    est <- g(est1, est2)
    Est <- c(est1, est2, est)
    if (estimate == "median") {
      names(Est) <- c("median1", "median2", "g")
      if (parmtype == "difference") {
        names(Est)[3] <- c("median2-median1")
      } else if (parmtype == "ratio") {
        names(Est)[3] <- c("median2/median1")
      } else if (parmtype == "oddsratio") {
        names(Est)[3] <- c("((med2)(1-med1))/(med1(1-med2))")
      }
    } else if (estimate == "mean") {
      names(Est) <- c("mean1", "mean2", "g")
      if (parmtype == "difference") {
        names(Est)[3] <- c("mean2-mean1")
      } else if (parmtype == "ratio") {
        names(Est)[3] <- c("mean2/mean1")
      } else if (parmtype == "oddsratio") {
        names(Est)[3] <- c("((mean2)(1-mean1))/(mean1(1-mean2))")
      }
    }
    
    attr(ci, "conf.level") <- conf.level
    method <- "Melding Test on Confidence Distributions"
    if (calcmethod == "mc") {
      method <- paste(method, "(Monte Carlo implementation, nmc=",
                      nmc, ")")
    }
    rval <-
      list(
        statistic = NULL,
        parameter = NULL,
        p.value = p.value,
        conf.int = ci,
        estimate = Est,
        null.value = nullparm,
        alternative = alternative,
        method = method,
        data.name = dname
      )
    class(rval) <- "htest"
    return(rval)
    
  }

