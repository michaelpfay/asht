metaNorm <-
  function(y,
           s2,
           method = c("PM", "DL", "fixed"),
           df = NULL,
           nullparm = 0,
           conf.level = 0.95,
           alternative = c("two.sided", "less", "greater"),
           niter = 100,
           epsilon = 1e-10) {
    ## see Desimonian and Kacker, 2007, Contemporary Clinical Trials, 28:105-114
    wi <- 1 / s2
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    wmean <- function(x, w) {
      sum(w * x) / sum(w)
    }
    yw <- wmean(y, wi)
    k <- length(y)
    
    
    
    
    
    
    
    getci <-
      function(u,
               stderr,
               DF = df,
               nparm = nullparm,
               clev = conf.level,
               alt = alternative) {
        alpha <- 1 - clev
        tstat <- (u - nparm) / stderr
        if (alt == "two.sided") {
          ci <- u + c(-1, 1) * stderr * qt(1 - alpha / 2, df = DF)
          p.value <- 2 * (1 - pt(abs(tstat), df = DF))
        } else if (alt == "less") {
          ci <- c(
            -Inf,
            u + stderr * qt(1 - alpha, df = DF) )
            p.value <- pt(tstat, df = DF)
        } else if (alt == "greater") {
          ci <- c(u - stderr * qt(1 - alpha, df = DF), Inf)
          p.value <- pt(tstat, df = DF, lower.tail = FALSE)
        }
        attr(ci, "conf.level")<- clev
        names(df) <- "df"
        if (df == Inf) {
          names(tstat) <- "Z statistic"
        } else {
          names(tstat) <- "t statistic"
        }
        
        names(u) <- names(nparm)<-"weighted mean"
  
        list(
          stderr = stderr,
          null.value = nparm,
          statistic = tstat,
          parameter = df,
          estimate = u,
          conf.int = ci,
          alternative=alt,
          p.value = p.value
        )
      }
    
    
    if (method == "fixed") {
      t2 <- 0
      se <-  1 / sqrt(sum(wi))
      # if degrees of freedom is unspecified, use df=Inf => normal
      METHOD <- "Fixed Effects Meta Analysis"
      if (is.null(df)) {
        df <- Inf
      } else {
        METHOD <- paste0(METHOD, " (using t-distribution with df specfied)")
      }
      result <-   getci(yw, se, DF = df)
    } else if (method == "DL" | method == "PM") {
      #Dersimonian-Laird method (no iteration) or Paule-Mandel method (iterate)
      ## see eq 8, Desimonian and Kacker, 2007
      WSSE <- sum(wi * (y - yw) ^ 2)
      
      t2DL <- (WSSE - (k - 1)) / (sum(wi) - sum(wi ^ 2) / sum(wi))
      t2DL <- max(0, t2DL)
      t2<- t2DL
      ## D-Laird estimator (see right under eq. 8)
      ai <- 1 / (t2DL + s2)
      mwDL <- wmean(y, ai)
      seDL <-  1 / sqrt(sum(ai))
      
      if (method == "DL") {
        METHOD <- "Random Effects Meta Analysis (Dersimonian-Laird Method)"
        if (is.null(df)) {
          df <- Inf
        } else {
          METHOD <- paste0(METHOD, " (using t-distribution with df specfied)")
        }
        result <- getci(mwDL, seDL, DF = df)
      } else if (method == "PM") {
        # Paule-Mandel method
        ## eq 12 of D and K, 2007
        eq12 <- function(wiD, yi, mwDL, si2) {
          max(0,
              (sum(wiD * (yi - mwDL) ^ 2) - (
                sum(wiD * si2) - sum(wiD ^ 2 * si2) / sum(wiD)
              )) /
                (sum(wiD) - sum(wiD ^ 2) / sum(wiD)))
        }
        
        ## eq 12
        t2DL2 <- eq12(ai, y, mwDL, s2)
        mwDLi <- c(mwDL, wmean(y, 1 / (t2DL2 + s2)), rep(NA, niter))
        t2DLi <- c(t2DL, t2DL2, rep(NA, niter))
        
        for (i in 3:(niter + 2)) {
          wiD <- 1 / (t2DLi[i - 1] + s2)
          t2DLi[i] <- eq12(wiD, y, mwDLi[i - 1], s2)
          mwDLi[i] <- wmean(y, 1 / (t2DLi[i] + s2))
          if (abs(mwDLi[i] - mwDLi[i - 1]) < epsilon)
            break()
        }
        mwDLi <- mwDLi[!is.na(mwDLi)]
        N <- length(mwDLi)
        METHOD <- "Random Effects Meta Analysis, Paule-Mandel Method"
        if (is.null(df)) {
          df <- k - 1
          METHOD <- paste0(METHOD, " (using t-distribution with df=k-1)")
        } else {
          METHOD <- paste0(METHOD, " (using prespecified df)")
        }
        statlist <-
          list(
            mu = mwDLi,
            tau2 = t2DLi[1:N],
            mu0 = mwDLi[1],
            tau2DL = t2DLi[1],
            mu.iter = mwDLi[N],
            tau2DL.iter = t2DLi[N]
          )
        t2<- t2DLi[N]
        ai <- 1 / (t2 + s2)
        mwPM <- wmean(y, ai)
        sePM <-  1 / sqrt(sum(ai))
        result <- getci(mwPM, sePM, DF = df)
      }
      
    }
    names(t2)<-"tau squared"
    stats<- c(t2,result$statistic)
    result$statistic<- stats
    stderr<- result$stderr
    names(stderr)<-"Std Err Wgt Mean"
    result$estimate<-c(result$estimate, stderr)
    result <- c(result,
                method = METHOD,
                data.name = deparse(substitute(y)))
    
    class(result) <- "htest"
    result
  }


# Table 3 of Hartung and Knapp 2001, Stat in Med 20:1771-1782
#absDiff<-c(0.2343,0.2541,0.1451,-0.1347,0.1566,0.0894,0.6669,0.1423)
#varAbsDiff<- 0.001*c(4.92,9.18,9.55,15.64,5.81,9.61,62.78,5.38)

#metaNorm(absDiff,varAbsDiff,method="fixed")
#metaNorm(absDiff,varAbsDiff,method="DL")
#metaNorm(absDiff,varAbsDiff,method="PM")

