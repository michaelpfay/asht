abcnonHtest<-function (x, tt, nullValue=NULL, conf.level=0.95, 
    alternative=c("two.sided", "less", "greater"), 
    epsilon = 0.001, minp=0.001) {
    ## This is a modification of the abcnon from the boot R package
    ## first modification, replace alpha with conf.level in call
    ## keep using alpha as in abcnon, so 95% central limit would 
    ## use alpha=c(0.025,0.975)
    ## use alternative also
    dname <- deparse(substitute(x))
    fname<- deparse(substitute(tt))
    alternative <- match.arg(alternative)
    if (alternative=="two.sided"){
        alpha<- (1-conf.level)/2
        alpha<-c(alpha,1-alpha)
    } else if (alternative=="less"){
        # for one-sided upper limit use conf.level as alpha 
        alpha<- conf.level
    } else if (alternative=="greater"){
        alpha<- 1-conf.level
    }
    if (!is.null(nullValue) & is.null(names(nullValue))) names(nullValue)<-"parameter"
    ## end of first modification
    call <- match.call()
    ## allow data.frames
    if (is.matrix(x) | is.data.frame(x)) {
        n <- nrow(x)
    }
    else {
        n <- length(x)
    }
    ep <- epsilon/n
    I <- diag(n)
    P0 <- rep(1/n, n)
    t0 <- tt(P0, x)
    t. <- t.. <- numeric(n)
    for (i in 1:n) {
        di <- I[i, ] - P0
        tp <- tt(P0 + ep * di, x)
        tm <- tt(P0 - ep * di, x)
        t.[i] <- (tp - tm)/(2 * ep)
        t..[i] <- (tp - 2 * t0 + tm)/ep^2
    }
    sighat <- sqrt(sum(t.^2))/n
    a <- (sum(t.^3))/(6 * n^3 * sighat^3)
    delta <- t./(n^2 * sighat)
    cq <- (tt(P0 + ep * delta, x) - 2 * t0 + tt(P0 - ep * delta, 
        x))/(2 * sighat * ep^2)
    bhat <- sum(t..)/(2 * n^2)
    curv <- bhat/sighat - cq
    z0 <- qnorm(2 * pnorm(a) * pnorm(-curv))
    Z <- z0 + qnorm(alpha)
    za <- Z/(1 - a * Z)^2
    ## mod: do not need stan anymore
    #stan <- t0 + sighat * qnorm(alpha)
    abc <- seq(alpha)
    ## mod: do not need pp for output
    #pp <- matrix(0, nrow = n, ncol = length(alpha))
    for (i in seq(alpha)) {
        abc[i] <- tt(P0 + za[i] * delta, x)
        #pp[, i] <- P0 + za[i] * delta
    }
    ## modification: find p-value=value alpha where conf.int limit
    ## equals the null hypothesis value
    if (!is.null(nullValue)){
        getABC<-function(alpha,target=nullValue){
            Z<- z0+qnorm(alpha)
            za<-Z/(1 - a * Z)^2
            tt(P0+za*delta, x) - target
        }
        ## find limits for search, if not between minp and 1-minp, 
        ## do not search
        if (abc[1]>nullValue){
             if (getABC(minp)>0){
                  intLimits<-NULL
             } else {
                  intLimits<-c(minp,alpha[1])
             }
        } else if (length(abc)==1){
             if (getABC(1-minp)<0){
                  intLimits<-NULL
             } else {
                  intLimits<- c(alpha[1],1-minp)
             }
        } else if (abc[1]<= nullValue & nullValue<=abc[2]){
             intLimits<- alpha
        } else {
             if (getABC(1-minp)<0){
                  intLimits<-NULL
             } else {
                  intLimits<- c(alpha[2],1-minp)
             }
        }
        if (is.null(intLimits)){
             p.value<-minp
        } else {
             p.value<-uniroot(getABC,interval=intLimits)$root
        }
        if (alternative=="two.sided"){
            p.value<- min(1,2*p.value,2*(1-p.value))
        } else if (alternative=="less"){
            p.value<- 1-p.value
        }
    } else {
        p.value<-NULL
    }
    # mod: change output, comment out old output
    #limits <- cbind(alpha, abc, stan)
    #dimnames(limits)[[2]] <- c("alpha", "abc", "stan")
    #return(list(limits = limits, stats = list(t0 = t0, sighat = sighat, 
    #    bhat = bhat), constants = list(a = a, z0 = z0, cq = cq), 
    #    tt.inf = t., pp = pp, call = call))

    ## mod: create hyptest class object
    if (alternative=="less"){
        abc<-c(-Inf,abc)
    } else if (alternative=="greater"){
        abc<-c(abc,Inf)
    }
    attr(abc, "conf.level") <- conf.level
    method<- "Nonparametric ABC confidence interval method"
    names(t0)<-"estimate"
    dname<-paste(dname,"  function:",fname)
    out <- list(statistic = NULL, parameter = NULL, p.value = p.value, 
        conf.int = abc, estimate = t0, null.value = nullValue, 
        alternative = alternative, method = method, data.name = dname)
    class(out) <- "htest"
    return(out)
}