wsrTest <-
function(x,y=NULL,conf.int=TRUE,conf.level=.95,mu=0,
    alternative=c("two.sided","less","greater"), 
    digits=NULL, tieDigits=8){
  
    if (is.null(y)){
        dname<-paste0("single sample=",deparse(substitute(x)))
        y<-rep(0,length(x))
    } else {
        dname<-paste0(deparse(substitute(x))," minus ",deparse(substitute(y)))
    }
    
    x<-round(x,tieDigits)
    y<-round(y,tieDigits)
    # if you use d<- x-y-mu, then the range searched (min(d) to max(d)) 
    # to look for the CIs is not wide enough
    # d<-x-y-mu
    d<-x-y
    Est<-NA
    
    alt<-match.arg(alternative)
    if (conf.int){
        ci<-c(-Inf,Inf)
        alpha<-1-conf.level
        # using uniroot.integer, we estimate theta to within
        # epsilon of its true value
        # we want to take maximum step sizes of 
        # 1/2^(STEP.POWER+1) <= epsilon/(max(d)-min(d))
        # where 
        Range<-max(d)-min(d)
        # if digits missing, round to the nearest 1 if Range=10,000
        #   or round to nearest 0.0001 if Range=1
        if (is.null(digits)) digits<- -(log10(Range)-4) 
        epsilon<- 10^(-digits)
        # so... 
        STEP.POWER<- - floor(log2(epsilon/Range))
        # theta takes the interval from min(d) to max(d)
        # and divides it up into 1/2^(STEP.POWER+1) equal parts
        # theta(0) = min(d)
        # and theta(2^(STEP.POWER+1))=max(d)
        # and theta(2^STEP.POWER) = min(d) + .5*(max(d)-min(d))  
        theta<-function(i){
             min(d) + (i/(2^(STEP.POWER+1)))*(max(d)-min(d)) 
        }
     
        
            EstRootFunc<-function(i){
              ytemp<- round(y+ theta(i),tieDigits)
              pvalue(wilcoxsign_test(x~ytemp,
                    zero.method="Pratt",
                    alternative="less",distribution=exact())) -
              pvalue(wilcoxsign_test(x~ytemp,
                    zero.method="Pratt",
                    alternative="greater",distribution=exact())) 
            } 
            irootPos<-uniroot.integer(EstRootFunc,c(0,2^(STEP.POWER+1)),
                step.power=STEP.POWER,step.up=TRUE,pos.side=TRUE,
                print.steps=FALSE)$root
            irootNeg<-uniroot.integer(EstRootFunc,c(0,2^(STEP.POWER+1)),
                step.power=STEP.POWER,step.up=TRUE,pos.side=FALSE,
                print.steps=FALSE)$root
            Est<- 0.5*theta(irootPos) + 0.5*theta(irootNeg)
            
 
        ## now to root function for CIs
        rootfunc<-function(i,Alt="less"){
            ytemp<- round(y+ theta(i),tieDigits)
            alpha - pvalue(wilcoxsign_test(x~ytemp,zero.method="Pratt",alternative=Alt,distribution=exact()))
        }
    
        if (alt=="two.sided") alpha<-alpha/2
        if (alt=="greater" | alt=="two.sided"){
            # using the Pratt method, values tied for the lowest 
            # difference will be zero if the shift is min(d)
            # so if m=number greater than lowest value
            # the one-sided p-value for theta0=min(d)
            # is 1/2^m
            # 
            m<-length(d[d>min(d)])
            minpval<-(1/2^m)
            # if alpha<minpval then ci[1] remains -Inf
            if (alpha==minpval){
               ci[1]<-min(d)
            } else if (alpha>minpval){
                #ci[1]<-uniroot(rootfunc,c(min(d)+epsilon,max(d)-epsilon))$root
                # We are stepping up from the minimum value
                # p values are getting larger as i increases
                # we start out with alpha-p >0      i.e., alpha>p
                # we want to get to inf(alpha-p)>0  i.e., alpha>sup(p)
                iroot<-uniroot.integer(rootfunc,c(0,2^(STEP.POWER+1)),
                       step.power=STEP.POWER,step.up=TRUE,pos.side=TRUE,
                       print.steps=FALSE,Alt="greater")$root
                ci[1]<-theta(iroot)
            }
        }
        if (alt=="less" | alt=="two.sided"){
            # using the Pratt method, values tied for the highest 
            # difference will be zero if the shift is max(d)
            # so if m=number less than highest value
            # the one-sided p-value for theta0=max(d)
            # is 1/2^m
            # 
            m<-length(d[d<max(d)])
            minpval<-(1/2^m)
            # if alpha<minpval then ci[2] remains Inf
            if (alpha==minpval){
               ci[2]<-max(d)
            } else if (alpha>minpval){
                #ci[2]<-uniroot(rootfunc,c(min(d)+epsilon,max(d)-epsilon),Alt="greater")$root
                # We are stepping down from the maximum value
                # p values are getting larger as i decreases
                # we start out with alpha-p >0      i.e., alpha>p
                # we want to get to inf(alpha-p)>0  i.e., alpha>sup(p)
                iroot<-uniroot.integer(rootfunc,c(0,2^(STEP.POWER+1)),
                       step.power=STEP.POWER,step.up=FALSE,pos.side=TRUE,
                       print.steps=FALSE,Alt="less")$root
                ci[2]<-theta(iroot)
            }
        }
        # round ci to the number of significant digits set by digits argument
        ci<-round(ci,digits)
    } else {
      # conf.int=FALSE
      ci<-c(NA,NA)
    }
    attr(ci,"conf.level")<-conf.level
    
    names(Est)<-"median (Hodges-Lehmann estimator)"
    xstar<- round(x-mu,tieDigits)
    ystar<- round(y,tieDigits)
    out<-wilcoxsign_test(xstar~ystar,zero.method="Pratt",alternative=alt,distribution=exact())
    names(mu)<-"median"
    S3out<-list(
        statistic=NULL,
        method="Exact Wilcoxon Signed-Rank Test (with Pratt modification if zeros)",
        estimate=Est,conf.int=ci,
        p.value=pvalue(out),
        null.value=mu,
        data.name=dname,
        alternative=alt
        )
     class(S3out)<-"htest"
     S3out
            
}
