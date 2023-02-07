signTest <-
  function(x,stat=c("cd","cpp","ud"),nullparm=NULL, alternative=c("two.sided","less","greater"), conf.level=0.95,...){
    c("pos-neg","prop pos")
    stat<-match.arg(stat,choices=c("cd","cpp","ud","pos-neg","prop pos"))
    if (stat=="pos-neg"){
      stat<-"cd"
      warning("stat='pos-neg' may not be used in future versions, you should use stat='cd' instead")
    } else if (stat=="prop pos"){
      stat<-"cpp"
      warning("stat='prop pos' may not be used in future versions, you should use stat='cpp' instead")
    }
    
    dname <- deparse1(substitute(x))
    # the null parameter in the binom.exact is called 'p', but we do not want to pass that to binom.exact because it 
    # may be the wrong parameter associated with the stat argument, instead we want to use the nullparm argument
    dots<- match.call(expand.dots=FALSE)$"..."
    
    if (any(names(dots)=="p")){
      if (stat=="cpp"){
        stop("use 'nullparm' instead of 'p' ")               
      } else {
        stop("the 'p' argument does not agree with the stat argument")
      }
    }
    if (any(names(dots)=="midp")){
      midp<- dots$midp
    } else {
      midp<-FALSE
    }
    
    # check that the nullparm argument makes sense
    if ((stat=="cd" | stat=="ud") & !is.null(nullparm)){
      if (nullparm>1 | nullparm< (-1)) stop("nullparm out of range of possible values")
    } else if (stat=="cpp" & !is.null(nullparm)){
      if (nullparm>1 | nullparm<0) stop("nullparm out of range of possible values")
    }
    
    
    if (is.null(nullparm)){
      if (stat=="cd" | stat=="ud"){
        # default null for difference in proportions is 0
        nullparm<-0
      } else {
        # default null for conditional proportion positive is 0.5
        nullparm<- 0.5
      }
    }
    
    n<-length(x)
    yg<-length(x[x>0])
    yless<-length(x[x<0])
    m<-yg+yless
    if (stat=="cpp"){
      statName<-"conditional proportion positive"
      b<-binom.exact(yg,m,p=nullparm,alternative=alternative, conf.level=conf.level,...)
      #b$statistic<- yg/m
      #names(b$statistic)<-"conditional proportion positive"
      #b$parameter<- yless/m
      #names(b$parameter)<-"conditional proportion negative"   
      names(b$null.value)<-statName
      names(b$estimate)<-"conditional proportion positive (n.pos/n.nonzero)"
      b$method<-"Exact Sign Test" 
      if (midp){ b$method<- "mid-P Sign Test" }
      # I thought about calling it a Generalized sign test when the null value is not the default, but decided against it
      #if (b$null.value !=0.5) b$method<- paste0("Generalized ",b$method)
    } else if (stat=="cd"){
      statName<-"conditional difference: (prop pos)-(prop neg)"
      tr<-function(p){ 2*p-1 }
      trinv<-function(diff){ (1+diff)/2 }
      b<-binom.exact(yg,m,p=trinv(nullparm),alternative=alternative,conf.level=conf.level,...)
      #b$statistic<- yg/m
      #names(b$statistic)<-"conditional proportion positive"
      #b$parameter<- yless/m
      #names(b$parameter)<-"conditional proportion negative"   
      b$conf.int<-tr(b$conf.int)
      #  null.value=nullparm=tr(p)= tr(trinv(nullparm))
      b$null.value<-tr(b$null.value)
      b$method<-"Exact Sign Test"                             
      if (midp){ b$method<- "mid-P Sign Test" }
      # I thought about calling it a Generalized sign test when the null value is not the default, but decided against it
      #if (b$null.value !=0) b$method<- paste0("Generalized ",b$method)
      names(b$null.value)<-statName
      b$estimate<-c(yg/m,yless/m,tr(b$estimate))
      names(b$estimate)<-c("n.pos/n.nonzero","n.neg/n.nonzero","cond diff: pos-neg")
      #names(b$estimate)<-statName
    }else if (stat=="ud"){
      statName<-"unconditional difference: (prop pos)-(prop neg)"
      # because there is no midp version of menemarExactDP give error if midp tried
      if (midp){ stop("no midp version for stat='ud' ") }
      b<- mcnemarExactDP(n=n,m=m,x=yg, nullparm=nullparm, alternative=alternative, conf.level=conf.level,...)
      #b$statistic<- c(yg/n, yless/n, (n-yg-yless)/n)
      #names(b$statistic)<-c("prop pos","prop neg","prop zero")
      #b$parameter<- n
      #names(b$parameter)<-"n"   
      b$method<-"Exact Sign Test" 
      # I thought about calling it a Generalized sign test when the null value is not the default, but decided against it
      #if (b$null.value !=0) b$method<- paste0("Generalized ",b$method)
      names(b$null.value)<-statName
      b$estimate<- b$estimate
      names(b$estimate)<-c("n.pos/n","n.neg/n","uncond diff: pos-neg")
    }
    # if  tsmethod or some other argument (except midp) is added using "..." put that in description
    if (!is.null(dots)){
      names.dots<- names(dots)
      if (!midp){
        dotsDescription<- paste(names(dots),dots,sep="=", collapse = ",")                   
        b$method<-paste0(b$method,", using (",dotsDescription,")")                
      } else if (midp & length(names.dots)>1){
        # if midp=TRUE, take out midp=TRUE from dotsDescription
        midp.index<- c(1:length(names.dots))[names.dots %in% "midp"]
        dotsDescription<- paste(names(dots)[-midp.index],dots[[-midp.index]],sep="=", collapse = ",")                   
        b$method<-paste0(b$method,", using (",dotsDescription,")")                
      } # else if (midp & length(names.dots)==1) then no need to add any extra description to method
      
      
    }       
    b$statistic<- c(yg, yless, (n-yg-yless), yg+yless)
    names(b$statistic)<-c("n.pos","n.neg","n.zero","n.nonzero")
    b$data.name<-dname
    
    
    b
  }

#library(exact2x2)
#library(exactci)
#set.seed(1)
#times<-rpois(100,5) - 4
#signTest(times, stat='cd')
#signTest(times, stat='cpp')
#signTest(times, stat='ud')
#signTest(x,stat="pos-neg")
#signTest(x,stat="prop pos")
#signTest(x,stat="cd", alternative = "two.sided", nullparm=0.5, midp=TRUE, tsmethod="minlike")
