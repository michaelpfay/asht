##  one-way anova function
anovaOneWay<-function(y,g,var.equal=TRUE,nullValue=0,parm=c("ICC","varb"),conf.level=0.90){
  dname <- paste("response=",deparse(substitute(y)), " group=", deparse(substitute(g)))
  parm <- match.arg(parm)
  g<-as.character(g)
  ug<-unique(g)
  k<-length(ug)
  n<-length(y)
  if (length(g)!=n) stop("length of y and g must be the same")
  ua<-na<-s2a<-rep(NA,k)
  for (i in 1:k){
    na[i]<- length(y[g==ug[i]])
    if (na[i]<2) warning("not sure if this works for na[i]=1")
    ua[i]<-mean(y[g==ug[i]])
    s2a[i]<-var(y[g==ug[i]])
  }
  wa<-na/n
  u<- sum(wa*ua)
  df1<- k-1
  
  if (var.equal){
    df2<- n-k    
    sigma2hat<- sum( (na-1)*s2a )/(n-k)
  } else {
    # see Brown and Forsythe, 1974 Biometrics, 719-724
    ci<- (1- na/n)* s2a
    ci<- ci/sum(ci)
    df2<- ( sum(ci^2/(na-1)))^(-1)
    sigma2hat<- sum( (1-wa)*s2a )/(k-1)
  }
  MSTR<- n*sum(wa*(ua-u)^2)/df1
  TF<- MSTR/sigma2hat
  
  p.value<- 1- pf(TF,df1,df2)
  
  rootfunc<-function(t2,q=0.05){
    pf(TF,df1,df2,n*t2) - q
  }
  alpha<- 1-conf.level
  
  if (1-p.value< alpha/2){ 
    # for TF very low, we can get upper limits equal to 0
    tau2.over.sig2.U<-0
  } else {
    tau2.over.sig2.U<-uniroot(rootfunc,c(10^-8, 10^8),q=alpha/2)$root
  }
  if (p.value>=alpha/2){ 
     tau2.over.sig2.L<-0
  } else {
    tau2.over.sig2.L<- uniroot(rootfunc,c(10^-8, 10^8),q=1-alpha/2)$root
  }
  tau2hat.over.sig2<- (1/n)*(TF*df1*(df2-2)/df2 - df1)
  tau2hat.over.sig2<- max(0,tau2hat.over.sig2)
  #tau2hat.over.sig2<- (df1/n) * TF
  names(TF)<-"F"
  df<-c(df1,df2)
  names(df)<-paste0("df",1:2)
  ci<-c(tau2.over.sig2.L,tau2.over.sig2.U)
  attr(ci,"conf.level")<-conf.level
  est<- tau2hat.over.sig2
  names(est)<-parm
  nullParm<- nullValue
  names(nullParm)<-parm
  if (nullValue!=0) stop("p-value does not calculate for nullValue!=0. To be programmed in the future")
  if (parm=="varb"){
    ci<- ci*sigma2hat
    est<- est*sigma2hat
  } else if (parm=="ICC"){
    est<- est/(est+1)
    ci<- ci/(ci+1)
  }
  
  
  out<-list(statistic=TF, parameter=df,p.value=p.value,
       conf.int=ci,
       estimate=est,
       null.value=nullParm,
       alternative="greater",
       method="One-way ANOVA",
       data.name=dname)
  class(out)<-"htest"
  out
}




#na<-c(400,900,1200,300,400)
#ua<- 1:5/40
#n<-sum(na)
#k<-length(na)

#set.seed(1)
#Ua<-rep(ua,na)
#y<-rnorm(n)+Ua
#g<- as.factor(rep(1:k,na))
#anovaOneWay(y,g,var.equal=TRUE,parm="ICC")
#anovaOneWay(y,g,var.equal=TRUE,parm="varb")
#t.test(y~g)

#lout<-lm(y~g)
#a<- anova(lout)
#a



