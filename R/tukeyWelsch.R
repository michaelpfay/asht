# stepdownOneWay
library(perm)


aovMeanVarDF<-function(y,g,var.equal=TRUE){
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
  list(na=na,ua=ua,sigma2hat=sigma2hat,df2=df2)
}
anovaOneWayFast<-function(na,ua,sigma2hat,df2){
  k<-length(na)
  n<-sum(na)
  wa<-na/n
  u<- sum(wa*ua)
  df1<- k-1
  MSTR<- n*sum(wa*(ua-u)^2)/df1
  TF<- MSTR/sigma2hat
  p.value<- 1- pf(TF,df1,df2)
  p.value
}
numberOfTests<-function(k){
  if (k<2) stop("number of groups must be at least 2")
  if (k==1){
    ntests<-0
  } else if (k==2){
    ntests<-1
  }  else {
    ntests<-1
    for (i in (k-1):2){
      ntests<-ntests+choose(k,i)
    }
  }
  ntests
}

#numberOfTests(14)

padjTW<- function(unadjP,k,subsetSize){
  i<-subsetSize
  if (i==k | i==k-1){
    adjP<- unadjP 
  } else {
    adjP<- 1- (1-unadjP)^(k/i)
  }
  adjP
}

tukeyFast<-function(na,ua,sigma2hat,df2){
  ## modify code from TukeyHSD
  center <- outer(ua, ua, "-")
  keep <- lower.tri(center)
  center <- center[keep]
  est <- center/(sqrt((sigma2hat/2) * outer(1/na, 1/na, "+"))[keep])
  pvals <- ptukey(abs(est), length(ua), df2, 
                  lower.tail = FALSE)
  min(pvals)
}



testLevel.aov.sr<-function(na,ua,sigma2hat,df2,subsetSize,padjfunc=padjTW,pvalFast=anovaOneWayFast){
  # Hochberg and Tamhane use p instead of j to mean 
  # number of elements in subset hypotheses 
  j<-subsetSize
  k<-length(na)
  chmat<-chooseMatrix(k,j)
  nlevelTests<-nrow(chmat)
  pvals<-rep(NA,nlevelTests)
  for (i in 1:nlevelTests){
    uasubset<- ua[chmat[i,]==1]
    nasubset<- na[chmat[i,]==1]
    pvals[i]<- pvalFast(nasubset,uasubset,sigma2hat,df2)
  }
  list(chmat=chmat,pvals=pvals,padj=padjfunc(pvals,k,j))
}

testLevel<-function(y,g,ug,k,subsetSize,pvalfunc,padjfunc=padjTW,...){
  # Hochberg and Tamhane use p instead of j to mean 
  # number of elements in subset hypotheses 
  j<-subsetSize
  chmat<-chooseMatrix(k,j)
  nlevelTests<-nrow(chmat)
  pvals<-rep(NA,nlevelTests)
  for (i in 1:nlevelTests){
    ugsubset<- ug[chmat[i,]==1]
    I<- g %in% ugsubset
    yi<-y[I]
    # to stop issues p-value function from giving error 
    # with levels with no elements, reduce the level set to the set in g[I]
    gi<-factor(g[I], exclude=NULL)
    pvals[i]<-pvalfunc(yi,gi,...)
  }
  list(chmat=chmat,pvals=pvals,padj=padjfunc(pvals,k,j))
}
#unique(g)
#testLevel(y,g,1:4,4,2,pvalfunc)

pairwiseSDpvals<-function(outlist,ug){
  ## take out list and get p-values for pairwise comparisons
  ## use Tukey-Welsch (p. 111) method, and get minimum alpha over
  ## all levels
  k<-length(ug)
  if (length(outlist)!=k-1) stop("outlist must have all subset size hypotheses from k,k-1,...,2")
  pairMat<- outlist[[k-1]]$chmat
  npairs<- nrow(pairMat)
  pairNames<-rep(NA,npairs)
  padj<-rep(NA,npairs)
  for (i in 1:npairs){
    # I picks out columns from the pairMat that defines each pair
    I<- pairMat[i,]==1
    pairNames[i]<- paste(ug[I],collapse = "-")
    # 
    
    padjLevel<-rep(NA,k-1)
    padjLevel[1]<-outlist[[1]]$padj
    for (j in 2:(k-1)){
      # outlist[[j]] gives results for subsetSize=k+1-j
      # so outlist[[1]] has subsetSize=k
      # and outlist[[k-1]] has subsetSize=2
      chMat<-outlist[[j]]$chmat
      # pick out all subsets that contain both elements of the pair
      # to do this
      # pick out 2 columns associated with the pair, and
      # if both are equal to 1, that means that the pair
      # is a subset of the hypothesis
      pairIn<- apply(chMat[,I],1,sum)==2
      padjLevel[j]<- max(outlist[[j]]$padj[pairIn])
    }
    padj[i]<-max(padjLevel) 
  }
  names(padj)<-pairNames
  padj
}
#pairwiseSDpvals(out,ug)



tukeyWelsch<-function(y,g,method=c("aov","kw","sr","user"),pvalfunc=NULL, 
                         padjfunc=padjTW,maxnTest=10^4,nTestMessage=FALSE,...){
  dname <- paste(deparse(substitute(y)), "and", deparse(substitute(g)))
  if (any(is.na(g)) | any(is.na(y))) stop("function is not defined for missing y or g values")
  # g needs to be a factor for the aov method to work right
  # Sept 16, 2022: change from 
  #     ###  if (class(g)!="factor"){ ###   
  #  to (as suggested by check --as-cran)...
  if (!inherits(g,"factor")){
    warning("g treated as a factor")
    g<-factor(g)
  }
  ug<- unique(as.character(g))
  k<-length(ug)
  
  method<-match.arg(method)
  if (method=="kw"){
    pvalfunc<- function(y,g){
      kruskal.test(y,g,...)$p.value
    }
  } else if (method=="aov" | method=="sr"){
    aout<-aovMeanVarDF(y,g,var.equal=TRUE)
    na<-aout$na
    ua<-aout$ua
    sigma2hat<-aout$sigma2hat
    df2<-aout$df2
  } else if (method=="user") {
    if (is.null(pvalfunc)){
      stop("must supply pvalfunc")
    }
  }  
  ntests<-numberOfTests(k)
  if (nTestMessage) message(paste0("Starting tests, requires testing ",ntests," tests..."))
  if (maxnTest<ntests) stop("number of groups uses too many tests, increase maxnTest")
  out<-list()
  for (j in k:2){
    if (method=="aov"){
      out[[length(out)+1]]<-testLevel.aov.sr(na,ua,sigma2hat,df2,j,
                                            padjfunc=padjTW,pvalFast=anovaOneWayFast)
    } else if (method=="sr") {
      out[[length(out)+1]]<-testLevel.aov.sr(na,ua,sigma2hat,df2,j,
                                             padjfunc=padjTW,pvalFast=tukeyFast)
    } else   {
      # for method=='kw' or 'user'
      out[[length(out)+1]]<-testLevel(y,g,ug,k,j,pvalfunc)
    }
  }
  
  methodlong<-switch(method,aov="one-way ANOVA",sr="studentized range test",
                     kw="Kruskal-Wallis test",user="user supplied function")
  
  METHOD<-paste0("Tukey-Welsch method with ",methodlong)
  pairPvals<-pairwiseSDpvals(out,ug)
  out<-list(fullResults=out,data.name=dname,method=METHOD,ksample.pvalue=out[[1]]$pvals, pairwise.pvalues=pairPvals)
  class(out)<-"tukeyWelsch"
  out
}



print.tukeyWelsch<-function (x, digits = max(1L, getOption("digits") - 5L), ...) 
{
  # mostly copied from print.pairwise.htest
  cat("\n\tPairwise comparisons using\n", x$method, "\n\n")
  cat("data: ", x$data.name, "\n\n")
  pp <- format.pval(x$pairwise.pvalues, digits = digits, na.form = "-")
  kp<-  format.pval(x$ksample.pvalue, digits = digits)
  attributes(pp) <- attributes(x$pairwise.pvalues)
  cat("\nPairwise p-values:\n")
  print(pp, quote = FALSE, ...)
  cat("\n K-sample p-value:", kp, "\n")
  invisible(x)
}




##################
# Test it
##################


#createData<-function(n,props,shifts,ry=rnorm){
#  k<-length(props)
#  if (round(sum(props),8)!=1) stop("sum of props must be 1")
#  props<- props/sum(props)
#  if (length(shifts)!=k) stop("length of shifts must equal length of props")
#  g<-rep(1:k,as.vector(rmultinom(1,n,props)))
#  y<-ry(n)
#  for (i in 1:k){
#    y[g==i]<-y[g==i]+shifts[i]
#    }
#  list(y=y,g=g)
#}
#set.seed(1)
#d<-createData(100,c(.2,.3,.2,.3),c(0,0,0,1))
#tukeyWelsch(d$y,factor(d$g),method="kw")
#tukeyWelsch(d$y,factor(d$g),method="aov")
#tukeyWelsch(d$y,factor(d$g),method="sr")
#TukeyHSD(aov(d$y~factor(d$g)))[[1]][,"p adj"]


