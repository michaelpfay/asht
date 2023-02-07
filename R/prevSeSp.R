prevSeSp <-
function(AP,nP,Se,nSe,Sp, nSp, conf.level=0.95, neg.to.zero=TRUE){
  # NOtation from Lang and Reiczigel, 2014, 
  # Preventative Veterinary Medicine 113:13-22.
  # Prevalence=P (eq 2)
  P<- (AP + Sp -1)/(Se+Sp-1)
  if (conf.level<0.90 | conf.level>0.99) warning("method tested mostly on 95% CIs")

  if (Se<0.50 | Sp<0.5)  warning("method developed for Se>=0.50 and Sp>=0.50") 
  
  # Use "_" for prime notation in paper
  nSe_ <- nSe+2
  nSp_ <- nSp+2
  Se_ <-  (nSe * Se + 1)/(nSe+2)
  Sp_ <- (nSp*Sp +1)/(nSp+2)
  Zcrit<- qnorm(1-(1-conf.level)/2)
  # see eq 6
  nP_  <- nP + Zcrit^2
  # APprime, see eq 7
  AP_ <- (nP*AP + Zcrit^2/2)/(nP+Zcrit^2) 
  P_ <- (AP_ + Sp_ -1)/(Se_+Sp_ -1)
  # eq 17
  dP<- 2*Zcrit^2*( P_*Se_*(1-Se_)/nSe_ - (1-P_)*Sp_*(1-Sp_)/nSp_)
  # VarP_, see eq 4 with primes replacing equivalent 
  # (see text after eq 15)
  # (I think there is a typo, it should be "(1-P_)^2" instead of "(1+P_)^2")
  VarP_ <-  (AP_*(1-AP_)/nP_ + 
            P_^2*Se_*(1-Se_)/nSe_ + 
            (1-P_)^2 * Sp_*(1-Sp_)/nSp_)/
            (Se_ + Sp_ -1)^2
  # CI defined in eq 16
  ci<- P_ + dP + c(-1,1)*Zcrit*sqrt(VarP_)
  if (ci[1]<0) ci[1]<-0
  if (ci[2]>1) ci[2]<-1
  attr(ci,"conf.level")<-conf.level
  if (neg.to.zero){
    if (P<0) P<-0
    if (ci[1]<0) ci[1]<-0
  }
  estimate<- P
  names(estimate)<-"adjusted prevalence"
  data<- paste0("Unadjusted prevalence=",round(AP,4), "(nP=",nP,")")
  statistic<-Se
  names(statistic)<-paste0("Sensitivity (using nSe=",nSe,")")
  parameter<-Sp
  names(parameter)<-paste0("Specificity (using nSp=",nSp,")")
  method<-"Prevalence Adjusted for Sensitivity and Specificity (CI by Lang-Reiczigel method)"
  if (neg.to.zero)  method<-"Prevalence Adjusted for Sensitivity and Specificity (CI by Lang-Reiczigel method, with negatives set to zero)"
  output<-list(estimate=estimate,statistic=statistic,parameter=parameter,conf.int=ci, data.name=data, method=method)
  class(output)<-"htest"
  output
}
