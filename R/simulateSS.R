simulateSS<-function(decFunc,dataGenFunc,nstart=100,numBatches=100,repsPerBatch=100,power=0.30, 
                     alpha=0.025,nrepeatSwitch=3,printSteps=TRUE){
  
  if (nrepeatSwitch>5 | nrepeatSwitch<2) stop("nrepeatSwitch must be one of 2,3,4, or 5")

  if (numBatches<5) stop("numBatches must be at least 5")
  # this restriction is somewhat arbitrary
  if (repsPerBatch<10) stop("repsPerBatch must be at least 10")
  # keep track of all the sample sizes tried, N, and proportion rejected for each, P
  N<-P<-rep(NA,numBatches)
  
  # define function to estimate power from one batch only
  doBatch<-function(Ni,reps=repsPerBatch){
    dec<-rep(NA,reps)
    for (i in 1:reps){
      d<-dataGenFunc(Ni)
      dec[i]<-decFunc(d)
    }
    powerHat<- sum(dec)/reps
    # cannot have simulated power equal to 0 or 1, then the normal approximation will not work
    # change it by adding or subtracting 1/2 to the count
    if (sum(dec)==0){ powerHat<- 0.5/reps
    } else if (sum(dec)==reps){ powerHat<- (reps-0.5)/reps }
    powerHat
  }
  
 # estimate power from first batch
  N[1]<- nstart
  P[1]<- doBatch(nstart)
  if (printSteps) print(paste("Step 1:", "n=",N[1],"simulated power(n)=",P[1]))
  # estimate standardized effect, (u1-u0)/tau, from  normal approximation
  # basically set P[1] = 1- pnorm( qnorm(1-alpha) + sqrt(N[1])* SEFF )
  #  or equivalently set N[1] = (qnorm(1-alpha) + qnorm(P[1]))^2/SEFF^2
  
  SEFF<- (qnorm(1-alpha)+qnorm(P[1]))/sqrt(N[1])
  if (SEFF<0) stop("function defined so that standardized effect should be positive")
  # Now use the normal aproximation to estimate the sample size given the estimated SEFF
  Nnorm<-  (qnorm(1-alpha)+qnorm(power))^2/SEFF^2
  # Next we do simulation for Nnorm/2, Nnorm, and 2*Nnorm
  N[2:4]<-c(Nnorm/2,Nnorm,2*Nnorm)
  for (i in 2:4){
    P[i]<-doBatch(N[i])
    if (printSteps) print(paste("Step ",i,":", "n=",N[i],"simulated power(n)=",P[i]))
  }
  # Now use isotonic regression to nonparametrically estimate the power
  estNpower<-function(Nc,Pc){
    # Nc= vector of N[i] values already calculated
    # Pc= vector of P[i] values already calculated
    irout<-isoreg(Nc,Pc)  
    # predict by linear interpolation the value of N that gives P=power
    # the fitted values from the isoreg list by the ordered Nc values
    # so we have to reorder the Nc values=irout$x so that they match the y fitted
    X<- irout$yf
    # correct error: Sept 15, 2022, irout$ord is NULL if already ordered
    #Y<- irout$x[irout$ord]
    if (irout$isOrd){
      Y<-irout$x
    } else {
      Y <- irout$x[irout$ord]
    } 
    # also since approx predicts y not x, we change the names
    imin<- max((1:length(X))[X< power])
    imax<- min((1:length(X))[X>power])
    X2<- X[c(imin,imax)]
    Y2<- Y[c(imin,imax)]
    
    Npower<- approx(X2,Y2,xout=power)$y
    Npower
  }  
  if (max(P[1:4])<power | min(P[1:4])>power)  stop("initial estimates based nstart of on double and half the 
                                                   estimated normal approx (based on nstart) did not give 
                                                   simulated power greater and less than 'power'. 
                                                   Perhaps try increasing 'repsPerBatch' ")
  
  # if we use estNpower and it gets a too high estimate early on, it may get stuck too high
  # so if we get 3 in a row with the same, then switch to the up-and-down method for the rest
  # i.e., subtract 1 to the current estimate of N if previous power estimate is greater than power, 
  # otherwise add 1
  upDown<-FALSE
  upDownFunc<-function(Ni,Pi){ Nnext<-ifelse(Pi>power,Ni-1,Ni+1); Nnext }
  for (i in 5:numBatches){
    if (upDown){
      N[i]<-upDownFunc(N[i-1],P[i-1])
    } else {
      # for the intermediate calculates use round, to get the best estimates (nearest integer)
      N[i]<- round(estNpower(N[1:(i-1)],P[1:(i-1)]))
    }
    P[i]<- doBatch(N[i])
    if (printSteps) print(paste("Step ",i,":", "n=",N[i],"simulated power(n)=",P[i]))
    # if nrepeatSwitch (e.g., 3)  N in a row are the same, switch to up-and-down method
    if (!upDown && all(diff(N[(i-nrepeatSwitch+1):i])==rep(0,nrepeatSwitch-1))) upDown<-TRUE
  }
  # for the final Nestimate use ceiling instead of round, since we usually want to be 
  # conservative in the sample size estimate
  Nstar<- estNpower(N,P)
  list(N=N,P=P,Nstar=Nstar,Nestimate=ceiling(Nstar))
}

