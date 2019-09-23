
simulateDRsplinedata.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){ #
  #This function simulate ns studies for dose of doserange with splines. The coefficients from the two dose transformations
  # are beta1 and beta2 and we assume a single tau for their RE
  # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
  # samplesize is the average sample size per study arm

  library(rms)

  #Create the dose and its spline transformations
  d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2]),2),nrow=ns))##
  d<-t(apply(d,1,sort))

  knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
  trans.d<-rcs(c(t(d)),knots)

  #nr of observations: I assume each study has 0 dose and 3 levels of dose
  nobs<-ns*3

  #the event rate in the zero dose has a maximum limit at p0<1/RR
  maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2])#
  maxRR<-exp(maxlogRR)
  p0<-0.5/maxRR#set p0 to be half the maximum allowed, just to be on the safe side!


  #create the dose-specific logRR, cases and controls
  # Underlying RR
  beta1<-c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
  beta2<-c(sapply(rnorm(ns,beta2.pooled,tau),rep,3)) #random effects of the slopes
  logRR<- beta1*trans.d[,1]+beta2*trans.d[,2]   #derive study-specific logRR using regression
  p1 <-exp(logRR)*p0

  # estimated RR
  uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
  ss<-c(sapply(uniquess,rep,3)) #sample size per study arm
  cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     #events per study at zero dose
  noncases<-matrix(c(ss-cases),nrow = 3)     #events per study at zero dose
hatRR <- cases[2:3,]/cases[1,]
hatlogRR <- log(rbind(rep(1,ns),hatRR))

  #a much easier way to calculate the SE
  SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
  selogRR<-c(rbind(NA,SEhatlogRR))

  Study_No<-rep(1:ns,each=3)

  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logRR=c(hatlogRR),dose1=c(trans.d[,1]),dose2=c(trans.d[,2]),cases=as.vector(cases),noncases=as.vector(noncases),
                                    selogRR =selogRR, type=rep('cc',3*ns))


  return(list(simulatedDRdata=simulatedDRdata,knots=knots, all.dose=trans.d))
}
