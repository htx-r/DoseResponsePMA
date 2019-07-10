#This function simulate ns studies for dose of doserange with splines. The coefficients from the two dose transformations
# are beta1 and beta2 and we assume a single tau for their RE
# p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
# samplesize is the average sample size per study arm


simulateDRsplinedataOR.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1){ #

  library(rms)

  #Create the dose and its spline transformations
  d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2]),2),nrow=ns))##
  d<-t(apply(d,1,sort))

  knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
  trans.d<-rcs(c(t(d)),knots)

  #nr of observations: I assume each study has 0 dose and 3 levels of dose
  nobs<-ns*3

  #the event rate in the zero dose has a maximum limit at p0<1/RR
  # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2])#
  # maxRR<-exp(maxlogRR)
  # p0<-0.5/maxRR#set p0 to be half the maximum allowed, just to be on the safe side!


  #create the dose-specific logRR, cases and controls

  beta1<-c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
  beta2<-c(sapply(rnorm(ns,beta2.pooled,tau),rep,3)) #random effects of the slopes
  uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
  cases0<-rbinom(ns,uniquess,p0)#events per study at zero dose
  ss<-c(sapply(uniquess,rep,3)) #sample size per study arm
  c0<-c(sapply(cases0,rep,3))
  logOR<- beta1*trans.d[,1]+beta2*trans.d[,2]   #derive study-specific logRR using regression
  OR<-exp(logOR)
  p0hat <- cases0/uniquess
  odds0 <- p0hat/(1-p0hat)
  p1hat <- OR*odds0/(1-(OR*odds0))
  cdose <- round(p1hat*uniquess)

  cases<-c0*(logOR==0)+cdose*(logOR!=0)  #merge events in zero and non-zero studies

  Study_No<-rep(1:ns,each=3)

  #a much easier way to calculate the SE
 cases <- matrix(as.vector(cases),ns,3,byrow = T)
 SE<-1/cases[,1]+1/cases[,c(2,3)]+1/(uniquess-cases[,1]) + 1/(uniquess-cases[,c(2,3)])
 selogOR<-c(t(cbind(NA,sqrt(SE))))

  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logOR=c(t(logOR)),dose1=c(trans.d[,1]),dose2=c(trans.d[,2]),cases=as.vector(cases),noncases=as.vector(ss-cases),
                                    selogOR =selogOR, type=rep('cc',3*ns))


  return(list(simulatedDRdata=simulatedDRdata,knots=knots, all.dose=trans.d))
}


