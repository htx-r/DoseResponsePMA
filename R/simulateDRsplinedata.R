simulateDRsplinedata.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,samplesize=200){ # ,doserange=c(1, 80)
  #simulate ns studies for dose of doserange with splines. The coefficients from the two dose transformations
  # are beta1 and beta2 and we assume a single tau for their RE
  # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
  # samplesize is the average sample size per study arm

library(rms)

#Create the dose and its spline transformations
d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2])),nrow=ns))
d<-t(apply(d,1,sort))
knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
trans.d<-rcs(c(t(d)),knots)

#nr of observations: I assume each study has 0 dose and 3 levels of dose
nobs<-ns*3

#the event rate in the zero dose has a maximum limit at p0<1/RR
maxlogRR<-(beta1.pooled+2*tau)*max(trans.d[,1])+(beta2.pooled+2*tau)*max(trans.d[,2]) #the maximum possible value of logRR
maxRR<-exp(maxlogRR)
p0<-0.5/maxRR#set p0 to be hlf the maximum allowed, just to be on the safe side!


#create the dose-specific logRR, cases and controls

beta1<-c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
beta2<-c(sapply(rnorm(ns,beta2.pooled,tau),rep,3)) #random effects of the slopes
uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
cases0<-rbinom(ns,uniquess,p0)#events per study at zero dose
ss<-c(sapply(uniquess,rep,3)) #sample size per study arm
c0<-c(sapply(cases0,rep,3))
logRR<-beta1*trans.d[,1]+beta2*trans.d[,2] #derive study-specific logRR using regression
RR<-exp(logRR)
pevent<-c0*RR/ss #calculate the event rate in non-zero doses using the dose- and study-specitic RR

cdose<-c()
for(i in 1:nobs){cdose[i]<-rbinom(1,ss[i],as.numeric(pevent[i]))} #calculate the number of events in non-zero doses

cases<-c0*(logRR==0)+cdose*(logRR!=0)  #merge events in zero and non-zero studies

Study_No<-c(sapply(1:ns,rep,3))

simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logRR=as.vector(logRR),dose=c(t(d)),cases=as.vector(cases),noncases=as.vector(ss-cases))

return(list(simulatedDRdata=simulatedDRdata,knots=knots))
}


