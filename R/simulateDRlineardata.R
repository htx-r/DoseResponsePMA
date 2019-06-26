# This function generate a dataset based on linear dose-response model.
# Arguments:
# beta.pooled: the underlying slope of the linear dose-response model
# tau: the heterogenity across studies.
# ns: number of studies
# doserange: the range of the doses in the generated dataset.
# sample size: it is not the actual sample size but it is used to draw a sample size for each dose from a uinform distribution around this value(200), i.e. Unif(samplesize-20,samplesize+20)
      #  Within each study, each dose assumed to have the same drawn sample size


simulateDRlineardata.fun=function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){ #
  #simulate ns studies for dose of doserange with linear. The linear coefficient is beta and we assume a tau for their RE.
  # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
  # samplesize is the average sample size per study arm


  #Create the dose
  d<-cbind(rep(0,ns),matrix(round(seq(doserange[1],doserange[2],l=2*ns)),nrow=ns))

  #nr of observations: I assume each study has 0 dose and 3 levels of dose
  nobs<-ns*3

  #the event rate in the zero dose has a maximum limit at p0<1/RR
  maxlogRR<-(beta.pooled+2*tau)*max(d)#the maximum possible value of logRR
  maxRR<-exp(maxlogRR)
  p0<-0.5/maxRR #set p0 to be half the maximum allowed, just to be on the safe side!

  #create the dose-specific logRR, cases and controls
  beta<-c(sapply(rnorm(ns,beta.pooled,tau),rep,3)) #random effects of the slopes
  uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
  cases0<-rbinom(ns,uniquess,p0)#events per study at zero dose
  ss<-c(sapply(uniquess,rep,3)) #sample size per study arm # rep(uniquess, each=3)
  c0<-c(sapply(cases0,rep,3))
  logRR<-beta*d #derive study-specific logRR using regression
  RR<-exp(logRR)
  pevent<-c0*RR/ss #calculate the event rate in non-zero doses using the dose- and study-specitic RR
  #pevent<-p0/RR
  cdose<-c()
  for(i in 1:nobs){cdose[i]<-rbinom(1,ss[i],as.numeric(pevent[i]))} #calculate the number of events in non-zero doses

  cases<-c0*(logRR==0)+cdose*(logRR!=0)  #merge events in zero and non-zero studies

  Study_No<-c(sapply(1:ns,rep,3))  # another option rep(1:ns,each=3)

  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logRR=c(t(logRR)),dose=c(t(d)),cases=as.vector(cases),noncases=as.vector(ss-cases))

  ###%%% TASNIM: I added the seLogRR to the simulated data
  simulatedDRdata$compSElogRR <-  simulatedDRdata$noncases/(simulatedDRdata$cases*ss)
  selogRR <- sapply(1:ns, function(i) c(NA,sqrt(simulatedDRdata[simulatedDRdata$Study_No==i,]$compSElogRR[1]+simulatedDRdata[simulatedDRdata$Study_No==i,]$compSElogRR[c(2,3)])),simplify = F)
  simulatedDRdata$selogRR <- unlist(selogRR)
  simulatedDRdata$type <- rep('cc',3*ns)

  return(simulatedDRdata=simulatedDRdata[,-6])
}





















#simulateDRlineardata.fun()


# simulateDRlineardata.fun=function(beta.pooled=0.2,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){ #
#   #simulate ns studies for dose of doserange with splines. The coefficients from the two dose transformations
#   # are beta1 and beta2 and we assume a single tau for their RE
#   # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
#   # samplesize is the average sample size per study arm
#
#
#   #Create the dose and its spline transformations
#   d<-cbind(rep(0,ns),matrix(round(seq(doserange[1],doserange[2],l=2*ns)),nrow=ns))
#   #d <- seq(doserange[1],doserange[2],l=ns)
#   #nr of observations: I assume each study has 0 dose and 3 levels of dose
#   nobs<-ns*3
#
#   # zero dose
#   #the event rate in the zero dose has a maximum limit at p0<1/RR
#   maxlogRR<-(beta.pooled+2*tau)*max(d)#the maximum possible value of logRR
#   maxRR<-exp(maxlogRR)
#   p0<-0.5/maxRR #set p0 to be half the maximum allowed, just to be on the safe side!
#   n0 <- n1 <- round(runif(ns,samplesize-20,samplesize+20))
#   c0 <- rbinom(ns,n0,p0)
#
#   # non-zero dose
#   beta<-cbind(rep(0,ns),matrix(rnorm(2*ns,beta.pooled,tau),ns,2,byrow=T))
#   logRR<-beta*d
#   RR<-exp(logRR)
#   pevent <- c0*RR/n0
#   Study_No<- rep(1:ns,each=3)  # another option
#   simulatedDRdata <- data.frame(Study_No=Study_No,logRR=c(t(logRR)),dose=c(t(d)),cases=NA)
#   simulatedDRdata$cases <- unlist(sapply(1:ns, function(i) simulatedDRdata$cases[simulatedDRdata$Study_No==i] <- c(c0[i],rbinom(1,n1[i],pevent[i,1]),rbinom(1,n1[i],pevent[i,2])),simplify = F))
#   simulatedDRdata$noncases <- rep(n1,each=3)-simulatedDRdata$cases
# simulatedDRdata[simulatedDRdata$dose!=0,]$cases/rep(simulatedDRdata[simulatedDRdata$dose==0,]$cases,each=2)
#   #create the dose-specific logRR, cases and controls
#
#   # beta<-c(sapply(rnorm(ns,beta.pooled,tau),rep,2)) #random effects of the slopes
#   # uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
#   # cases0<-rbinom(ns,uniquess,p0)#events per study at zero dose
#   # #ss<-c(sapply(uniquess,rep,3)) #sample size per study arm # rep(uniquess, each=3)
#   # c0<-cases0#c(sapply(cases0,rep,3))
#   # logRR<-beta*d #derive study-specific logRR using regression
#   # RR<-exp(logRR)
#   # pevent<-c0*RR/ss #calculate the event rate in non-zero doses using the dose- and study-specitic RR
#   # #pevent<-p0/RR
#   # cdose<-c()
#   # for(i in 1:nobs){cdose[i]<-rbinom(1,ss[i],as.numeric(pevent[i]))} #calculate the number of events in non-zero doses
#   #
#   # cases<-c0*(logRR==0)+cdose*(logRR!=0)  #merge events in zero and non-zero studies
#   #
#   Study_No<- rep(1:ns,each=3)  # another option
#
#   simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logRR=as.vector(logRR),dose=c(t(d)),cases=as.vector(cases),noncases=as.vector(ss-cases))
#
#   ### TASNIM: I also changed dose
#   simulatedDRdata$compSElogRR <-  simulatedDRdata$noncases/(simulatedDRdata$cases*ss)
#   selogRR <- sapply(1:ns, function(i) c(NA,sqrt(simulatedDRdata[simulatedDRdata$Study_No==i,]$compSElogRR[1]+simulatedDRdata[simulatedDRdata$Study_No==i,]$compSElogRR[c(2,3)])),simplify = F)
#   simulatedDRdata$selogRR <- unlist(selogRR)
#   simulatedDRdata$type <- rep('cc',3*ns)
#
#   return(simulatedDRdata=simulatedDRdata[,-6])
# }
# #simulateDRlineardata.fun()
#
