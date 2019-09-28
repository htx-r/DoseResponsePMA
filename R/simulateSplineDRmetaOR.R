simulateSplineDRmetaOR.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1){ #

  # This function simulates ns studies for dose of doserange with splines. The coefficients from the two dose transformations
  # are beta1 and beta2 and we assume a single tau for their random effects
  # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
  # samplesize is the average sample size per study arm
  # The function assumes 3 dose levels per study, the lowest of which is zero.

  library(rms)

  #Create the dose and its spline transformations
  d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2]),2),nrow=ns))##
  d<-t(apply(d,1,sort))

  knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
  trans.d<-rcs(c(t(d)),knots)

  #nr of observations: I assume each study has 0 dose and 3 levels of dose
  nobs<-ns*3


  #create the dose-specific logOR, cases and controls
  beta1<-c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
  beta2<-c(sapply(rnorm(ns,beta2.pooled,tau),rep,3)) #random effects of the slopes
  logOR <- beta1*trans.d[,1]+beta2*trans.d[,2]   #derive study-specific logOR using regression
  odds0 <- p0/(1-p0)
  odds1 <- exp(logOR)*odds0
  p1<- odds1/(1+odds1)#estimate p1 the probability of a case at any dose level

  uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size of each study arm
  ss<-c(sapply(uniquess,rep,3))#sample size for all study arms
  cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) #drow the cases for all levels
  noncases<-matrix(c(ss-cases),nrow=3)
  hatodds<-cases/noncases
  hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #sample logOR
  SEhatlogOR<-1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),]) #SE of the sample logOR
  SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR)

  Study_No<-rep(1:ns,each=3)


  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logOR=c(hatlogOR),dose1=c(trans.d[,1]),dose2=c(trans.d[,2]),cases=c(cases),noncases=c(noncases),
                                    selogOR =c(SEhatlogOR), type=rep('cc',3*ns))


  return(list(simulatedDRdata=simulatedDRdata,knots=knots, all.dose=trans.d))
}


