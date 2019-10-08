
simulateDRmeta.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

  # This function generate a dataset based on linear dose-response model.
  # Arguments:
  # beta.pooled: (vector) the underlying slope of the linear dose-response model and rcs trnasformation
  # tau: the heterogenity across studies.
  # ns: number of studies
  # doserange: the range of the doses in the generated dataset.
  # sample size: it is not the actual sample size but it is used to draw a sample size for each dose from a uinform distribution around this value(200), i.e. Unif(samplesize-20,samplesize+20)
  #  Within each study, each dose assumed to have the same drawn sample size
library(rms)
  # 1. Create the doses and its transformation
  d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2]),2),nrow=ns))##
  d<-t(apply(d,1,sort))
  dose <- c(t(d))

  #nr of observations: I assume each study has 3 levels of doses
  nobs<-ns*3


  # 2. Create the dose-specific logOR,
  if(splines){

    knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
    trans.d<-rcs(c(t(d)),knots)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])

    maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2])#
    beta1<-c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
    beta2<-c(sapply(rnorm(ns,beta2.pooled,tau),rep,3)) #random effects of the slopes
    logrr<- beta1*trans.d[,1]+beta2*trans.d[,2]   #derive study-specific logOR using regression

      }else{
  maxlogRR<- (beta1.pooled+2*tau)*max(dose) #the maximum possible value of logRR
  beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes
  logrr <- beta*dose  #derive study-specific logOR using regression
  dose1 <- dose2 <- dose
      }


  # 3.
  if(OR){ ## odds ratio
    odds0 <- p0/(1-p0)
    odds1 <- exp(logrr)*odds0
    p1<- odds1/(1+odds1)#estimate p1 the probability of a case at any dose level

    uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size of each study arm
    ss<-c(sapply(uniquess,rep,3))#sample size for all study arms
    cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) #drow the cases for all levels
    noncases<-matrix(c(ss-cases),nrow=3)
    hatodds<-cases/noncases
    hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #sample logOR
    SEhatlogOR<-1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),]) #SE of the sample logOR
    SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR)

  hatlogrr <- hatlogOR
  SEhatlogrr <- SEhatlogOR
  type=rep('cc',3*ns)
  }else{ ## Risk ratio
    maxRR<-exp(maxlogRR)
    p0<-0.5/maxRR#set p0 to be half the maximum allowed, just to be on the safe side!
    p1 <-exp(logrr)*p0

    uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) #sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     #events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     #events per study at zero dose
    hatRR <- cases[2:3,]/cases[1,]
    hatlogRR <- log(rbind(rep(1,ns),hatRR))

    #a much easier way to calculate the SE
    SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
    selogRR<-c(rbind(NA,SEhatlogRR))

    hatlogrr <- hatlogRR
    SEhatlogrr <- selogRR
    type=rep('ci',3*ns)
  }


  Study_No<-rep(1:ns,each=3)

  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
                                    selogrr =c(SEhatlogrr), type=type)

  return(simulatedDRdata=simulatedDRdata)
  }
#







































# This function generate a dataset based on linear dose-response model.
# Arguments:
# beta.pooled: the underlying slope of the linear dose-response model
# tau: the heterogenity across studies.
# ns: number of studies
# doserange: the range of the doses in the generated dataset.
# sample size: it is not the actual sample size but it is used to draw a sample size for each dose from a uinform distribution around this value(200), i.e. Unif(samplesize-20,samplesize+20)
#  Within each study, each dose assumed to have the same drawn sample size


# simulateDRlineardataOR.fun=function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1){ #
#   #simulate ns studies for dose of doserange with linear. The linear coefficient is beta and we assume a tau for their RE.
#   # p0 is the event rate at dose 0 and it is set by the simulations according to the maximum allowed value
#   # samplesize is the average sample size per study arm
#
#
#   #Create the dose
#   d<-cbind(rep(0,ns),matrix(round(runif(2*ns,doserange[1],doserange[2]),2),nrow=ns))##
#   d<-t(apply(d,1,sort))
#
#   #nr of observations: I assume each study has 0 dose and 3 levels of dose
#   nobs<-ns*3
#
#   # the event rate in the zero dose has a maximum limit at p0<1/RR
#   # maxlogOR<-(beta.pooled+2*tau)*max(d)#the maximum possible value of logRR
#   # maxOR<-exp(maxlogOR)
#   # p0<-0.5/maxRR #set p0 to be half the maximum allowed, just to be on the safe side!
#
#   #create the dose-specific logRR, cases and controls
#   beta<-c(sapply(rnorm(ns,beta.pooled,tau),rep,3)) #random effects of the slopes
#   uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
#   cases0<-rbinom(ns,uniquess,p0)#events per study at zero dose
#   ss<-c(sapply(uniquess,rep,3)) #sample size per study arm # rep(uniquess, each=3)
#   c0<-c(sapply(cases0,rep,3))
#   logOR<-beta*d #derive study-specific logRR using regression
#   OR<-exp(logOR)
# p0hat <- cases0/uniquess
# odds0 <- p0hat/(1-p0hat)
# odds1 <- OR*odds0
# p1hat <- odds1/(1+odds1)
# #p1hat <- OR*odds0/(1-(OR*odds0))
# cdose <- round(p1hat*uniquess)
# cases<-c0*(logOR==0)+cdose*(logOR!=0)  #merge events in zero and non-zero studies
#
#   Study_No<-rep(1:ns,each=3)
#
#   #a much easier way to calculate the SE
#   SE<-sqrt(1/cases[,1]+1/cases[,c(2,3)]+1/(uniquess-cases[,1]) + 1/(uniquess-cases[,c(2,3)]))
#   selogOR<-c(t(cbind(NA,sqrt(SE))))
#
#   simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logOR=c(t(logOR)),dose=c(t(d)),cases=as.vector(t(cases)),noncases=as.vector(t(ss-cases))
#                                     ,selogOR=selogOR, type=rep('cc',3*ns))
#
#   return(simulatedDRdata=simulatedDRdata)
# }

