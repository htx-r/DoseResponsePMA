
simulateDRmeta.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

  # This function generate a dataset with odds ratio OR or risk ratio RR based on dose-response meta-analysis model.

  # Arguments:
  # beta1.pooled: (numeric) the underlying slope of the linear dose-response model
  # beta2.pooled: (numeric) the underlying coeffiecient of cubic spline transformation of dose-response model
  # tau: the commom heterogenity for both regression coeffiecients across studies.
  # ns: number of studies
  # doserange: the range of the generated dosages.
  # sample size: the mean value of the uinform distribution that we generate the sample size on each dose i.e. Unif(samplesize-20,samplesize+20)
  # p0: is the probability of the event in the zero dose only for OR.
  # OR: (logical) indicate whether for OR (TRUE) or RR (FALSE) the data need to be generated.
  # splines: (logical) indicate whether linear or restricted cubic spline dose-response models.

  # load libraries
  require(Hmisc)

  # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
   d<-cbind(rep(0,ns),matrix(round(c(runif(2*ns,doserange[1],doserange[2])),2),nrow=ns))##
   d<-t(apply(d,1,sort))
   dose <- c(t(d))


  # 2. Create the study-specific logOR and logRR either for splines or linear
  if(splines==TRUE){ ## for splines
  # find the dose cubic transformation
    knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
  # implement the dose-response model
    maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    logrr<- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model

      }else{ ## for linear
  maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
  beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
  logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
  dose1 <- dose2 <- dose
      }


  # 3. Generate the dose-specific logOR and logRR
  if(OR==TRUE){ ## odds ratio (OR)
  # find the probabilities of the event to occur in zero dose;p0 and in the nonzero dose; p1
    odds0 <- p0/(1-p0)  # as p0 is provided as an argument we compute the odds in zero dose
    odds1 <- exp(logrr)*odds0 # from above we get the study-specific OR or RR
    p1<- odds1/(1+odds1) #compute p1 the probability of a case at any dose level

  # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size of each study arm
    ss<-c(sapply(uniquess,rep,3)) # sample size for all study arms
    cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) # generate the cases for all dose levels
    noncases<-matrix(c(ss-cases),nrow=3) # calculate the noncases = n - cases
    #%%%!!!!! To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1

  # compute the estimated odds ratio (hatOR) and its standard error for each dose level
    hatodds<-cases/noncases # hat odds
    hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #  log hat OR
    SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of log hat OR
    SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR) # add NA for SE in the zero dose1

  # the final results go to the the returned data
    hatlogrr <- hatlogOR
    SEhatlogrr <- SEhatlogOR
    type=rep('cc',3*ns) # type of the data is cc= case-control for OR. This is needed in dosresmeta function to choose how to compute the SE
  }else{ ## Risk ratio (RR)

  # restrictions to ensure that we will not get very large or small value for maxRR and so p0 and p1
maxRR<-exp(maxlogRR)
if(maxRR>10){
  p0 <- 0.05
  p1 <-ifelse(exp(logrr)*p0>1,0.95,exp(logrr)*p0)
}else if(maxRR<0.5){
p0 <- 0.95
p1 <-ifelse(exp(logrr)*p0<1,0.05,exp(logrr)*p0)
}else{
  p0 <- 0.5/maxRR
  p1 <-exp(logrr)*p0
}

   # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose

    # compute the estimated risk ratio (hatRR) and its standard error for each dose level
    hatRR <- cases[2:3,]/cases[1,]
    hatlogRR <- log(rbind(rep(1,ns),hatRR))
    SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
    selogRR<-c(rbind(NA,SEhatlogRR))

    # the final results go to the the returned data object
    hatlogrr <- hatlogRR
    SEhatlogrr <- selogRR
    type=rep('ci',3*ns) # type of the data is ci= cumulative incidance for RR. This is needed in dosresmeta function to choose how to compute the SE
  }


  Study_No<-rep(1:ns,each=3) # label the studies in numbers

  # a data frame of the simulated data that we need to get from this function
  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
                                    selogrr =c(SEhatlogrr), type=type)

  return(simulatedDRdata=simulatedDRdata)
  }
# End

































# d<-cbind(rep(0,ns),matrix(round(c(runif(ns,doserange[1],doserange[2]/2),runif(ns,doserange[2]/2,doserange[2])),2),nrow=ns))##
# d<-t(apply(d,1,sort))
# dose <- c(t(d))
# d1 <- round(c(rchisq(2*ns,5)),2)
# d1[d1>10] <- 10
# d <- cbind(rep(0,ns),matrix(d1,nrow=ns))
#
# d<-t(apply(d,1,sort))
# dose <- c(t(d))

#    cbind(d[,2:3],d2[,2:3])
#
#   par(mfrow=c(1,2))
# plot(1:40,d[,3],ylim = c(0,10))
# points(1:40,d[,2],col=2)
# abline(v=1:40,lty=1,col=1:2)
#
# plot(1:40,d2[,3],ylim = c(0,10))
# points(1:40,d2[,2],col=2)
# abline(v=1:40,lty=1,col=1:2)

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

