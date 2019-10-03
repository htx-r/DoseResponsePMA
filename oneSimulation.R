OneSimulation <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
  sim.data <- simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines)
  if(splines==FALSE){
  # 1. Freq: dosresmeta
  linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                               se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

  # 2. Bayes Normal: jags
  jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelNorLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  # 3. Bayes Binomial:jags
if(OR){
  linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
  }else{
    linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','beta','tau'),model.file = modelBinLinearDRmetaRR,
                                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
  }

  # Results
  f <-coef(linearDRmetaFreq)[1]
  bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
  bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
  tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
  tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau

  sig.testNor <- ifelse(linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
  sig.testBin <- ifelse(linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
  sig.testF <- ifelse(summary(linearDRmetaFreq)$coefficients[,'Pr(>|z|)']<0.05,1,0)
rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),tauN=tn,tauB=tb,sig.testBin=sig.testBin,sig.testNor=sig.testNor,sig.testF =sig.testF)
  }else{#
    # 1. Freq: dosresmeta
    rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                                     se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

    # 2.Bayes Normal: jags
    jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=T,new.dose.range = c(1,10))

    rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelNorSplineDRmeta,
                                             n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
    # 3.Bayes Binomial: jags
    if(OR==TRUE){
    splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                              n.chains=2,n.iter = 10000,n.burnin =500,DIC=F,n.thin = 1)
    }else{
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

    }
    # beta1.pooled
  f1 <-coef(rcsplineDRmetaFreq)[1]
  b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
  b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled

  t1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau1
  t1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau1

  sig.testNor1 <- ifelse(rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','2.5%']*rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','97.5%']>0, 1,0)
  sig.testBin1 <- ifelse(splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','2.5%']*splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','97.5%']>0, 1,0)
  sig.testF1 <- ifelse(summary(rcsplineDRmetaFreq)$coefficients['dose1','Pr(>|z|)']<0.05,1,0)

  # beta2.pooled
  f2 <-coef(rcsplineDRmetaFreq)[2]
  b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
  b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled

  t2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau2
  t2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau2

  sig.testNor2 <- ifelse(rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','2.5%']*rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','97.5%']>0, 1,0)
  sig.testBin2 <- ifelse(splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','2.5%']*splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','97.5%']>0, 1,0)
  sig.testF2 <- ifelse(summary(rcsplineDRmetaFreq)$coefficients['dose2','Pr(>|z|)']<0.05,1,0)

  rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),tauN1=t1n,tauB1=t1b,sig.testBin1=sig.testBin1,sig.testNor1=sig.testNor1,sig.testF1 =sig.testF1,
            BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),tauN2=t2n,tauB2=t2b,sig.testBin2=sig.testBin2,sig.testNor2=sig.testNor2,sig.testF2 =sig.testF2)
  }
  return(rval)

}
