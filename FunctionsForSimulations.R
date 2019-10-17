

OneSimulation <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
  sim.data <- simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines)
  if(splines==FALSE){
    # 1. Freq: dosresmeta
    linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                                 se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')

    # 2. Bayes Normal: jags
    jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))

    linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
                                           n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

    # 3. Bayes Binomial:jags
    if(OR){
      linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                                n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
    }else{
      linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
                                                n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
    }

    # Results
    f <-coef(linearDRmetaFreq)[1]
    bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
    bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled

    sdF <- sqrt(linearDRmetaFreq$vcov)
    sdBin <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','sd']
    sdNor <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','sd']

    tf <- sqrt(linearDRmetaFreq$Psi)
    tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau

    sig.testNor <- ifelse(linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
    sig.testBin <- ifelse(linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
    sig.testF <- ifelse(summary(linearDRmetaFreq)$coefficients[,'Pr(>|z|)']<0.05,1,0)

    RhatN <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','Rhat']
    RhatB <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','Rhat']

    rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),sdF=sdF,sdNor=sdNor,sdBin=sdBin
              ,tauN=tn,tauB=tb,tauF=tf,RhatN=RhatN,RhatB=RhatB)
  }else{#
    # 1. Freq: dosresmeta
    rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                                     se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

    # 2.Bayes Normal: jags
    jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=T,new.dose.range = c(1,10))

    rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
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

    sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
    sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
    sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']

    tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
    tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau

    RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
    RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']

    sig.testNor1 <- ifelse(rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','2.5%']*rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','97.5%']>0, 1,0)
    sig.testBin1 <- ifelse(splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','2.5%']*splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','97.5%']>0, 1,0)
    sig.testF1 <- ifelse(summary(rcsplineDRmetaFreq)$coefficients['dose1','Pr(>|z|)']<0.05,1,0)

    # beta2.pooled
    f2 <-coef(rcsplineDRmetaFreq)[2]
    b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
    b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled

    sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
    sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
    sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']

    tf2 <- sqrt(rcsplineDRmetaFreq$Psi[2,2])

    RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
    RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']

    sig.testNor2 <- ifelse(rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','2.5%']*rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','97.5%']>0, 1,0)
    sig.testBin2 <- ifelse(splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','2.5%']*splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','97.5%']>0, 1,0)
    sig.testF2 <- ifelse(summary(rcsplineDRmetaFreq)$coefficients['dose2','Pr(>|z|)']<0.05,1,0)

    rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf1,RhatN1=RhatN1,RhatB1=RhatB1,
              BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,tauF2=tf2,RhatN2=RhatN2,RhatB2=RhatB2)
  }
  return(rval)

}

library('rsimsum')
library(tidyr)


simpower <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
  res <- replicate(nsim,OneSimulation(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
  df <- data.frame(beta=c(t(res)[,c('BayesB','BayesN','Freq')]),se=c(t(res)[,c('sdBin','sdNor','sdF')]),par=rep(c('BayesB','BayesN','Freq'),each=nsim))
  ms <- multisimsum(
    data = df,
    par = "par", true = c(BayesB=beta1.pooled,BayesN=beta1.pooled,Freq=beta1.pooled),
    estvarname = "beta", se = "se",x=TRUE)
  sms <- summary(ms)
  rval <- sms$summ[sms$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]

  m <- spread(rval,par,est)
  rownames(m) <- m[,'stat']
  m <- m[,-1]
  mm <- cbind(m[1,],m[2,],m[3,],m[4,],m[5,])
  colnames(mm) <- paste0(rep(c('BayesB','BayesN','Freq'),5),rep(rownames(m),each=3))
  rownames(mm) <-NULL


  mm$mcse.biasB <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesB','mcse']
  mm$mcse.biasN <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesN','mcse']
  mm$mcse.biasF <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='Freq','mcse']

  mdf <- colMeans(t(res))
  mm$RhatB<-mdf['RhatB']
  mm$RhatN<-mdf['RhatN']
  dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])

  result <- cbind.data.frame(true.beta=beta1.pooled,mm,dftau)
  rownames(result) <- NULL

  return(result)
  #


  }


simpower()





