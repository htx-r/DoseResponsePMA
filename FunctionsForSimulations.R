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

simpower <- function(nrep=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneSimulation(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
if(splines==FALSE){
  # Biases
  biasBbin <- colMeans(t(res))['BayesB']-beta1.pooled
  biasBnor <- colMeans(t(res))['BayesN']-beta1.pooled
  biasF <- colMeans(t(res))['Freq']- beta1.pooled

  # heterogenity
  tauN.hat <- colMeans(t(res))['tauN']
  tauB.hat <- colMeans(t(res))['tauB']

  # MSE: Mean square error
  mseBbin <- mean((t(res)[,'BayesB']-beta1.pooled)^2)
  mseBnor <- mean((t(res)[,'BayesN']-beta1.pooled)^2)
  mseF <- mean((t(res)[,'Freq']-beta1.pooled)^2)

  # Type 1 error
  alphaNor <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testNor',]==1),NA)
  alphaBin <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testBin',]==1),NA)
  alphaF <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testF',]==1),NA)

  # Type 2 error
  betaNor <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testNor',]==0),NA)
  betaBin <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testBin',]==0),NA)
  betaF <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testF',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin <- sqrt(sum((t(res)[,'BayesB']-colMeans(t(res))['BayesB'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor <- sqrt(sum((t(res)[,'BayesN']-colMeans(t(res))['BayesN'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF <- sqrt(sum((t(res)[,'Freq']-colMeans(t(res))['Freq'])^2)/(ncol(res)*(ncol(res)-1)))

  rval <- cbind(beta.pooled=beta1.pooled,tau=tau, # true values
                   tauN.hat=tauN.hat,tauB.hat=tauB.hat, # estimation of tau
                   biasBnor=biasBnor,biasBbin=biasBbin,biasF=biasF, # bias of beta
                   mseBnor=mseBnor,mseBbin=mseBbin,mseF=mseF, # mean squared error for beta
                   alphaBin=alphaBin,alphaNor=alphaNor,alphaF=alphaF, # type 1 error (alpha)
                   betaBin=betaBin,betaNor=betaNor,betaF=betaF, # type 2 error (beta)
                   MCseBin=MCseBin,MCseNor=MCseNor,MCseF=MCseF) # monte carlo standard error
  row.names(rval) <- NULL
}else{
  biasBbin1 <- colMeans(t(res))['BayesB1']-beta1.pooled
  biasBnor1 <- colMeans(t(res))['BayesN1']-beta1.pooled
  biasF1 <- colMeans(t(res))['Freq1']- beta1.pooled

  # heterogenity
  tauN.hat1 <- colMeans(t(res))['tauN1']
  tauB.hat1 <- colMeans(t(res))['tauB1']

  # MSE: Mean square error
  mseBbin1 <- mean((t(res)[,'BayesB1']-beta1.pooled)^2)
  mseBnor1 <- mean((t(res)[,'BayesN1']-beta1.pooled)^2)
  mseF1 <- mean((t(res)[,'Freq1']-beta1.pooled)^2)

  # Type 1 error
  alphaNor1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testNor1',]==1),NA)
  alphaBin1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testBin1',]==1),NA)
  alphaF1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testF1',]==1),NA)

  # Type 2 error
  betaNor1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testNor1',]==0),NA)
  betaBin1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testBin1',]==0),NA)
  betaF1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testF1',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin1 <- sqrt(sum((t(res)[,'BayesB1']-colMeans(t(res))['BayesB1'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor1 <- sqrt(sum((t(res)[,'BayesN1']-colMeans(t(res))['BayesN1'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF1 <- sqrt(sum((t(res)[,'Freq1']-colMeans(t(res))['Freq1'])^2)/(ncol(res)*(ncol(res)-1)))

  #### #### #### ####
  #### beta2.pooled
  #### #### #### ####
  # Biases
  biasBbin2 <- colMeans(t(res))['BayesB2']-beta2.pooled
  biasBnor2 <- colMeans(t(res))['BayesN2']-beta2.pooled
  biasF2 <- colMeans(t(res))['Freq2']- beta2.pooled

  # heterogenity
  tauN.hat2 <- colMeans(t(res))['tauN2']
  tauB.hat2 <- colMeans(t(res))['tauB2']

  # MSE: Mean square error
  mseBbin2 <- mean((t(res)[,'BayesB2']-beta2.pooled)^2)
  mseBnor2 <- mean((t(res)[,'BayesN2']-beta2.pooled)^2)
  mseF2 <- mean((t(res)[,'Freq2']-beta2.pooled)^2)

  # Type 1 error
  alphaNor2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testNor2',]==1),NA)
  alphaBin2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testBin2',]==1),NA)
  alphaF2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testF2',]==1),NA)

  # Type 2 error
  betaNor2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testNor2',]==0),NA)
  betaBin2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testBin2',]==0),NA)
  betaF2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testF2',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin2 <- sqrt(sum((t(res)[,'BayesB2']-colMeans(t(res))['BayesB2'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor2 <- sqrt(sum((t(res)[,'BayesN2']-colMeans(t(res))['BayesN2'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF2 <- sqrt(sum((t(res)[,'Freq2']-colMeans(t(res))['Freq2'])^2)/(ncol(res)*(ncol(res)-1)))


  # Return
  rval <- c( # beta1.pooled
    beta1.pooled=beta1.pooled,tau=tau, # true values
    tauN.hat1=tauN.hat1,tauB.hat1=tauB.hat1, # estimation of tau
    biasBnor1=biasBnor1,biasBbin1=biasBbin1,biasF1=biasF1, # bias of beta
    mseBnor1=mseBnor1,mseBbin1=mseBbin1,mseF1=mseF1, # mean squared error for beta
    alphaBin1=alphaBin1,alphaNor1=alphaNor1,alphaF1=alphaF1, # type 1 error (alpha)
    betaBin1=betaBin1,betaNor1=betaNor1,betaF1=betaF1, # type 2 error (beta)
    MCseBin1=MCseBin1,MCseNor1=MCseNor1,MCseF1=MCseF1,
    # beta2.pooled
    beta2.pooled=beta2.pooled,tau=tau, # true values
    tauN.hat2=tauN.hat2,tauB.hat2=tauB.hat2, # estimation of tau
    biasBnor2=biasBnor2,biasBbin2=biasBbin2,biasF2=biasF2, # bias of beta
    mseBnor2=mseBnor2,mseBbin2=mseBbin2,mseF2=mseF2, # mean squared error for beta
    alphaBin2=alphaBin2,alphaNor2=alphaNor2,alphaF2=alphaF2, # type 1 error (alpha)
    betaBin2=betaBin2,betaNor2=betaNor2,betaF2=betaF2, # type 2 error (beta)
    MCseBin2=MCseBin2,MCseNor2=MCseNor2,MCseF2=MCseF2) # monte carlo standard error
  row.names(rval) <- NULL
}
  return(rval)
}
