

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

simpower <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

# linear model
if(splines==FALSE){
  # Porduce the simulated dataset nsim times
  res <- replicate(nsim,OneSimulation(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

  # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
  df <- data.frame(beta=c(t(res)[,c('BayesB','BayesN','Freq')]),se=c(t(res)[,c('sdBin','sdNor','sdF')]),par=rep(c('BayesB','BayesN','Freq'),each=nsim))
  ms <- multisimsum(data = df,par = "par", true = c(BayesB=beta1.pooled,BayesN=beta1.pooled,Freq=beta1.pooled),estvarname = "beta", se = "se")
  sms <- summary(ms)

  # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
  rval <- sms$summ[sms$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
  m <- spread(rval,par,est)
  rownames(m) <- m[,'stat']
  m <- m[,-1]

  # Convert the matrix above of PM to a single row of data.frame
  dfbeta <- cbind(m[1,],m[2,],m[3,],m[4,],m[5,])
  colnames(dfbeta) <- paste0(rep(c('BayesB','BayesN','Freq'),5),rep(rownames(m),each=3))
  rownames(dfbeta) <-NULL

  # Add Monte carlo standard error of bias to the row before (as dataframe)
  dfbeta$mcse.biasB <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesB','mcse']
  dfbeta$mcse.biasN <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesN','mcse']
  dfbeta$mcse.biasF <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='Freq','mcse']

  # Add the convergence measure values for Rhat
  mdf <- colMeans(t(res))
  dfbeta$RhatB<-mdf['RhatB']
  dfbeta$RhatN<-mdf['RhatN']

  # Calculate as a single row dataframe the true tau with its three estimated tau
  dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])

  # End: bind the true beta with the two row dataframe dfbeta and dftau
  result <- cbind.data.frame(true.beta=beta1.pooled,dfbeta,dftau)
  rownames(result) <- NULL

  return(result)
  #
}else{
  # Porduce the simulated dataset nsim times
  res <- replicate(n=nsim,OneSimulation(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

# 1. beta1
  # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
  df1 <- data.frame(beta1=c(t(res)[,c('BayesB1','BayesN1','Freq1')]),se1=c(t(res)[,c('sdBin1','sdNor1','sdF1')]),par1=rep(c('BayesB1','BayesN1','Freq1'),each=nsim))
  ms1 <- multisimsum(data = df1,par = "par1", true = c(BayesB1=beta1.pooled,BayesN1=beta1.pooled,Freq1=beta1.pooled),estvarname = "beta1", se = "se1")
  sms1 <- summary(ms1)

  # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
  rval1 <- sms1$summ[sms1$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
  m1 <- spread(rval1,par1,est)
  rownames(m1) <- m1[,'stat']
  m1 <- m1[,-1]

  # Convert the matrix above of PM to a single row of data.frame
  dfbeta1 <- cbind(m1[1,],m1[2,],m1[3,],m1[4,],m1[5,])
  colnames(dfbeta1) <- paste0(rep(c('BayesB1','BayesN1','Freq1'),5),rep(rownames(m1),each=3))
  rownames(dfbeta1) <-NULL

  # Add Monte carlo standard error of bias to the row before (as dataframe)
  dfbeta1$mcse.biasB1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesB1','mcse']
  dfbeta1$mcse.biasN1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesN1','mcse']
  dfbeta1$mcse.biasF1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='Freq1','mcse']

  # Add the convergence measure values for Rhat
  mdf <- colMeans(t(res))
  dfbeta1$RhatB1<-mdf['RhatB1']
  dfbeta1$RhatN1<-mdf['RhatN1']


  # End for beta1



  # 2. beta2
  # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
  df2 <- data.frame(beta2=c(t(res)[,c('BayesB2','BayesN2','Freq2')]),se2=c(t(res)[,c('sdBin2','sdNor2','sdF2')]),par2=rep(c('BayesB2','BayesN2','Freq2'),each=nsim))
  ms2 <- multisimsum(data = df2,par = "par2", true = c(BayesB2=beta2.pooled,BayesN2=beta2.pooled,Freq2=beta2.pooled),estvarname = "beta2", se = "se2")
  sms2 <- summary(ms2)

  # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
  rval2 <- sms2$summ[sms2$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
  m2 <- spread(rval2,par2,est)
  rownames(m2) <- m2[,'stat']
  m2 <- m2[,-1]

  # Convert the matrix above of PM to a single row of data.frame
  dfbeta2 <- cbind(m2[1,],m2[2,],m2[3,],m2[4,],m2[5,])
  colnames(dfbeta2) <- paste0(rep(c('BayesB2','BayesN2','Freq2'),5),rep(rownames(m2),each=3))
  rownames(dfbeta2) <-NULL

  # Add Monte carlo standard error of bias to the row before (as dataframe)
  dfbeta2$mcse.biasB2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesB2','mcse']
  dfbeta2$mcse.biasN2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesN2','mcse']
  dfbeta2$mcse.biasF2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='Freq2','mcse']

  # Add the convergence measure values for Rhat
  dfbeta2$RhatB2<-mdf['RhatB2']
  dfbeta2$RhatN2<-mdf['RhatN2']

  # Calculate as a single row dataframe the true tau with its three estimated tau
  dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF1'],tau.hatF=mdf['tauF2'])

  # End: bind the true beta with the two row dataframe dfbeta and dftau
  result <- cbind.data.frame(true.beta1=beta1.pooled,dfbeta1,true.beta2=beta2.pooled,dfbeta2,dftau)
  rownames(result) <- NULL

  return(result)
}


  }


simpower()





