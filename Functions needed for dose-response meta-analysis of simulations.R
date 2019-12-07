# This file contains two functions, I will use them to analyze the simulated data
# 1. OneSimulation() is used to obtain the results from only one simulation with different specification: OR or RR, linear or spline
# 2. simpower() is used to repeat the OneSimulation() function many times (nsim) and then use rsimsum package to combine
          # all 'nsim' results and give the final performance measure: bias, MSE, MCerror, coverage, power, ....

OneSimulation <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
        #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
        #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
        #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
    v<- class(sim.data)
    }

# if(class(sim.data)=='try-error'){
# rval <- rep(NA,22)
# }else{

  if(splines==FALSE){
    #** 2l. linear inferences based on the three approaches: freq, normal bayes and binomial bayes
    # Freq: dosresmeta
    linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                                 se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')

    # Bayes Normal: jags
    jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,splines=FALSE)

    linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
                                           n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)

    # Bayes Binomial:jags
    if(OR){
      linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                                n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
    }else{
      linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
                                                n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
    }

    #** 3l. linear results

    # beta
     # mean
    f <-coef(linearDRmetaFreq)[1]
    bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
    bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled

     # standard error
    sdF <- sqrt(linearDRmetaFreq$vcov)
    sdBin <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','sd']
    sdNor <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','sd']

     # measure to check the convergence: Rhat gelamn statistic
    RhatN <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','Rhat']
    RhatB <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','Rhat']

    # heterogenity tau
     # mean
    tf <- sqrt(linearDRmetaFreq$Psi)
    tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau

     # standard error
    sdBintau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
    sdNortau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']

     # measure to check the convergence: Rhat gelamn statistic
    RhatBtau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
    RhatNtau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']


# the return object is vector that combine all results: linear
    rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),sdF=sdF,sdNor=sdNor,sdBin=sdBin
              ,tauN=tn,tauB=tb,tauF=tf,RhatN=RhatN,RhatB=RhatB,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
  }else{#
    #** 2s. spline inferences based on the three approaches: freq, normal bayes and binomial bayes

    # Freq: dosresmeta
    rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                                     se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

    # Bayes Normal: jags
    jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,splines=T)

    rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                             n.chains=3,n.iter = 200000,n.burnin = 20000,DIC=F,n.thin = 5)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 5)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 5)
    }

    #** 3s. spline results

    # beta1
     # mean
    f1 <-coef(rcsplineDRmetaFreq)[1]
    b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
    b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled

     # standard error
    sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
    sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
    sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']

     # measure to check the convergence: Rhat gelamn statistic
    RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
    RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']

    # heterogenity tau
     # mean
    tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
    tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau

    # standard error
    sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
    sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
    sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']


    # beta2
     # mean
    f2 <-coef(rcsplineDRmetaFreq)[2]
    b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
    b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled

     # standard error
    sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
    sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
    sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']

     # measure to check the convergence: Rhat gelamn statistic
    RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
    RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']

    # commomn heterogenity tau in both beta1 and beta2
     # mean
    tf <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
    tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau

     # standard error
    sdBintau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
    sdNortau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']

     # measure to check the convergence: Rhat gelamn statistic
    RhatNtau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
    RhatBtau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']

     # the return object is vector that combine all results: spline
    rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf,RhatN1=RhatN1,RhatB1=RhatB1,
              BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,RhatN2=RhatN2,RhatB2=RhatB2
              ,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
  }

  return(rval)

}
# End of OneSimulation()

simpower <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

# linear model
if(splines==FALSE){
  # repeat the simulation nsim times
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
  mdf <- colMeans(t(res),na.rm = TRUE)
  dfbeta$RhatB<-mdf['RhatB']
  dfbeta$RhatN<-mdf['RhatN']

  # Calculate as a single row dataframe the true tau with its three estimated tau
  #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])

  # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
  dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
  mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
  smstau <- summary(mstau)

  # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
  rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
  mtau <- spread(rvaltau,par,est)
  rownames(mtau) <- mtau[,'stat']
  mtau <- mtau[,-1]

  # Convert the matrix above of PM to a single row of data.frame
  dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
  colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
  rownames(dftau) <-NULL

  # Add Monte carlo standard error of bias to the row before (as dataframe)
  dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
  dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']

  # Add the convergence measure values for Rhat
  mdf <- colMeans(t(res),na.rm = TRUE)
  dftau$RhatBtau <-mdf['RhatBtau']
  dftau$RhatNtau <-mdf['RhatNtau']

  # End: bind the true beta with the two row dataframe dfbeta and dftau
  result <- cbind.data.frame(true.beta=beta1.pooled,dfbeta,true.tau=tau,dftau)
  rownames(result) <- NULL

  # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'

  return(list(res1=result,res2=res))
  #
}else{
  # Repeat the simulation nsim times
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
  mdf <- colMeans(t(res),na.rm = TRUE)
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
  #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])

  # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
  dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
  mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
  smstau <- summary(mstau)

  # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
  rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
  mtau <- spread(rvaltau,par,est)
  rownames(mtau) <- mtau[,'stat']
  mtau <- mtau[,-1]

  # Convert the matrix above of PM to a single row of data.frame
  dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
  colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
  rownames(dftau) <-NULL

  # Add Monte carlo standard error of bias to the row before (as dataframe)
  dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
  dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']

  # Add the convergence measure values for Rhat
  mdf <- colMeans(t(res),na.rm = TRUE)
  dftau$RhatBtau <-mdf['RhatBtau']
  dftau$RhatNtau <-mdf['RhatNtau']

  # End: bind the true beta with the two row dataframe dfbeta and dftau
  result <- cbind.data.frame(true.beta1=beta1.pooled,dfbeta1,true.beta2=beta2.pooled,dfbeta2,true.tau=tau,dftau)
  rownames(result) <- NULL

  # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
  return(list(res1=result,res2=res))
}

  }
# end of simpower()

