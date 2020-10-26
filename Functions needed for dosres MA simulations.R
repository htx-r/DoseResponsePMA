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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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
}else{ #
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

# 1. Sample
simulateDRmeta.funSample=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

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
    logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
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
    uniquess<-round(runif(ns,20,100)) # sample size of each study arm
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationSample <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funSample(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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

simpowerSample <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # linear model
  if(splines==FALSE){
    # repeat the simulation nsim times
    res <- replicate(nsim,OneSimulationSample(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
  }else{ #
    # Repeat the simulation nsim times
    res <- replicate(n=nsim,OneSimulationSample(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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


# 2. Dose
simulateDRmeta.funDose=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

  # This function generate a dataset with odds ratio OR or risk ratio RR based on dose-response meta-analysis model.

  # Arguments:
  # beta1.pooled: (numeric) the underlying slope of the linear dose-response model
  # beta2.pooled: (numeric) the underlying coeffiecient of cubic spline transformation of dose-response model
  # tau: the commom heterogenity for both regression coeffiecients across studies.
  # ns: number of studies (assumed to be even)
  # doserange: the range of the generated dosages.
  # sample size: the mean value of the uinform distribution that we generate the sample size on each dose i.e. Unif(samplesize-20,samplesize+20)
  # p0: is the probability of the event in the zero dose only for OR.
  # OR: (logical) indicate whether for OR (TRUE) or RR (FALSE) the data need to be generated.
  # splines: (logical) indicate whether linear or restricted cubic spline dose-response models.

  # load libraries
  require(Hmisc)

  # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
  d1<-cbind(rep(0,ns/2),matrix(round(c(runif(ns,1,6)),2),nrow=ns/2))##
  d2<-cbind(rep(0,ns/2),matrix(round(c(runif(ns,4,10)),2),nrow=ns/2))##

  d <- rbind(d1,d2)
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
    logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationDose <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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

simpowerDose <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # linear model
  if(splines==FALSE){
    # repeat the simulation nsim times
    res <- replicate(nsim,OneSimulationDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
  }else{ #
    # Repeat the simulation nsim times
    res <- replicate(n=nsim,OneSimulationDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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

# 3. half-sigmoid & uniform
simulateDRmeta.funSigm=function(ns=20,doserange=c(1, 10),samplesize=200,OR=TRUE,splines=TRUE,p0=0.1){ #

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
    # knots<-unlist(round(quantile(d[,2:3],c(0.1,0.5,0.9))))
    knots <- c(0,1,3)
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
    #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
    # implement the dose-response model
    # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
    #logrr <- ifelse(dose==0,0,log(log(dose)+1))
    logrr <- dose/(sqrt(1+dose^2))
  }else{ ## for linear
    # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
    # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
    # dose1 <- dose2 <- dose
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationSigm <- function(ns=20,doserange=c(1, 10),samplesize=200,OR=TRUE,splines = TRUE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funSigm(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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

simpowerSigm <- function(nsim=3,beta1.pooled=1,beta2.pooled=1,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # linear model
  if(splines==FALSE){
    # repeat the simulation nsim times
    res <- replicate(nsim,OneSimulationSigm(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
  }else{ #
    # Repeat the simulation nsim times
    res <- replicate(n=nsim,OneSimulationSigm(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
    #res <- replicate(n=nsim,OneSimulationLog(splines = T))
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

# 3. half-sigmoid & chi2
simulateDRmeta.funSigmChi=function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

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
  d<-cbind(rep(0,ns),matrix(round(c(rchisq(2*ns,2)),2),nrow=ns))##
  d<-t(apply(d,1,sort))
  dose <- c(t(d))


  # 2. Create the study-specific logOR and logRR either for splines or linear
  if(splines==TRUE){ ## for splines
    # find the dose cubic transformation
    #knots<-unlist(round(quantile(d[,2:3],c(0.1,0.5,0.9))))
    knots <- c(0,1,3)
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
    #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
    # implement the dose-response model
    # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
    #logrr <- ifelse(dose==0,0,log(log(dose)+1))
    logrr <- dose/(sqrt(1+dose^2))
  }else{ ## for linear
    # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
    # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
    # dose1 <- dose2 <- dose
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationSigmChi <- function(ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funSigmChi(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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

simpowerSigmChi <- function(nsim=3,beta1.pooled=1,beta2.pooled=1,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # linear model
  if(splines==FALSE){
    # repeat the simulation nsim times
    res <- replicate(nsim,OneSimulationSigmChi(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
  }else{ #
    # Repeat the simulation nsim times
    res <- replicate(n=nsim,OneSimulationSigmChi(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
    #res <- replicate(n=nsim,OneSimulationLog(splines = T))
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

# Log & uniform
simulateDRmeta.funLog=function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

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
    #knots<-unlist(round(quantile(dose,c(0.1,0.5,0.9))))
    knots <- c(0,1,3)
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
    #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
    # implement the dose-response model
    # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
    # logrr <- log(dose+1)
    # logrr <- ifelse(dose==0,0,log(log(dose)+1))
    # logrr <- log(log(dose+1))
    logrr <- log(log(dose+1)+1)
  }else{ ## for linear
    # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
    # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
    # dose1 <- dose2 <- dose
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationLog <- function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines = TRUE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funLog(ns=ns,doserange = doserange,samplesize = samplesize,p0=p0,OR=OR,splines=splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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


# Log & chi2
simulateDRmeta.funLogChi=function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

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
  d<-cbind(rep(0,ns),matrix(round(c(rchisq(2*ns,2)),2),nrow=ns))##
  d<-t(apply(d,1,sort))
  dose <- c(t(d))


  # 2. Create the study-specific logOR and logRR either for splines or linear
  if(splines==TRUE){ ## for splines
    # find the dose cubic transformation
    #knots<-unlist(round(quantile(dose,c(0.1,0.5,0.9))))
    knots <- c(0,1,3)
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
    #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
    # implement the dose-response model
    # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
    #logrr <- ifelse(dose==0,0,log(log(dose)+1))
    logrr <- log(log(dose+1)+1)
  }else{ ## for linear
    # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
    # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
    # dose1 <- dose2 <- dose
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationLogChi <- function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines = TRUE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funLogChi(ns=ns,doserange = doserange,samplesize = samplesize,p0=p0,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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



# 4. DiscDose
simulateDRmeta.funDiscDose=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

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
  d<-cbind(rep(0,ns),matrix(round(c(sample(doserange[1]:doserange[2],2*ns,replace = T)),2),nrow=ns))##
  # d1 <- runif(ns,1,7)
  # d2 <- runif(ns,2,10)
  # d <- cbind(rep(0,ns),d1,d2)
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
    logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
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
    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
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

OneSimulationDiscDose <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funDiscDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
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
                                             n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
    if(OR==TRUE){ ## model for OR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                                n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
    }else{ ## model for RR
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                                n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
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

simpowerDiscDose <- function(nsim=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # linear model
  if(splines==FALSE){
    # repeat the simulation nsim times
    res <- replicate(nsim,OneSimulationDiscDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
  }else{ #
    # Repeat the simulation nsim times
    res <- replicate(n=nsim,OneSimulationDiscDose(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)

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
























# # shape 1 & chi-square
# simulateDRmeta.funShape1Chi=function(ns=20,samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #
#
#   # This function generate a dataset with odds ratio OR or risk ratio RR based on dose-response meta-analysis model.
#
#   # Arguments:
#   # beta1.pooled: (numeric) the underlying slope of the linear dose-response model
#   # beta2.pooled: (numeric) the underlying coeffiecient of cubic spline transformation of dose-response model
#   # tau: the commom heterogenity for both regression coeffiecients across studies.
#   # ns: number of studies
#   # doserange: the range of the generated dosages.
#   # sample size: the mean value of the uinform distribution that we generate the sample size on each dose i.e. Unif(samplesize-20,samplesize+20)
#   # p0: is the probability of the event in the zero dose only for OR.
#   # OR: (logical) indicate whether for OR (TRUE) or RR (FALSE) the data need to be generated.
#   # splines: (logical) indicate whether linear or restricted cubic spline dose-response models.
#
#   # load libraries
#   require(Hmisc)
#
#   # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
#   d<-cbind(rep(0,ns),matrix(round(c(rchisq(2*ns,2)),2),nrow=ns))##
#   d<-t(apply(d,1,sort))
#   dose <- c(t(d))
#
#
#   # 2. Create the study-specific logOR and logRR either for splines or linear
#   if(splines==TRUE){ ## for splines
#     # find the dose cubic transformation
#     #knots<-unlist(round(quantile(d[,2:3],c(0.1,0.5,0.9))))
#     knots <- c(0,0.5,3)
#     trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     dose1 <- c(trans.d[,1])
#     dose2 <- c(trans.d[,2])
#     #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
#     # implement the dose-response model
#     # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
#     # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
#     # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
#     #logrr <- ifelse(dose==0,0,log(log(dose)+1))
#     logrr <- sqrt(dose)/(1+sqrt(dose))
#   }else{ ## for linear
#     # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
#     # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
#     # dose1 <- dose2 <- dose
#   }
#
#
#   # 3. Generate the dose-specific logOR and logRR
#   if(OR==TRUE){ ## odds ratio (OR)
#     # find the probabilities of the event to occur in zero dose;p0 and in the nonzero dose; p1
#     odds0 <- p0/(1-p0)  # as p0 is provided as an argument we compute the odds in zero dose
#     odds1 <- exp(logrr)*odds0 # from above we get the study-specific OR or RR
#     p1<- odds1/(1+odds1) #compute p1 the probability of a case at any dose level
#
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size of each study arm
#     ss<-c(sapply(uniquess,rep,3)) # sample size for all study arms
#     cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) # generate the cases for all dose levels
#     noncases<-matrix(c(ss-cases),nrow=3) # calculate the noncases = n - cases
#     #%%%!!!!! To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#
#     # compute the estimated odds ratio (hatOR) and its standard error for each dose level
#     hatodds<-cases/noncases # hat odds
#     hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #  log hat OR
#     SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of log hat OR
#     SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR) # add NA for SE in the zero dose1
#
#     # the final results go to the the returned data
#     hatlogrr <- hatlogOR
#     SEhatlogrr <- SEhatlogOR
#     type=rep('cc',3*ns) # type of the data is cc= case-control for OR. This is needed in dosresmeta function to choose how to compute the SE
#   }else{ ## Risk ratio (RR)
#
#     # restrictions to ensure that we will not get very large or small value for maxRR and so p0 and p1
#     maxRR<-exp(maxlogRR)
#     if(maxRR>10){
#       p0 <- 0.05
#       p1 <-ifelse(exp(logrr)*p0>1,0.95,exp(logrr)*p0)
#     }else if(maxRR<0.5){
#       p0 <- 0.95
#       p1 <-ifelse(exp(logrr)*p0<1,0.05,exp(logrr)*p0)
#     }else{
#       p0 <- 0.5/maxRR
#       p1 <-exp(logrr)*p0
#     }
#     p1 <- ifelse(p1>1, 0.97,p1)
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
#     ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
#     cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
#     noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#     # compute the estimated risk ratio (hatRR) and its standard error for each dose level
#     hatRR <- cases[2:3,]/cases[1,]
#     hatlogRR <- log(rbind(rep(1,ns),hatRR))
#     SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
#     selogRR<-c(rbind(NA,SEhatlogRR))
#
#     # the final results go to the the returned data object
#     hatlogrr <- hatlogRR
#     SEhatlogrr <- selogRR
#     type=rep('ci',3*ns) # type of the data is ci= cumulative incidance for RR. This is needed in dosresmeta function to choose how to compute the SE
#   }
#
#
#   Study_No<-rep(1:ns,each=3) # label the studies in numbers
#
#   # a data frame of the simulated data that we need to get from this function
#   simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
#                                     selogrr =c(SEhatlogrr), type=type)
#
#   return(simulatedDRdata=simulatedDRdata)
#
#
# }
# # End
#
# OneSimulationShape1Chi <- function(ns=20,samplesize=200,OR=FALSE,splines = FALSE){
#
#   #** 1. simulate the data;
#   #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
#   #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
#   #!!! for specific doses espacially for narrow dose range.
#   v <- 'try-error'
#   while (v=='try-error') {
#     sim.data <- try(simulateDRmeta.funShape1Chi(ns=ns,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
#     v<- class(sim.data)
#   }
#
#   # if(class(sim.data)=='try-error'){
#   # rval <- rep(NA,22)
#   # }else{
#
#   if(splines==FALSE){
#     #** 2l. linear inferences based on the three approaches: freq, normal bayes and binomial bayes
#     # Freq: dosresmeta
#     linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
#                                  se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,splines=FALSE)
#
#     linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
#                                            n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#
#     # Bayes Binomial:jags
#     if(OR){
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }else{
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }
#
#     #** 3l. linear results
#
#     # beta
#     # mean
#     f <-coef(linearDRmetaFreq)[1]
#     bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
#     bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
#
#     # standard error
#     sdF <- sqrt(linearDRmetaFreq$vcov)
#     sdBin <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','sd']
#     sdNor <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','Rhat']
#     RhatB <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf <- sqrt(linearDRmetaFreq$Psi)
#     tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatBtau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#     RhatNtau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#
#
#     # the return object is vector that combine all results: linear
#     rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),sdF=sdF,sdNor=sdNor,sdBin=sdBin
#               ,tauN=tn,tauB=tb,tauF=tf,RhatN=RhatN,RhatB=RhatB,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }else{#
#     #** 2s. spline inferences based on the three approaches: freq, normal bayes and binomial bayes
#
#     # Freq: dosresmeta
#     rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
#                                      se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,splines=T)
#
#     rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
#                                              n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     # Bayes Binomial: jags
#     if(OR==TRUE){ ## model for OR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
#                                                 n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
#     }else{ ## model for RR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
#                                                 n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     }
#
#     #** 3s. spline results
#
#     # beta1
#     # mean
#     f1 <-coef(rcsplineDRmetaFreq)[1]
#     b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
#     b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
#     RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#
#     # beta2
#     # mean
#     f2 <-coef(rcsplineDRmetaFreq)[2]
#     b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
#     b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled
#
#     # standard error
#     sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
#     sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
#     sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
#     RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']
#
#     # commomn heterogenity tau in both beta1 and beta2
#     # mean
#     tf <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatNtau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#     RhatBtau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#
#     # the return object is vector that combine all results: spline
#     rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf,RhatN1=RhatN1,RhatB1=RhatB1,
#               BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,RhatN2=RhatN2,RhatB2=RhatB2
#               ,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }
#
#   return(rval)
#
# }
# # End of OneSimulation()
#
# simpowerShape1Chi <- function(nsim=3,beta1.pooled=1,beta2.pooled=1,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
#
#   # linear model
#   if(splines==FALSE){
#     # repeat the simulation nsim times
#     res <- replicate(nsim,OneSimulationShape1Chi(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df <- data.frame(beta=c(t(res)[,c('BayesB','BayesN','Freq')]),se=c(t(res)[,c('sdBin','sdNor','sdF')]),par=rep(c('BayesB','BayesN','Freq'),each=nsim))
#     ms <- multisimsum(data = df,par = "par", true = c(BayesB=beta1.pooled,BayesN=beta1.pooled,Freq=beta1.pooled),estvarname = "beta", se = "se")
#     sms <- summary(ms)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval <- sms$summ[sms$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m <- spread(rval,par,est)
#     rownames(m) <- m[,'stat']
#     m <- m[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta <- cbind(m[1,],m[2,],m[3,],m[4,],m[5,])
#     colnames(dfbeta) <- paste0(rep(c('BayesB','BayesN','Freq'),5),rep(rownames(m),each=3))
#     rownames(dfbeta) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta$mcse.biasB <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesB','mcse']
#     dfbeta$mcse.biasN <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesN','mcse']
#     dfbeta$mcse.biasF <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='Freq','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta$RhatB<-mdf['RhatB']
#     dfbeta$RhatN<-mdf['RhatN']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta=beta1.pooled,dfbeta,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#
#     return(list(res1=result,res2=res))
#     #
#   }else{ #
#     # Repeat the simulation nsim times
#     res <- replicate(n=nsim,OneSimulationShape1Chi(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#     #res <- replicate(n=nsim,OneSimulationLog(splines = T))
#     # 1. beta1
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df1 <- data.frame(beta1=c(t(res)[,c('BayesB1','BayesN1','Freq1')]),se1=c(t(res)[,c('sdBin1','sdNor1','sdF1')]),par1=rep(c('BayesB1','BayesN1','Freq1'),each=nsim))
#     ms1 <- multisimsum(data = df1,par = "par1", true = c(BayesB1=beta1.pooled,BayesN1=beta1.pooled,Freq1=beta1.pooled),estvarname = "beta1", se = "se1")
#     sms1 <- summary(ms1)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval1 <- sms1$summ[sms1$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m1 <- spread(rval1,par1,est)
#     rownames(m1) <- m1[,'stat']
#     m1 <- m1[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta1 <- cbind(m1[1,],m1[2,],m1[3,],m1[4,],m1[5,])
#     colnames(dfbeta1) <- paste0(rep(c('BayesB1','BayesN1','Freq1'),5),rep(rownames(m1),each=3))
#     rownames(dfbeta1) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta1$mcse.biasB1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesB1','mcse']
#     dfbeta1$mcse.biasN1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesN1','mcse']
#     dfbeta1$mcse.biasF1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='Freq1','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta1$RhatB1<-mdf['RhatB1']
#     dfbeta1$RhatN1<-mdf['RhatN1']
#
#
#     # End for beta1
#
#
#
#     # 2. beta2
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df2 <- data.frame(beta2=c(t(res)[,c('BayesB2','BayesN2','Freq2')]),se2=c(t(res)[,c('sdBin2','sdNor2','sdF2')]),par2=rep(c('BayesB2','BayesN2','Freq2'),each=nsim))
#     ms2 <- multisimsum(data = df2,par = "par2", true = c(BayesB2=beta2.pooled,BayesN2=beta2.pooled,Freq2=beta2.pooled),estvarname = "beta2", se = "se2")
#     sms2 <- summary(ms2)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval2 <- sms2$summ[sms2$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m2 <- spread(rval2,par2,est)
#     rownames(m2) <- m2[,'stat']
#     m2 <- m2[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta2 <- cbind(m2[1,],m2[2,],m2[3,],m2[4,],m2[5,])
#     colnames(dfbeta2) <- paste0(rep(c('BayesB2','BayesN2','Freq2'),5),rep(rownames(m2),each=3))
#     rownames(dfbeta2) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta2$mcse.biasB2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesB2','mcse']
#     dfbeta2$mcse.biasN2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesN2','mcse']
#     dfbeta2$mcse.biasF2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='Freq2','mcse']
#
#     # Add the convergence measure values for Rhat
#     dfbeta2$RhatB2<-mdf['RhatB2']
#     dfbeta2$RhatN2<-mdf['RhatN2']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta1=beta1.pooled,dfbeta1,true.beta2=beta2.pooled,dfbeta2,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#     return(list(res1=result,res2=res))
#   }
#
# }
#
#
# # shape 2
# simulateDRmeta.funShape2=function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #
#
#   # This function generate a dataset with odds ratio OR or risk ratio RR based on dose-response meta-analysis model.
#
#   # Arguments:
#   # beta1.pooled: (numeric) the underlying slope of the linear dose-response model
#   # beta2.pooled: (numeric) the underlying coeffiecient of cubic spline transformation of dose-response model
#   # tau: the commom heterogenity for both regression coeffiecients across studies.
#   # ns: number of studies
#   # doserange: the range of the generated dosages.
#   # sample size: the mean value of the uinform distribution that we generate the sample size on each dose i.e. Unif(samplesize-20,samplesize+20)
#   # p0: is the probability of the event in the zero dose only for OR.
#   # OR: (logical) indicate whether for OR (TRUE) or RR (FALSE) the data need to be generated.
#   # splines: (logical) indicate whether linear or restricted cubic spline dose-response models.
#
#   # load libraries
#   require(Hmisc)
#
#   # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
#   d<-cbind(rep(0,ns),matrix(round(c(runif(2*ns,doserange[1],doserange[2])),2),nrow=ns))##
#   d<-t(apply(d,1,sort))
#   dose <- c(t(d))
#
#
#   # 2. Create the study-specific logOR and logRR either for splines or linear
#   if(splines==TRUE){ ## for splines
#     # find the dose cubic transformation
#     #knots<-unlist(round(quantile(d[,2:3],c(0.1,0.5,0.9))))
#     knots <- c(0,1,3)
#     trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     dose1 <- c(trans.d[,1])
#     dose2 <- c(trans.d[,2])
#     #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
#     # implement the dose-response model
#     # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
#     # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
#     # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
#     #logrr <- ifelse(dose==0,0,log(log(dose)+1))
#     logrr <- (dose^2)/(0.5+(0.5*dose^3))
#   }else{ ## for linear
#     # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
#     # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
#     # dose1 <- dose2 <- dose
#   }
#
#
#   # 3. Generate the dose-specific logOR and logRR
#   if(OR==TRUE){ ## odds ratio (OR)
#     # find the probabilities of the event to occur in zero dose;p0 and in the nonzero dose; p1
#     odds0 <- p0/(1-p0)  # as p0 is provided as an argument we compute the odds in zero dose
#     odds1 <- exp(logrr)*odds0 # from above we get the study-specific OR or RR
#     p1<- odds1/(1+odds1) #compute p1 the probability of a case at any dose level
#
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size of each study arm
#     ss<-c(sapply(uniquess,rep,3)) # sample size for all study arms
#     cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) # generate the cases for all dose levels
#     noncases<-matrix(c(ss-cases),nrow=3) # calculate the noncases = n - cases
#     #%%%!!!!! To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#
#     # compute the estimated odds ratio (hatOR) and its standard error for each dose level
#     hatodds<-cases/noncases # hat odds
#     hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #  log hat OR
#     SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of log hat OR
#     SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR) # add NA for SE in the zero dose1
#
#     # the final results go to the the returned data
#     hatlogrr <- hatlogOR
#     SEhatlogrr <- SEhatlogOR
#     type=rep('cc',3*ns) # type of the data is cc= case-control for OR. This is needed in dosresmeta function to choose how to compute the SE
#   }else{ ## Risk ratio (RR)
#
#     # restrictions to ensure that we will not get very large or small value for maxRR and so p0 and p1
#     maxRR<-exp(maxlogRR)
#     if(maxRR>10){
#       p0 <- 0.05
#       p1 <-ifelse(exp(logrr)*p0>1,0.95,exp(logrr)*p0)
#     }else if(maxRR<0.5){
#       p0 <- 0.95
#       p1 <-ifelse(exp(logrr)*p0<1,0.05,exp(logrr)*p0)
#     }else{
#       p0 <- 0.5/maxRR
#       p1 <-exp(logrr)*p0
#     }
#     p1 <- ifelse(p1>1, 0.97,p1)
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
#     ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
#     cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
#     noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#     # compute the estimated risk ratio (hatRR) and its standard error for each dose level
#     hatRR <- cases[2:3,]/cases[1,]
#     hatlogRR <- log(rbind(rep(1,ns),hatRR))
#     SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
#     selogRR<-c(rbind(NA,SEhatlogRR))
#
#     # the final results go to the the returned data object
#     hatlogrr <- hatlogRR
#     SEhatlogrr <- selogRR
#     type=rep('ci',3*ns) # type of the data is ci= cumulative incidance for RR. This is needed in dosresmeta function to choose how to compute the SE
#   }
#
#
#   Study_No<-rep(1:ns,each=3) # label the studies in numbers
#
#   # a data frame of the simulated data that we need to get from this function
#   simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
#                                     selogrr =c(SEhatlogrr), type=type)
#
#   return(simulatedDRdata=simulatedDRdata)
#
#
# }
# # End
#
# OneSimulationShape2 <- function(ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
#
#   #** 1. simulate the data;
#   #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
#   #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
#   #!!! for specific doses espacially for narrow dose range.
#   v <- 'try-error'
#   while (v=='try-error') {
#     sim.data <- try(simulateDRmeta.funShape2(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
#     v<- class(sim.data)
#   }
#
#   # if(class(sim.data)=='try-error'){
#   # rval <- rep(NA,22)
#   # }else{
#
#   if(splines==FALSE){
#     #** 2l. linear inferences based on the three approaches: freq, normal bayes and binomial bayes
#     # Freq: dosresmeta
#     linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
#                                  se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,splines=FALSE)
#
#     linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
#                                            n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#
#     # Bayes Binomial:jags
#     if(OR){
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }else{
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }
#
#     #** 3l. linear results
#
#     # beta
#     # mean
#     f <-coef(linearDRmetaFreq)[1]
#     bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
#     bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
#
#     # standard error
#     sdF <- sqrt(linearDRmetaFreq$vcov)
#     sdBin <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','sd']
#     sdNor <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','Rhat']
#     RhatB <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf <- sqrt(linearDRmetaFreq$Psi)
#     tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatBtau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#     RhatNtau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#
#
#     # the return object is vector that combine all results: linear
#     rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),sdF=sdF,sdNor=sdNor,sdBin=sdBin
#               ,tauN=tn,tauB=tb,tauF=tf,RhatN=RhatN,RhatB=RhatB,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }else{#
#     #** 2s. spline inferences based on the three approaches: freq, normal bayes and binomial bayes
#
#     # Freq: dosresmeta
#     rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
#                                      se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,splines=T)
#
#     rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
#                                              n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     # Bayes Binomial: jags
#     if(OR==TRUE){ ## model for OR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
#                                                 n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
#     }else{ ## model for RR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
#                                                 n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     }
#
#     #** 3s. spline results
#
#     # beta1
#     # mean
#     f1 <-coef(rcsplineDRmetaFreq)[1]
#     b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
#     b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
#     RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#
#     # beta2
#     # mean
#     f2 <-coef(rcsplineDRmetaFreq)[2]
#     b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
#     b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled
#
#     # standard error
#     sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
#     sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
#     sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
#     RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']
#
#     # commomn heterogenity tau in both beta1 and beta2
#     # mean
#     tf <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatNtau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#     RhatBtau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#
#     # the return object is vector that combine all results: spline
#     rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf,RhatN1=RhatN1,RhatB1=RhatB1,
#               BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,RhatN2=RhatN2,RhatB2=RhatB2
#               ,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }
#
#   return(rval)
#
# }
# # End of OneSimulation()
#
# simpowerShape2 <- function(nsim=3,beta1.pooled=1,beta2.pooled=1,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
#
#   # linear model
#   if(splines==FALSE){
#     # repeat the simulation nsim times
#     res <- replicate(nsim,OneSimulationShape2(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df <- data.frame(beta=c(t(res)[,c('BayesB','BayesN','Freq')]),se=c(t(res)[,c('sdBin','sdNor','sdF')]),par=rep(c('BayesB','BayesN','Freq'),each=nsim))
#     ms <- multisimsum(data = df,par = "par", true = c(BayesB=beta1.pooled,BayesN=beta1.pooled,Freq=beta1.pooled),estvarname = "beta", se = "se")
#     sms <- summary(ms)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval <- sms$summ[sms$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m <- spread(rval,par,est)
#     rownames(m) <- m[,'stat']
#     m <- m[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta <- cbind(m[1,],m[2,],m[3,],m[4,],m[5,])
#     colnames(dfbeta) <- paste0(rep(c('BayesB','BayesN','Freq'),5),rep(rownames(m),each=3))
#     rownames(dfbeta) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta$mcse.biasB <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesB','mcse']
#     dfbeta$mcse.biasN <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesN','mcse']
#     dfbeta$mcse.biasF <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='Freq','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta$RhatB<-mdf['RhatB']
#     dfbeta$RhatN<-mdf['RhatN']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta=beta1.pooled,dfbeta,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#
#     return(list(res1=result,res2=res))
#     #
#   }else{ #
#     # Repeat the simulation nsim times
#     res <- replicate(n=nsim,OneSimulationShape2(ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#     #res <- replicate(n=nsim,OneSimulationLog(splines = T))
#     # 1. beta1
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df1 <- data.frame(beta1=c(t(res)[,c('BayesB1','BayesN1','Freq1')]),se1=c(t(res)[,c('sdBin1','sdNor1','sdF1')]),par1=rep(c('BayesB1','BayesN1','Freq1'),each=nsim))
#     ms1 <- multisimsum(data = df1,par = "par1", true = c(BayesB1=beta1.pooled,BayesN1=beta1.pooled,Freq1=beta1.pooled),estvarname = "beta1", se = "se1")
#     sms1 <- summary(ms1)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval1 <- sms1$summ[sms1$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m1 <- spread(rval1,par1,est)
#     rownames(m1) <- m1[,'stat']
#     m1 <- m1[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta1 <- cbind(m1[1,],m1[2,],m1[3,],m1[4,],m1[5,])
#     colnames(dfbeta1) <- paste0(rep(c('BayesB1','BayesN1','Freq1'),5),rep(rownames(m1),each=3))
#     rownames(dfbeta1) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta1$mcse.biasB1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesB1','mcse']
#     dfbeta1$mcse.biasN1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesN1','mcse']
#     dfbeta1$mcse.biasF1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='Freq1','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta1$RhatB1<-mdf['RhatB1']
#     dfbeta1$RhatN1<-mdf['RhatN1']
#
#
#     # End for beta1
#
#
#
#     # 2. beta2
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df2 <- data.frame(beta2=c(t(res)[,c('BayesB2','BayesN2','Freq2')]),se2=c(t(res)[,c('sdBin2','sdNor2','sdF2')]),par2=rep(c('BayesB2','BayesN2','Freq2'),each=nsim))
#     ms2 <- multisimsum(data = df2,par = "par2", true = c(BayesB2=beta2.pooled,BayesN2=beta2.pooled,Freq2=beta2.pooled),estvarname = "beta2", se = "se2")
#     sms2 <- summary(ms2)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval2 <- sms2$summ[sms2$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m2 <- spread(rval2,par2,est)
#     rownames(m2) <- m2[,'stat']
#     m2 <- m2[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta2 <- cbind(m2[1,],m2[2,],m2[3,],m2[4,],m2[5,])
#     colnames(dfbeta2) <- paste0(rep(c('BayesB2','BayesN2','Freq2'),5),rep(rownames(m2),each=3))
#     rownames(dfbeta2) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta2$mcse.biasB2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesB2','mcse']
#     dfbeta2$mcse.biasN2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesN2','mcse']
#     dfbeta2$mcse.biasF2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='Freq2','mcse']
#
#     # Add the convergence measure values for Rhat
#     dfbeta2$RhatB2<-mdf['RhatB2']
#     dfbeta2$RhatN2<-mdf['RhatN2']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta1=beta1.pooled,dfbeta1,true.beta2=beta2.pooled,dfbeta2,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#     return(list(res1=result,res2=res))
#   }
#
# }
#
# simulateDRmeta.funShape2Chi=function(ns=20,samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #
#
#   # This function generate a dataset with odds ratio OR or risk ratio RR based on dose-response meta-analysis model.
#
#   # Arguments:
#   # beta1.pooled: (numeric) the underlying slope of the linear dose-response model
#   # beta2.pooled: (numeric) the underlying coeffiecient of cubic spline transformation of dose-response model
#   # tau: the commom heterogenity for both regression coeffiecients across studies.
#   # ns: number of studies
#   # doserange: the range of the generated dosages.
#   # sample size: the mean value of the uinform distribution that we generate the sample size on each dose i.e. Unif(samplesize-20,samplesize+20)
#   # p0: is the probability of the event in the zero dose only for OR.
#   # OR: (logical) indicate whether for OR (TRUE) or RR (FALSE) the data need to be generated.
#   # splines: (logical) indicate whether linear or restricted cubic spline dose-response models.
#
#   # load libraries
#   require(Hmisc)
#
#   # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
#   d<-cbind(rep(0,ns),matrix(round(c(rchisq(2*ns,2)),2),nrow=ns))##
#   d<-t(apply(d,1,sort))
#   dose <- c(t(d))
#
#
#   # 2. Create the study-specific logOR and logRR either for splines or linear
#   if(splines==TRUE){ ## for splines
#     # find the dose cubic transformation
#     #knots<-unlist(round(quantile(d[,2:3],c(0.1,0.5,0.9))))
#     knots <- c(0,1,3)
#     trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     dose1 <- c(trans.d[,1])
#     dose2 <- c(trans.d[,2])
#     #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
#     # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
#     # implement the dose-response model
#     # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
#     # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
#     # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
#     #logrr <- ifelse(dose==0,0,log(log(dose)+1))
#     logrr <- (dose^2)/(0.5+(0.5*dose^3))
#   }else{ ## for linear
#     # maxlogRR<- (beta1.pooled+2*tau)*max(dose) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
#     # beta  <- c(sapply(rnorm(ns,beta1.pooled,tau),rep,3)) #random effects of the slopes: generate the study-specific linear coeff from a normal distribution
#     # logrr <- beta*dose  # derive study-specific logOR by implementing the spline dose-response model
#     # dose1 <- dose2 <- dose
#   }
#
#
#   # 3. Generate the dose-specific logOR and logRR
#   if(OR==TRUE){ ## odds ratio (OR)
#     # find the probabilities of the event to occur in zero dose;p0 and in the nonzero dose; p1
#     odds0 <- p0/(1-p0)  # as p0 is provided as an argument we compute the odds in zero dose
#     odds1 <- exp(logrr)*odds0 # from above we get the study-specific OR or RR
#     p1<- odds1/(1+odds1) #compute p1 the probability of a case at any dose level
#
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size of each study arm
#     ss<-c(sapply(uniquess,rep,3)) # sample size for all study arms
#     cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) # generate the cases for all dose levels
#     noncases<-matrix(c(ss-cases),nrow=3) # calculate the noncases = n - cases
#     #%%%!!!!! To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#
#     # compute the estimated odds ratio (hatOR) and its standard error for each dose level
#     hatodds<-cases/noncases # hat odds
#     hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #  log hat OR
#     SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of log hat OR
#     SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR) # add NA for SE in the zero dose1
#
#     # the final results go to the the returned data
#     hatlogrr <- hatlogOR
#     SEhatlogrr <- SEhatlogOR
#     type=rep('cc',3*ns) # type of the data is cc= case-control for OR. This is needed in dosresmeta function to choose how to compute the SE
#   }else{ ## Risk ratio (RR)
#
#     # restrictions to ensure that we will not get very large or small value for maxRR and so p0 and p1
#     maxRR<-exp(maxlogRR)
#     if(maxRR>10){
#       p0 <- 0.05
#       p1 <-ifelse(exp(logrr)*p0>1,0.95,exp(logrr)*p0)
#     }else if(maxRR<0.5){
#       p0 <- 0.95
#       p1 <-ifelse(exp(logrr)*p0<1,0.05,exp(logrr)*p0)
#     }else{
#       p0 <- 0.5/maxRR
#       p1 <-exp(logrr)*p0
#     }
#     p1 <- ifelse(p1>1, 0.97,p1)
#     # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
#     uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
#     ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
#     cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
#     noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
#     cases[cases==0] <- 1
#     noncases[noncases==0] <- 1
#     # compute the estimated risk ratio (hatRR) and its standard error for each dose level
#     hatRR <- cases[2:3,]/cases[1,]
#     hatlogRR <- log(rbind(rep(1,ns),hatRR))
#     SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
#     selogRR<-c(rbind(NA,SEhatlogRR))
#
#     # the final results go to the the returned data object
#     hatlogrr <- hatlogRR
#     SEhatlogrr <- selogRR
#     type=rep('ci',3*ns) # type of the data is ci= cumulative incidance for RR. This is needed in dosresmeta function to choose how to compute the SE
#   }
#
#
#   Study_No<-rep(1:ns,each=3) # label the studies in numbers
#
#   # a data frame of the simulated data that we need to get from this function
#   simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
#                                     selogrr =c(SEhatlogrr), type=type)
#
#   return(simulatedDRdata=simulatedDRdata)
#
#
# }
# # End
#
# OneSimulationShape2Chi <- function(ns=20,samplesize=200,OR=FALSE,splines = FALSE){
#
#   #** 1. simulate the data;
#   #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
#   #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
#   #!!! for specific doses espacially for narrow dose range.
#   v <- 'try-error'
#   while (v=='try-error') {
#     sim.data <- try(simulateDRmeta.funShape2Chi(ns=ns,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
#     v<- class(sim.data)
#   }
#
#   # if(class(sim.data)=='try-error'){
#   # rval <- rep(NA,22)
#   # }else{
#
#   if(splines==FALSE){
#     #** 2l. linear inferences based on the three approaches: freq, normal bayes and binomial bayes
#     # Freq: dosresmeta
#     linearDRmetaFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
#                                  se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,splines=FALSE)
#
#     linearDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
#                                            n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#
#     # Bayes Binomial:jags
#     if(OR){
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }else{
#       linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
#                                                 n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)
#     }
#
#     #** 3l. linear results
#
#     # beta
#     # mean
#     f <-coef(linearDRmetaFreq)[1]
#     bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
#     bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
#
#     # standard error
#     sdF <- sqrt(linearDRmetaFreq$vcov)
#     sdBin <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','sd']
#     sdNor <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN <- linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','Rhat']
#     RhatB <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf <- sqrt(linearDRmetaFreq$Psi)
#     tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatBtau <- linearDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#     RhatNtau <- linearDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#
#
#     # the return object is vector that combine all results: linear
#     rval <- c(BayesB=bBin,BayesN=bNor,Freq=unname(f),sdF=sdF,sdNor=sdNor,sdBin=sdBin
#               ,tauN=tn,tauB=tb,tauF=tf,RhatN=RhatN,RhatB=RhatB,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }else{#
#     #** 2s. spline inferences based on the three approaches: freq, normal bayes and binomial bayes
#
#     # Freq: dosresmeta
#     rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
#                                      se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')
#
#     # Bayes Normal: jags
#     jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,splines=T)
#
#     rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
#                                              n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     # Bayes Binomial: jags
#     if(OR==TRUE){ ## model for OR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
#                                                 n.chains=3,n.iter = 100000,n.burnin =10000,DIC=F,n.thin = 1)
#     }else{ ## model for RR
#       splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
#                                                 n.chains=3,n.iter = 100000,n.burnin = 10000,DIC=F,n.thin = 1)
#     }
#
#     #** 3s. spline results
#
#     # beta1
#     # mean
#     f1 <-coef(rcsplineDRmetaFreq)[1]
#     b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
#     b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
#     RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']
#
#     # heterogenity tau
#     # mean
#     tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
#     sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
#     sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']
#
#
#     # beta2
#     # mean
#     f2 <-coef(rcsplineDRmetaFreq)[2]
#     b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
#     b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled
#
#     # standard error
#     sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
#     sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
#     sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
#     RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']
#
#     # commomn heterogenity tau in both beta1 and beta2
#     # mean
#     tf <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
#     tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
#     tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
#     # standard error
#     sdBintau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
#     sdNortau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']
#
#     # measure to check the convergence: Rhat gelamn statistic
#     RhatNtau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
#     RhatBtau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']
#
#     # the return object is vector that combine all results: spline
#     rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf,RhatN1=RhatN1,RhatB1=RhatB1,
#               BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,RhatN2=RhatN2,RhatB2=RhatB2
#               ,sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)
#   }
#
#   return(rval)
#
# }
# # End of OneSimulation()
#
# simpowerShape2Chi <- function(nsim=3,beta1.pooled=1,beta2.pooled=1,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){
#
#   # linear model
#   if(splines==FALSE){
#     # repeat the simulation nsim times
#     res <- replicate(nsim,OneSimulationShape2Chi(ns=ns,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df <- data.frame(beta=c(t(res)[,c('BayesB','BayesN','Freq')]),se=c(t(res)[,c('sdBin','sdNor','sdF')]),par=rep(c('BayesB','BayesN','Freq'),each=nsim))
#     ms <- multisimsum(data = df,par = "par", true = c(BayesB=beta1.pooled,BayesN=beta1.pooled,Freq=beta1.pooled),estvarname = "beta", se = "se")
#     sms <- summary(ms)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval <- sms$summ[sms$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m <- spread(rval,par,est)
#     rownames(m) <- m[,'stat']
#     m <- m[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta <- cbind(m[1,],m[2,],m[3,],m[4,],m[5,])
#     colnames(dfbeta) <- paste0(rep(c('BayesB','BayesN','Freq'),5),rep(rownames(m),each=3))
#     rownames(dfbeta) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta$mcse.biasB <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesB','mcse']
#     dfbeta$mcse.biasN <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='BayesN','mcse']
#     dfbeta$mcse.biasF <- sms$summ[sms$summ$stat=='bias'&sms$summ$par=='Freq','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta$RhatB<-mdf['RhatB']
#     dfbeta$RhatN<-mdf['RhatN']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta=beta1.pooled,dfbeta,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#
#     return(list(res1=result,res2=res))
#     #
#   }else{ #
#     # Repeat the simulation nsim times
#     res <- replicate(n=nsim,OneSimulationShape2Chi(ns=ns,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
#     #res <- replicate(n=nsim,OneSimulationLog(splines = T))
#     # 1. beta1
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df1 <- data.frame(beta1=c(t(res)[,c('BayesB1','BayesN1','Freq1')]),se1=c(t(res)[,c('sdBin1','sdNor1','sdF1')]),par1=rep(c('BayesB1','BayesN1','Freq1'),each=nsim))
#     ms1 <- multisimsum(data = df1,par = "par1", true = c(BayesB1=beta1.pooled,BayesN1=beta1.pooled,Freq1=beta1.pooled),estvarname = "beta1", se = "se1")
#     sms1 <- summary(ms1)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval1 <- sms1$summ[sms1$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m1 <- spread(rval1,par1,est)
#     rownames(m1) <- m1[,'stat']
#     m1 <- m1[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta1 <- cbind(m1[1,],m1[2,],m1[3,],m1[4,],m1[5,])
#     colnames(dfbeta1) <- paste0(rep(c('BayesB1','BayesN1','Freq1'),5),rep(rownames(m1),each=3))
#     rownames(dfbeta1) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta1$mcse.biasB1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesB1','mcse']
#     dfbeta1$mcse.biasN1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='BayesN1','mcse']
#     dfbeta1$mcse.biasF1 <- sms1$summ[sms1$summ$stat=='bias'&sms1$summ$par=='Freq1','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dfbeta1$RhatB1<-mdf['RhatB1']
#     dfbeta1$RhatN1<-mdf['RhatN1']
#
#
#     # End for beta1
#
#
#
#     # 2. beta2
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     df2 <- data.frame(beta2=c(t(res)[,c('BayesB2','BayesN2','Freq2')]),se2=c(t(res)[,c('sdBin2','sdNor2','sdF2')]),par2=rep(c('BayesB2','BayesN2','Freq2'),each=nsim))
#     ms2 <- multisimsum(data = df2,par = "par2", true = c(BayesB2=beta2.pooled,BayesN2=beta2.pooled,Freq2=beta2.pooled),estvarname = "beta2", se = "se2")
#     sms2 <- summary(ms2)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rval2 <- sms2$summ[sms2$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     m2 <- spread(rval2,par2,est)
#     rownames(m2) <- m2[,'stat']
#     m2 <- m2[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dfbeta2 <- cbind(m2[1,],m2[2,],m2[3,],m2[4,],m2[5,])
#     colnames(dfbeta2) <- paste0(rep(c('BayesB2','BayesN2','Freq2'),5),rep(rownames(m2),each=3))
#     rownames(dfbeta2) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dfbeta2$mcse.biasB2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesB2','mcse']
#     dfbeta2$mcse.biasN2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='BayesN2','mcse']
#     dfbeta2$mcse.biasF2 <- sms2$summ[sms2$summ$stat=='bias'&sms2$summ$par=='Freq2','mcse']
#
#     # Add the convergence measure values for Rhat
#     dfbeta2$RhatB2<-mdf['RhatB2']
#     dfbeta2$RhatN2<-mdf['RhatN2']
#
#     # Calculate as a single row dataframe the true tau with its three estimated tau
#     #dftau <- data.frame(true.tau=tau,tau.hatB=mdf['tauB'],tau.hatN=mdf['tauN'],tau.hatF=mdf['tauF'])
#
#     # Calculate the performance measure (PM) using multisimsum() for beta and display them in one row (as dataframe)
#     dftau <- data.frame(tau=c(t(res)[,c('tauN','tauB')]),seTau=c(t(res)[,c('sdNortau','sdBintau')]),par=rep(c('BayesN','BayesB'),each=nsim))
#     mstau <- multisimsum(data = dftau,par = "par", true = c(BayesN=tau,BayesB=tau),estvarname = "tau", se = "seTau")
#     smstau <- summary(mstau)
#
#     # Extract only some of the PM: 'bias','se2mean','mse','cover','power'
#     rvaltau <- smstau$summ[smstau$summ$stat%in%c('bias','se2mean','mse','cover','power'),c(1,2,4)]
#     mtau <- spread(rvaltau,par,est)
#     rownames(mtau) <- mtau[,'stat']
#     mtau <- mtau[,-1]
#
#     # Convert the matrix above of PM to a single row of data.frame
#     dftau <- cbind(mtau[1,],mtau[2,],mtau[3,],mtau[4,],mtau[5,])
#     colnames(dftau) <- paste0(rep(c('BayesBtau','BayesNtau'),5),rep(rownames(mtau),each=2))
#     rownames(dftau) <-NULL
#
#     # Add Monte carlo standard error of bias to the row before (as dataframe)
#     dftau$mcse.biastauB <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesB','mcse']
#     dftau$mcse.biastauN <- smstau$summ[smstau$summ$stat=='bias'&smstau$summ$par=='BayesN','mcse']
#
#     # Add the convergence measure values for Rhat
#     mdf <- colMeans(t(res),na.rm = TRUE)
#     dftau$RhatBtau <-mdf['RhatBtau']
#     dftau$RhatNtau <-mdf['RhatNtau']
#
#     # End: bind the true beta with the two row dataframe dfbeta and dftau
#     result <- cbind.data.frame(true.beta1=beta1.pooled,dfbeta1,true.beta2=beta2.pooled,dfbeta2,true.tau=tau,dftau)
#     rownames(result) <- NULL
#
#     # the return object is list of the summarized results 'res1' and the results for all simulations 'res2'
#     return(list(res1=result,res2=res))
#   }
#
# }

