
# load libraries
source('Functions needed for dosres MA simulations.R')
library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
library('rsimsum')
library(tidyr)
devAskNewPage(ask=F)

#
# load and exclude single arm studies
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]

#
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders

# function to compute the relative odds ratio as: odds of non-reference dose / odds of refernce dose using metabin
createORreference.fun=function(r,n)
{

  logOR=c(0)
  selogOR=c(NA)

  for(i in 2:c(length(n)))
  {
    calculate=metabin(r[i],n[i],r[1],n[1],sm="OR")
    logOR=c(logOR,calculate$TE)
    selogOR=c(selogOR,calculate$seTE)

  }
  return(cbind(logOR=logOR,selogOR=selogOR))
}

# apply the function above to all studies
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# restricted cubic spline transformation doses
knots = c(10,20,50)
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

## I modified the origional Onesimulation function by adding the part where we remove the
    # a one dos-effect from 48 of the studies out of 60 and keep the rest (12) with 3 arm doses

OneSimulation <- function(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    antidepsimulation <- try(simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),silent = TRUE)
    v<- class(antidepsimulation)
  }

  antidepsimulation2arm1 <-sapply(1:48, function(i) antidepsimulation[antidepsimulation$Study_No==i,][-sample(c(2,3),1),],simplify = FALSE)
  antidepsimulation2arm<- do.call(rbind,antidepsimulation2arm1)
  sim.data <- rbind(antidepsimulation2arm,antidepsimulation[145:180,])
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### simulate the antidepressant freq estimated coefficients of the OR spline

load('antidepORsplineFINAL')

beta1.pooled <- unname(coef(doseresORsplineFreq)[1])
beta2.pooled <- unname(coef(doseresORsplineFreq)[2])
tau <- doseresORsplineBin$BUGSoutput$mean$tau
ns <- 60
doserange <- c(1,80)
nsim <- 1000


start <- Sys.time()
antidepSim <- simpower(nsim=nsim,beta1.pooled = beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,OR=TRUE,doserange=doserange,ns=ns,splines = TRUE)
end=Sys.time()
(end-start)

save(antidepSim ,file = 'antidepSimFINAL')

## table the results to display it in the table
beta1est <- as.numeric(antidepSim$res1['true.beta1'])+as.numeric(antidepSim$res1[c('BayesB1bias','BayesN1bias','Freq1bias')])

beta1estse <-sqrt(as.numeric(antidepSim$res1[c( 'BayesB1mse','BayesN1mse','Freq1mse')])
              - as.numeric(antidepSim$res1[c('BayesB1bias','BayesN1bias','Freq1bias')]^2))


beta2est <- as.numeric(antidepSim$res1['true.beta2'])+as.numeric(antidepSim$res1[c('BayesB2bias','BayesN2bias','Freq2bias')])

beta2estse <-sqrt(as.numeric(antidepSim$res1[c( 'BayesB2mse','BayesN2mse','Freq2mse')])
                  - as.numeric(antidepSim$res1[c('BayesB2bias','BayesN2bias','Freq2bias')]^2))

rvalSimantidep <- round(c(beta1.pooled,beta1est,beta1estse,beta2.pooled,beta2est,beta2estse),4)
names(rvalSimantidep) <- c('true.beta1','BayesB1','BayesN1','Freq1','BayesB1se','BayesN1se','Freq1se'
                           ,'true.beta2','BayesB2','BayesN2','Freq2','BayesB2se','BayesN2se','Freq2se')

tauest <- as.numeric(antidepSim$res1['true.tau'])+as.numeric(antidepSim$res1[c('BayesBtaubias','BayesNtaubias')])

tauestse <-sqrt(as.numeric(antidepSim$res1[c( 'BayesBtaumse','BayesNtaumse')])
                  - as.numeric(antidepSim$res1[c('BayesBtaubias','BayesNtaubias')]^2))

rvalSimantideptau <- round(c(tau,tauest,tauestse),4)
names(rvalSimantideptau) <- c('true.tau','BayesBtau','BayesNtau',
                             'BayesBtause','BayesNtause')
#




































































#
# antidepsimulation <- simulateDRmeta.fun(beta1.pooled = beta1,beta2.pooled = beta2,tau=tau,ns=ns,doserange = c(1,80))
# antidepsimulation2arm1 <-sapply(1:48, function(i) antidepsimulation[antidepsimulation$Study_No==i,][-sample(c(2,3),1),],simplify = FALSE)
# antidepsimulation2arm<- do.call(rbind,m)
# antidepsimulation2and3arms <- rbind(mm,antidepsimulation[145:180,])
#
#
#
# ORsplineFreqAntidepSim <- dosresmeta(formula=logrr~dose1+dose2, proc="1stage",id=Study_No, type=type,cases=cases,
#                                      n=cases+noncases,se=selogrr,data=antidepsimulation2and3arms,method = 'reml')
#
# jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=antidepsimulation2and3arms,splines=T)
#
# ORsplineNorBayesAntidepSim <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
#                                             n.chains=3,n.iter = 2000,n.burnin = 200,DIC=F,n.thin = 1)
# # Bayes Binomial: jags
# ORsplineBinBayesAntidepSim <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
#                                             n.chains=3,n.iter = 1000,n.burnin =100,DIC=F,n.thin = 1)
# #
