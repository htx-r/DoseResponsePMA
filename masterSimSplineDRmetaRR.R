### RR (relative risk): In this file I run jags normal, jags binomial and dosresmeta for the simulated data for cubic spline transformation
### OR (odds ratio): In this file I run jags normal, jags binomial and dosresmeta (freq) for the simulated data for linear transformation

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

OneSimSplineRR <- function(beta1.pooled=0.03,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  sim.data <- simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=FALSE,splines=TRUE)
  # 1. Freq: dosresmeta
  rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                                   se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

  # 2.Bayes Normal: jags
  jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=T,new.dose.range = c(1,10))

  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelNorSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
  # 3.Bayes Binomial: jags


  splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
  # Results

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

  return(c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),tauN1=t1n,tauB1=t1b,sig.testBin1=sig.testBin1,sig.testNor1=sig.testNor1,sig.testF1 =sig.testF1,
           BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),tauN2=t2n,tauB2=t2b,sig.testBin2=sig.testBin2,sig.testNor2=sig.testNor2,sig.testF2 =sig.testF2))


}


simpowerSplineRR <- function(nrep=3,beta1.pooled=0.02,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneSimSplineRR(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = T)

  #### #### #### ####
  #### beta1.pooled
  #### #### #### ####
  # Biases
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
  ret.obj <- c( # beta1.pooled
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
  return(ret.obj)
}



nrep <- 10
beta1.pooled <- c(0,0.02,0.04,0.06,0.1,0.2,0.3)
beta2.pooled <- c(0,0,0.02,0.05,0.03,0.1,0.2 )
tau <- c(0.001,0.05)

## %% smaller tau
# Scenario 1
set.seed('197')
S1RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1])

# Scenario 2
set.seed('297')

S2RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1])

# Scenario 3
set.seed('397')

S3RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1])

# Scenario 4
set.seed('497')

S4RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1])

# Scenario 5
set.seed('597')
S5RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[1])

# Scenario 6
set.seed('697')

S6RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[1])

# Scenario 7
set.seed('797')

S7RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[7],beta2.pooled = beta2.pooled[7],tau=tau[1])

## %% Larger tau

# Scenario 8:
set.seed('897')

S8RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2])

# Scenario 9
set.seed('997')

S9RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2])


# Scenario 10
set.seed('1097')

S10RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2])


#Scenario 11
set.seed('1197')

S11RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2])


# Scenario 12
set.seed('1279')

S12RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[2]-0.04)


# Scenario 13
set.seed('1397')

S13RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[2]-0.04)

# Scenario 14
set.seed('1497')
S14RRspline <- simpowerSplineRR(nrep=nrep,beta1.pooled = beta1.pooled[7],beta2.pooled = beta2.pooled[7],tau=tau[2]-0.04)



# Save the results in a file

resRRspline <- rbind(S1RRspline,S2RRspline,S3RRspline,S4RRspline,S5RRspline,S6RRspline,S7RRspline,S8RRspline,S9RRspline,S10RRspline, S11RRspline, S12RRspline,S13RRspline,S14RRspline)
write.csv(resRRspline,file=paste0(Sys.Date(),"resRRspline.csv")) # keeps the rownames















# resRRspline1_df <- as.data.frame(resRRspline1)
# plot(resRRspline1_df$beta1.pooled,abs(resRRspline1_df$biasB1n),ylim = c(-0.001,0.03),pch=19,las=1,xlab='true.beta1',ylab='bias')#,col=as.numeric(as.factor(resRRspline1_df$tau)))
# points(resRRspline1_df$beta1.pooled,abs(resRRspline1_df$biasB1b),col=2,pch=19)
# points(resRRspline1_df$beta1.pooled,abs(resRRspline1_df$biasF1),col=3,pch=19)
# legend('topright',legend=c('Normal','Binomial','Freq.'),col=1:3,pch=19,bty='n')
# title(' RR spline, beta1')
#
# resRRspline2_df <- as.data.frame(resRRspline2)
# plot(resRRspline2_df$beta2.pooled,abs(resRRspline2_df$biasB2n),ylim = c(-0.001,0.03),pch=19,las=1,xlab='true.beta2',ylab='bias')#,col=as.numeric(as.factor(resRRspline1_df$tau)))
# points(resRRspline2_df$beta2.pooled,abs(resRRspline2_df$biasB2b),col=2,pch=19)
# points(resRRspline2_df$beta2.pooled,abs(resRRspline2_df$biasF2),col=3,pch=19)
# legend('topright',legend=c('Normal','Binomial','Freq.'),col=1:3,pch=19,bty='n')
# title(' RR spline, beta2')









# Arguments
# beta1.pooled <- rep(c(0,0.02,0.03,0.05,0.03),2)
# beta2.pooled <- rep(c(0,0,0.02,0.03,0.05),2)
# tau <- rep(c(0.001,0.01),each=6)
# ns <- 20
# doserange <- c(0,10)
# samplesize<- 200
# nrep <- 100
#
# # Results
# #MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
# res <- mapply(MultiRunSimulateDRSplineRR, beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
# res.mat1 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef1[1,]))
# res.mat2 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef2[1,]))
# rownames(res.mat1) <- paste0('Scenario ',1:nrow(res.mat1))
# rownames(res.mat2) <- paste0('Scenario ',1:nrow(res.mat2))
# res.mat1
# res.mat2
#return(list(sum.coef1=res.mat1,sum.coef2=res.mat2))
#}








n.sim.data <- 100
ns <- 20
start <- Sys.time()
beta1.pooled=0.03
beta2.pooled=0.01
bayesCoef1RRnor<-bayesCoef2RRnor<-bayesCoef1RRbin<-bayesCoef2RRbin<-freqCoef1RR<-freqCoef2RR<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Freq
  DRfreqmodelRCS <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  #rcsplineDRmetaFreq2 <- dosresmeta(formula = logRR~rcs(dose1,sim.data$knots), id = Study_No,type=type,
  #                                 se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  #coef2<-cbind(coef2,c(coef(rcsplineDRmetaFreq2)[1],coef(rcsplineDRmetaFreq2)[2]))

  # Bayes normal
  jagsdataRCSnor <- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(1,10))
  ## ----!!!!!---- PLEASE PUT THE COMMAND OF THE PRECISION IN THE makeJAGSDRmeta FUNCTION. And as we agreed you need to find out why you calculated it wrong-------------
  jagsdataRCSnor$prec <-  matrix(unlist(sapply(DRfreqmodelRCS$Slist,solve,simplify = F)),40,2,byrow = T)
  DRjagsmodelRCSnor <- jags.parallel(data = jagsdataRCSnor,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
  # traceplot(DRjagsmodelRCSnor$BUGSoutput)

  #Bayes binomial
  #jagsdataRCSbin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,noncases=noncases,data=sim.data$simulatedDRdata,Splines=T)

  #DRjagsmodelRCSbin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmeta,
  #                                         n.chains=1,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)

  bayesCoef1RRnor <- c(bayesCoef1RRnor,DRjagsmodelRCSnor$BUGSoutput$mean$beta1.pooled)
  bayesCoef2RRnor <- c(bayesCoef2RRnor,DRjagsmodelRCSnor$BUGSoutput$mean$beta2.pooled)

  bayesCoef1RRbin <- c(bayesCoef1RRbin,DRjagsmodelRCSbin$BUGSoutput$mean$beta1.pooled)
  bayesCoef2RRbin <- c(bayesCoef2RRbin,DRjagsmodelRCSbin$BUGSoutput$mean$beta2.pooled)

  freqCoef1RR<-c(freqCoef1RR,coef(DRfreqmodelRCS)[1])
  freqCoef2RR<-c(freqCoef2RR,coef(DRfreqmodelRCS)[2])
}

(mean(bayesCoef1RRnor)-beta1.pooled) ## 20/07 bias = -0.0002586919  22/8: 0.001274669
(mean(bayesCoef1RRnor)-beta2.pooled) ## 20/07 bias = 0.0004295099   22/8: -0.004455597


(mean(freqCoef1RR)-beta1.pooled) # 20/07 bias = -7.182011e-06  22/08: 0.001327889
(mean(freqCoef2RR)-beta2.pooled) # 20/07 bias = 2.164394e-05   22/08: -0.00396104


cbind(bayes1=quantile(bayesCoef1RR), freq1=quantile(freqCoef1RR))
