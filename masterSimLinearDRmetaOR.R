### OR (odds ratio): In this file I run jags normal, jags binomial and dosresmeta (freq) for the simulated data for linear transformation

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

#covar.logrr( cases = cases, n = cases+noncases, y=logOR,v=selogOR^2,type = 'cc',data = sim.data[sim.data$Study_No==20,])
#Slist <- sapply(unique(sim.data$Study_No), function(i) covar.logrr( cases = cases, n = cases+noncases, y=logOR,v=selogOR^2,type = type,data = sim.data[sim.data$Study_No==i,]) ,simplify = F)

OneSimLinearOR <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  sim.data <- simulateLinearDRmetaOR.fun(beta.pooled,tau = tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # 1. Freq: dosresmeta
  linearDRmetaFreq <- dosresmeta(formula = logOR~dose, id = Study_No,type=type,
                                 se = selogOR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  # 2.Bayes Normal: jags
  jagsdatalinear<- makejagsNorDRmeta(Study_No,logOR,dose,dose2=NULL,cases,noncases,se=selogOR,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))
  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelNorLinearDRmeta,
                                         n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
  # 3.Bayes Binomial: jags
  jagsdataLinearBin <- makejagsBinDRmeta(studyid=Study_No,dose1 = dose,dose2=NULL,cases=cases,noncases=noncases,data=sim.data,Splines=F)
  linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

  # Results
  f <-coef(linearDRmetaFreq)[1]
  bNor <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
  bBin <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
  tn <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
  tb <- linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau

  sig.testNor <- ifelse(linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodel$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
  sig.testBin <- ifelse(linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','2.5%']*linearDRmetaJAGSmodelBin$BUGSoutput$summary['beta.pooled','97.5%']>0, 1,0)
  sig.testF <- ifelse(summary(linearDRmetaFreq)$coefficients[,'Pr(>|z|)']<0.05,1,0)
  return(c(BayesB=bBin,BayesN=bNor,Freq=unname(f),tauN=tn,tauB=tb,sig.testBin=sig.testBin,sig.testNor=sig.testNor,sig.testF =sig.testF))

}


simpowerLinearOR <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneSimLinearOR(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = T)
  # res.mat <- do.call(rbind,res)
  # res.m <- apply(res, 1, mean)

  # Biases
  biasBbin <- colMeans(t(res))['BayesB']-beta.pooled
  biasBnor <- colMeans(t(res))['BayesN']-beta.pooled
  biasF <- colMeans(t(res))['Freq']- beta.pooled

  # heterogenity
  tauN.hat <- colMeans(t(res))['tauN']
  tauB.hat <- colMeans(t(res))['tauB']

  # MSE: Mean square error
  mseBbin <- mean((t(res)[,'BayesB']-beta.pooled)^2)
  mseBnor <- mean((t(res)[,'BayesN']-beta.pooled)^2)
  mseF <- mean((t(res)[,'Freq']-beta.pooled)^2)

  # Type 1 error
  alphaNor <- ifelse(beta.pooled==0,mean(beta.pooled==0&res['sig.testNor',]==1),NA)
  alphaBin <- ifelse(beta.pooled==0,mean(beta.pooled==0&res['sig.testBin',]==1),NA)
  alphaF <- ifelse(beta.pooled==0,mean(beta.pooled==0&res['sig.testF',]==1),NA)

  # Type 2 error
  betaNor <- ifelse(beta.pooled!=0,mean(beta.pooled!=0 & res['sig.testNor',]==0),NA)
  betaBin <- ifelse(beta.pooled!=0,mean(beta.pooled!=0 & res['sig.testBin',]==0),NA)
  betaF <- ifelse(beta.pooled!=0,mean(beta.pooled!=0 & res['sig.testF',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin <- sqrt(sum((t(res)[,'BayesB']-colMeans(t(res))['BayesB'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor <- sqrt(sum((t(res)[,'BayesN']-colMeans(t(res))['BayesN'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF <- sqrt(sum((t(res)[,'Freq']-colMeans(t(res))['Freq'])^2)/(ncol(res)*(ncol(res)-1)))

  ret.obj <- c(beta.pooled=beta.pooled,tau=tau, # true values
                   tauN.hat=tauN.hat,tauB.hat=tauB.hat, # estimation of tau
                   biasBnor=biasBnor,biasBbin=biasBbin,biasF=biasF, # bias of beta
                   mseBnor=mseBnor,mseBbin=mseBbin,mseF=mseF, # mean squared error for beta
                   alphaBin=alphaBin,alphaNor=alphaNor,alphaF=alphaF, # type 1 error (alpha)
                   betaBin=betaBin,betaNor=betaNor,betaF=betaF, # type 2 error (beta)
                   MCseBin=MCseBin,MCseNor=MCseNor,MCseF=MCseF) # monte carlo standard error
  #,beta.pooled.hatBnor=colMeans(res.mat)[1],beta.pooled.hatBbin=colMeans(res.mat)[2],beta.pooled.hatF=colMeans(res.mat)[3],
  return(ret.obj)
}



#MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrep <- 3
beta.pooled <- c(0,0.02,0.04,0.06,0.1,0.2,0.3)
tau <- c(0.001,0.05)

## %% smaller tau
# Scenario 1
set.seed('145')
S1ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[1],tau=tau[1])


# Scenario 2
set.seed('245')
S2ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[2],tau=tau[1])

# Scenario 3
set.seed('345')

S3ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[3],tau=tau[1])


# Scenario 4
set.seed('445')

S4ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[4],tau=tau[1])


# Scenario 5
set.seed('545')
S5ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[5],tau=tau[1])


# Scenario 6:
set.seed('645')

S6ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[6],tau=tau[1])

# Scenario 7:
set.seed('745')

S7ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[7],tau=tau[1])

## %% larger tau
# Scenario 8:
set.seed('845')

S8ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[1],tau=tau[2])

# Scenario 9:
set.seed('945')

S9ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[2],tau=tau[2])

# Scenario 10:
set.seed('1045')

S10ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[3],tau=tau[2])

# Scenario 11:
set.seed('1145')

S11ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[4],tau=tau[2])

# Scenario 12:
set.seed('1245')

S12ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[5],tau=tau[2])

# Scenario 13:
set.seed('1345')

S13ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[6],tau=tau[2])
# Scenario 14:
set.seed('1445')
S14ORlinear <- simpowerLinearOR(nrep,beta.pooled=beta.pooled[7],tau=tau[2])

# Save the results in a file
resORlinear <- rbind(S1ORlinear,S2ORlinear,S3ORlinear,S4ORlinear,S5ORlinear,S6ORlinear,S7ORlinear,S8ORlinear,S9ORlinear,S10ORlinear,S11ORlinear,S12ORlinear,S13ORlinear,S14ORlinear)
write.csv(resORlinear,file=paste0(Sys.Date(),'resORlinear.csv')) # keeps the rownames
















# par(mfrow=c(1,2))
# resORlinear_df <- as.data.frame(resORlinear)
# plot(resORlinear_df$beta.pooled,abs(resORlinear_df$biasBnor),ylim = c(-0.001,0.03),pch=19,las=1,xlab='true.beta',ylab='bias')#,col=as.numeric(as.factor(resORlinear_df$tau)))
# points(resORlinear_df$beta.pooled,abs(resORlinear_df$biasBbin),col=2,pch=19)
# points(resORlinear_df$beta.pooled,abs(resORlinear_df$biasF),col=3,pch=19)
# legend('topright',legend=c('Normal','Binomial','Freq.'),col=1:3,pch=19,bty='n')
# title(' OR linear')



# beta.pooled <- rep(c(0.02)) ## 0,0.05
# tau <- rep(c(0.001)) # ,0.01
# ns <- 20
# doserange <- c(0,10)
# samplesize<- 200
# nrep <- 100
# res <- mapply(simpowerLinearOR, beta.pooled=beta.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
# result <- do.call(rbind,res)
# rownames(result) <- paste0('Scenario ',1:nrow(result))
# result

#}
## Run all scenarios
