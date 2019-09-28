### RR (relative risk): In this file I run jags normal, jags binomial and dosresmeta (freq) for the simulated data for linear transformation

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

OneSimLinearRR <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateLinearDRmetaRR.fun(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # 1. Freq: dosresmeta
  linearDRmetaFreq<-dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                               se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')

  # 2. Bayes Normal: jags
  jagsdatalinear<- makejagsNorDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,se=selogRR,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelNorLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
# jpeg('plotn.jpg')
# devAskNewPage(ask=TRUE)
# par(ask=TRUE)
# traceplot(linearDRmetaJAGSmodel$BUGSoutput,varname='tau')
# dev.off()
  # 3. Bayes Binomial:jags
  jagsdataLinearBin <- makejagsBinDRmeta(studyid=Study_No,dose1 = dose,dose2=NULL,cases=cases,noncases=noncases,data=sim.data,Splines=F)

  linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','beta','tau'),model.file = modelBinLinearDRmetaRR,
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


simpowerLinearRR <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneSimLinearRR(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = T)

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

  ret.obj <- cbind(beta.pooled=beta.pooled,tau=tau, # true values
                   tauN.hat=tauN.hat,tauB.hat=tauB.hat, # estimation of tau
                   biasBnor=biasBnor,biasBbin=biasBbin,biasF=biasF, # bias of beta
                   mseBnor=mseBnor,mseBbin=mseBbin,mseF=mseF, # mean squared error for beta
                   alphaBin=alphaBin,alphaNor=alphaNor,alphaF=alphaF, # type 1 error (alpha)
                   betaBin=betaBin,betaNor=betaNor,betaF=betaF, # type 2 error (beta)
                   MCseBin=MCseBin,MCseNor=MCseNor,MCseF=MCseF) # monte carlo standard error
  row.names(ret.obj) <- NULL
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
set.seed('123')
S1RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[1],tau=tau[1])


# Scenario 2
set.seed('223')
S2RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[2],tau=tau[1])

# Scenario 3
set.seed('323')

S3RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[3],tau=tau[1])


# Scenario 4
set.seed('423')

S4RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[4],tau=tau[1])


# Scenario 5
set.seed('523')
S5RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[5],tau=tau[1])


# Scenario 6:
set.seed('623')

S6RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[6],tau=tau[1])

# Scenario 7:
set.seed('723')

S7RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[7],tau=tau[1])

## %% larger tau
# Scenario 8:
set.seed('823')

S8RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[1],tau=tau[2])

# Scenario 9:
set.seed('923')

S9RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[2],tau=tau[2])

# Scenario 10:
set.seed('1023')

S10RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[3],tau=tau[2])

# Scenario 11:
set.seed('1123')

S11RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[4],tau=tau[2])

# Scenario 12:
set.seed('1223')

S12RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[5],tau=tau[2])

# Scenario 13:
set.seed('1323')

S13RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[6],tau=tau[2])
# Scenario 14:
set.seed('1423')
S14RRlinear <- simpowerLinearRR(nrep,beta.pooled=beta.pooled[7],tau=tau[2])

# Save the results in a file
#eval(parse(text=paste0('S',1:2,'RRlinear')))
resRRlinear <- rbind(S1RRlinear,S2RRlinear,S3RRlinear,S4RRlinear,S5RRlinear,S6RRlinear,S7RRlinear,S8RRlinear,S9RRlinear,S10RRlinear,S11RRlinear,S12RRlinear,S13RRlinear,S14RRlinear)
write.csv(resRRlinear,file=paste0(Sys.Date(),"RRlinear.csv")) # keeps the rownames














# resRRlinear_df <- as.data.frame(resRRlinear)
# plot(resRRlinear_df$beta.pooled,abs(resRRlinear_df$biasBnor),ylim = c(-0.001,0.03),pch=19,las=1,xlab='true.beta',ylab='bias')#,col=as.numeric(as.factor(resRRlinear_df$tau)))
# points(resRRlinear_df$beta.pooled,abs(resRRlinear_df$biasBbin),col=2,pch=19)
# points(resRRlinear_df$beta.pooled,abs(resRRlinear_df$biasF),col=3,pch=19)
# legend('topright',legend=c('Normal','Binomial','Freq.'),col=1:3,pch=19,bty='n')
# title(' RR linear')









# beta.pooled <- rep(c(0.02)) ## 0,0.05
# tau <- rep(c(0.001)) # ,0.01
# ns <- 20
# doserange <- c(0,10)
# samplesize<- 200
# nrep <- 100
# simulation <- function(nrep=3,beta.pooled=c(0.02,0.01,0.2),tau=c(0.001,0.01),ns=20,doserange=c(1, 10),samplesize=200){
# # set random generator
#   # seed <- .Random.seed
#   res <- mapply(simpowerLinearRR, beta.pooled=beta.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
#   rval <- do.call(rbind,res)
#    #rval <-cbind(rval,seed=seed)
#   rownames(rval) <- paste0('Scenario ',1:nrow(rval))
#   return(rval)
# }
# simulaton(beta.pooled,tau)
#}
## Run all scenarios
