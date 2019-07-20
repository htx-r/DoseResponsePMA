# Libraries
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

# Analysis:
# a. Bayes: JAGS model
# jagsdataLinearBin$new.dose <- c(5,10,15)
# jagsdataLinearBin$new.n <- length(jagsdataLinearBin$new.dose)
bayesCoef1RR <- bayesCoef2RR<- freqCoef1RR <- freqCoef2RR<- c()
beta1.pooled <- 0.03
beta2.pooled <- 0.05
n.sim.data <- 100
for (i in 1:n.sim.data) {
  # Simulated data
  sim <- simulateDRsplinedata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots

  # Bayes
  jagsdataSplineBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,noncases=noncases,data=sim.data,Splines=T)

  splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmeta,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

  bayesCoef1RR <- c(bayesCoef1RR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled)
  bayesCoef2RR <- c(bayesCoef2RR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled)
  # traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varnames='beta1.pooled')
  # traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varnames='beta2.pooled')

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~rcs(sim.data$dose1,knots), id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',covariance = 'gl')

  freqCoef1RR<-c(freqCoef1RR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2RR<-c(freqCoef2RR,coef(rcsplineDRmetaFreq)[2])

}

(mean(bayesCoef1RR)-beta1.pooled) ## 19/07 bias = 0.01143533
(mean(bayesCoef2RR)-beta2.pooled) ## bias = 0.03684286


(mean(freqCoef1RR)-beta1.pooled) # 3.515168e-05
(mean(freqCoef2RR)-beta2.pooled) # 3.136834e-05


cbind(bayes1=quantile(bayesCoef1RR), freq1=quantile(freqCoef1RR))
# bayes1      freq1
# 0%   0.02129746 0.02905141
# 25%  0.03896093 0.02984452
# 50%  0.04167763 0.03006621
# 75%  0.04390764 0.03025678
# 100% 0.05358901 0.03068128
cbind(bayes2=quantile(bayesCoef2RR), freq2=quantile(freqCoef2RR))
# bayes2      freq2
# 0%   0.06813511 0.04811490
# 25%  0.08198005 0.04948280
# 50%  0.08600072 0.04998800
# 75%  0.09159239 0.05068297
# 100% 0.12627021 0.05265776



# the Bayes is still biased and the estimates is largely varied
