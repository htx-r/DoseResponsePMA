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

(mean(bayesCoef1RR)-beta1.pooled) ## bias = 0.0001831384
(mean(bayesCoef2RR)-beta2.pooled) ## bias = -0.0006625379


(mean(freqCoef1RR)-beta1.pooled) # 5.354848e-05
(mean(freqCoef2RR)-beta2.pooled) # -2.350605e-05


cbind(quantile(bayesCoef1RR), quantile(freqCoef1RR))
cbind(quantile(bayesCoef2RR), quantile(freqCoef2RR))

# 0%   0.02755065 0.02865196
# 25%  0.02955543 0.02985764
# 50%  0.03027903 0.03004675
# 75%  0.03087558 0.03029513
# 100% 0.03281189 0.03090807


# 0%   0.03936567 0.04805174
# 25%  0.04819887 0.04934958
# 50%  0.04928875 0.04982736
# 75%  0.05045750 0.05068508
# 100% 0.05622783 0.05186370

# the Bayes is still biased and the estimates is largely varied
