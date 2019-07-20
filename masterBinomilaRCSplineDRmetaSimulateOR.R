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
bayesCoef1OR <- bayesCoef2OR<- freqCoef1OR <- freqCoef2OR<- c()
beta1.pooled <- 0.03
beta2.pooled <- 0.05
n.sim.data <- 100
for (i in 1:n.sim.data) {
  # Simulated data
  sim <- simulateDRsplinedataOR.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots

  # Bayes
  jagsdataSplineBinOR <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,noncases=noncases,data=sim.data,Splines=T)

  splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBinOR,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmetaOR,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

  bayesCoef1OR <- c(bayesCoef1OR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled)
  bayesCoef2OR <- c(bayesCoef2OR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled)
   #traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varnames='beta1.pooled')
  # traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varnames='beta2.pooled')

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logOR~rcs(sim.data$dose1,knots), id = Study_No,type=type,
                                   se = selogOR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',covariance = 'gl')

  freqCoef1OR<-c(freqCoef1OR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2OR<-c(freqCoef2OR,coef(rcsplineDRmetaFreq)[2])

}

(mean(bayesCoef1OR)-beta1.pooled) ## bias = 0.0001831384
(mean(bayesCoef2OR)-beta2.pooled) ## bias = -0.0006625379


(mean(freqCoef1OR)-beta1.pooled) # 5.354848e-05
(mean(freqCoef2OR)-beta2.pooled) # -2.350605e-05


cbind(quantile(bayesCoef1OR), quantile(freqCoef1OR))
cbind(quantile(bayesCoef2OR), quantile(freqCoef2OR))

