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
bayesCoef <- c()
n.sim.data <- 100
for (i in 1:n.sim.data) {
  sim <- simulateDRsplinedata.fun(beta1.pooled=0.02,beta2.pooled=0.03,tau=0.001,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots
  jagsdataSplineBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose = dose,cases=cases,controls=noncases,data=sim.data,Splines=T,knots=knots)

  splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmeta,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)
  bayesCoef <- c(bayesCoef,linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled)
}
mean(bayesCoef)-0.01
quantile(bayesCoef)
#traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varnames='beta.pooled')

### SIMULATIONS FOR DOSERESMETA
freqCoef<-c()
for(i in 1:100){
  sim <- simulateDRsplinedata.fun(beta1.pooled=0.02,beta2.pooled=0.03,tau=0.001,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~rcs(sim.data$dose,knots), id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',covariance = 'gl')

  freqCoef<-c(freqCoef,coef(linearDRmetaFreq))
}
mean(freqCoef)-0.01
quantile(freqCoef) #median and mean are equal to beta.pooled 0.01
#


