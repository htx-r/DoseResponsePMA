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
sim.data <- simulateDRlineardata.fun(beta.pooled = 0.01,tau = 0.001)
jagsdataLinearBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose = dose,cases=cases,controls=noncases,data=sim.data,Splines=F)

linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinomialLinearDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)
bayesCoef <- c(bayesCoef,linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled)
}
mean(bayesCoef)-0.01
quantile(bayesCoef)
#traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varnames='beta.pooled')

### SIMULATIONS FOR DOSERESMETA
freqCoef<-c()
for(i in 1:100){
sim.data <- simulateDRlineardata.fun(beta.pooled = 0.01,tau = 0.001)
linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ dose, type = type, id = Study_No,
                                           se = selogRR, cases = cases, n = cases+noncases  , data = sim.data,covariance = 'gl',proc = '2stage',method = 'fixed')#!!!!!!!!!!!!!!
freqCoef<-c(freqCoef,coef(linearDRmetaFreq))
}
mean(freqCoef)-0.01
quantile(freqCoef) #median and mean are equal to beta.pooled 0.01
#










# lm(logRR~dose-1,data=sim.data)
# linearDRmetaJAGSmodelBin$BUGSoutput$mean
#
#
#
# exp(0.01*sim.data$dose)
#
# linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
# linearDRmetaJAGSmodelBin$BUGSoutput$mean$tau
#
# traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varname='beta.pooled')
#
# lm(simulatedDRdata[simulatedDRdata$dose!=0,]$cases/simulatedDRdata[simulatedDRdata$dose==0,]$cases
# ~simulatedDRdata$dose[simulatedDRdata$dose!=0,]-1)
#



