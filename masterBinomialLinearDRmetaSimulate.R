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
freqCoef<-c()
n.sim.data <- 100
beta.pooled = 0.01
tau <- 0.001
for (i in 1:n.sim.data) {
sim.data <- simulateDRlineardata.fun(beta.pooled ,tau)
#sim.data$logRR[sim.data$dose!=0] <- log(sim.data$cases[sim.data$dose!=0]/sim.data$cases[sim.data$dose==0])
jagsdataLinearBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose,dose2=NULL,cases=cases,controls=noncases,data=sim.data,Splines=F)
#jagsdataLinearBin$prec.beta <- 1/(0.001)^2
linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','beta'),model.file = modelBinomialLinearDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)
linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
#traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varname='beta.pooled')
bayesCoef <- c(bayesCoef,linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled)

linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ dose, type = type, id = Study_No,
                                           se = selogRR, cases = cases, n = cases+noncases  , data = sim.data,covariance = 'gl',proc = '2stage',method = 'fixed')#!!!!!!!!!!!!!!
freqCoef<-c(freqCoef,coef(linearDRmetaFreq))
}
mean(bayesCoef)-0.01 ## notfixtau: bias=0.001623046, true=0.01, n.sim.data =100  ##fix tau: bias=-0.003465108, true=0.01, n.sim.data =100
quantile(bayesCoef)
mean(freqCoef)-0.01 # bias = 3.414024e-05, true =0.01
quantile(freqCoef) #median and mean are equal to beta.pooled 0.01
#  Bayes is more biased compared to Freq.

# Based on sim.data, It should be that RR = beta*X ??=?? pevent/p0= cases_nonreferent/cases_referent (since n0=n1)

cbind(RR1=sim.data$cases[sim.data$dose!=0]/sim.data$cases[sim.data$dose==0],RR2=exp(sim.data$logRR[sim.data$dose!=0]))

# But actually they are not exactly equal because we draw cases from random binomial distibution.
# Binomila Bayes model is builded up using cases and noncases rather than logRR directly
# While Freq make the estimates based on the logRR which is by definition has true value for beta
# Based on that I modify the simulated data












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



