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
bayesCoefOR <- c()
freqCoefOR<-c()
n.sim.data <- 100
beta.pooled = 0.01
tau <- 0.001
for (i in 1:n.sim.data) {
  sim.data <- simulateDRlineardataOR.fun(beta.pooled ,tau)
  #sim.data$logRR[sim.data$dose!=0] <- log(sim.data$cases[sim.data$dose!=0]/sim.data$cases[sim.data$dose==0])
  jagsdataLinearBinOR <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose,dose2=NULL,cases=cases,noncases=noncases,data=sim.data,Splines=F)
  #jagsdataLinearBin$prec.beta <- 1/(0.001)^2
  linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBinOR,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinomialLinearDRmetaOR,
                                            n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)
  linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
  #traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varname='beta.pooled')
  bayesCoefOR <- c(bayesCoefOR,linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled)

  linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logOR ~ dose, type = type, id = Study_No,
                                             se = selogOR, cases = cases, n = cases+noncases  , data = sim.data,covariance = 'gl',proc = '2stage',method = 'fixed')#!!!!!!!!!!!!!!
  freqCoefOR<-c(freqCoefOR,coef(linearDRmetaFreq))
}
mean(bayesCoefOR)-beta.pooled # 19/07: bias=-0.0004581789
mean(freqCoefOR)-beta.pooled  # 19/07: bias=1.331602e-05

#  Bayes is unbiased but the values are largely varied compared to Freq.
cbind(bayes=quantile(bayesCoefOR),freq=quantile(freqCoefOR))


# 19/07
# bayes        freq
# 0%   -0.011492108 0.009062892
# 25%   0.001889695 0.009783084
# 50%   0.009374487 0.010019618
# 75%   0.014383270 0.010287524
# 100%  0.037490445 0.010706817





















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



