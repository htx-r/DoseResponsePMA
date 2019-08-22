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
bayesCoefRR <- c()
freqCoefRR<-c()
n.sim.data <- 100
beta.pooled = 0.01
tau <- 0.001
for (i in 1:n.sim.data) {
sim.data <- simulateDRlineardata.fun(beta.pooled ,tau)
#sim.data$logRR[sim.data$dose!=0] <- log(sim.data$cases[sim.data$dose!=0]/sim.data$cases[sim.data$dose==0])
jagsdataLinearBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose,dose2=NULL,cases=cases,noncases=noncases,data=sim.data,Splines=F)
#jagsdataLinearBin$prec.beta <- 1/(0.001)^2
linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','beta','tau'),model.file = modelBinomialLinearDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled
#traceplot(linearDRmetaJAGSmodelBin$BUGSoutput,varname='beta.pooled')
bayesCoefRR <- c(bayesCoefRR,linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled)

linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ dose, type = type, id = Study_No,
                                           se = selogRR, cases = cases, n = cases+noncases  , data = sim.data,covariance = 'gl',proc = '2stage',method = 'fixed')#!!!!!!!!!!!!!!
freqCoefRR<-c(freqCoefRR,coef(linearDRmetaFreq))
}
mean(bayesCoefRR)-beta.pooled ## 19/07: bias=0.008004524, true=0.01,  21/08: 0.0001429092
mean(freqCoefRR)-beta.pooled # 19/07: bias = 2.697325e-05, true =0.01, 21/08: 0.0002989553
cbind(bayes=quantile(bayesCoefRR), freq=quantile(freqCoefRR))

#19/07     bayes        freq
# 0%   -0.008135956 0.009072720
# 25%   0.010400212 0.009766708
# 50%   0.018660350 0.009995822
# 75%   0.025347307 0.010282267
# 100%  0.040481854 0.011272515

# 21/08
# bayes         freq
# 0%   0.001019433 -0.002693500
# 25%  0.007889885  0.007547016
# 50%  0.010408897  0.010039148
# 75%  0.011995450  0.013721754
# 100% 0.020547748  0.021663729








# Based on sim.data, It should be that RR = beta*X ??=?? pevent/p0= cases_nonreferent/cases_referent (since n0=n1)

cbind(RR1=sim.data$cases[sim.data$dose!=0]/rep(sim.data$cases[sim.data$dose==0],each=2),RR2=exp(sim.data$logRR[sim.data$dose!=0]))
## Because we round cdose we dont get the exactly the same results

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



