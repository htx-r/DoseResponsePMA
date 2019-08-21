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
beta2.pooled <- 0.02
tau <- 0.001
n.sim.data <- 100
for (i in 1:n.sim.data) {
  # Simulated data
  sim <- simulateDRsplinedataOR.fun(beta1.pooled,beta2.pooled,tau=tau,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots

  # Bayes
  jagsdataSplineBinOR <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,noncases=noncases,data=sim.data,Splines=T)
#inits <- function(){list('beta1.pooled'=0.01,'beta2.pooled'=0.03)}#,'tau'=0.00001)}
#jagsdataSplineBinOR$prec.beta <- 1/tau^2
  splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBinOR,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmetaOR,
                                            n.chains=2,n.iter = 100000,n.burnin =5000,DIC=F,n.thin = 1)

  bayesCoef1OR <- c(bayesCoef1OR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled)
  bayesCoef2OR <- c(bayesCoef2OR,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled)
   traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varname='beta1.pooled')
    traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varname='beta2.pooled')
    traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varname='tau')

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logOR~rcs(sim.data$dose1,knots), id = Study_No,type=type,
                                   se = selogOR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',covariance = 'gl')

  freqCoef1OR<-c(freqCoef1OR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2OR<-c(freqCoef2OR,coef(rcsplineDRmetaFreq)[2])

}

#COMPARING BIAS
cbind(biasBayesCoef1=(mean(bayesCoef1OR-beta1.pooled)),
      biasFreqCoef1=(mean(freqCoef1OR)-beta1.pooled),
      biasBayesCoef2=(mean(bayesCoef2OR-beta2.pooled)),
      biasFreqCoef2=(mean(freqCoef2OR)-beta2.pooled))
## I DO NOT SEE THEM BEING BIASED, AT LEAST NOT MUCH

# %%% TASNIM (21/08/2019): Do the biases need to be closer to zero or that reasonable values for biases?
    # However I checked the traceplot for some of the simulations and it does not look good!

#biasBayesCoef1 biasFreqCoef1 biasBayesCoef2 biasFreqCoef2
#[1,]    0.002934424   0.003401241   -0.006126321   0.009240196


#COMPARING RELATIVE BIAS
cbind(biasBayesCoef1=(mean(bayesCoef1OR-beta1.pooled))/beta1.pooled,
      biasFreqCoef1=(mean(freqCoef1OR)-beta1.pooled)/beta1.pooled,
      biasBayesCoef2=(mean(bayesCoef2OR-beta2.pooled))/beta2.pooled,
      biasFreqCoef2=(mean(freqCoef2OR)-beta2.pooled)/beta2.pooled)

#biasBayesCoef1 biasFreqCoef1 biasBayesCoef2 biasFreqCoef2
#[1,]     0.09781414     0.1133747     -0.1225264     0.1848039

cbind(biasBayesCoef1=(mean(bayesCoef1OR-beta1.pooled)),
      biasFreqCoef1=(mean(freqCoef1OR)-beta1.pooled),
      biasBayesCoef2=(mean(bayesCoef2OR-beta2.pooled)),
      biasFreqCoef2=(mean(freqCoef2OR)-beta2.pooled))
## I DO NOT SEE THEM BEING BIASED, AT LEAST NOT MUCH

cbind(bayes1=quantile(bayesCoef1OR), freq1=quantile(freqCoef1OR))
# bayes1      freq1
# 0%   -0.052715682 0.02930667
# 25%  -0.015105731 0.02982458
# 50%   0.003872969 0.02998341
# 75%   0.018825339 0.03016611
# 100%  0.062772456 0.03074494
cbind(bayes2=quantile(bayesCoef2OR), freq2=quantile(freqCoef2OR))
# bayes2      freq2
# 0%   -0.229376824 0.04811173
# 25%  -0.030164388 0.04962202
# 50%  -0.003988628 0.05000004
# 75%   0.036190734 0.05057417
# 100%  0.140007732 0.05241626

# ------
## Bayes is highly biased. I think the problem in the convergence not in the model itself;
  # the traceplotfor beta1.pooled, beta2.pooled do not converge (bad mixing +  not stationary)!
  # I tried also to fix tau -> it does not get better
  # I tried to increase iterations and burn in -> no change
  # I initialize with true values ->
