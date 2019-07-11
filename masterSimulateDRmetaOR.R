### In this file I run jags and dosresmeta for the simulated data assuming: cubic spline, linear and quadratic.

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines
# compare dosresmeta vs Bayes


n.sim.data <- 100
ns <- 20
start <- Sys.time()
beta1.pooled=0.03
beta2.pooled=0.05
coef1<-coef2<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedataOR.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logOR~dose1+dose2, id = Study_No,type=type,
                                   se = selogOR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  #rcsplineDRmetaFreq2 <- dosresmeta(formula = logOR~rcs(dose1,sim.data$knots), id = Study_No,type=type,
   #                                 se = selogOR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  coef1<-cbind(coef1,c(coef(rcsplineDRmetaFreq)[1],coef(rcsplineDRmetaFreq)[2]))
  #coef2<-cbind(coef2,c(coef(rcsplineDRmetaFreq2)[1],coef(rcsplineDRmetaFreq2)[2]))

  # Bayes normal
  jagsdataRCS<- makeJAGSDRmeta(Study_No,logOR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(1,10))
  ## ----!!!!!---- PLEASE PUT THE COMMAND OF THE PRECISION IN THE makeJAGSDRmeta FUNCTION. And as we agreed you need to find out why you calculated it wrong-------------
  jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
  coef2<-cbind(coef2, c(rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled,rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled))

#   #Bayes binomial
#   jagsdataSplineBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,controls=noncases,data=sim.data$simulatedDRdata,Splines=T)
#
#   splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmeta,
#                                             n.chains=1,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)
#   coef4<-cbind(coef4, c(splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled,splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled))

}
end <- Sys.time()
(end-start)


biasF1<- (c(beta1.pooled, beta2.pooled)-coef1)
biasBNorm<- (c(beta1.pooled, beta2.pooled)-coef2)
BIAS=rbind(Freq=apply(biasF1,1,mean),
           BNorm=apply(biasBNorm,1,mean))
relativeBIAS<-round(BIAS/c(beta1.pooled, beta2.pooled)*100,2)
# the bayesian binomial model has the largest bias, which is weird... we need to investigate this better


# the results from 500 simulations are
#> relativeBIAS
# dose1 dose2
# Freq  -0.01 -0.30
# BNorm  0.65 -0.35
#> BIAS
#         dose1         dose2
# Freq  -3.692444e-06 -9.089506e-05
# BNorm  3.225625e-04 -1.727547e-04






# Linear
# compare dosresmeta vs. Bayes

b1 <-f1 <- vector()
n.sim.data <- 100
beta.pooled = 0.05
start <- Sys.time()
for (j in 1:n.sim.data) {
  sim.data <- simulateDRlineardataOR.fun(beta.pooled)

  # Freq
  linearDRmetaFreq <- dosresmeta(formula = logOR~dose, id = Study_No,type=type,
                              se = selogOR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logOR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))

  #####!!!!!!!! I THOUGHT YOU SAID YOU HAVE TO USE THE CORRECT MATRIX FROM THE FREQ MODEL!!!################

  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  f1[j] <-coef(linearDRmetaFreq)[1]
  b1[j] <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled


}
end <- Sys.time()
(end-start)
biasb<- quantile(beta.pooled-b1) #the Bayesian is biased- wrong?
biasf<- quantile(beta.pooled-f1) #the frequentist is unbiased
Linearres.100 <- cbind(freq=biasf,bayes=biasb)


#            freq         bayes
# 0%   -9.624503e-04 -0.0007143604
# 25%  -3.277294e-04 -0.0001535071
# 50%  -8.560753e-05  0.0001277164
# 75%   2.022103e-04  0.0004545837
# 100%  1.019967e-03  0.0011979086
