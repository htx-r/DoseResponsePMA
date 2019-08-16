### In this file I run jags and dosresmeta for the simulated data where the dose-response is on the OR scale assuming: cubic spline, linear and quadratic.

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines
# compare dosresmeta vs Bayes


n.sim.data <- 100
ns <- 20
beta1.pooled=0.03
beta2.pooled=0.05
bayesCoef1OR<-bayesCoef2OR <- freqCoef1OR<- freqCoef2OR<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedataOR.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logOR~dose1+dose2, id = Study_No,type=type,
                                   se = selogOR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')


  # Bayes normal
  jagsdataRCS<- makeJAGSDRmeta(Study_No,logOR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(1,10))
  ## ----!!!!!---- PLEASE PUT THE COMMAND OF THE PRECISION IN THE makeJAGSDRmeta FUNCTION. And as we agreed you need to find out why you calculated it wrong-------------
  jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)

  bayesCoef1OR <- c(bayesCoef1OR,rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled)
  bayesCoef2OR <- c(bayesCoef2OR,rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled)
  freqCoef1OR<-c(freqCoef1OR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2OR<-c(freqCoef2OR,coef(rcsplineDRmetaFreq)[2])
}


(mean(bayesCoef1OR)-beta1.pooled)/beta1.pooled ## 20/07 bias = -0.0003760008
(mean(bayesCoef2OR)-beta2.pooled)/beta2.pooled ## 20/07 bias = 7.249017e-05

(mean(freqCoef1OR)-beta1.pooled)/beta1.pooled # 20/07 bias = -2.985656e-06
(mean(freqCoef2OR)-beta2.pooled)/beta2.pooled # 20/07  bias = 3.118053e-05


cbind(bayes1=quantile(bayesCoef1OR), freq1=quantile(freqCoef1OR))

#20/07  bayes1      freq1
# 0%   0.02178371 0.02918613
# 25%  0.02780766 0.02981357
# 50%  0.02974700 0.02999414
# 75%  0.03118109 0.03015892
# 100% 0.03587717 0.03101574
cbind(bayes2=quantile(bayesCoef2OR), freq2=quantile(freqCoef2OR))
# 20/07  bayes2      freq2
# 0%   0.03398417 0.04831166
# 25%  0.04614443 0.04957128
# 50%  0.04991107 0.04989353
# 75%  0.05470165 0.05053147
# 100% 0.06738744 0.05183385











# Linear
# compare dosresmeta vs. Bayes

bayesCoefOR <-freqCoefOR <- c()
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

  bayesCoefOR <- c(bayesCoefOR,linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled)
  freqCoefOR<-c(freqCoefOR,coef(linearDRmetaFreq))


}

mean(bayesCoefOR)-beta.pooled # 20/07: bias=-0.0002572297
mean(freqCoefOR)-beta.pooled  # 20/07: bias=9.161834e-05

cbind(bayes=quantile(bayesCoefOR),freq=quantile(freqCoefOR))

#20/07  bayes       freq
# 0%   0.04830245 0.04901056
# 25%  0.04941113 0.04985690
# 50%  0.04975684 0.05009685
# 75%  0.05008100 0.05031505
# 100% 0.05118960 0.05110579














# biasF1<- (c(beta1.pooled, beta2.pooled)-coef1)
# biasBNorm<- (c(beta1.pooled, beta2.pooled)-coef2)
# BIAS=rbind(Freq=apply(biasF1,1,mean),
#            BNorm=apply(biasBNorm,1,mean))
# relativeBIAS<-round(BIAS/c(beta1.pooled, beta2.pooled)*100,2)
