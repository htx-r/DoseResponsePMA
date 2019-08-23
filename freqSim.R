library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

################################
# Relative Risk (RR)
################################


## RR: Cubic Splines

n.sim.data <- 500
ns <- 20
start <- Sys.time()
beta1.pooled=0.02
beta2.pooled=0.01
freqCoef1RR<-freqCoef2RR<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,12))

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!

  freqCoef1RR<-c(freqCoef1RR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2RR<-c(freqCoef2RR,coef(rcsplineDRmetaFreq)[2])
}


(mean(freqCoef1RR)-beta1.pooled) # 23/08 bias =0.0001154705  -0.0004600168 (n.sim=500)
(mean(freqCoef2RR)-beta2.pooled) # 23/08 bias =-0.0001067644  0.001292083 (n.sim=500)


cbind(freq1=quantile(freqCoef1RR), freq2=quantile(freqCoef2RR))

# 23/8 freq1        freq2
# 0%   -0.007075189 -0.066064181
# 25%   0.015518117 -0.002794734
# 50%   0.020226890  0.009464867
# 75%   0.024572409  0.022630662
# 100%  0.045103164  0.094466002





## RR: Linear


freqCoefRR <- vector()
n.sim.data <- 500
beta.pooled = 0.02
tau <- 0.01
for (j in 1:n.sim.data) {
  sim.data <- simulateDRlineardata.fun(beta.pooled,tau)

  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  freqCoefRR<-c(freqCoefRR,coef(linearDRmetaFreq))

}

mean(freqCoefRR)-beta.pooled # 23/08: bias = -0.0001926179  -9.631433e-05 true =0.02
quantile(freqCoefRR)
sd(freqCoefRR)^2
library(MASS)
truehist(freqCoefRR)
# 0%         25%         50%         75%        100%
# 0.003362988 0.015926584 0.019900525 0.023634805 0.034700581

# 0%         25%         50%         75%        100%
# 0.004374584 0.015977389 0.019984619 0.023866652 0.036444412


################################
# Odds ratio (OR)
################################


## OR: Cubic Spline

n.sim.data <- 500
ns <- 20
beta1.pooled=0.02
beta2.pooled=0.01
freqCoef1OR<- freqCoef2OR<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedataOR.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,12))

  rcsplineDRmetaFreq <- dosresmeta(formula = logOR~dose1+dose2, id = Study_No,type=type,
                                   se = selogOR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')

  freqCoef1OR<-c(freqCoef1OR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2OR<-c(freqCoef2OR,coef(rcsplineDRmetaFreq)[2])
}


(mean(freqCoef1OR)-beta1.pooled)# 23/08 bias = 0.009246271 (n.sim=100) 0.005380327 (n.sim=500)
(mean(freqCoef2OR)-beta2.pooled) # 23/08  bias = -0.01659826 (n.sim=100) -0.01336922 (n.sim=500)

cbind(freq2=quantile(freqCoef2OR), freq1=quantile(freqCoef1OR))

# freq2       freq1
# 0%   -0.164943217 -0.02230373
# 25%  -0.041447550  0.01172301
# 50%  -0.006205263  0.02650169
# 75%   0.037352265  0.04390756
# 100%  0.152434301  0.10246233

# freq2       freq1
# 0%   -0.393633904 -0.04335618
# 25%  -0.035959631  0.01087476
# 50%  -0.001183256  0.02529229
# 75%   0.033633162  0.03808721
# 100%  0.167295848  0.09950265


## OR: Linear

freqCoefOR <- c()
n.sim.data <- 500
beta.pooled = 0.05
for (j in 1:n.sim.data) {
  sim.data <- simulateDRlineardataOR.fun(beta.pooled)

  # Freq
  linearDRmetaFreq <- dosresmeta(formula = logOR~dose, id = Study_No,type=type,
                                 se = selogOR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  freqCoefOR<-c(freqCoefOR,coef(linearDRmetaFreq))


}

mean(freqCoefOR)-beta.pooled  # 23/08: bias=0.003339908,  0.003351409

quantile(freqCoefOR)
# 0%        25%        50%        75%       100%
# 0.02315905 0.04607257 0.05262326 0.06102271 0.08836876

# 0%        25%        50%        75%       100%
# 0.01909019 0.04513631 0.05228773 0.06102233 0.09037502

# End
