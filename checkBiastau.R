library(R2jags)
library(dosresmeta)
library(mixmeta)

library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)


## To check the estimate of tau in freq settings for many simulatuions

freqfun <- function(tau,OR=FALSE){

sim.data <- simulateDRmeta.fun(beta1.pooled=0.01,beta2.pooled=NULL,tau=tau,ns=20,doserange=c(1, 10),samplesize=200,OR=OR,splines = FALSE)

tauHatFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                             se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')$Psi

 # tauHat.mixmeta <- mixmeta(logrr ~ 0 + dose1, S=selogrr^2, random= ~ 0 + dose1|Study_No, data=sim.data,
 #                         subset=!is.na(selogrr))$Psi
 return(sqrt(tauHatFreq))
 }
nrep <- 100
tau <- c(0.06,0.05,0.02,0.01,0.005,0.001)
tauhat1 <- mean(replicate(nrep,freqfun(tau = tau[1])))
tauhat2 <- mean(replicate(nrep,freqfun(tau = tau[2])))
tauhat3 <- mean(replicate(nrep,freqfun(tau = tau[3])))
tauhat4 <- mean(replicate(nrep,freqfun(tau = tau[4])))
tauhat5 <- mean(replicate(nrep,freqfun(tau = tau[5])))
tauhat6 <- mean(replicate(nrep,freqfun(tau = tau[6])))

tauhat1OR <- mean(replicate(nrep,freqfun(tau = tau[1],OR=TRUE)))
tauhat2OR <- mean(replicate(nrep,freqfun(tau = tau[2],OR=TRUE)))
tauhat3OR <- mean(replicate(nrep,freqfun(tau = tau[3],OR=TRUE)))
tauhat4OR <- mean(replicate(nrep,freqfun(tau = tau[4],OR=TRUE)))
tauhat5OR <- mean(replicate(nrep,freqfun(tau = tau[5],OR=TRUE)))
tauhat6OR <- mean(replicate(nrep,freqfun(tau = tau[6],OR=TRUE)))
cbind(tau,tauhat=c(tauhat1,tauhat2,tauhat3,tauhat4,tauhat5,tauhat6),tauhat=c(tauhat1OR,tauhat2OR,tauhat3OR,tauhat4OR,tauhat5OR,tauhat6OR))






#######################
# check bias in tau
#######################
source('FunctionsForSimulations.R')

nrep <- 100
tau <- c(0.06,0.03,0.01,0.005,0.001)
beta.pooled <- c(0.03, 0.1)

# small beta
S1RRlinear <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[1],OR=FALSE,splines=FALSE)

S2RRlinear <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[2],OR=FALSE,splines=FALSE)

S3RRlinear <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[3],OR=FALSE,splines=FALSE)

S4RRlinear <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[4],OR=FALSE,splines=FALSE)

S5RRlinear <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[5],OR=FALSE,splines=FALSE)

# larger beta
S1RRlinear2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[1],OR=FALSE,splines=FALSE)

S2RRlinear2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[2],OR=FALSE,splines=FALSE)

S3RRlinear2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[3],OR=FALSE,splines=FALSE)

S4RRlinear2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[4],OR=FALSE,splines=FALSE)

S5RRlinear2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[5],OR=FALSE,splines=FALSE)

resRRlinearInfprior <- rbind(S1RRlinear,S2RRlinear,S3RRlinear,S4RRlinear,S5RRlinear,
                     S1RRlinear2,S2RRlinear2,S3RRlinear2,S4RRlinear2,S5RRlinear2)
write.csv(resRRlinearInfprior,file=paste0(Sys.Date(),'resRRlinearInfprior.csv')) # keeps the rownames




# small beta
S1RRlinearN <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[1],OR=FALSE,splines=FALSE)

S2RRlinearN <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[2],OR=FALSE,splines=FALSE)

S3RRlinearN <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[3],OR=FALSE,splines=FALSE)

S4RRlinearN <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[4],OR=FALSE,splines=FALSE)

S5RRlinearN <- simpower(nrep,beta1.pooled=beta.pooled[1],tau=tau[5],OR=FALSE,splines=FALSE)

# larger beta
S1RRlinearN2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[1],OR=FALSE,splines=FALSE)

S2RRlinearN2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[2],OR=FALSE,splines=FALSE)

S3RRlinearN2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[3],OR=FALSE,splines=FALSE)

S4RRlinearN2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[4],OR=FALSE,splines=FALSE)

S5RRlinearN2 <- simpower(nrep,beta1.pooled=beta.pooled[2],tau=tau[5],OR=FALSE,splines=FALSE)

resRRlinearNonInfprior <- rbind(S1RRlinearN,S2RRlinearN,S3RRlinearN,S4RRlinearN,S5RRlinearN,
                     S1RRlinearN2,S2RRlinearN2,S3RRlinearN2,S4RRlinearN2,S5RRlinearN2)
write.csv(resRRlinearNonInfprior,file=paste0(Sys.Date(),'resRRlinearNonInfprior.csv')) # keeps the rownames
# end for linear RR model


















######

beta1.pooled=0.02
beta2.pooled=0
tau=0.001
ns=20
doserange=c(1, 10)
samplesize=200
OR=TRUE
splines = TRUE
simulateDRmeta.fun()
sim.data <- simulateDRmeta.fun(beta1.pooled=0.3,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE)


library(lme4)


fm1 <- lmer(logrr ~ 0 + dose1+ (dose1+0 | Study_No), sim.data)
sqrt(summary(fm1)$vcov)


