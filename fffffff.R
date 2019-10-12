library(R2jags)
library(dosresmeta)
library(mixmeta)

library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)


## To check the estimate of tau in freq settings for many simulatuions

freqfun <- function(tau,OR=FALSE){

sim.data <- simulateDRmeta.fun(beta1.pooled=0.01,beta2.pooled=NULL,tau=tau,ns=20,doserange=c(1, 10),samplesize=200,OR=OR,splines = FALSE)

tauHatFreq<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                             se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')$Psi

 # tauHat.mixmeta <- mixmeta(logrr ~ 0 + dose1, S=selogrr^2, random= ~ 0 + dose1|Study_No, data=sim.data,
 #                         subset=!is.na(selogrr))$Psi
 return(sqrt(tauHatFreq))
 }
nrep <- 1000
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


beta1.pooled=0.3
beta2.pooled=NULL
tau=0.001
ns=20
doserange=c(1, 10)
samplesize=200
OR=FALSE
splines = FALSE

sim.data <- simulateDRmeta.fun(beta1.pooled=0.3,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE)


library(lme4)


fm1 <- lmer(logrr ~ 0 + dose1+ (dose1+0 | Study_No), sim.data)
sqrt(summary(fm1)$vcov)


################3
# check if we fix tau, did we get less biased beta?
##############3

jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))





