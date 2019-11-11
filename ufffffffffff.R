source('FunctionsForSimulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)

nsim <- 1000
beta1.pooled <- c(0,0.04,0.1,0.2,0.2)
beta2.pooled <- c(0,0,0.03,-0.2,-0.3 )
tau <- c(0.001,0.01)
ns <- 40


beta1.pooled=0.1
beta2.pooled=0.03
tau=0.001
ns=40
doserange=c(1, 10)
samplesize=200
OR=TRUE
splines = TRUE

traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varname='beta1.pooled')

traceplot(splineDRmetaJAGSmodelBin$BUGSoutput,varname='beta2.pooled')
splineDRmetaJAGSmodelBin$BUGSoutput$mean

### 1. odds ratio (OR)
## %% smaller tau
# Scenario 1
set.seed('122')
S1ORspline <- OneSimulation(beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
S2linear <- OneSimulation(beta1.pooled = beta1.pooled[1],beta2.pooled = NULL,tau=tau[1],OR=TRUE,ns=ns,splines = FALSE)

#
rval <- list()

for (i in 1:1000) {
myd <- simulateDRmeta.fun(beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled,tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
myd$p <- myd$cases/(myd$noncases+myd$cases)
myd$t <- ifelse(myd$dose1==0,0,1)
rval[[i]] <- coef(glm(p~dose1+dose2-1+t,data = myd,family = binomial(link = 'logit')))[1:2]

}

colMeans( do.call(rbind,rval))


rval <- list()

for (i in 1:100) {
  sim.data <- simulateDRmeta.fun(beta1.pooled = beta1.pooled[1],beta2.pooled = NULL,tau=tau[1],OR=TRUE,ns=ns,splines = FALSE)
jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2=NULL,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=F,new.dose.range = c(5,10))
rval[[i]] <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                          n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 1)$BUGSoutput$mean
}

(rval[[100]])
colMeans(matrix(unlist(rval),100,2,byrow = T))



