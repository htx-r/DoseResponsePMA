
## To check the estimate of tau in freq settings
sim.data <- list()
linearDRmetaFreq <- vector()
mixmetaModel <- vector()
for (i in 1:100) {


sim.data[[i]] <- simulateDRmeta.fun(beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE)

linearDRmetaFreq[i]<-dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
                             se = selogrr, cases = cases, n = cases+noncases, data = sim.data[[i]], proc='2stage',method = 'fixed',covariance = 'gl')$Psi
mixmetaModel[i] <- mixmeta(logrr ~ 0 + dose1, S=selogrr^2, random= ~ 0 + dose1|Study_No, data=sim.data,
                        subset=!is.na(selogrr))$Psi
}
0.001^2-mean(linearDRmetaFreq)
0.001-mean(sqrt(linearDRmetaFreq))


beta1.pooled=0.3
beta2.pooled=NULL
tau=0.001
ns=20
doserange=c(1, 10)
samplesize=200
OR=FALSE
splines = FALSE

sim.data <- simulateDRmeta.fun(beta1.pooled=0.3,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE)
#sim.data$cases[is.na(sim.data$selogrr)]
dosresModel <- dosresmeta(formula = logrr~dose1, id = Study_No,type=type,
           se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',covariance = 'gl')
summary(dosresModel)
sqrt(summary(dosresModel)$Psi)

mixmetaModel <- mixmeta(logrr ~ 0 + dose1, S=selogrr^2, random= ~ 0 + dose1|Study_No, data=sim.data,
                   subset=!is.na(selogrr))
summary(mixmetaModel)
sqrt(summary(mixmetaModel)$Psi)

library(mixmeta)


fm1 <- lmer(logrr ~ 0 + dose1+ (dose1+0 | Study_No), sim.data)
sqrt(summary(fm1)$vcov)







