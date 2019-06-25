library(R2jags)
sim.data <- simulateDRlineardata.fun(beta.pooled = 0.01)
jagsdataLinear2 <- make2JAGSDRmeta(Study_No,logRR,dose,cases,noncases,data=sim.data,Splines=F)

## If the simulated data is correct then cases/noncases should equals exp(logRR) because the sample size is the same
# But they are not equal, because of that I got wrong

# Anaylsis:
# a. Bayes: JAGS model
linearDRmetaJAGSmodel <- jags.parallel(data = jagsdataLinear2,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = model2LinearDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)

lm(logRR~dose-1,data=sim.data)
linearDRmetaJAGSmodel$BUGSoutput$mean



exp(0.01*sim.data$dose)

linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
linearDRmetaJAGSmodel$BUGSoutput$mean$tau

traceplot(linearDRmetaJAGSmodel$BUGSoutput,varname='beta.pooled')

lm(simulatedDRdata[simulatedDRdata$dose!=0,]$cases/simulatedDRdata[simulatedDRdata$dose==0,]$cases
~simulatedDRdata$dose[simulatedDRdata$dose!=0,]-1)


linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ dose, type = type, id = Study_No,
                                           se = selogRR, cases = cases, n = cases+noncases  , data = sim.data,covariance = 'gl',proc = '2stage',method = 'fixed')#!!!!!!!!!!!!!!

coef(linearDRmetaFreq)

