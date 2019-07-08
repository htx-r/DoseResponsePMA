OneRunSimulateDRlinear <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRlineardata.fun(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))
  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  f <-coef(linearDRmetaFreq)[1]
  b <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
return(cbind(Bayes=b,Freq=f))

}

OneRunSimulateDRlinear <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRlineardata.fun(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))
  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  f <-coef(linearDRmetaFreq)[1]
  b <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
  return(cbind(Bayes=b,Freq=f))

}

MultiRunSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

# Repeat the simulation nrep times
res <- replicate(nrep,OneRunSimulateDRlinear(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = F)
res.mat <- do.call(rbind,res)

# Biases
biasB <- beta.pooled-colMeans(res.mat)[1]
biasF <- beta.pooled-colMeans(res.mat)[2]

# standard error
se <- tau/nrep

# Mean square error
mseB <- se^2 + biasB^2
mseF <- se^2 + biasF^2

ret.obj <- cbind(biasB=biasB,biasF=biasF,mseB=mseB,mseF=mseF,beta.pooled=beta.pooled,tau=tau,est.mean.beta.pooledB=colMeans(res.mat)[1],est.mean.beta.pooledF=colMeans(res.mat)[2])
row.names(ret.obj) <- NULL
return(ret.obj)
  }
