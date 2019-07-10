OneRunSimulateDRlinear <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRlineardata.fun(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  linearDRmetaFreq<-dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))
  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  f <-coef(linearDRmetaFreq)[1]
  b <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
  t <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
return(cbind(Bayes=b,Freq=f,tau=t))

}


MultiRunSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

# Repeat the simulation nrep times
res <- replicate(nrep,OneRunSimulateDRlinear(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = F)
res.mat <- do.call(rbind,res)

# Biases
biasB <- beta.pooled-colMeans(res.mat)[1]
biasF <- beta.pooled-colMeans(res.mat)[2]

# standard error
tau.hat <- colMeans(res.mat)[3]
se <- tau.hat/nrep

# Mean square error
mseB <- se^2 + biasB^2
mseF <- se^2 + biasF^2

ret.obj <- cbind(beta.pooled=beta.pooled,tau=tau,beta.pooled.hatB=colMeans(res.mat)[1],beta.pooled.hatF=colMeans(res.mat)[2],tau.hat=tau.hat,biasB=biasB,biasF=biasF,mseB=mseB,mseF=mseF)
row.names(ret.obj) <- NULL
return(ret.obj)
  }


#MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta.pooled <- rep(c(0,0.02,0.03,0.05),2)
tau <- rep(c(0.001,0.05),each=4)
ns <- 20
doserange <- c(0,10)
samplesize<- 200
nrep <- 100
res <- mapply(MultiRunSimulateDRlinear, beta.pooled=beta.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
result <- do.call(rbind,res)
rownames(result) <- paste0('Scenario ',1:nrow(result))
result

#}
## Run all scenarios
