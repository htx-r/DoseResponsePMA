###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Restricted Cubic Spline
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# There are two functions:
# OneRunSimulateDRspline: run the simulation for one set of arguments and for one simulation
# MultiRunSimulateDRspline: run the simulation for one set of arguments and for 'nrep' simulations

# And to run the simulation for different arguments, I apply the MultiRunSimulateDRspline() for the different scenarios
     # using mapply (we can make it as a function if Georgia would like)
OneRunSimulateDRspline <- function(beta1.pooled=0.03,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRsplinedata.fun(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!

  # Bayes
  jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(5,10))
  #jagsdataRCS$prec.beta <- 1/(0.001)^2 I tried a fixed precision but it gives similar result
  jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)

  f1 <-coef(rcsplineDRmetaFreq)[1]
  b1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

  f2 <-coef(rcsplineDRmetaFreq)[2]
  b2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled

  t1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau1
  t2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau2

  return(cbind(Bayes=c(b1,b2),Freq=c(f1,f2),tau.hat=c(t1,t2)))

}



MultiRunSimulateDRspline <- function(nrep=3,beta1.pooled=0.02,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneRunSimulateDRspline(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = F)
  res.mat1 <- t(sapply(1:nrep, function(i) res[[i]][1,]))
  res.mat2 <- t(sapply(1:nrep, function(i) res[[i]][2,]))


  # Biases
  biasB1 <- beta1.pooled-colMeans(res.mat1)[1]
  biasF1 <- beta1.pooled-colMeans(res.mat1)[2]

  biasB2 <- beta2.pooled-colMeans(res.mat2)[1]
  biasF2 <- beta2.pooled-colMeans(res.mat2)[2]


  # standard error
  tau1.hat <- colMeans(res.mat1)[3]
  se1 <- tau1.hat/nrep

  tau2.hat <- colMeans(res.mat2)[3]
  se2 <- tau2.hat/nrep
  # Mean square error
  mseB1 <- se1^2 + biasB1^2
  mseF1 <- se1^2 + biasF1^2

  mseB2 <- se2^2 + biasB2^2
  mseF2 <- se2^2 + biasF2^2

  ret.obj <- list(sum.coef1=cbind(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,beta1.pooled.hatB=colMeans(res.mat1)[1],beta1.pooled.hatF=colMeans(res.mat1)[2],tau1.hat=tau1.hat,biasB1=biasB1,biasF1=biasF1,mseB1=mseB1,mseF1=mseF1),
                  sum.coef2=cbind(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,beta2.pooled.hatB=colMeans(res.mat2)[1],beta2.pooled.hatF=colMeans(res.mat2)[2],tau2.hat=tau2.hat,biasB2=biasB2,biasF2=biasF2,mseB2=mseB2,mseF2=mseF2))
  row.names(ret.obj[[1]]) <-row.names(ret.obj[[2]]) <- NULL
  return(ret.obj)
}

# Arguments
beta1.pooled <- rep(c(0,0.02,0.03,0.05,0.03),2)
beta2.pooled <- rep(c(0,0,0.02,0.03,0.05),2)
tau <- rep(c(0.001,0.01),each=6)
ns <- 20
doserange <- c(0,10)
samplesize<- 200
nrep <- 100

# Results
#MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
res <- mapply(MultiRunSimulateDRspline, beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
res.mat1 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef1[1,]))
res.mat2 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef2[1,]))
rownames(res.mat1) <- paste0('Scenario ',1:nrow(res.mat1))
rownames(res.mat2) <- paste0('Scenario ',1:nrow(res.mat2))
res.mat1
res.mat2
#return(list(sum.coef1=res.mat1,sum.coef2=res.mat2))
#}
